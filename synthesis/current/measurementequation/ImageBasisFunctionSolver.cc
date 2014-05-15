/// @copyright (c) 2007 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
#include <boost/shared_ptr.hpp>

ASKAP_LOGGER(logger, ".measurementequation.imagebasisfunctionsolver");

#include <askap/AskapError.h>

// need it just for null deleter
#include <askap/AskapUtil.h>

#include <utils/MultiDimArrayPlaneIter.h>
#include <profile/AskapProfiler.h>

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/Vector.h>

#include <casa/BasicSL/Complex.h>

#include <askap_synthesis.h>
#include <measurementequation/ImageBasisFunctionSolver.h>
#include <deconvolution/DeconvolverBasisFunction.h>
#include <lattices/Lattices/ArrayLattice.h>

using namespace casa;
using namespace askap;
using namespace askap::scimath;

#include <iostream>

#include <cmath>
using std::abs;

#include <map>
#include <vector>
#include <string>

using std::map;
using std::vector;
using std::string;

namespace askap
{
  namespace synthesis
  {
    ImageBasisFunctionSolver::ImageBasisFunctionSolver()
    {
      // Now set up controller
      itsControl = boost::shared_ptr<DeconvolverControl<Float> >(new DeconvolverControl<Float>());
      // Now set up monitor
      itsMonitor = boost::shared_ptr<DeconvolverMonitor<Float> >(new DeconvolverMonitor<Float>());

      // Make the basis function
      std::vector<float> defaultScales(3);
      defaultScales[0]=0.0;
      defaultScales[1]=10.0;
      defaultScales[2]=30.0;
      itsBasisFunction=BasisFunction<Float>::ShPtr(new MultiScaleBasisFunction<Float>(defaultScales));
    }
    
    ImageBasisFunctionSolver::ImageBasisFunctionSolver(casa::Vector<float>& scales)
    {
      // Now set up controller
      itsControl = boost::shared_ptr<DeconvolverControl<Float> >(new DeconvolverControl<Float>());
      // Now set up monitor
      itsMonitor = boost::shared_ptr<DeconvolverMonitor<Float> >(new DeconvolverMonitor<Float>());

      itsBasisFunction=BasisFunction<Float>::ShPtr(new MultiScaleBasisFunction<Float>(scales));
    }
    
    void ImageBasisFunctionSolver::setBasisFunction(BasisFunction<Float>::ShPtr bf) {
      itsBasisFunction=bf;
    }

    BasisFunction<Float>::ShPtr ImageBasisFunctionSolver::basisFunction() {
      return itsBasisFunction;
    }

    void ImageBasisFunctionSolver::configure(const LOFAR::ParameterSet &parset) {

      ImageSolver::configure(parset);

      ASKAPASSERT(this->itsMonitor);
      this->itsMonitor->configure(parset);
      ASKAPASSERT(this->itsControl);
      this->itsControl->configure(parset);
    }
    
    void ImageBasisFunctionSolver::init()
    {
      resetNormalEquations();
    }
    
    /// @brief Solve for parameters, updating the values kept internally
    /// The solution is constructed from the normal equations. The parameters named 
    /// image* are interpreted as images and solved for.
    /// @param[in] ip current model (to be updated)
    /// @param[in] quality Solution quality information
    bool ImageBasisFunctionSolver::solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality)
    {
      ASKAPTRACE("ImageBasisFunctionSolver::solveNormalEquations");

      // Solving A^T Q^-1 V = (A^T Q^-1 A) P
      uint nParameters=0;
      
      // Find all the free parameters beginning with image
      vector<string> names(ip.completions("image"));
      map<string, uint> indices;
      
      for (vector<string>::const_iterator  it=names.begin();it!=names.end();it++)
	{
	  const std::string name="image"+*it;
	  if(ip.isFree(name)) {
	    indices[name]=nParameters;
	    nParameters+=ip.value(name).nelements();
	  }
	}
      ASKAPCHECK(nParameters>0, "No free parameters in ImageBasisFunctionSolver");
      
      for (map<string, uint>::const_iterator indit=indices.begin();indit!=indices.end();++indit)
	{
	  // Axes are dof, dof for each parameter
	  //const casa::IPosition vecShape(1, ip.value(indit->first).nelements());
	  for (scimath::MultiDimArrayPlaneIter planeIter(ip.value(indit->first).shape());
	       planeIter.hasMore(); planeIter.next()) {
	    
	    ASKAPCHECK(normalEquations().normalMatrixDiagonal().count(indit->first)>0, "Diagonal not present for "<<
		       indit->first);
	    casa::Vector<double> diag(normalEquations().normalMatrixDiagonal().find(indit->first)->second);
	    ASKAPCHECK(normalEquations().dataVector(indit->first).size()>0, "Data vector not present for "<<
		       indit->first);
	    casa::Vector<double> dv = normalEquations().dataVector(indit->first);
	    ASKAPCHECK(normalEquations().normalMatrixSlice().count(indit->first)>0, "PSF Slice not present for "<<
		       indit->first);
	    casa::Vector<double> slice(normalEquations().normalMatrixSlice().find(indit->first)->second);
	    
	    if (planeIter.tag()!="") {
	      // it is not a single plane case, there is something to report
	      ASKAPLOG_INFO_STR(logger, "Processing plane "<<planeIter.sequenceNumber()<<
				" tagged as "<<planeIter.tag());
	    }
	    
	    casa::Array<float> dirtyArray = padImage(planeIter.getPlane(dv));
	    casa::Array<float> psfArray = padImage(planeIter.getPlane(slice));
	    casa::Array<float> basisFunctionArray = padImage(planeIter.getPlane(ip.value(indit->first)));
	    casa::Array<float> maskArray(dirtyArray.shape());
	    ASKAPLOG_INFO_STR(logger, "Plane shape "<<planeIter.planeShape()<<" becomes "<<
			      dirtyArray.shape()<<" after padding");

	    if(doPreconditioning(psfArray,dirtyArray)) {
	      // Normalize
	      doNormalization(padDiagonal(planeIter.getPlane(diag)),tol(),psfArray,dirtyArray, 
			      boost::shared_ptr<casa::Array<float> >(&maskArray, utility::NullDeleter()));
	      // Store the new PSF in parameter class to be saved to disk later
	      saveArrayIntoParameter(ip, indit->first, planeIter.shape(), "psf.image", unpadImage(psfArray),
				     planeIter.position());
	    }
	    else {
	      // Normalize
	      doNormalization(padDiagonal(planeIter.getPlane(diag)),tol(),psfArray,dirtyArray, 
			      boost::shared_ptr<casa::Array<float> >(&maskArray, utility::NullDeleter()));
	    }
	    
	    // optionally clip the image and psf if there was padding
	    ASKAPLOG_INFO_STR(logger, "Peak data vector flux (derivative) before clipping "<<max(dirtyArray));
	    clipImage(dirtyArray);
	    clipImage(psfArray);
	    ASKAPLOG_INFO_STR(logger, "Peak data vector flux (derivative) after clipping "<<max(dirtyArray));
	    
	    // This takes up some memory and we have to ship the residual image out inside
	    // the parameter class. Therefore, we may not need this functionality in the 
	    // production version (or may need to implement it in a different way).
	    saveArrayIntoParameter(ip, indit->first, planeIter.shape(), "residual",
				   unpadImage(dirtyArray), planeIter.position());
	    
	    // uncomment the code below to save the mask
	    saveArrayIntoParameter(ip, indit->first, planeIter.shape(), "mask", unpadImage(maskArray),
				   planeIter.position());
	    
	    // Startup now costs so little it's better to create a new
	    // deconvolver each time we need it	    
	    boost::shared_ptr<DeconvolverBasisFunction<float, casa::Complex> >
	      basisFunctionDec(new DeconvolverBasisFunction<float,
		    casa::Complex>(dirtyArray, psfArray));
	    ASKAPASSERT(basisFunctionDec);     
	    basisFunctionDec->setMonitor(itsMonitor);
	    basisFunctionDec->setControl(itsControl);
	    basisFunctionDec->setWeight(maskArray);

        casa::Array<float> cleanArray(planeIter.planeShape());
        casa::convertArray<float, double>(cleanArray, planeIter.getPlane(ip.value(indit->first)));
        basisFunctionDec->setModel(cleanArray);

	    itsBasisFunction->initialise(dirtyArray.shape());
	    basisFunctionDec->setBasisFunction(itsBasisFunction);

	    // We have to reset the initial objective function
	    // so that the fractional threshold mechanism will work.
	    basisFunctionDec->state()->resetInitialObjectiveFunction();
	    // By convention, iterations are counted from scratch each
	    // major cycle
	    basisFunctionDec->state()->setCurrentIter(0);
	    
	    basisFunctionDec->control()->setTargetObjectiveFunction(threshold().getValue("Jy"));
	    basisFunctionDec->control()->setFractionalThreshold(fractionalThreshold());

	    ASKAPLOG_INFO_STR(logger, "Starting basis function deconvolution");
	    basisFunctionDec->deconvolve();
	    ASKAPLOG_INFO_STR(logger, "Peak flux of the Basis function image "
			      << max(basisFunctionDec->model()));
	    ASKAPLOG_INFO_STR(logger, "Peak residual of Basis function image "
			      << max(abs(basisFunctionDec->dirty())));
	    
	    const std::string peakResParam = std::string("peak_residual.") + indit->first;
	    if (ip.has(peakResParam)) {
	      ip.update(peakResParam, basisFunctionDec->state()->peakResidual());
	    } else {
	      ip.add(peakResParam, basisFunctionDec->state()->peakResidual());
	    }
	    ip.fix(peakResParam);	    
            planeIter.getPlane(ip.value(indit->first)).nonDegenerate()=unpadImage(basisFunctionDec->model());
	  } // loop over all planes of the image cube
	} // loop over map of indices
      
      quality.setDOF(nParameters);
      quality.setRank(0);
      quality.setCond(0.0);
      quality.setInfo("BasisFunction deconvolver");
      
      /// Save the PSF and Weight
      saveWeights(ip);      
      savePSF(ip);
      
      return true;
    };
    
    Solver::ShPtr ImageBasisFunctionSolver::clone() const
    {
      return Solver::ShPtr(new ImageBasisFunctionSolver(*this));
    }
    
  }
}
