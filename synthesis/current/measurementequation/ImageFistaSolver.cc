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

ASKAP_LOGGER(logger, ".measurementequation.imagefistasolver");

#include <askap/AskapError.h>

// need it just for null deleter
#include <askap/AskapUtil.h>

#include <utils/MultiDimArrayPlaneIter.h>

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/Vector.h>

#include <casa/BasicSL/Complex.h>

#include <askap_synthesis.h>
#include <measurementequation/ImageFistaSolver.h>
#include <deconvolution/DeconvolverFista.h>
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
    ImageFistaSolver::ImageFistaSolver()
    {
      // Now set up controller
      itsControl = boost::shared_ptr<DeconvolverControl<Float> >(new DeconvolverControl<Float>());
      // Now set up monitor
      itsMonitor = boost::shared_ptr<DeconvolverMonitor<Float> >(new DeconvolverMonitor<Float>());
    }
    
    void ImageFistaSolver::configure(const LOFAR::ParameterSet &parset) {
      ASKAPASSERT(this->itsMonitor);
      this->itsMonitor->configure(parset);
      ASKAPASSERT(this->itsControl);
      this->itsControl->configure(parset);

      if(parset.isDefined("scales")) {
	std::vector<float> defaultScales(6);
	defaultScales[0]=0.0;
	defaultScales[1]=2.0;
	defaultScales[2]=4.0;
	defaultScales[2]=8.0;
	defaultScales[2]=16.0;
	defaultScales[2]=32.0;
	std::vector<float> scales=parset.getFloatVector("scales", defaultScales);
	
	ASKAPLOG_INFO_STR(decfistalogger, "Constructing Multiscale basis function with scales " << scales);
        Bool orthogonal=parset.getBool("orthogonal", "false");

	itsBasisFunction = BasisFunction<Float>::ShPtr(new MultiScaleBasisFunction<Float>(scales,
                                                                                          orthogonal));
      }
    }
    
    void ImageFistaSolver::init()
    {
      resetNormalEquations();
    }
    
    void ImageFistaSolver::setBasisFunction(BasisFunction<Float>::ShPtr bf) {
      itsBasisFunction=bf;
    }

    BasisFunction<Float>::ShPtr ImageFistaSolver::basisFunction() {
      return itsBasisFunction;
    }

    /// @brief Solve for parameters, updating the values kept internally
    /// The solution is constructed from the normal equations. The parameters named 
    /// image* are interpreted as images and solved for.
    /// @param[in] ip current model (to be updated)
    /// @param[in] quality Solution quality information
    bool ImageFistaSolver::solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality)
    {
      
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
      ASKAPCHECK(nParameters>0, "No free parameters in ImageFistaSolver");
      
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
	    casa::Array<float> fistaArray = padImage(planeIter.getPlane(ip.value(indit->first)));
	    casa::Array<float> maskArray(dirtyArray.shape());
	    ASKAPLOG_INFO_STR(logger, "Plane shape "<<planeIter.planeShape()<<" becomes "<<
			      dirtyArray.shape()<<" after padding");
	    
	    // Precondition the PSF and DIRTY images before solving.
	    if(doPreconditioning(psfArray,dirtyArray)) {
	      // Normalize	         
	      doNormalization(padDiagonal(planeIter.getPlane(diag)),tol(),psfArray,dirtyArray, 
			      boost::shared_ptr<casa::Array<float> >(&maskArray, utility::NullDeleter()));
	      // Store the new PSF in parameter class to be saved to disk later
	      saveArrayIntoParameter(ip, indit->first, planeIter.shape(), "psf.image", unpadImage(psfArray),
				     planeIter.position());
	    } // if there was preconditioning
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
	    
	    
	    
	    // Startup costs so little it's better to create a new
	    // deconvolver each time we need it	    
	    boost::shared_ptr<DeconvolverFista<float, casa::Complex> >
	      fistaDec(new DeconvolverFista<float, casa::Complex>(dirtyArray, psfArray));
	    ASKAPDEBUGASSERT(fistaDec);     
	    fistaDec->setMonitor(itsMonitor);
	    fistaDec->setControl(itsControl);
	    fistaDec->setWeight(maskArray);
	    
	    if(itsBasisFunction) {
	      itsBasisFunction->initialise(dirtyArray.shape());
	      fistaDec->setBasisFunction(itsBasisFunction);
	    }

	    // We have to reset the initial objective function
	    // so that the fractional threshold mechanism will work.
	    fistaDec->state()->resetInitialObjectiveFunction();
	    // By convention, iterations are counted from scratch each
	    // major cycle
	    fistaDec->state()->setCurrentIter(0);
	    
            fistaDec->setModel(fistaArray);

	    ASKAPLOG_INFO_STR(logger, "Starting FISTA deconvolution");
	    // FISTA is not incremental so we need to set the
	    // background image which remains fixed during one 
	    // deconvolve step
	    fistaDec->deconvolve();
	    ASKAPLOG_INFO_STR(logger, "Peak flux of the FISTA image "
			      << max(fistaDec->model()));
	    ASKAPLOG_INFO_STR(logger, "Peak residual of FISTA image "
			      << max(abs(fistaDec->dirty())));
	    
	    const std::string deconvolverKey = indit->first + planeIter.tag();
	    const std::string peakResParam = std::string("peak_residual.") + deconvolverKey;
	    if (ip.has(peakResParam)) {
	      ip.update(peakResParam, fistaDec->state()->peakResidual());
	    } else {
	      ip.add(peakResParam, fistaDec->state()->peakResidual());
	    }
	    ip.fix(peakResParam);
            planeIter.getPlane(ip.value(indit->first)).nonDegenerate() = unpadImage(fistaDec->model());
	  } // loop over all planes of the image cube
	} // loop over map of indices
      
      quality.setDOF(nParameters);
      quality.setRank(0);
      quality.setCond(0.0);
      quality.setInfo("Fista deconvolver");
      
      /// Save the PSF and Weight
      saveWeights(ip);      
      savePSF(ip);
      
      return true;
    };
    
    Solver::ShPtr ImageFistaSolver::clone() const
    {
      return Solver::ShPtr(new ImageFistaSolver(*this));
    }
    
  }
}
