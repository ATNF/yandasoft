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

#include <measurementequation/ImageMultiScaleSolver.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
#include <boost/shared_ptr.hpp>

ASKAP_LOGGER(logger, ".measurementequation.imagemultiscalesolver");

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

#include <lattices/Lattices/LatticeCleaner.h>
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
    // note, constructor parameter for itsCleaners hard codes the cache size.
    // ideally it should be the number of spectral channels times number of facets
    // processed by every thread.
    
    ImageMultiScaleSolver::ImageMultiScaleSolver() : 
      itsCleaners(8), itsDoSpeedUp(false), itsSpeedUpFactor(1.)
    {
      itsScales.resize(3);
      itsScales(0)=0;
      itsScales(1)=10;
      itsScales(2)=30;
    }
    
    ImageMultiScaleSolver::ImageMultiScaleSolver(const casa::Vector<float>& scales) : 
      itsCleaners(8), itsDoSpeedUp(false), itsSpeedUpFactor(1.)
    {
      itsScales.resize(scales.size());
      itsScales=scales;
    }
    
    void ImageMultiScaleSolver::init()
    {
      resetNormalEquations();
    }
    
    /// @brief switch the speed up on
    /// @param[in] factor speed up factor
    void ImageMultiScaleSolver::setSpeedUp(float factor)
    {
      itsDoSpeedUp = true;
      itsSpeedUpFactor = factor;
    }
    
    
    /// @brief Solve for parameters, updating the values kept internally
    /// The solution is constructed from the normal equations. The parameters named 
    /// image* are interpreted as images and solved for.
    /// @param[in] ip current model (to be updated)
    /// @param[in] quality Solution quality information
    bool ImageMultiScaleSolver::solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality)
    {
      ASKAPTRACE("ImageMultiScaleSolver::solveNormalEquations");
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
      ASKAPCHECK(nParameters>0, "No free parameters in ImageMultiScaleSolver");
      
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
	    casa::Array<float> cleanArray = padImage(planeIter.getPlane(ip.value(indit->first)));
	    casa::Array<float> maskArray(dirtyArray.shape());
	    ASKAPLOG_INFO_STR(logger, "Plane shape "<<planeIter.planeShape()<<" becomes "<<
			      dirtyArray.shape()<<" after padding");

	    // Precondition the PSF and DIRTY images before solving.
	    const bool wasPreconditioning = doPreconditioning(psfArray,dirtyArray);
        doNormalization(padDiagonal(planeIter.getPlane(diag)),tol(),psfArray,dirtyArray, 
                        boost::shared_ptr<casa::Array<float> >(&maskArray, utility::NullDeleter()));
	    if (wasPreconditioning) {
	      // Store the new PSF in parameter class to be saved to disk later
	      saveArrayIntoParameter(ip, indit->first, planeIter.shape(), "psf.image", unpadImage(psfArray),
				     planeIter.position());
	    } // if there was preconditioning
	    
	    // optionally clip the image and psf if there was padding
	    ASKAPLOG_INFO_STR(logger, "Peak data vector flux (derivative) before clipping "<<max(dirtyArray));
	    clipImage(dirtyArray);
	    clipImage(psfArray);
	    ASKAPLOG_INFO_STR(logger, "Peak data vector flux (derivative) after clipping "<<max(dirtyArray));
	    
	    // We need lattice equivalents. We can use ArrayLattice which involves
	    // no copying
	    casa::ArrayLattice<float> dirty(dirtyArray);
	    casa::ArrayLattice<float> psf(psfArray);
	    casa::ArrayLattice<float> clean(cleanArray);
	    casa::ArrayLattice<float> mask(maskArray);
	    
	    // uncomment the code below to save the residual image
	    // This takes up some memory and we have to ship the residual image out inside
	    // the parameter class. Therefore, we may not need this functionality in the 
	    // production version (or may need to implement it in a different way).
	    saveArrayIntoParameter(ip, indit->first, planeIter.shape(), "residual", unpadImage(dirtyArray),
				   planeIter.position());
	    
            /*
	    // uncomment the code below to save the mask
	    saveArrayIntoParameter(ip, indit->first, planeIter.shape(), "mask", unpadImage(maskArray),
				   planeIter.position());
            */
	    
	    
	    // Create a lattice cleaner to do the dirty work :)
	    /// @todo More checks on reuse of LatticeCleaner
	    // every plane should have its own LatticeCleaner, therefore we should ammend the 
	    // key somehow to make it individual for each plane. Adding tag seems to be a good idea
	    const std::string cleanerKey = indit->first + planeIter.tag();
	    
	    itsCleaners.find(cleanerKey);                     
	    boost::shared_ptr<casa::LatticeCleaner<float> > lc = itsCleaners.cachedItem();
	    if(!itsCleaners.notFound()) {
	      ASKAPDEBUGASSERT(lc);
	      lc->update(dirty);
	    } else {
	      itsCleaners.cachedItem().reset(new casa::LatticeCleaner<float>(psf, dirty));
	      lc = itsCleaners.cachedItem();
	      ASKAPDEBUGASSERT(lc);     
	      if (itsDoSpeedUp) {
		lc->speedup(itsSpeedUpFactor);
	      }          
	      lc->setMask(mask,maskingThreshold());	  
	      
	      if(algorithm()=="Hogbom") {
		casa::Vector<float> scales(1);
		scales(0)=0.0;
		lc->setscales(scales);
		lc->setcontrol(casa::CleanEnums::HOGBOM, niter(), gain(), threshold(),
			       fractionalThreshold(), false);
	      } else {
		lc->setscales(itsScales);
		lc->setcontrol(casa::CleanEnums::MULTISCALE, niter(), gain(), threshold(), 
			       fractionalThreshold(),false);
	      } // if algorithm == Hogbom, else case (other algorithm)
	      lc->ignoreCenterBox(true);
	    } // if cleaner found in the cache, else case - new cleaner needed
	    lc->clean(clean);
	    ASKAPLOG_INFO_STR(logger, "Peak flux of the clean image "<<max(cleanArray));
	    
	    const std::string peakResParam = std::string("peak_residual.") + cleanerKey;
	    if (ip.has(peakResParam)) {
	      ip.update(peakResParam, lc->strengthOptimum());
	    } else {
	      ip.add(peakResParam, lc->strengthOptimum());
	    }
	    ip.fix(peakResParam);	    
	    planeIter.getPlane(ip.value(indit->first)) = unpadImage(cleanArray);
	  } // loop over all planes of the image cube
	} // loop over map of indices
      
      quality.setDOF(nParameters);
      quality.setRank(0);
      quality.setCond(0.0);
      quality.setInfo("Multiscale Clean");
      
      /// Save the PSF and Weight
      saveWeights(ip);      
      savePSF(ip);
      
      return true;
    };
    
    Solver::ShPtr ImageMultiScaleSolver::clone() const
    {
      return Solver::ShPtr(new ImageMultiScaleSolver(*this));
    }
    
  }
}


