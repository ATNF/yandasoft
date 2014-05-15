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
ASKAP_LOGGER(logger, ".measurementequation.imageamsmfsolver");

#include <askap/AskapError.h>
#include <profile/AskapProfiler.h>

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/Vector.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <measurementequation/ImageParamsHelper.h>
#include <utils/MultiDimArrayPlaneIter.h>

#include <deconvolution/DeconvolverMultiTermBasisFunction.h>

#include <lattices/Lattices/LatticeCleaner.h>
#include <lattices/Lattices/MultiTermLatticeCleaner.h>
#include <lattices/Lattices/ArrayLattice.h>

#include <measurementequation/ImageAMSMFSolver.h>

using namespace casa;
using namespace askap;
using namespace askap::scimath;

#include <iostream>

#include <cmath>
using std::abs;

#include <map>
#include <vector>
#include <string>
#include <set>

using std::map;
using std::vector;
using std::string;

namespace askap
{
  namespace synthesis
  {
    
    
    ImageAMSMFSolver::ImageAMSMFSolver() : itsScales(3,0.), itsNumberTaylor(0),
					   itsSolutionType("MINCHISQ"), itsOrthogonal(False)
    {
      ASKAPDEBUGASSERT(itsScales.size() == 3);
      itsScales(1)=10;
      itsScales(2)=30;
      // Now set up controller
      itsControl.reset(new DeconvolverControl<Float>());
      // Now set up monitor
      itsMonitor.reset(new DeconvolverMonitor<Float>());
      const BasisFunction<Float>::ShPtr bfPtr(new MultiScaleBasisFunction<Float>(itsScales,
										      itsOrthogonal));
      itsBasisFunction = bfPtr;
    }
    
    ImageAMSMFSolver::ImageAMSMFSolver(const casa::Vector<float>& scales) : 
      itsScales(scales), itsNumberTaylor(0), itsSolutionType("MINCHISQ"), itsOrthogonal(False)
    {
      // Now set up controller
      itsControl.reset(new DeconvolverControl<Float>());
      // Now set up monitor
      itsMonitor.reset(new DeconvolverMonitor<Float>());

      const BasisFunction<Float>::ShPtr bfPtr(new MultiScaleBasisFunction<Float>(scales,
										      itsOrthogonal));
      itsBasisFunction = bfPtr;
    }
    
    Solver::ShPtr ImageAMSMFSolver::clone() const
    {
      return Solver::ShPtr(new ImageAMSMFSolver(*this));
    }
    
    void ImageAMSMFSolver::init()
    {
      resetNormalEquations();
    }
    
    /// @brief Solve for parameters
    /// The solution is constructed from the normal equations
    /// The solution is constructed from the normal equations. The parameters named 
    /// image* are interpreted as images and solved for.
    /// @param[in] ip current model (to be updated)        
    /// @param[in] quality Solution quality information
    bool ImageAMSMFSolver::solveNormalEquations(askap::scimath::Params& ip,askap::scimath::Quality& quality)
    {
      ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations");

      // Solving A^T Q^-1 V = (A^T Q^-1 A) P
      
      // Find all the free parameters beginning with image
      vector<string> names(ip.completions("image"));
      for (vector<string>::iterator it = names.begin(); it!=names.end(); ++it) {
           *it = "image" + *it;
      }
      // this should work for faceting as well, taylorMap would contain one element
      // per facet in this case
      std::map<std::string, int> taylorMap;
      SynthesisParamsHelper::listTaylor(names, taylorMap);
      
      uint nParameters=0;
      ASKAPCHECK(taylorMap.size() != 0, "Solver doesn't have any images to solve for");
      for (std::map<std::string, int>::const_iterator tmIt = taylorMap.begin(); 
           tmIt!=taylorMap.end(); ++tmIt) {
	
	// The MSMF Solver expects 2xNTaylor-1 image parameter for each Stokes parameter.
	// 
	// Desired loop structure for multiple stokes and Taylor terms.
	// For image.i loop over Taylor 0,1,2,...
	// For image.q loop over Taylor 0,1,2,...
	// ...
	//
	// stokeslist = ['i','q']
	// ntaylor = 3
	//
	// for ( stokes in stokeslist )
	//    latticecleaner[stokes]->setup();
	//    for ( order in [0:(2*ntaylor-1)-1] )
	//        latticecleaner[stokes]->setpsf(order,psf[order]);
	//        if ( order < ntaylor )
	//           latticecleaner[stokes]->setresidual(order,residual[order]);
	//           latticecleaner[stokes]->setmodel(order,model[order]);
	//    latticecleaner[stokes]->mtclean();
	//    for ( order in [0:ntaylor-1] )
	//        latticecleaner[stokes]->getmodel(order,model[order]);
	//
	
	ASKAPLOG_INFO_STR(logger, "AMSMFS minor cycle, processing image "<<tmIt->first);
	// Determine the number of stokes planes and ensuring that all Taylor terms
	// have the same number of polarisations
	ASKAPDEBUGASSERT(tmIt->second != 0);
	// this can be a facet, hence create a helper
	ImageParamsHelper iph(tmIt->first);
	
	// nOrders is the number of free parameters
	const uInt nOrders = tmIt->second;
	
	// nOrders is the total number of free parameters. Initially this is 2 * nTaylor - 1.
	// This will not work correctly if the number of terms differs between images!
	if(itsNumberTaylor==0) {
	  itsNumberTaylor=(nOrders+1)/2;
	  ASKAPLOG_INFO_STR(logger, "There are " << itsNumberTaylor << " Taylor terms");
	  ASKAPLOG_INFO_STR(logger, "There are " << nOrders << " PSFs calculated for this first pass");
	}
	else {
	  ASKAPLOG_INFO_STR(decmtbflogger, "There are " << itsNumberTaylor << " Taylor terms");
	}
	
	if(this->itsNumberTaylor>1) {
	  iph.makeTaylorTerm(0);
	}
	else {
	  ASKAPLOG_INFO_STR(logger, "No Taylor terms will be solved");
	}
	const casa::IPosition imageShape = ip.value(iph.paramName()).shape();               
	const uint nPol = imageShape.nelements()>=3 ? uint(imageShape(2)) : 1;
	ASKAPLOG_INFO_STR(logger, "There are " << nPol << " polarisation planes to solve" );
	nParameters += imageShape.product(); // add up the number of pixels for zero order
	// check consistency
	for (uInt order=1;order<uInt(tmIt->second);++order) {
	  // make the helper a Taylor term of the given order
	  iph.makeTaylorTerm(order);
	  const casa::IPosition thisShape = ip.value(iph.paramName()).shape();               
	  const uint thisNPol = thisShape.nelements()>=3 ? uint(thisShape(2)) : 1;
	  ASKAPCHECK(thisNPol == nPol, "Number of polarisations are supposed to be consistent for all Taylor terms, order="<<
		     order<<" has "<<thisNPol<<" polarisation planes");
	  nParameters += thisShape.product(); // add up the number of pixels for this order
	}
	
	/*
	// this check is temporary, to avoid unnecessary surprises while further developing the code
	if (imageShape.nelements()>=4) {
	  ASKAPCHECK(imageShape(3) == 1, "Output cube for MSMFS solver should have just one spectral plane, shape="<<
		     imageShape<<" nPol="<<nPol);
	}
	//
	*/
	
	// as polarisations are not necessarily represented by a different parameter
	// we have to build a set of parameters which are going to be fixed inside the loop
	// (or alternatively fix them multiple times, which is also a reasonable solution)
	std::set<std::string> parametersToBeFixed;
	
	// Iterate through Polarisations
	for (scimath::MultiDimArrayPlaneIter planeIter(imageShape); planeIter.hasMore(); planeIter.next()) {
	  const uint plane = planeIter.sequenceNumber();
	  std::string tagLogString(planeIter.tag());
	  if (tagLogString.size()) {
	    tagLogString = "tagged as " + tagLogString;
	  } else {
	    tagLogString = "not tagged";
	  }
	  
	  ASKAPLOG_INFO_STR(logger, "Preparing iteration for plane " 
			    << plane<<" ("<<tagLogString<<") in image "<<tmIt->first);
	  // make the helper a 0-order Taylor term
	  if(this->itsNumberTaylor>1) {
	    ASKAPLOG_INFO_STR(logger, "Solving for " << this->itsNumberTaylor << " Taylor terms");
	    iph.makeTaylorTerm(0);
	  }
	  else {
	    ASKAPLOG_INFO_STR(logger, "No Taylor terms will be solved");
	  }
	  const std::string zeroOrderParam = iph.paramName();
	  
	  // Setup the normalization vector	  
	  ASKAPLOG_INFO_STR(logger, "Reading the normalization vector from : " << zeroOrderParam);
	  ASKAPCHECK(normalEquations().normalMatrixDiagonal().count(zeroOrderParam)>0, "Diagonal not present");
	  casa::Vector<double> normdiag(normalEquations().normalMatrixDiagonal().find(zeroOrderParam)->second);
	  
	  ASKAPDEBUGASSERT(planeIter.planeShape().nelements()>=2);
	  
	  const double maxDiag = casa::max(planeIter.getPlaneVector(normdiag));
	  ASKAPLOG_INFO_STR(logger, "Maximum of weights = " << maxDiag );
	  
	  // a unique string for every Taylor decomposition (unique for every facet for faceting)
	  const std::string imageTag = tmIt->first + planeIter.tag();
	  const bool firstcycle = !SynthesisParamsHelper::hasValue(itsCleaners,imageTag);          
	  
	  Vector<Array<Float> > cleanVec(itsNumberTaylor);
	  Vector<Array<Float> > dirtyVec(itsNumberTaylor);
	  Vector<Array<Float> > dirtyLongVec(2*itsNumberTaylor-1);
	  Vector<Array<Float> > psfVec(itsNumberTaylor);
	  Vector<Array<Float> > psfLongVec(2*itsNumberTaylor-1);
	  
	  // Setup the PSFs - all ( 2 x ntaylor - 1 ) of them for the first time. We keep a copy of
	  // the first since we will need it for preconditioning the others

	  const uInt limit = firstcycle ? 2*this->itsNumberTaylor-1 : this->itsNumberTaylor;
	  
	  // Now precondition the residual images using the zeroth order psf. We need to 
	  // keep a copy of the zeroth PSF to avoid having it overwritten each time.
	  Array<float> psfWorkArray;
	  for( uInt order=0; order < limit; ++order) {
	       ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._initcopy");
	    
	    // make helper to represent the given order
	    if(this->itsNumberTaylor>1) {
	      ASKAPLOG_INFO_STR(logger, "Solving for Taylor term " << order);
	      iph.makeTaylorTerm(order);
	    }
	    else {
	      ASKAPLOG_INFO_STR(logger, "No Taylor terms will be solved");
	    }
	    const std::string thisOrderParam = iph.paramName();
	    ASKAPLOG_INFO_STR(logger, "AMSMFS solver: processing order "
			      << order << " (" << itsNumberTaylor <<
			      " Taylor terms + " << itsNumberTaylor-1 << " cross-terms), parameter name: " << thisOrderParam);
	    
	    // Always need the actual PSFs, and on first cycle we also need the
	    // higher order terms
	    if ((order<this->itsNumberTaylor)||(firstcycle)) {
	      ASKAPCHECK(normalEquations().normalMatrixSlice().count(thisOrderParam)>0,
			 "PSF Slice for plane="<< plane<<" and order="<<order<<" is not present");
	      casa::Vector<double> slice(normalEquations().normalMatrixSlice().find(thisOrderParam)->second);
	      psfLongVec(order).resize(planeIter.planeShape());
	      casa::convertArray<float, double>(psfLongVec(order), planeIter.getPlane(slice));
	    }
	    
	    // For the dirty images, we need all terms
	    ASKAPCHECK(normalEquations().dataVector(thisOrderParam).size()>0, "Data vector not present for cube plane="<<
		       plane<<" and order="<<order);
	    dirtyLongVec(order).resize(planeIter.planeShape());
	    casa::Vector<double> dv = normalEquations().dataVector(thisOrderParam);
	    casa::convertArray<float, double>(dirtyLongVec(order), planeIter.getPlane(dv));
	    
	    // For the clean images, we need only the first nTaylor
	    if(order<this->itsNumberTaylor) {
	      cleanVec(order).resize(planeIter.planeShape());
	      casa::convertArray<float, double>(cleanVec(order), 
						planeIter.getPlane(ip.value(thisOrderParam)));
	    }
	  } // Loop over order
	  
	  casa::Array<float> maskArray(planeIter.planeShape());
	  
	  uInt nx(planeIter.planeShape()(0));
	  uInt ny(planeIter.planeShape()(1));
	  IPosition centre(2, nx/2, ny/2);

	  itsPSFZeroArray=psfLongVec(0).copy();
	  if (firstcycle) {
          ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._fc_norm+precnd");
	    
	    ASKAPLOG_DEBUG_STR(logger, "Deriving scale from PSF(0) centre value " << psfLongVec(0).nonDegenerate()(centre));
	    // For the first cycle we need to precondition and normalise all PSFs and all dirty images
	    itsPSFZeroCentre=-1;
	    for(uInt order=0; order < 2 * itsNumberTaylor - 1; ++order) {
	      // We need to work with the original preconditioning PSF since it gets overridden
	      psfWorkArray = itsPSFZeroArray.copy();
	      ASKAPLOG_DEBUG_STR(logger, "Initial PSF(" << order << ") centre value " << psfLongVec(order).nonDegenerate()(centre));
              // Precondition this order PSF using PSF(0)
	      if(doPreconditioning(psfWorkArray, psfLongVec(order))) {
             ASKAPLOG_DEBUG_STR(logger, "After preconditioning PSF(" << order << ") centre value " << psfLongVec(order).nonDegenerate()(centre));
             // Now we can precondition the dirty (residual) array using PSF(0)
             psfWorkArray = itsPSFZeroArray.copy();
             ASKAPLOG_INFO_STR(logger, "Preconditioning dirty image for plane=" << plane<<
                                " ("<<tagLogString<< ") and order=" << order);
             doPreconditioning(psfWorkArray,dirtyLongVec(order));
	      }
	      // Normalise. 
	      ASKAPLOG_DEBUG_STR(logger, "Normalising PSF and Dirty image for order " << order);
	      ASKAPLOG_DEBUG_STR(logger, "Before normalisation PSF(" << order << ") centre value " << psfLongVec(order).nonDegenerate()(centre));
	      // First call the scaling is via the psf and the value is returned. Thereafter we use that value for the normalisation
	      // of all the PSFs. Thus for MFS, the first PSF should have centre value 1.0 and the others lower values
	      if(order==0) {
             itsPSFZeroCentre = doNormalization(planeIter.getPlaneVector(normdiag),tol(),psfLongVec(order),dirtyLongVec(order),
						 boost::shared_ptr<casa::Array<float> >(&maskArray, utility::NullDeleter()));
	      }  else {
             doNormalization(planeIter.getPlaneVector(normdiag),tol(),psfLongVec(order),itsPSFZeroCentre,dirtyLongVec(order),
				boost::shared_ptr<casa::Array<float> >(&maskArray, utility::NullDeleter()));
	      }
	      ASKAPLOG_DEBUG_STR(logger, "After  normalisation PSF(" << order << ") centre value " << psfLongVec(order).nonDegenerate()(centre));
	      if (order<itsNumberTaylor) {
              psfVec(order)=psfLongVec(order);
              dirtyVec(order)=dirtyLongVec(order);
	      }
	    }// Loop over order
	  }
	  else {
        ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._norm+precnd");
	    // For the subsequent cycles cycle we need to precondition and normalise the updated dirty images
            // Precondition the dirty (residual) array
	    for(uInt order=0; order < itsNumberTaylor; ++order) {
	      psfWorkArray = itsPSFZeroArray.copy();
	      if(doPreconditioning(psfWorkArray,dirtyLongVec(order))) {
		ASKAPLOG_INFO_STR(logger, "Preconditioning dirty image for plane=" << plane<< " ("<<tagLogString<< ") and order=" << order);
	      }
	      // Normalise. 
	      psfWorkArray = itsPSFZeroArray.copy();
	      doNormalization(planeIter.getPlaneVector(normdiag),tol(),psfWorkArray,itsPSFZeroCentre,dirtyLongVec(order),
			      boost::shared_ptr<casa::Array<float> >(&maskArray, utility::NullDeleter()));
	      if(order<itsNumberTaylor) {
		dirtyVec(order)=dirtyLongVec(order);
	      }
	    }// Loop over order
	  }
	  
	  // Now that we have all the required images, we can initialise the deconvolver
	  if (firstcycle)  {// Initialize everything only once.
          ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._fc_initdeconvolver");
	      
	      ASKAPLOG_INFO_STR(logger, "Creating solver for plane " << plane <<" tag "<<imageTag);
	      itsCleaners[imageTag].reset(new DeconvolverMultiTermBasisFunction<Float, Complex>(dirtyVec, psfVec, psfLongVec));
	      ASKAPDEBUGASSERT(itsCleaners[imageTag]);
	      
	      itsCleaners[imageTag]->setMonitor(itsMonitor);
	      itsControl->setTargetObjectiveFunction(threshold().getValue("Jy"));
	      itsControl->setFractionalThreshold(fractionalThreshold());
	      
	      itsCleaners[imageTag]->setControl(itsControl);
	      
	      ASKAPCHECK(itsBasisFunction, "Basis function not initialised");

              ASKAPDEBUGASSERT(dirtyVec.nelements() > 0);
	      itsBasisFunction->initialise(dirtyVec(0).shape());
	      itsCleaners[imageTag]->setBasisFunction(itsBasisFunction);
	      itsCleaners[imageTag]->setSolutionType(itsSolutionType);
	      if (maskArray.nelements()) {
                  ASKAPLOG_INFO_STR(logger, "Defining mask as weight image");
		  itsCleaners[imageTag]->setWeight(maskArray);
	      }
	  } else {
          ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._updatedeconvolver");
	      ASKAPCHECK(itsCleaners[imageTag], "Deconvolver not yet defined");
	      // Update the dirty images
	      ASKAPLOG_INFO_STR(logger, "Multi-Term Basis Function deconvolver already exists - update dirty images");
	      itsCleaners[imageTag]->updateDirty(dirtyVec);
	      ASKAPLOG_INFO_STR(logger, "Successfully updated dirty images");
	  }
	    
	  // We have to reset the initial objective function
	  // so that the fractional threshold mechanism will work.
	  itsCleaners[imageTag]->state()->resetInitialObjectiveFunction();
	  // By convention, iterations are counted from scratch each
	  // major cycle
	  itsCleaners[imageTag]->state()->setCurrentIter(0);

	  for (uInt order=0; order < itsNumberTaylor; ++order) {
	    if (this->itsNumberTaylor>1) {
	      ASKAPLOG_INFO_STR(logger, "Solving for Taylor term " << order);
	      iph.makeTaylorTerm(order);
	    }
	    else {
	      ASKAPLOG_INFO_STR(logger, "No Taylor terms will be solved");
	    }
	    const std::string thisOrderParam = iph.paramName();
	    
	    if(saveIntermediate()) {
              ASKAPLOG_DEBUG_STR(logger, "Dirty(" << order << ") shape = " << dirtyVec(order).shape());
	      saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "residual",
				     unpadImage(dirtyVec(order)), planeIter.position());
              if(firstcycle) {
                ASKAPLOG_DEBUG_STR(logger, "PSF(" << order << ") shape = " << psfVec(order).shape());
                saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "psf.image",
                                       unpadImage(psfVec(order)), planeIter.position());
                if(order==0&&maskArray.nelements()) {
                  saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "mask", unpadImage(maskArray),
                                         planeIter.position());
                }
	      }
	    }
	    
            casa::Array<float> cleanArray(planeIter.planeShape());
            casa::convertArray<float, double>(cleanArray, planeIter.getPlane(ip.value(thisOrderParam)));
            itsCleaners[imageTag]->setModel(cleanArray, order);
	  } // end of 'order' loop
	  
	  {
         ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._calldeconvolver");	  
         ASKAPLOG_INFO_STR(logger, "Starting Minor Cycles ("<<imageTag<<").");
         itsCleaners[imageTag]->deconvolve();
         ASKAPLOG_INFO_STR(logger, "Finished Minor Cycles ("<<imageTag<<").");
      }
	  
          // Now update the stored peak residual
          const std::string peakResParam = std::string("peak_residual.") + imageTag;
          if (ip.has(peakResParam)) {
            ip.update(peakResParam, itsCleaners[imageTag]->state()->peakResidual());
          } else {
            ip.add(peakResParam, itsCleaners[imageTag]->state()->peakResidual());
          }
          ip.fix(peakResParam);	    

	  // Write the final vector of clean model images into parameters
	  for( uInt order=0; order < itsNumberTaylor; ++order) {
	    // make the helper to correspond to the given order
	    if(this->itsNumberTaylor>1) {
	      ASKAPLOG_INFO_STR(logger, "Solved for Taylor term " << order);
	      iph.makeTaylorTerm(order);
	    }
	    else {
	      ASKAPLOG_INFO_STR(logger, "No Taylor terms were solved");
	    }
	    const std::string thisOrderParam = iph.paramName();
	    ASKAPLOG_INFO_STR(logger, "About to get model for plane="<<plane<<" Taylor order="<<order<<
			      " for image "<<tmIt->first);
	    planeIter.getPlane(ip.value(thisOrderParam)).nonDegenerate()=unpadImage(itsCleaners[imageTag]->model(order));
	  }
	  // add extra parameters (cross-terms) to the to-be-fixed list
	  for (uInt order = itsNumberTaylor; order<uInt(tmIt->second); ++order) {
	    // make the helper to correspond to the given order
	    if(this->itsNumberTaylor>1) {
	      iph.makeTaylorTerm(order);
	    }
	    const std::string thisOrderParam = iph.paramName();
	    parametersToBeFixed.insert(thisOrderParam);
	  }
	} // end of polarisation (i.e. plane) loop 
	
	// Make sure that the next set of minor cycles does not redo unnecessary things.
	// Also "fix" parameters for order >= itsNumberTaylor. so that the gridding doesn't get done
	// for these extra terms.
	
	// Fix the params corresponding to extra Taylor terms.
	// (MV) probably this part needs another careful look
	for (std::set<std::string>::const_iterator ci = parametersToBeFixed.begin(); 
	     ci != parametersToBeFixed.end(); ++ci) {
	  if (ip.isFree(*ci)) {
	    ip.fix(*ci);
	  }
	}
      } // loop: tmIt
      
      ASKAPCHECK(nParameters>0, "No free parameters in ImageAMSMFSolver");
      
      quality.setDOF(nParameters);
      quality.setRank(0);
      quality.setCond(0.0);
      quality.setInfo("Multi-Scale Multi-Frequency Clean");
      
      /// Save PSFs and Weights into parameter class (to be exported later)
      saveWeights(ip);
      savePSF(ip);
      
      return true;
    };
    
    void ImageAMSMFSolver::setBasisFunction(BasisFunction<Float>::ShPtr bf) {
      itsBasisFunction=bf;
    }
    
    BasisFunction<Float>::ShPtr ImageAMSMFSolver::basisFunction() {
      return itsBasisFunction;
    }
    
    void ImageAMSMFSolver::configure(const LOFAR::ParameterSet &parset) {
      ImageSolver::configure(parset);

      ASKAPASSERT(this->itsMonitor);
      this->itsMonitor->configure(parset);
      ASKAPASSERT(this->itsControl);
      this->itsControl->configure(parset);
      
      String solutionType=parset.getString("solutiontype", "MAXCHISQ");
      if(solutionType=="MAXBASE") {
      }
      else if(solutionType=="MAXTERM0") {
      }
      else {
	solutionType="MAXCHISQ";
      }
      ASKAPLOG_INFO_STR(decmtbflogger, "Solution type = " << solutionType);
      this->itsSolutionType=solutionType;

      this->itsOrthogonal=parset.getBool("orthogonal", "false");
      if (this->itsOrthogonal) {
        ASKAPLOG_DEBUG_STR(decmtbflogger, "Multiscale basis functions will be orthogonalised");
      }

    }
  }
}
