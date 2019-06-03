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

// System includes
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <string>

// ASKAPsoft includes
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/scimath/utils/MultiDimArrayPlaneIter.h>
#include <askap/askap_synthesis.h>
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <profile/AskapProfiler.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/scimath/Mathematics/VectorKernel.h>
#include <casacore/images/Images/TempImage.h>

// Local package includes
#include <askap/measurementequation/ImageRestoreSolver.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/ImageParamsHelper.h>
#include <askap/measurementequation/ImageSolverFactory.h>
#include <askap/measurementequation/IImagePreconditioner.h>
#include <askap/measurementequation/WienerPreconditioner.h>
#include <askap/measurementequation/GaussianTaperPreconditioner.h>
#include <askap/measurementequation/Image2DConvolver.h>

ASKAP_LOGGER(logger, ".measurementequation.imagerestoresolver");

// Using
using namespace casa;
using namespace askap;
using namespace askap::scimath;
using std::abs;
using std::map;
using std::vector;
using std::string;

namespace askap
{
  namespace synthesis
  {
    ImageRestoreSolver::ImageRestoreSolver(const RestoringBeamHelper &beamHelper) :
	    itsBeamHelper(beamHelper), itsEqualiseNoise(false), itsModelNeedsConvolving(true), itsResidualNeedsUpdating(true)
    {
        setIsRestoreSolver();
    }

    void ImageRestoreSolver::init()
    {
	    resetNormalEquations();
    }

    /// @brief Solve for parameters, updating the values kept internally
    /// The solution is constructed from the normal equations. The parameters named
    /// image* are interpreted as images and solved for.
    /// @param[in] ip current model (to be updated)
    /// @param[in] quality Solution quality information
    bool ImageRestoreSolver::solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality)
    {
        ASKAPTRACE("ImageRestoreSolver::solveNormalEquations");
	// Solving A^T Q^-1 V = (A^T Q^-1 A) P

	// Find all the free parameters beginning with image
	vector<string> names(ip.completions("image"));
	uint nParameters=0;
	for (vector<string>::iterator it=names.begin(); it!=names.end(); ++it)
	{
		const string name="image"+*it;
		// completions should return only free parameters according to its code in Code/Base
		ASKAPDEBUGASSERT(ip.isFree(name));
		*it = name; // append the common part to the front of the parameter name
		nParameters+=ip.value(name).nelements();
	}

	ASKAPCHECK(nParameters>0, "No free parameters in ImageRestoreSolver");
	// get restoring beam
	ASKAPCHECK(itsBeamHelper.valid(), "Unable to obtain beam estimate as itsBeamHelper is not initialised");
	std::string psfName = SynthesisParamsHelper::findPSF(ip);
	if (psfName == "") {
	    // it doesn't matter that we don't have access to preconditioned PSF at this stage - only zero model will be
	    // convolved with this PSF
            ASKAPLOG_INFO_STR(logger, "ImageRestoreSolver::solveNormalEquations: Extracting PSF from normal equations");
            savePSF(ip);
            psfName = SynthesisParamsHelper::findPSF(ip);
            ASKAPCHECK(psfName != "", "Failed to find a PSF parameter");
    }

	// determine which images are faceted and, if need be, setup parameters
    // representing the result of a merge.
	map<string,int> facetmap;
	SynthesisParamsHelper::listFacets(names, facetmap);
	for (map<string,int>::const_iterator ci=facetmap.begin();ci!=facetmap.end();++ci) {
	     if ((ci->second != 1) && !ip.has(ci->first)) {
	         // this is a multi-facet image, add a fixed parameter representing the whole image
	         ASKAPLOG_INFO_STR(logger, "Adding a fixed parameter " << ci->first<<
                           " representing faceted image with "<<ci->second<<" facets");
	         SynthesisParamsHelper::add(ip,ci->first,ci->second);
	         ip.fix(ci->first);
	     }
	}
	//

	// iterate over all free parameters (i.e. parts of the image for faceted case)
	for (vector<string>::const_iterator ci=names.begin(); ci !=names.end(); ++ci) {
      ImageParamsHelper iph(*ci);
      // obtain name with just the taylor suffix, if present
      const std::string name = iph.taylorName();

	  if (facetmap[name] == 1) {
	      // this is not a faceting case, restore the image in situ and add residuals
	      ASKAPLOG_INFO_STR(logger, "Restoring " << *ci );

          // convolve with restoring beam before adding residuals, but wait until after
          // preconditioning to calculate restoring beam size.
          itsModelNeedsConvolving = true;

	      // add residuals
	      addResiduals(ip, *ci);

	      SynthesisParamsHelper::setBeam(ip, *ci, itsBeamHelper.value());
      } else {
          // this is a single facet of a larger image, just fill in the bigger image with the model
          ASKAPLOG_INFO_STR(logger, "Inserting facet " << iph.paramName()<<" into merged image "<<name);
          casacore::Array<double> patch = SynthesisParamsHelper::getFacet(ip,iph.paramName());
          const casacore::Array<double> model = scimath::PaddingUtils::centeredSubArray(ip.value(iph.paramName()),
                                          patch.shape());
          patch = model;
      }
	}

	// restore faceted images
    for (map<string,int>::const_iterator ci=facetmap.begin();ci!=facetmap.end();++ci) {
	     if (ci->second != 1) {
	         // this is a multi-facet image
	         ASKAPLOG_INFO_STR(logger, "Restoring faceted image " << ci->first );

             // convolve with restoring beam before adding residuals, but wait until after
             // preconditioning to calculate restoring beam size. This is done once for
             // the set of facets.
             itsModelNeedsConvolving = true;

	         // add residuals
	         for (int xFacet = 0; xFacet<ci->second; ++xFacet) {
	              for (int yFacet = 0; yFacet<ci->second; ++yFacet) {
	                   ASKAPLOG_INFO_STR(logger, "Adding residuals for facet ("<<xFacet<<","
	                        <<yFacet<<")");
	                   // ci->first may have taylor suffix defined, load it first and then add facet indices
	                   ImageParamsHelper iph(ci->first);
	                   iph.makeFacet(xFacet,yFacet);
	                   addResiduals(ip, ci->first, iph.paramName());

	              }
	         }

	         SynthesisParamsHelper::setBeam(ip, ci->first, itsBeamHelper.value());

	     }
	}

	quality.setDOF(nParameters);
	quality.setRank(0);
	quality.setCond(0.0);
	quality.setInfo("Restored image calculated");

	return true;
    };

    /// @brief solves for and adds residuals
    /// @details Restore solver convolves the current model with the beam and adds the
    /// residual image. The latter has to be "solved for" with a proper preconditioning and
    /// normalisation using the normal equations stored in the base class. All operations
    /// required to extract residuals from normal equations and fill an array with them
    /// are encapsulated in this method. Faceting needs a subimage only, hence the array
    /// to fill may not have exactly the same shape as the dirty (residual) image corresponding
    /// to the given parameter. This method assumes that the centres of both images are the same
    /// and extracts only data required.
    /// @param[in,out] ip param object containing all of the parameters
    /// @param[in] imagename name of the parameter to work with
    /// @param[in] facetname name of the current facet, if using facets
    void ImageRestoreSolver::addResiduals(askap::scimath::Params& ip,
                                          const std::string &imagename,
                                          const std::string facetname)

    {
       ASKAPTRACE("ImageRestoreSolver::addResiduals");

       // when convolving model images with the restoring beam and
       // adding residuals, the full, merged model needs to be used.
       std::string name;
       casacore::Array<double> out;
       if (facetname == "") {
         name = imagename;
         out.reference(ip.value(name));
       }
       else {
         name = facetname;
         out.reference(SynthesisParamsHelper::getFacet(ip,name));
       }
       const casacore::IPosition shape = ip.value(name).shape();
       const scimath::Axes axes = ip.axes(name);

	   // Axes are dof, dof for each parameter
	   //casacore::IPosition vecShape(1, out.shape().product());
       bool saveNewPSFRequired = true;

	   for (scimath::MultiDimArrayPlaneIter planeIter(shape); planeIter.hasMore(); planeIter.next()) {

	        ASKAPCHECK(normalEquations().normalMatrixDiagonal().count(name)>0,
                "Diagonal not present " << name);
	        casacore::Vector<double> diag(normalEquations().normalMatrixDiagonal().find(name)->second);
	        ASKAPCHECK(normalEquations().dataVector(name).size()>0,
                "Data vector not present " << name);
	        casacore::Vector<double> dv = normalEquations().dataVector(name);
	        ASKAPCHECK(normalEquations().normalMatrixSlice().count(name)>0,
                "PSF Slice not present " << name);

            casacore::Vector<double> slice(normalEquations().normalMatrixSlice().find(name)->second);
	        ASKAPCHECK(normalEquations().preconditionerSlice().count(name)>0,
                "Preconditioner fuction Slice not present for " << name);
	        casacore::Vector<double> pcf(normalEquations().preconditionerSlice().find(name)->second);

            if (planeIter.tag()!="") {
                // it is not a single plane case, there is something to report
                ASKAPLOG_INFO_STR(logger, "Processing plane "<<planeIter.sequenceNumber()<<
                                       " tagged as "<<planeIter.tag());

            }

	        ASKAPLOG_INFO_STR(logger, "Maximum of data vector corresponding to "<<name<<" is "<<casacore::max(dv));

            casacore::Array<float> dirtyArray(planeIter.planeShape());
	        casacore::convertArray<float, double>(dirtyArray,planeIter.getPlane(dv));

	        ASKAPLOG_INFO_STR(logger, "Maximum of data vector corresponding to "<<name<<" and plane "<<
	                 planeIter.sequenceNumber()<<" is "<<casacore::max(dirtyArray));

            casacore::Array<float> psfArray(planeIter.planeShape());
            casacore::convertArray<float, double>(psfArray, planeIter.getPlane(slice));

	        // send an anternative preconditioner function, if it isn't empty.
	        casacore::Array<float> pcfArray;
            if (pcf.shape() > 0) {
	          ASKAPDEBUGASSERT(pcf.shape() == slice.shape());
              pcfArray.resize(planeIter.planeShape());
	          casacore::convertArray<float, double>(pcfArray, planeIter.getPlane(pcf));
            }

            // uninitialised mask shared pointer means that we don't need it (i.e. no weight equalising)
            boost::shared_ptr<casacore::Array<float> > mask;
            if (itsEqualiseNoise) {
                ASKAPLOG_INFO_STR(logger, "Residual will be multiplied by sqrt(normalised weight) during restoration");
                // mask will have a noramised sqrt(weight) pattern after doNormalization
                mask.reset(new casacore::Array<float>(dirtyArray.shape()));
            } else {
                ASKAPLOG_INFO_STR(logger, "Restored image will have primary beam corrected noise (no equalisation)");
            }

            // Do the preconditioning
            doPreconditioning(psfArray,dirtyArray,pcfArray);

            // Normalize by the diagonal
            doNormalization(planeIter.getPlaneVector(diag),tol(),psfArray,dirtyArray,mask);

	        // we have to do noise equalisation for final residuals after preconditioning
	        if (itsEqualiseNoise) {
	            const casacore::IPosition vecShape(1,dirtyArray.nelements());
                casacore::Vector<float> dirtyVector(dirtyArray.reform(vecShape));
                ASKAPDEBUGASSERT(mask);
	            const casacore::Vector<float> maskVector(mask->reform(vecShape));
	            for (int i = 0; i<vecShape[0]; ++i) {
	                 dirtyVector[i] *= maskVector[i];
	            }
	        }
            if (itsResidualNeedsUpdating) {
            // Store the current dirtyImage parameter class to be saved to disk later
                ASKAPLOG_INFO_STR(logger, "Saving current residual image to model parameter");
                ASKAPLOG_INFO_STR(logger, "Shape is " << dirtyArray.shape() << " position is " << planeIter.position());
                saveArrayIntoParameter(ip,name,dirtyArray.shape(),"residual",
                dirtyArray,planeIter.position());
            }
            // Store the new PSF in parameter class to be saved to disk later


            if (saveNewPSFRequired == true) {
                ASKAPLOG_INFO_STR(logger, "Saving new PSF parameter as model NEW parameter -- needs full shape");
                saveArrayIntoParameter(ip, name, shape, "psf.image", psfArray,
      				     planeIter.position());
                saveNewPSFRequired = false;
            } else {
                ASKAPLOG_INFO_STR(logger, "Saving new PSF parameter as model EXISTING parameter -- using plane shape");
  	            saveArrayIntoParameter(ip, name, psfArray.shape(), "psf.image", psfArray,
  				           planeIter.position());
            }

	        // Add the residual image
            // First, convolve the model image to the resolution of the synthesised beam if not already done.
            casacore::Vector<casacore::Quantum<double> > restoringBeam;
            if (itsModelNeedsConvolving) {

                if (itsBeamHelper.fitRequired()) {
                    ASKAPLOG_INFO_STR(logger, "Fitting of Restoring beam required");
                    casacore::Array<double> psfDArray(psfArray.shape());
                    casacore::convertArray<double, float>(psfDArray, psfArray);
                    ASKAPLOG_INFO_STR(logger, "Fitting restoring beam");
                    restoringBeam = SynthesisParamsHelper::fitBeam(psfDArray, axes, itsBeamHelper.cutoff());
                    ASKAPDEBUGASSERT(restoringBeam.size() == 3);
                    ASKAPLOG_INFO_STR(logger, "Restore solver will convolve with the 2D gaussian: " <<
                        restoringBeam[0].getValue("arcsec") << " x "<<restoringBeam[1].getValue("arcsec") <<
                        " arcsec at position angle "<<restoringBeam[2].getValue("deg")<<" deg");
                    itsBeamHelper.assign(restoringBeam);
                } else {
                    restoringBeam = itsBeamHelper.value();
                }
                ASKAPLOG_INFO_STR(logger, "Convolving the model image to the resolution of the synthesised beam");
	            // Create a temporary image
                boost::shared_ptr<casacore::TempImage<float> >
                    image(SynthesisParamsHelper::tempImage(ip, imagename));
                askap::synthesis::Image2DConvolver<float> convolver;
                const casacore::IPosition pixelAxes(2, 0, 1);
                convolver.convolve(*image, *image, casacore::VectorKernel::GAUSSIAN,
                                   pixelAxes, restoringBeam, true, 1.0, false);
                SynthesisParamsHelper::update(ip, imagename, *image);
                // for some reason update makes the parameter free as well
                ip.fix(imagename);
                // this only needs to be done once per parameter
                itsModelNeedsConvolving = false;
            }

	        // The code below involves an extra copying. We can replace it later with a copyless version
	        // doing element by element adding explicitly.
	        const casacore::IPosition outSliceShape = planeIter.planeShape(out.shape());
	        // convertedResidual contains just one plane of residuals
	        casacore::Array<double> convertedResidual(outSliceShape);
	        convertArray(convertedResidual, scimath::PaddingUtils::centeredSubArray(dirtyArray,
	                     outSliceShape));
	        // figure out where to put the slice to (can't use planeIter functionality directly because out
	        // array can have a different shape
	        casacore::IPosition blc(out.shape().nelements(),0);
	        casacore::IPosition trc(out.shape());
	        const casacore::IPosition curPos(planeIter.position());
	        for (casacore::uInt dim = 0; dim<trc.nelements(); ++dim) {
	             trc[dim] -= 1;
	             ASKAPDEBUGASSERT(trc[dim]<out.shape()[dim]);
	             if ( (dim>=2) && (dim<curPos.nelements()) ) {
	                 blc[dim] = curPos[dim];
	                 trc[dim] = curPos[dim];
	             }
	        }
	        // copy stuff
	        Array<double> outSlice = out(blc,trc);
	        outSlice += convertedResidual;
	   }
    }

    /// @brief obtain an estimate of the restoring beam
    /// @details This method fits a 2D Gaussian into the central area of the PSF
    /// (a support is searched assuming 50% cutoff) if the appropriate option
    /// is set. Otherwise, it just returns the beam parameters passed in the constructor
    /// (i.e. user override).
    /// @param[in] name name of the parameter to work with
    casacore::Vector<casacore::Quantum<double> > ImageRestoreSolver::getBeam(const std::string &) const
    {
        return itsBeamHelper.value();
    }

    Solver::ShPtr ImageRestoreSolver::clone() const
    {
	    return Solver::ShPtr(new ImageRestoreSolver(*this));
    }

    /// @brief static method to create solver
    /// @details Each solver should have a static factory method, which is
    /// able to create a particular type of the solver and initialise it with
    /// the parameters taken from the given parset. It is assumed that the method
    /// receives a subset of parameters where the solver name, if it was present in
    /// the parset, is already taken out
    /// @param[in] parset input parset file
    /// @return a shared pointer to the solver instance
    boost::shared_ptr<ImageRestoreSolver> ImageRestoreSolver::createSolver(const LOFAR::ParameterSet &parset)
    {
       RestoringBeamHelper rbh;
       const vector<string> beam = parset.getStringVector("beam");
       if (beam.size() == 1) {
           ASKAPCHECK(beam[0] == "fit",
               "beam parameter should be either equal to 'fit' or contain 3 elements defining the beam size. You have "
               <<beam[0]);
           rbh.configureFit(parset.getDouble("beam.cutoff",0.05));
       } else {
          ASKAPCHECK(beam.size() == 3, "Need three elements for beam or a single word 'fit'. You have "<<beam);
          casacore::Vector<casacore::Quantum<double> > qBeam(3);
          for (int i=0; i<3; ++i) {
               casacore::Quantity::read(qBeam(i), beam[i]);
          }
          rbh.assign(qBeam);
       }
       //
       boost::shared_ptr<ImageRestoreSolver> result(new ImageRestoreSolver(rbh));
       const bool equalise = parset.getBool("equalise",false);
       result->equaliseNoise(equalise);
       const bool update = parset.getBool("updateresiduals",true);
       result->updateResiduals(update);

       return result;
    }

    /// @brief configure basic parameters of the restore solver
    /// @details This method configures basic parameters of this restore solver the same way as
    /// they are configured for normal imaging solver. We want to share the same parameters between
    /// these two types of solvers (e.g. weight cutoff tolerance, preconditioning, etc), but the
    /// appropriate parameters are given in a number of places of the parset, sometimes with
    /// solver-specific prefies, so parsing a parset in createSolver is not a good idea. This method
    /// does the job and encapsulates all related code.
    /// @param[in] ts template solver (to take parameters from)
    void ImageRestoreSolver::configureSolver(const ImageSolver &ts)
    {
      setThreshold(ts.threshold());
      setVerbose(ts.verbose());
      setTol(ts.tol());

      // behavior in the weight cutoff area
      zeroWeightCutoffMask(ts.zeroWeightCutoffMask());
      zeroWeightCutoffArea(ts.zeroWeightCutoffArea());
    }


  } // namespace synthesis
} // namespace askap