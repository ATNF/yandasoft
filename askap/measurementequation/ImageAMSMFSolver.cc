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

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.imageamsmfsolver");

#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/scimath/Mathematics/MatrixMathLA.h>
#include <casacore/casa/Arrays/Vector.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/ImageParamsHelper.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>
#include <askap/scimath/utils/PaddingUtils.h>

#include <askap/deconvolution/DeconvolverMultiTermBasisFunction.h>

#include <casacore/lattices/LatticeMath/LatticeCleaner.h>
#include <casacore/lattices/LatticeMath/MultiTermLatticeCleaner.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>

#include <askap/measurementequation/ImageAMSMFSolver.h>

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


    ImageAMSMFSolver::ImageAMSMFSolver() : itsScales(3,0.),itsNumberTaylor(0),
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

    ImageAMSMFSolver::ImageAMSMFSolver(const casacore::Vector<float>& scales) :
      itsScales(scales), itsNumberTaylor(0), itsSolutionType("MINCHISQ"), itsOrthogonal(False)
    {
      // Now set up controller
      itsControl.reset(new DeconvolverControl<Float>());
      // Now set up monitor
      itsMonitor.reset(new DeconvolverMonitor<Float>());

      const BasisFunction<Float>::ShPtr bfPtr(new MultiScaleBasisFunction<Float>(scales, itsOrthogonal));
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
    bool ImageAMSMFSolver::solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality)
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
      for (std::map<std::string, int>::const_iterator tmIt = taylorMap.begin(); tmIt!=taylorMap.end(); ++tmIt) {

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
            ASKAPLOG_INFO_STR(logger, "There are " << itsNumberTaylor << " Taylor terms");
        }

        casa::Matrix<casa::Double> couplingMatrix(itsNumberTaylor,itsNumberTaylor);

        if(this->itsNumberTaylor>1) {
            iph.makeTaylorTerm(0);
        }
        else {
            ASKAPLOG_INFO_STR(logger, "No Taylor terms will be solved");
        }
        const casacore::IPosition imageShape = ip.shape(iph.paramName());
        const uint nPol = imageShape.nelements()>=3 ? uint(imageShape(2)) : 1;
        ASKAPLOG_INFO_STR(logger, "There are " << nPol << " polarisation planes to solve" );
        nParameters += imageShape.product(); // add up the number of pixels for zero order
        // check consistency
        for (uInt order=1; order<nOrders; ++order) {
            // make the helper a Taylor term of the given order
            iph.makeTaylorTerm(order);
            const casacore::IPosition thisShape = ip.shape(iph.paramName());
            const uint thisNPol = thisShape.nelements()>=3 ? uint(thisShape(2)) : 1;
            ASKAPCHECK(thisNPol == nPol,
                "Number of polarisations are supposed to be consistent for all Taylor terms, order="<<
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

        // Support for cleaning at higher resolution than that used in calcNE
        // Check if a separate full-resolution clean model has been saved
        // It will be in a param with the starting "image" swapped to "fullres"):
        bool importModelFromNE = true;
       if (itsExtraOversamplingFactor) {
            ASKAPDEBUGASSERT(*itsExtraOversamplingFactor > 1.);
            // if the 0-order param exists, assume they all do. If not, assume they all need to be initialised
            if (this->itsNumberTaylor>1) {
                iph.makeTaylorTerm(0);
            }
            string fullResName = iph.paramName();
            const size_t index = fullResName.find("image");
            ASKAPCHECK(index == 0, "Trying to swap to full-resolution param name but something is wrong");
            fullResName.replace(index,5,"fullres");
            if (ip.has(fullResName)) {
                // the full-resolution 0-order param already exists, so update those
                importModelFromNE = false;
            }
            else {
                // the full-resolution 0-order param does not exist
                // use the NE models, but first set up a new params for storing the results
                // set the shape of a full-resolution model and also an iterator.
                casacore::IPosition fullResShape(scimath::PaddingUtils::paddedShape(ip.shape(iph.paramName()),
                                                                                    *itsExtraOversamplingFactor));
                casacore::IPosition fullResIterShape(ip.shape(iph.paramName()));
                fullResIterShape(0) = fullResShape(0);
                fullResIterShape(1) = fullResShape(1);
                // copy the low-resolution axes from the 0-order param
                scimath::Axes axes(ip.axes(iph.paramName()));
                // update for the new resolution
                casacore::DirectionCoordinate radec = axes.directionAxis();
                /// @todo double check that the rounding is correct for the ref pixel
                radec.setReferencePixel(radec.referencePixel()*double(*itsExtraOversamplingFactor));
                radec.setIncrement(radec.increment()/double(*itsExtraOversamplingFactor));
                axes.addDirectionAxis(radec);
                // I should not need to save full-resolution models for order >= nTaylor (i.e. upto 2*nTaylor-1)
                for( uInt order=0; order < this->itsNumberTaylor; ++order) {
                    if(this->itsNumberTaylor>1) {
                        iph.makeTaylorTerm(order);
                    }
                    fullResName = iph.paramName();
                    const size_t index = fullResName.find("image");
                    ASKAPCHECK(index == 0, "Trying to swap to full-resolution param name but something is wrong");
                    fullResName.replace(index,5,"fullres");
                    ASKAPLOG_INFO_STR(logger, "Creating empty parameter "<<fullResName<<
                        " with shape "<<fullResIterShape);
                    // add the new param
                    ip.add(fullResName, fullResIterShape, axes);
                }
            }
        }

        // as polarisations are not necessarily represented by a different parameter
        // we have to build a set of parameters which are going to be fixed inside the loop
        // (or alternatively fix them multiple times, which is also a reasonable solution)
        std::set<std::string> parametersToBeFixed;

        bool firstcycle;

        // Iterate through Polarisations
        for (imagemath::MultiDimArrayPlaneIter planeIter(imageShape); planeIter.hasMore(); planeIter.next()) {
            const uint plane = planeIter.sequenceNumber();
            std::string tagLogString(planeIter.tag());
            if (tagLogString.size()) {
                tagLogString = "tagged as " + tagLogString;
            } else {
                tagLogString = "not tagged";
            }

            ASKAPLOG_INFO_STR(logger, "Preparing iteration for plane "<<plane<<
                " ("<<tagLogString<<") in image "<<tmIt->first);
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
            ASKAPCHECK(normalEquations().normalMatrixDiagonal().count(zeroOrderParam)>0,
                "Diagonal not present " << zeroOrderParam);
            casacore::Vector<imtype> normdiag(normalEquations().normalMatrixDiagonal().find(zeroOrderParam)->second);
            // Setup the preconditioner function vector, if it has been defined
            ASKAPCHECK(normalEquations().preconditionerSlice().count(zeroOrderParam)>0,
                "Preconditioner function Slice not present for " << zeroOrderParam);
            casacore::Vector<imtype> pcf(normalEquations().preconditionerSlice().find(zeroOrderParam)->second);

            ASKAPDEBUGASSERT(planeIter.planeShape().nelements()>=2);

            const double maxDiag = casacore::max(planeIter.getPlaneVector(normdiag));
            ASKAPLOG_INFO_STR(logger, "Maximum of weights = " << maxDiag );

            // a unique string for every Taylor decomposition (unique for every facet for faceting)
            const std::string imageTag = tmIt->first + planeIter.tag();
            firstcycle = !SynthesisParamsHelper::hasValue(itsCleaners,imageTag);

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
            // also keep a copy of an alternative preconditioner function, if it isn't empty.
            // should only need the zeroth order for preconditioning, so set that now.
            Array<float> pcfZeroArray;
            if (pcf.shape() > 0) {
                #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                pcfZeroArray.reference(planeIter.getPlane(pcf));
                #else
                pcfZeroArray.resize(planeIter.planeShape());
                casacore::convertArray<float, double>(pcfZeroArray, planeIter.getPlane(pcf));
                #endif
            }
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
                ASKAPLOG_INFO_STR(logger, "AMSMFS solver: processing order "<<order<<
                                  " ("<<itsNumberTaylor<<" Taylor terms + " << itsNumberTaylor-1 <<
                                  " cross-terms), parameter name: " << thisOrderParam);

                // Always need the actual PSFs, and on first cycle we also need the
                // higher order terms
                if ((order<this->itsNumberTaylor)||(firstcycle)) {
                    ASKAPCHECK(normalEquations().normalMatrixSlice().count(thisOrderParam)>0,
                        "PSF Slice for plane="<< plane<<" and order="<<order<<" is not present");
                    casacore::Vector<imtype> slice(normalEquations().normalMatrixSlice().find(thisOrderParam)->second);
                    #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                    if (itsExtraOversamplingFactor) {
                        // copy for oversampling later
                        psfLongVec(order) = planeIter.getPlane(slice);
                    } else {
                        // use a reference
                        psfLongVec(order).reference(planeIter.getPlane(slice));
                    }
                    #else
                    psfLongVec(order).resize(planeIter.planeShape());
                    casacore::convertArray<float, double>(psfLongVec(order), planeIter.getPlane(slice));
                    #endif
                }

                // For the dirty images, we need all terms
                ASKAPCHECK(normalEquations().dataVectorT(thisOrderParam).size()>0,
                    "Data vector not present for cube plane="<<plane<<" and order="<<order);
                casacore::Vector<imtype> dv = normalEquations().dataVectorT(thisOrderParam);
                #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                if (itsExtraOversamplingFactor) {
                    // copy for oversampling later
                    dirtyLongVec(order) = planeIter.getPlane(dv);
                } else {
                    // use a reference
                    dirtyLongVec(order).reference(planeIter.getPlane(dv));
                }
                #else
                dirtyLongVec(order).resize(planeIter.planeShape());
                casacore::convertArray<float, double>(dirtyLongVec(order), planeIter.getPlane(dv));
                #endif

                // For the clean images, we need only the first nTaylor
                if(importModelFromNE && (order < this->itsNumberTaylor)) {
                    //cleanVec(order).assign(planeIter.getPlane(ip.valueF(thisOrderParam)));
                    // if extraOS>1 should only get here if the model is empty. But do the following for consistency
                    if (itsExtraOversamplingFactor) {
                        // copy for oversampling later
                        cleanVec(order) = planeIter.getPlane(ip.valueF(thisOrderParam));
                    } else {
                        // use a reference
                        cleanVec(order).reference(planeIter.getPlane(ip.valueF(thisOrderParam)));
                    }
                }
            } // Loop over order

            casacore::Array<float> maskArray(planeIter.planeShape());

            uInt nx(planeIter.planeShape()(0));
            uInt ny(planeIter.planeShape()(1));
            IPosition centre(2, nx/2, ny/2);

            itsPSFZeroArray=psfLongVec(0);
            if (firstcycle) {
                ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._fc_norm+precnd");

                ASKAPLOG_DEBUG_STR(logger, "Deriving scale from PSF(0) centre value " <<
                    psfLongVec(0).nonDegenerate()(centre));
                // For the first cycle we need to precondition and normalise all PSFs and all dirty images
                itsPSFZeroCentre=-1;
                for(uInt order=0; order < 2 * itsNumberTaylor - 1; ++order) {
                    if (itsNumberTaylor == 1) {
                        // shortcut for the case with no Taylor terms
                        ASKAPLOG_INFO_STR(logger, "Preconditioning dirty image for plane=" << plane<<
                            " ("<<tagLogString<< ")");
                        doPreconditioning(psfLongVec(0), dirtyLongVec(0), pcfZeroArray);
                    } else {
                        // We need to work with the original PSF since it gets overridden
                        psfWorkArray = itsPSFZeroArray;
                        ASKAPLOG_DEBUG_STR(logger, "Initial PSF(" << order <<
                            ") centre value " << psfLongVec(order).nonDegenerate()(centre));
                        // Precondition this order PSF using PSF(0)
                        if(doPreconditioning(psfWorkArray, psfLongVec(order), pcfZeroArray)) {
                            ASKAPLOG_DEBUG_STR(logger, "After preconditioning PSF(" << order <<
                                ") centre value " << psfLongVec(order).nonDegenerate()(centre));
                            // Now we can precondition the dirty (residual) array using PSF(0)
                            psfWorkArray = itsPSFZeroArray;
                            ASKAPLOG_INFO_STR(logger, "Preconditioning dirty image for plane=" << plane<<
                                " ("<<tagLogString<< ") and order=" << order);
                            doPreconditioning(psfWorkArray,dirtyLongVec(order),pcfZeroArray);
                        }
                    }
                    // Normalise.
                    ASKAPLOG_DEBUG_STR(logger, "Normalising PSF and Dirty image for order " << order);
                    ASKAPLOG_DEBUG_STR(logger, "Before normalisation PSF(" << order <<
                        ") centre value " << psfLongVec(order).nonDegenerate()(centre));
                    // First call the scaling is via the psf and the value is returned.
                    // Thereafter we use that value for the normalisation of all the PSFs.
                    // Thus for MFS, the first PSF should have centre value 1.0 and the others lower values
                    if(order==0) {
                        itsPSFZeroCentre = doNormalization(planeIter.getPlaneVector(normdiag), tol(),
                            psfLongVec(0), dirtyLongVec(0),
                            boost::shared_ptr<casacore::Array<float> >(&maskArray, utility::NullDeleter()));
                    }  else {
                        doNormalization(planeIter.getPlaneVector(normdiag), tol(),
                            psfLongVec(order), itsPSFZeroCentre, dirtyLongVec(order),
                            boost::shared_ptr<casacore::Array<float> >(&maskArray, utility::NullDeleter()));
                    }
                        for (uInt t1=0; t1 < itsNumberTaylor; ++t1) {
                            for (uInt t2=0; t2 < itsNumberTaylor; ++t2) {
                                if (t1+t2 == order) {
                                    couplingMatrix(t1,t2) = psfLongVec(order).nonDegenerate()(centre);
                                }
                            }
                        }

                    ASKAPLOG_DEBUG_STR(logger, "After  normalisation PSF(" << order <<
                        ") centre value " << psfLongVec(order).nonDegenerate()(centre));
                    if (order<itsNumberTaylor) {
                        psfVec(order).reference(psfLongVec(order));
                        dirtyVec(order).reference(dirtyLongVec(order));
                    }
                }// Loop over order
            }
            else {
                ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._norm+precnd");
                // For the subsequent cycles cycle we need to precondition and normalise the updated dirty images
                // Precondition the dirty (residual) array
                for(uInt order=0; order < itsNumberTaylor; ++order) {
                    psfWorkArray = itsPSFZeroArray;
                    if(doPreconditioning(psfWorkArray,dirtyLongVec(order),pcfZeroArray)) {
                        ASKAPLOG_INFO_STR(logger, "Preconditioning dirty image for plane="<<plane<<
                            " ("<<tagLogString<< ") and order=" << order);
                    }
                    // Normalise.
                    psfWorkArray = itsPSFZeroArray;
                    doNormalization(planeIter.getPlaneVector(normdiag), tol(),
                        psfWorkArray, itsPSFZeroCentre, dirtyLongVec(order),
                        boost::shared_ptr<casacore::Array<float> >(&maskArray, utility::NullDeleter()));
                    if(order<itsNumberTaylor) {
                        dirtyVec(order)=dirtyLongVec(order);
                    }
                }// Loop over order
            }

            // sinc interpolate via Fourier padding if cleaning requires higher resolution
            if (itsExtraOversamplingFactor) {
                ASKAPLOG_INFO_STR(logger,
                    "Oversampling by an extra factor of "<<*itsExtraOversamplingFactor<<" before cleaning");
                SynthesisParamsHelper::oversample(maskArray,*itsExtraOversamplingFactor);
                for (uInt order=0; order < limit; ++order) {
                    SynthesisParamsHelper::oversample(psfLongVec(order),*itsExtraOversamplingFactor);
                    if(order < this->itsNumberTaylor) {
                        SynthesisParamsHelper::oversample(dirtyVec(order),*itsExtraOversamplingFactor);
                        psfVec(order).reference(psfLongVec(order));
                        if (importModelFromNE) {
                            SynthesisParamsHelper::oversample(cleanVec(order),*itsExtraOversamplingFactor,false);
                        } else {
                            if (this->itsNumberTaylor>1) {
                                iph.makeTaylorTerm(order);
                            }
                            string fullResName = iph.paramName();
                            const size_t index = fullResName.find("image");
                            ASKAPCHECK(index == 0, "Swapping to full-resolution param name but something is wrong");
                            fullResName.replace(index,5,"fullres");
                            imagemath::MultiDimArrayPlaneIter fullResPlaneIter(ip.shape(fullResName));
                            cleanVec(order).reference(
                                fullResPlaneIter.getPlane( ip.valueT(fullResName), planeIter.position() ) );
                        }
                    }

                }
            }

            // Now that we have all the required images, we can initialise the deconvolver
            if (firstcycle) {// Initialize everything only once.
                ASKAPTRACE("ImageAMSMFSolver::solveNormalEquations._fc_initdeconvolver");

                ASKAPLOG_INFO_STR(logger, "Creating solver for plane " << plane <<" tag "<<imageTag);
                itsCleaners[imageTag].reset(
                    new DeconvolverMultiTermBasisFunction<Float, Complex>(dirtyVec, psfVec, psfLongVec));
                ASKAPDEBUGASSERT(itsCleaners[imageTag]);

                itsCleaners[imageTag]->setMonitor(itsMonitor);
                itsControl->setTargetObjectiveFunction(threshold().getValue("Jy"));
                itsControl->setTargetObjectiveFunction2(deepThreshold());
                itsControl->setFractionalThreshold(fractionalThreshold());

                itsCleaners[imageTag]->setControl(itsControl);

                ASKAPCHECK(itsBasisFunction, "Basis function not initialised");

                ASKAPDEBUGASSERT(dirtyVec.nelements() > 0);
                itsBasisFunction->initialise(dirtyVec(0).shape());
                itsCleaners[imageTag]->setBasisFunction(itsBasisFunction);
                itsCleaners[imageTag]->setSolutionType(itsSolutionType);
                itsCleaners[imageTag]->setDecoupled(itsDecoupled);
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
                    // Save create/update parameters after preconditioning.
                    // Will need downsampling if params are stored with lower resolution
                    if (itsExtraOversamplingFactor) {
                        Array<Float> tmpImg = dirtyVec(order);
                        SynthesisParamsHelper::downsample(tmpImg,*itsExtraOversamplingFactor);
                        ASKAPLOG_DEBUG_STR(logger, "Dirty(" << order << ") shape = " << dirtyVec(order).shape());
                        saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "residual",
                            unpadImage(tmpImg), planeIter.position());
                    } else {
                        saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "residual",
                            unpadImage(dirtyVec(order)), planeIter.position());
                    }
                    if(firstcycle) {
                        if (itsExtraOversamplingFactor) {
                            ASKAPLOG_DEBUG_STR(logger, "PSF(" << order << ") shape = " << psfVec(order).shape());
                            Array<Float> tmpImg = psfVec(order);
                            SynthesisParamsHelper::downsample(tmpImg,*itsExtraOversamplingFactor);
                            saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "psf.image",
                                unpadImage(tmpImg), planeIter.position());
                        } else {
                            saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "psf.image",
                                unpadImage(psfVec(order)), planeIter.position());
                        }
                        if(order==0&&maskArray.nelements()) {
                            if (itsExtraOversamplingFactor) {
                                Array<Float> tmpImg = maskArray;
                                SynthesisParamsHelper::downsample(tmpImg,*itsExtraOversamplingFactor);
                                saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "mask",
                                    unpadImage(tmpImg), planeIter.position());
                            } else {
                                saveArrayIntoParameter(ip, thisOrderParam, planeIter.shape(), "mask",
                                    unpadImage(maskArray), planeIter.position());
                            }
                        }
                    }
                }
                // cleanVec already references valueF, except in Nyquist node where it has been copied & oversampled
                //itsCleaners[imageTag]->setModel(planeIter.getPlane(ip.valueF(thisOrderParam)), order);
                itsCleaners[imageTag]->setModel(cleanVec(order), order);
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
                // if we cleaned at higher resolution, need to update the high-res params then downsample back
                if (itsExtraOversamplingFactor) {
                    // copy the new model back to the working vector of arrays
                    Array<Float> tmpImg = itsCleaners[imageTag]->model(order);
                    // store full-res model in a slice of the new fullres parameter
                    ASKAPCHECK(planeIter.position()(0)==0 && planeIter.position()(1)==0,
                        "Image offsets not supported with variable image resolution");
                    string fullResName = thisOrderParam;
                    const size_t index = fullResName.find("image");
                    ASKAPCHECK(index == 0, "Swapping to full-resolution param name but something is wrong");
                    fullResName.replace(index,5,"fullres");
                    const uInt numDegenerate = ip.shape(fullResName).size() - tmpImg.ndim();
                    saveArrayIntoParameter(ip, thisOrderParam, ip.shape(fullResName),
                        "fullres", tmpImg.addDegenerate(numDegenerate), planeIter.position());
                    // remove Fourier padding before returning to degridders
                    SynthesisParamsHelper::downsample(tmpImg,*itsExtraOversamplingFactor);
                    planeIter.getPlane(ip.valueT(thisOrderParam)).nonDegenerate() =
                        unpadImage(tmpImg);
                } else {
                    planeIter.getPlane(ip.valueT(thisOrderParam)).nonDegenerate() =
                        unpadImage(itsCleaners[imageTag]->model(order));
                }
            }
            // add extra parameters (cross-terms) to the to-be-fixed list
            for (uInt order = itsNumberTaylor; order<nOrders; ++order) {
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

        if (firstcycle) {
            ASKAPLOG_DEBUG_STR(logger, "Inverting and cacheing coupling matrix");
            casa::Matrix<casa::Double> inverseMatrix(itsNumberTaylor,itsNumberTaylor);
            casa::Double det;
            invertSymPosDef(inverseMatrix, det, couplingMatrix);
            setInverseCouplingMatrix(inverseMatrix);
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

      // if not using taylor terms we can save memory by throwing away the cleaners here
      // the extra work of recreating them is minimal in that case
      // Disabled for now - this messes up deep cleaning because scale masks are lost
      //if (itsNumberTaylor==1) {
      //    ASKAPLOG_INFO_STR(logger,"Clearing out the cleaners");
      //    itsCleaners.clear();
      //}

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
      ASKAPLOG_INFO_STR(logger, "Solution type = " << solutionType);
      this->itsSolutionType=solutionType;

      this->itsOrthogonal=parset.getBool("orthogonal", false);
      if (this->itsOrthogonal) {
          ASKAPLOG_DEBUG_STR(logger, "Multiscale basis functions will be orthogonalised");
          // need to reset the basisFunctions
          const BasisFunction<Float>::ShPtr bfPtr(new MultiScaleBasisFunction<Float>(itsScales, itsOrthogonal));
          itsBasisFunction = bfPtr;
      }
      this->itsDecoupled = parset.getBool("decoupled", false);
      if (this->itsDecoupled) {
          ASKAPLOG_DEBUG_STR(logger, "Using decoupled residuals");
      }

    }
  }
}
