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
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>
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
            itsBeamHelper(beamHelper), itsEqualiseNoise(false), itsModelNeedsConvolving(true),
            itsResidualNeedsUpdating(true), itsSaveRawPsf(false), itsSavePsfImage(false)
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
                nParameters+=ip.valueT(name).nelements();
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


            // need to make sure this works with facetting
            // this should work for faceting as well, taylorMap would contain one element
            // per facet in this case
            std::map<std::string, int> taylorMap;
            SynthesisParamsHelper::listTaylor(names, taylorMap);
            ASKAPCHECK(taylorMap.size() != 0, "Solver doesn't have any images to solve for");

            uInt nOrders = 1;

            for (std::map<std::string, int>::const_iterator tmIt = taylorMap.begin();
                                                            tmIt != taylorMap.end(); ++tmIt) {
                // create a helper
                ImageParamsHelper iph(tmIt->first);
                // nOrders is the number of separate image cubes (e.g. Taylor terms)
                nOrders = tmIt->second;

                // obtain name with just the taylor suffix, if present (i.e. no facet suffixes)
                if(nOrders>1) {
                    // facetmap has been setup with the full names, so need to use those here
                    iph.makeTaylorTerm(0);
                }
                const std::string name = iph.taylorName();

                if (facetmap[name] == 1) {

                    // this is not a faceting case, restore the image in situ and add residuals
                    if(nOrders == 1) {
                        ASKAPLOG_INFO_STR(logger, "Restoring " << tmIt->first );
                    } else {
                        ASKAPLOG_INFO_STR(logger, "Restoring " << nOrders << " terms of " << tmIt->first );
                    }

                    // convolve with restoring beam before adding residuals, but wait until after
                    // preconditioning to calculate restoring beam size.
                    // needs to be set before addResiduals to support facets
                    itsModelNeedsConvolving = true;

                    // add residuals
                    addResiduals(ip, tmIt->first, nOrders);

                    for( uInt order=0; order < nOrders; ++order) {
                        if(nOrders>1) {
                            iph.makeTaylorTerm(order);
                        }
                        SynthesisParamsHelper::setBeam(ip, iph.paramName(), itsBeamHelper.value());
                        if (itsExtraOversamplingFactor) {
                            // also setBeam for the fullres param
                            std::string fullresname = iph.paramName();
                            const size_t index = fullresname.find("image");
                            ASKAPCHECK(index == 0, "Swapping to full-resolution param name but something is wrong");
                            fullresname.replace(index,5,"fullres");
                            ASKAPASSERT(ip.has(fullresname));
                            SynthesisParamsHelper::setBeam(ip, fullresname, itsBeamHelper.value());
                        }
                    }

                } else {

                    for( uInt order=0; order < nOrders; ++order) {
                        if(nOrders>1) {
                            iph.makeTaylorTerm(order);
                        }

                        // this is a single facet of a larger image, just fill in the bigger image with the model
                        ASKAPLOG_INFO_STR(logger, "Inserting facet " << iph.paramName()<<" into merged image "<<name);
                        casa::Array<imtype> patch = SynthesisParamsHelper::getFacet(ip,iph.paramName());
                        const casa::Array<imtype> model = scimath::PaddingUtils::centeredSubArray(ip.valueT(iph.paramName()),
                        patch.shape());
                        patch = model;

                    }

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
                            ASKAPLOG_INFO_STR(logger,
                                "Adding residuals for facet ("<<xFacet<<","<<yFacet<<")");
                            // ci->first may have taylor suffix defined, load it first and then add facet indices
                            ImageParamsHelper iph(ci->first);
                            iph.makeFacet(xFacet,yFacet);
                            addResiduals(ip, ci->first, nOrders, iph.paramName());

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
            const std::string &basename,
            const casa::uInt nOrders,
            const std::string facetname)

        {
            ASKAPTRACE("ImageRestoreSolver::addResiduals");

            // name with the taylor-related suffixes removed (and facet suffixes preserved, if present)
            ImageParamsHelper iph(basename);

            casa::Array<float> pcfArray;
            casa::Array<float> psfWorkArray;
            casa::Array<float> psfZeroArray;

            // initialise some arrays
            bool decouplingRequired = true;
            casa::Matrix<double> inverseCouplingMatrix;
            if(nOrders>1) {
                iph.makeTaylorTerm(0);
                if (getInverseCouplingMatrix().shape() > 0) {
                    inverseCouplingMatrix = getInverseCouplingMatrix();
                    ASKAPCHECK((inverseCouplingMatrix.nrow() == nOrders) &&
                               (inverseCouplingMatrix.ncolumn() == nOrders),
                        "Inconsistent inv coupling matrix. " << inverseCouplingMatrix.shape() <<
                        " != [" << nOrders << ", " << nOrders << "]");
                    ASKAPLOG_INFO_STR(logger, "Inverse coupling matrix = " << inverseCouplingMatrix.row(0));
                    for (uInt order = 1; order < nOrders; order++) {
                        ASKAPLOG_INFO_STR(logger, "                          " << inverseCouplingMatrix.row(order));
                    }
                } else {
                    decouplingRequired = false;
                    ASKAPLOG_WARN_STR(logger, "No decoupling matrix available. Restoring Taylor terms as they come.");
                }
            } else {
                decouplingRequired = false;
            }
            std::string imagename = iph.taylorName();
            std::string name;
            if (facetname == "") {
                name = imagename;
            } else {
                name = facetname;
            }
            const casa::IPosition shape = ip.shape(name);

            // The following need to be done once and only once for each order (i.e. once for all planeIter)
            // So set to true here then to false at the end of the first planeIter loop.
            // Someone should check this...
            // Axes are dof, dof for each parameter
            //casa::IPosition vecShape(1, out.shape().product());
            bool saveNewPSFRequired = true;

            for (imagemath::MultiDimArrayPlaneIter planeIter(shape); planeIter.hasMore(); planeIter.next()) {

                float itsPSFZeroCentre;

                // Set up an accumulation array for Taylor term decoupling
                casa::Vector<casa::Array<float> > dirtyVector(nOrders);
                casa::Vector<casa::Array<float> > psfVector(nOrders);
                if (decouplingRequired) {
                    for( uInt order=0; order < nOrders; ++order) {
                        dirtyVector(order).resize(planeIter.planeShape());
                        dirtyVector(order).set(0.0);
                    }
                }

                for( uInt order=0; order < nOrders; ++order) {
                    if(nOrders>1) {
                        iph.makeTaylorTerm(order);
                    }
                    std::string imagename = iph.taylorName();
                    ASKAPLOG_INFO_STR(logger, "Preparing to restore "<<imagename);

                    // when convolving model images with the restoring beam and
                    // adding residuals, the full, merged model needs to be used.
                    std::string name;
                    if (facetname == "") {
                        name = imagename;
                    } else {
                        name = facetname;
                    }

                    ASKAPCHECK(normalEquations().normalMatrixDiagonal().count(name)>0,
                        "Diagonal not present " << name);
                    casa::Vector<imtype> diag(normalEquations().normalMatrixDiagonal().find(name)->second);
                    ASKAPCHECK(normalEquations().dataVectorT(name).size()>0,
                        "Data vector not present " << name);
                    casa::Vector<imtype> dv = normalEquations().dataVectorT(name);
                    ASKAPCHECK(normalEquations().normalMatrixSlice().count(name)>0,
                        "PSF Slice not present " << name);
                    casa::Vector<imtype> slice(normalEquations().normalMatrixSlice().find(name)->second);
                    ASKAPCHECK(normalEquations().preconditionerSlice().count(name)>0,
                        "Preconditioner fuction Slice not present for " << name);
                    casa::Vector<imtype> pcf(normalEquations().preconditionerSlice().find(name)->second);

                    if (planeIter.tag()!="") {
                        // it is not a single plane case, there is something to report
                        ASKAPLOG_INFO_STR(logger, "Processing plane "<<planeIter.sequenceNumber()<<
                            " tagged as "<<planeIter.tag());

                    }

                    ASKAPLOG_INFO_STR(logger, "Maximum of data vector corresponding to "<<name<<" is "<<casa::max(dv));

                    #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                    casa::Array<float> dirtyArray(planeIter.getPlane(dv));
                    #else
                    casa::Array<float> dirtyArray(planeIter.planeShape());
                    casa::convertArray<float, double>(dirtyArray,planeIter.getPlane(dv));
                    #endif

                    ASKAPLOG_INFO_STR(logger, "Maximum of data vector corresponding to "<<name<<" and plane "<<
                        planeIter.sequenceNumber()<<" is "<<casa::max(dirtyArray));

                    #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                    casa::Array<float> psfArray(planeIter.getPlane(slice));
                    #else
                    casa::Array<float> psfArray(planeIter.planeShape());
                    casa::convertArray<float, double>(psfArray, planeIter.getPlane(slice));
                    #endif

                    // send an alternative preconditioner function, if it isn't empty.
                    // whether the psf or pcf is used, preconditioning needs to be the same for all Taylor terms
                    if (order == 0) {
                        if (nOrders > 0) {
                            psfZeroArray = psfArray;
                        }
                        if (pcf.shape() > 0) {
                            ASKAPDEBUGASSERT(pcf.shape() == slice.shape());
                            #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                            pcfArray.reference(planeIter.getPlane(pcf));
                            #else
                            pcfArray.resize(planeIter.planeShape());
                            casa::convertArray<float, double>(pcfArray, planeIter.getPlane(pcf));
                            #endif
                        }
                    }

                    // uninitialised mask shared pointer means that we don't need it (i.e. no weight equalising)
                    boost::shared_ptr<casa::Array<float> > mask;
                    if (itsEqualiseNoise) {
                        ASKAPLOG_INFO_STR(logger,
                            "Residual will be multiplied by sqrt(normalised weight) during restoration");
                        // mask will have a normalised sqrt(weight) pattern after doNormalization
                        mask.reset(new casa::Array<float>(dirtyArray.shape()));
                    } else {
                        ASKAPLOG_INFO_STR(logger,
                            "Restored image will have primary beam corrected noise (no equalisation)");
                    }
                    // Save unnormalised PSF
                    if (itsSaveRawPsf) {
                        saveArrayIntoParameter(ip, name, shape, "psf.raw", psfArray,
                        planeIter.position());
                    }
                    // Do the preconditioning
                    if (nOrders == 1) {
                        // avoid the extra work & memory if we only have 1 taylor term
                        doPreconditioning(psfArray,dirtyArray,pcfArray);
                    } else {
                        psfWorkArray = psfZeroArray;
                        doPreconditioning(psfWorkArray,dirtyArray,pcfArray);
                        psfWorkArray = psfZeroArray;
                        doPreconditioning(psfWorkArray,psfArray,pcfArray);
                    }

                    // Normalize by the diagonal of Taylor term 0
                    if (order == 0) {
                        itsPSFZeroCentre = doNormalization(planeIter.getPlaneVector(diag),tol(),psfArray,dirtyArray,mask);
                    } else {
                        doNormalization(planeIter.getPlaneVector(diag),tol(),psfArray,itsPSFZeroCentre,dirtyArray,mask);
                    }

                    // we have to do noise equalisation for final residuals after preconditioning
                    if (itsEqualiseNoise) {
                        const casa::IPosition vecShape(1,dirtyArray.nelements());
                        casa::Vector<float> tmpVector(dirtyArray.reform(vecShape));
                        ASKAPDEBUGASSERT(mask);
                        const casa::Vector<float> maskVector(mask->reform(vecShape));
                        for (int i = 0; i<vecShape[0]; ++i) {
                            tmpVector[i] *= maskVector[i];
                        }
                    }

                    // If doing Taylor imaging, loop over decoupled orders and accumulate residuals
                    if (decouplingRequired) {
                        for( uInt dorder=0; dorder < nOrders; ++dorder) {
                            //dirtyVector(dorder) += inverseCouplingMatrix(dorder,order) * dirtyArray;
                            const casa::IPosition vecShape(1,dirtyArray.nelements());
                            casa::Vector<float> inVector(dirtyArray.reform(vecShape));
                            casa::Vector<float> outVector(dirtyVector(dorder).reform(vecShape));
                            for (int i = 0; i<vecShape[0]; ++i) {
                                outVector[i] += inVector[i] * inverseCouplingMatrix(dorder,order);
                            }
                        }
                    } else {
                        // otherwise just point at the dirtyArray.
                        dirtyVector(order).reference(dirtyArray);
                    }

                    // don't decouple PSF images
                    psfVector(order).reference(psfArray);

                }

                for( uInt order=0; order < nOrders; ++order) {

                    casa::Array<float> dirtyArray(dirtyVector(order));
                    casa::Array<float> psfArray(psfVector(order));

                    if(nOrders>1) {
                        iph.makeTaylorTerm(order);
                    }
                    std::string imagename = iph.taylorName();
                    ASKAPLOG_INFO_STR(logger, "Restoring "<<imagename);

                    // when convolving model images with the restoring beam and
                    // adding residuals, the full, merged model needs to be used.
                    std::string name;
                    casa::Array<imtype> out;
                    if (facetname == "") {
                        name = imagename;
                        // Check if a separate full-resolution clean model has been saved
                        // It will be in a param with the starting "image" swapped to "fullres"):
                        if (itsExtraOversamplingFactor) {
                            const size_t index = imagename.find("image");
                            ASKAPCHECK(index == 0, "Swapping to full-resolution param name but something is wrong");
                            imagename.replace(index,5,"fullres");
                            ASKAPASSERT(ip.has(imagename));
                            out.reference(ip.valueT(imagename));
                        } else {
                            out.reference(ip.valueT(name));
                        }
                    } else {
                        name = facetname;
                        out.reference(SynthesisParamsHelper::getFacet(ip,name));
                    }

                    if (itsResidualNeedsUpdating) {
                        // Store the current dirtyImage parameter class to be saved to disk later
                        ASKAPLOG_INFO_STR(logger, "Saving current residual image to model parameter");
                        ASKAPLOG_INFO_STR(logger, "Shape is " << dirtyArray.shape() << " position is " <<
                            planeIter.position());
                        saveArrayIntoParameter(ip,name,dirtyArray.shape(),"residual",
                            dirtyArray,planeIter.position());
                    }

                    if (itsSavePsfImage) {
                        if (saveNewPSFRequired) {
                            // Store the new PSF in parameter class to be saved to disk later
                            ASKAPLOG_INFO_STR(logger,
                                "Saving new PSF parameter as model NEW parameter -- needs full shape");
                            saveArrayIntoParameter(ip, name, shape, "psf.image", psfArray, planeIter.position());
                        } else {
                            ASKAPLOG_INFO_STR(logger,
                                "Saving new PSF parameter as model EXISTING parameter -- using plane shape");
                                saveArrayIntoParameter(ip, name, psfArray.shape(), "psf.image", psfArray,
                                    planeIter.position());
                        }
                    }

                    // sinc interpolate via Fourier padding if restoring requires higher resolution
                    if (itsExtraOversamplingFactor) {
                        ASKAPLOG_INFO_STR(logger,
                            "Oversampling dirty image and PSF by an extra factor of "<<*itsExtraOversamplingFactor);
                        SynthesisParamsHelper::oversample(dirtyArray,*itsExtraOversamplingFactor);
                        SynthesisParamsHelper::oversample(psfArray,*itsExtraOversamplingFactor);
                    }

                    // Add the residual image
                    // First, convolve the model image to the resolution of the synthesised beam if not already done.
                    casa::Vector<casa::Quantum<double> > restoringBeam;
                    if (itsModelNeedsConvolving) {

                        if (itsBeamHelper.fitRequired()) {
                            ASKAPLOG_INFO_STR(logger, "Fitting of Restoring beam required");
                            #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                            casa::Array<float>& psfDArray(psfArray);
                            #else
                            casa::Array<double> psfDArray(psfArray.shape());
                            casa::convertArray<double, float>(psfDArray, psfArray);
                            #endif
                            ASKAPLOG_INFO_STR(logger, "Fitting restoring beam");
                            // this was set by "name" for facets, but surely it should be equal for all facets
                            const scimath::Axes axes = ip.axes(imagename);
                            restoringBeam = SynthesisParamsHelper::fitBeam(psfDArray, axes, itsBeamHelper.cutoff());
                            ASKAPDEBUGASSERT(restoringBeam.size() == 3);
                            ASKAPLOG_INFO_STR(logger, "Restore solver will convolve with the 2D gaussian: " <<
                                restoringBeam[0].getValue("arcsec") << " x "<<restoringBeam[1].getValue("arcsec") <<
                                " arcsec at position angle "<<restoringBeam[2].getValue("deg")<<" deg");
                            itsBeamHelper.assign(restoringBeam);
                        } else {
                            restoringBeam = itsBeamHelper.value();
                        }
                        ASKAPLOG_INFO_STR(logger,
                            "Convolving the model image to the resolution of the synthesised beam");
                        // Create a temporary image
                        boost::shared_ptr<casa::TempImage<float> >
                            image(SynthesisParamsHelper::tempImage(ip, imagename));

                        askap::synthesis::Image2DConvolver<float> convolver;
                        const casa::IPosition pixelAxes(2, 0, 1);
                        convolver.convolve(*image, *image, casa::VectorKernel::GAUSSIAN,
                                           pixelAxes, restoringBeam, true, 1.0, false);
                        SynthesisParamsHelper::update(ip, imagename, *image);
                        // for some reason update makes the parameter free as well
                        ip.fix(imagename);
                    }

                    // The following can probably be done once.
                    // Even for facets, if addResiduals is called separately for each facet

                    // The code below involves an extra copying. We can replace it later with a copyless version
                    // doing element by element adding explicitly.
                    const casa::IPosition outSliceShape = planeIter.planeShape(out.shape());
                    // convertedResidual contains just one plane of residuals
                    #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                    casa::Array<float> convertedResidual(scimath::PaddingUtils::centeredSubArray(dirtyArray,
                        outSliceShape));
                    #else
                    casa::Array<double> convertedResidual(outSliceShape);
                    convertArray(convertedResidual, scimath::PaddingUtils::centeredSubArray(dirtyArray,
                        outSliceShape));
                    #endif
                    // figure out where to put the slice to (can't use planeIter functionality directly because out
                    // array can have a different shape
                    casa::IPosition blc(out.shape().nelements(),0);
                    casa::IPosition trc(out.shape());
                    const casa::IPosition curPos(planeIter.position());
                    for (casa::uInt dim = 0; dim<trc.nelements(); ++dim) {
                        trc[dim] -= 1;
                        ASKAPDEBUGASSERT(trc[dim]<out.shape()[dim]);
                        if ( (dim>=2) && (dim<curPos.nelements()) ) {
                            blc[dim] = curPos[dim];
                            trc[dim] = curPos[dim];
                        }
                    }
                    // copy stuff
                    Array<imtype> outSlice = out(blc,trc);
                    outSlice += convertedResidual;

                } // order loop

                // these only need to be done once per parameter
                saveNewPSFRequired = false;
                itsModelNeedsConvolving = false;

            } // planeIter loop
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
                    "beam parameter should be either equal to 'fit' or contain 3 elements defining the beam size."<<
                    " You have "<<beam[0]);
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
            const bool saverawpsf = parset.getBool("saverawpsf",false);
            result->saveRawPsf(saverawpsf);
            const bool savepsfimage = parset.getBool("savepsfimage",false);
            result->savePsfImage(savepsfimage);

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
        /// @param[in] doPrecon, copy the PreConditioners from other solver
        void ImageRestoreSolver::configureSolver(const ImageSolver &ts, bool doPrecon)
        {
            setThreshold(ts.threshold());
            setVerbose(ts.verbose());
            setTol(ts.tol());

            // behavior in the weight cutoff area
            zeroWeightCutoffMask(ts.zeroWeightCutoffMask());
            zeroWeightCutoffArea(ts.zeroWeightCutoffArea());

            //preconditioners
            if (doPrecon) {
                setPreconditioners(ts);
            }
        }

    } // namespace synthesis
} // namespace askap
