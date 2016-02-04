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

// Package level header file
#include <askap_synthesis.h>

// ASKAPsoft includes
#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.awprojectvisgridder");
#include <gridding/AWProjectVisGridder.h>
#include <fft/FFTWrapper.h>
#include <utils/PaddingUtils.h>
#include <casa/Arrays/ArrayIter.h>
#include <casa/BasicSL/Complex.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>
#include <measures/Measures/MDirection.h>
#include <casa/Quanta/MVDirection.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Quanta/MVTime.h>

// Local package includes
#include <measurementequation/SynthesisParamsHelper.h>
#include <gridding/SupportSearcher.h>
#include <profile/AskapProfiler.h>


namespace askap {
namespace synthesis {

AWProjectVisGridder::AWProjectVisGridder(const boost::shared_ptr<IBasicIllumination const> &illum,
        const double wmax, const int nwplanes,
        const double cutoff, const int overSample,
        const int maxSupport, const int limitSupport,
        const int maxFeeds, const int maxFields, const double pointingTol,
        const double paTol, const double freqTol,
        const bool frequencyDependent, const std::string& name) :
        AProjectGridderBase(maxFeeds, maxFields, pointingTol, paTol, freqTol),
        WProjectVisGridder(wmax, nwplanes, cutoff, overSample, maxSupport, limitSupport, name),
        itsReferenceFrequency(0.0), itsIllumination(illum),
        itsFreqDep(frequencyDependent),
        itsMaxFeeds(maxFeeds), itsMaxFields(maxFields)
{
    ASKAPDEBUGASSERT(itsIllumination);
    ASKAPCHECK(maxFeeds > 0, "Maximum number of feeds must be one or more");
    ASKAPCHECK(maxFields > 0, "Maximum number of fields must be one or more");
    ASKAPCHECK(overSample > 0, "Oversampling must be greater than 0");
    ASKAPCHECK(maxSupport > 0, "Maximum support must be greater than 0")
    setTableName(name);
}

/// @brief copy constructor
/// @details It is required to decouple internal array arrays, otherwise
/// those arrays are shared between all cloned gridders of this type
/// @note illumination model is copied as a pointer, so the same model is referenced
/// @param[in] other input object
AWProjectVisGridder::AWProjectVisGridder(const AWProjectVisGridder &other) :
        IVisGridder(other), AProjectGridderBase(other), WProjectVisGridder(other),
        itsReferenceFrequency(other.itsReferenceFrequency),
        itsIllumination(other.itsIllumination), itsFreqDep(other.itsFreqDep),
        itsMaxFeeds(other.itsMaxFeeds), itsMaxFields(other.itsMaxFields) {}

/// Clone a copy of this Gridder
IVisGridder::ShPtr AWProjectVisGridder::clone()
{
    return IVisGridder::ShPtr(new AWProjectVisGridder(*this));
}

/// Initialize the indices into the cube.
void AWProjectVisGridder::initIndices(const accessors::IConstDataAccessor& acc)
{
    ASKAPTRACE("AWProjectVisGridder::initIndices");    
    // calculate currentField
    indexField(acc);

    /// We have to calculate the lookup function converting from
    /// row and channel to plane of the w-dependent convolution
    /// function
    const int nSamples = acc.nRow();
    const int nChan = acc.nChannel();

    const int nPol = acc.nPol();
    itsCMap.resize(nSamples, nPol, nChan);
    itsCMap.set(0);

    const casa::Vector<casa::RigidVector<double, 3> > &rotatedUVW = acc.rotatedUVW(getTangentPoint());

    for (int i = 0; i < nSamples; ++i) {
        const int feed = acc.feed1()(i);
        ASKAPCHECK(feed < itsMaxFeeds, "Exceeded specified maximum number of feeds");
        ASKAPCHECK(feed > -1, "Illegal negative feed number");

        const double w = (rotatedUVW(i)(2)) / (casa::C::c);

        for (int chan = 0; chan < nChan; ++chan) {
            const double freq = acc.frequency()[chan];
            const int iw = getWPlane(w * freq);

            for (int pol = 0; pol < nPol; ++pol) {
                /// Order is (iw, chan, feed)
                if (itsFreqDep) {
                    itsCMap(i, pol, chan) = iw + nWPlanes() * (chan + nChan * (feed + itsMaxFeeds * currentField()));
                    ASKAPCHECK(itsCMap(i, pol, chan) < nWPlanes()*itsMaxFeeds*itsMaxFields*nChan,
                               "CMap index too large");
                    ASKAPCHECK(itsCMap(i, pol, chan) > -1, "CMap index less than zero");
                } else {
                    itsCMap(i, pol, chan) = iw + nWPlanes() * (feed + itsMaxFeeds * currentField());
                    ASKAPCHECK(itsCMap(i, pol, chan) < nWPlanes()*itsMaxFeeds*itsMaxFields,
                               "CMap index too large");
                    ASKAPCHECK(itsCMap(i, pol, chan) > -1, "CMap index less than zero");
                }
            }
        }
    }
}
/// @brief initialise sum of weights
/// @details We keep track the number of times each convolution function is used per
/// channel and polarisation (sum of weights). This method is made virtual to be able
/// to do gridder specific initialisation without overriding initialiseGrid.
/// This method accepts no parameters as itsShape, itsNWPlanes, etc should have already
/// been initialised by the time this method is called.
void AWProjectVisGridder::initialiseSumOfWeights()
{
    // this method is hopefully just a temporary stub until we figure out a better way of
    // managing a cache of convolution functions. It skips initialisation if itsSupport is
    // not zero, which means that some initialisation has been done before.
    // Note, it is not a very good way of doing things!
    if (itsSupport == 0) {
        WProjectVisGridder::initialiseSumOfWeights();
    }

    // Reset the weights
    zeroSumOfWeights();
}

/// @brief Initialise the gridding
/// @param axes axes specifications
/// @param shape Shape of output image: u,v,pol,chan
/// @param dopsf Make the psf?
void AWProjectVisGridder::initialiseGrid(const scimath::Axes& axes,  const casa::IPosition& shape,
        const bool dopsf, const bool dopcf)
{
    ASKAPTRACE("AWProjectVisGridder::initialiseGrid");
    WProjectVisGridder::initialiseGrid(axes, shape, dopsf, dopcf);

    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casa::uInt nx = std::min(maxSupport(), int(itsShape(0)));
    const casa::uInt ny = std::min(maxSupport(), int(itsShape(1)));

    ASKAPLOG_DEBUG_STR(logger, "Shape for calculating gridding convolution function = "
                           << nx << " by " << ny << " pixels");

    // this is just a buffer in the uv-space, oversampling is
    // taken into account inside the UVPattern object (in the past we handled
    // oversampling explicitly by using qnx and qny instead of nx and ny and
    // passing 1 instead of itsOverSample, but it caused scaling problems for
    // offset feeds).
    initUVPattern(nx, ny, itsUVCellSize(0), itsUVCellSize(1), itsOverSample);

    // this is a buffer for full-sized convolution function (nx by ny) before
    // a support is cut out. We initialise it here to put intensive operation
    // out of the loop.
    initCFBuffer(nx, ny);
}

/// @brief Initialise the degridding
/// @param axes axes specifications
/// @param image Input image: cube: u,v,pol,chan
void AWProjectVisGridder::initialiseDegrid(const scimath::Axes& axes,
        const casa::Array<double>& image)
{
    ASKAPTRACE("AWProjectVisGridder::initialiseDegrid");
    WProjectVisGridder::initialiseDegrid(axes, image);
    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casa::uInt nx = std::min(maxSupport(), int(itsShape(0)));
    const casa::uInt ny = std::min(maxSupport(), int(itsShape(1)));

    ASKAPLOG_DEBUG_STR(logger, "Shape for calculating degridding convolution function = "
                           << nx << " by " << ny << " pixels");

    // this is just a buffer in the uv-space, oversampling is
    // taken into account inside the UVPattern object (in the past we handled
    // oversampling explicitly by using qnx and qny instead of nx and ny and
    // passing 1 instead of itsOverSample, but it caused scaling problems for
    // offset feeds).
    initUVPattern(nx, ny, itsUVCellSize(0), itsUVCellSize(1), itsOverSample);

    // this is a buffer for full-sized convolution function (nx by ny) before
    // a support is cut out. We initialise it here to put intensive operation
    // out of the loop.
    initCFBuffer(nx, ny);
}


/// Initialize the convolution function into the cube. If necessary this
/// could be optimized by using symmetries.
/// @todo Make initConvolutionFunction more robust
void AWProjectVisGridder::initConvolutionFunction(const accessors::IConstDataAccessor& acc)
{
    ASKAPTRACE("AWProjectVisGridder::initConvolutionFunction");    
    casa::MVDirection out = getImageCentre();
    const int nSamples = acc.nRow();

    ASKAPDEBUGASSERT(itsIllumination);
    // just to avoid a repeated call to a virtual function from inside the loop
    const bool hasSymmetricIllumination = itsIllumination->isSymmetric();

    validateCFCache(acc, hasSymmetricIllumination);

    /// We have to calculate the lookup function converting from
    /// row and channel to plane of the w-dependent convolution
    /// function
    const int nChan = itsFreqDep ? acc.nChannel() : 1;

    if (itsSupport == 0) {
        itsConvFunc.resize(itsOverSample*itsOverSample*nWPlanes()*itsMaxFeeds*itsMaxFields*nChan);
        resizeSumOfWeights(nWPlanes()*itsMaxFeeds*itsMaxFields*nChan);
        zeroSumOfWeights();

        if (isOffsetSupportAllowed()) {
            initConvFuncOffsets(nWPlanes()*itsMaxFeeds*itsMaxFields*nChan);
        }
    }

    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casa::uInt nx = std::min(maxSupport(), int(itsShape(0)));
    const casa::uInt ny = std::min(maxSupport(), int(itsShape(1)));

    const casa::uInt qnx = nx / itsOverSample;
    const casa::uInt qny = ny / itsOverSample;

    // Find the actual cellsizes in x and y (radians)
    // corresponding to the limited support
    const double ccellx = 1.0 / (double(qnx) * itsUVCellSize(0));
    const double ccelly = 1.0 / (double(qny) * itsUVCellSize(1));

    /*
    casa::Vector<double> ccfx(nx);
    casa::Vector<double> ccfy(ny);
    for (casa::uInt ix = 0; ix < nx; ++ix) {
        const double nux = std::abs(double(ix) - double(nx / 2)) / double(nx / 2);
        ccfx(ix) = grdsf(nux); // /double(qnx);
    }
    for (casa::uInt iy = 0; iy < ny; ++iy) {
        const double nuy = std::abs(double(iy) - double(ny / 2)) / double(ny / 2);
        ccfy(iy) = grdsf(nuy); // /double(qny);
    }
    */

    UVPattern &pattern = uvPattern();
    casa::Matrix<casa::DComplex> thisPlane = getCFBuffer();
    ASKAPDEBUGASSERT(thisPlane.nrow() == nx);
    ASKAPDEBUGASSERT(thisPlane.ncolumn() == ny);

    int nDone = 0;

    for (int row = 0; row < nSamples; ++row) {
        const int feed = acc.feed1()(row);

        if (!isCFValid(feed, currentField())) {
            makeCFValid(feed, currentField());
            nDone++;
            casa::MVDirection offset(acc.pointingDir1()(row).getAngle());
            rwSlopes()(0, feed, currentField()) = isPSFGridder() || isPCFGridder() ? 0. : sin(offset.getLong()
                                                  - out.getLong()) * cos(offset.getLat());
            rwSlopes()(1, feed, currentField()) = isPSFGridder() || isPCFGridder() ? 0. : sin(offset.getLat())
                                                  * cos(out.getLat()) - cos(offset.getLat()) * sin(out.getLat())
                                                  * cos(offset.getLong() - out.getLong());

            const double parallacticAngle = hasSymmetricIllumination ? 0. : acc.feed1PA()(row);

            for (int chan = 0; chan < nChan; ++chan) {

                /// Extract illumination pattern for this channel
                itsIllumination->getPattern(acc.frequency()[chan], pattern,
                                            rwSlopes()(0, feed, currentField()),
                                            rwSlopes()(1, feed, currentField()), parallacticAngle);

                scimath::fft2d(pattern.pattern(), false);


                /// Calculate the total convolution function including
                /// the w term and the antenna convolution function

                for (int iw = 0; iw < nWPlanes(); ++iw) {
                    thisPlane.set(0.0);


                    // Loop over the central nx, ny region, setting it to the product
                    // of the phase screen and the spheroidal function
                    double maxCF = 0.0;
                    const double w = 2.0f * casa::C::pi * getWTerm(iw);
                    //std::cout<<"plane "<<iw<<" w="<<w<<std::endl;

                    for (int iy = 0; iy < int(ny); ++iy) {
                        const double y2 = casa::square((double(iy) - double(ny) / 2) * ccelly);

                        for (int ix = 0; ix < int(nx); ++ix) {
                            const double x2 = casa::square((double(ix) - double(nx) / 2) * ccellx);
                            const double r2 = x2 + y2;

                            if (r2 < 1.0) {
                                const double phase = w * (1.0 - sqrt(1.0 - r2));
                                // grid correction is temporary disabled as otherwise the fluxes are overestimated
                                const casa::DComplex wt = pattern(ix, iy) * conj(pattern(ix, iy));
                                //*casa::DComplex(ccfx(ix)*ccfy(iy));
                                // this ensures the oversampling is done
                                thisPlane(ix, iy) = wt * casa::DComplex(cos(phase), -sin(phase));
                                //thisPlane(ix, iy)=wt*casa::DComplex(cos(phase));
                                maxCF += casa::abs(wt);
                            }
                        }
                    }

                    ASKAPCHECK(maxCF > 0.0, "Convolution function is empty");


                    // At this point, we have the phase screen multiplied by the spheroidal
                    // function, sampled on larger cellsize (itsOverSample larger) in image
                    // space. Only the inner qnx, qny pixels have a non-zero value

                    // Now we have to calculate the Fourier transform to get the
                    // convolution function in uv space
                    scimath::fft2d(thisPlane, true);

                    // Now correct for normalization of FFT
                    thisPlane *= casa::DComplex(1.0 / (double(nx) * double(ny)));
                    // use this norm later on during normalisation
                    const double thisPlaneNorm = sum(real(thisPlane));
                    ASKAPDEBUGASSERT(thisPlaneNorm > 0.);

                    const int zIndex = iw + nWPlanes() * (chan + nChan * (feed + itsMaxFeeds * currentField()));

                    // If the support is not yet set, find it and size the
                    // convolution function appropriately

                    // by default the common support without offset is used
                    CFSupport cfSupport(itsSupport);

                    if (isSupportPlaneDependent() || (itsSupport == 0)) {
                        //  SynthesisParamsHelper::saveAsCasaImage("dbg.img", amplitude(thisPlane));
                        cfSupport = extractSupport(thisPlane);
                        const int support = cfSupport.itsSize;

                        ASKAPCHECK(support*itsOverSample < int(nx) / 2,
                                   "Overflowing convolution function - increase maxSupport or decrease overSample. " <<
                                   "Current support size = " << support << " oversampling factor=" << itsOverSample <<
                                   " image size nx=" << nx)

                        cfSupport.itsSize = limitSupportIfNecessary(support);

                        if (itsSupport == 0) {
                            itsSupport = cfSupport.itsSize;
                            ASKAPLOG_DEBUG_STR(logger, "Number of planes in convolution function = "
                                                   << itsConvFunc.size() << " or " << itsConvFunc.size() / itsOverSample / itsOverSample <<
                                               " before oversampling with factor " << itsOverSample);
                        }

                        if (isOffsetSupportAllowed()) {
                            setConvFuncOffset(zIndex, cfSupport.itsOffsetU, cfSupport.itsOffsetV);
                        }

                        // just for log output
                        const double cell = std::abs(itsUVCellSize(0) * (casa::C::c
                                                     / acc.frequency()[chan]));
                        ASKAPLOG_DEBUG_STR(logger, "CF cache w-plane=" << iw << " feed=" << feed << " field=" << currentField() <<
                                           ": maximum extent = " << support*cell << " (m) sampled at " << cell / itsOverSample << " (m)" <<
                                           " offset (m): " << cfSupport.itsOffsetU*cell << " " << cfSupport.itsOffsetV*cell);
                    }

                    // use either support determined for this particular plane or a generic one,
                    // determined from the first plane (largest support as we have the largest w-term)
                    const int support = isSupportPlaneDependent() ? cfSupport.itsSize : itsSupport;

                    // Since we are decimating, we need to rescale by the
                    // decimation factor
                    const double rescale = double(itsOverSample * itsOverSample);
                    const int cSize = 2 * support + 1;

                    for (int fracu = 0; fracu < itsOverSample; fracu++) {
                        for (int fracv = 0; fracv < itsOverSample; fracv++) {
                            const int plane = fracu + itsOverSample * (fracv + itsOverSample
                                              * zIndex);
                            ASKAPDEBUGASSERT(plane >= 0 && plane < int(itsConvFunc.size()));
                            itsConvFunc[plane].resize(cSize, cSize);
                            itsConvFunc[plane].set(0.0);

                            // Now cut out the inner part of the convolution function and
                            // insert it into the convolution function
                            for (int iy = -support; iy < support; iy++) {
                                for (int ix = -support; ix < support; ix++) {
                                    ASKAPDEBUGASSERT((ix + support >= 0) && (iy + support >= 0));
                                    ASKAPDEBUGASSERT(ix + support < int(itsConvFunc[plane].nrow()));
                                    ASKAPDEBUGASSERT(iy + support < int(itsConvFunc[plane].ncolumn()));
                                    ASKAPDEBUGASSERT((ix + cfSupport.itsOffsetU)*itsOverSample + fracu + int(nx) / 2 >= 0);
                                    ASKAPDEBUGASSERT((iy + cfSupport.itsOffsetV)*itsOverSample + fracv + int(ny) / 2 >= 0);
                                    ASKAPDEBUGASSERT((ix + cfSupport.itsOffsetU)*itsOverSample + fracu + int(nx) / 2 < int(thisPlane.nrow()));
                                    ASKAPDEBUGASSERT((iy + cfSupport.itsOffsetV)*itsOverSample + fracv + int(ny) / 2 < int(thisPlane.ncolumn()));

                                    itsConvFunc[plane](ix + support, iy + support)
                                    = rescale * thisPlane((ix + cfSupport.itsOffsetU) * itsOverSample + fracu + nx / 2,
                                                          (iy + cfSupport.itsOffsetV) * itsOverSample + fracv + ny / 2);
                                } // for ix
                            } // for iy


                            /*
                            // force normalization for all fractional offsets (or planes)
                            const double norm = sum(real(itsConvFunc[plane]));
                            //    ASKAPLOG_INFO_STR(logger, "Sum of convolution function = " << norm<<" for plane "<<plane<<
                            //       " full buffer has sum="<<thisPlaneNorm<<" ratio="<<norm/thisPlaneNorm);
                            ASKAPDEBUGASSERT(norm>0.);
                            if (norm>0.) {
                                //itsConvFunc[plane]*=casa::Complex(norm/thisPlaneNorm);
                            }
                            */

                        } // for fracv
                    } // for fracu

                } // w loop
            } // chan loop

        } // row of the accessor
    }


    if (nDone == itsMaxFeeds*itsMaxFields*nWPlanes()) {
        if (isSupportPlaneDependent()) {
            ASKAPLOG_DEBUG_STR(logger, "Convolution function cache has " << itsConvFunc.size() << " planes");
            ASKAPLOG_DEBUG_STR(logger, "Variable support size is used:");
            const size_t step = casa::max(itsConvFunc.size() / itsOverSample / itsOverSample / 10, 1);

            for (size_t plane = 0; plane < itsConvFunc.size(); plane += step * itsOverSample * itsOverSample) {
                ASKAPLOG_DEBUG_STR(logger, "CF cache plane " << plane << " (" << plane / itsOverSample / itsOverSample <<
                                   " prior to oversampling) shape is " << itsConvFunc[plane].shape());
            }
        } else {
            ASKAPLOG_INFO_STR(logger, "Shape of convolution function = "
                                  << itsConvFunc[0].shape() << " by " << itsConvFunc.size() << " planes");
        }
    }

    ASKAPCHECK(itsSupport > 0, "Support not calculated correctly");
    updateStats(nDone);
}


/// To finalize the transform of the weights, we use the following steps:
/// 1. For each plane of the convolution function, transform to image plane
/// and multiply by conjugate to get abs value squared.
/// 2. Sum all planes weighted by the weight for that convolution function.
void AWProjectVisGridder::finaliseWeights(casa::Array<double>& out)
{
    ASKAPTRACE("AWProjectVisGridder::finaliseWeights");
    ASKAPLOG_DEBUG_STR(logger, "Calculating sum of weights image");
    ASKAPDEBUGASSERT(itsShape.nelements() >= 3);

    const int nx = itsShape(0);
    const int ny = itsShape(1);
    const int nPol = itsShape(2);
    const int nChan = itsShape(3);

    const int nZ = sumOfWeights().nrow();

    /// We must pad the convolution function to full size, reverse transform
    /// square, and sum multiplied by the corresponding weight
    const int cnx = std::min(maxSupport(), nx);
    const int cny = std::min(maxSupport(), ny);
    const int ccenx = cnx / 2;
    const int cceny = cny / 2;

    /*
    // the following code is for grid-correction of the weight
    casa::Vector<double> ccfx(cnx);
    casa::Vector<double> ccfy(cny);
    for (int ix=0; ix<cnx; ++ix) {
         const double nux = std::abs(double(ix)-double(ccenx))/double(ccenx);
         const double val = grdsf(nux);
         ccfx(ix) = val; //casa::abs(val) > 1e-10 ? 1./val : 0.;
    }
    for (int iy=0; iy<cny; ++iy) {
         const double nuy = std::abs(double(iy)-double(cceny))/double(cceny);
         const double val = grdsf(nuy);
         ccfx(iy) = val; //casa::abs(val) > 1e-10 ? 1./val : 0.;
    }
    //
    */

    /// This is the output array before sinc padding
    casa::Array<double> cOut(casa::IPosition(4, cnx, cny, nPol, nChan));
    cOut.set(0.0);

    // for debugging
    double totSumWt = 0.;

    for (int iz = 0; iz < nZ; ++iz) {
        const int plane = cfIndexFromSumOfWeightsRow(iz);

        bool hasData = false;

        for (int chan = 0; chan < nChan; chan++) {
            for (int pol = 0; pol < nPol; pol++) {
                const double wt = sumOfWeights()(iz, pol, chan);
                ASKAPCHECK(!std::isnan(wt), "sumOfWeights returns NaN for row=" << iz <<
                           " pol=" << pol << " chan=" << chan);

                if (wt > 0.0) {
                    hasData = true;
                    totSumWt += wt;
                    //          break;
                }
            }
        }

        if (hasData) {

            // Now fill the inner part of the uv plane with the convolution function
            // and transform to obtain the image. The uv sampling is fixed here
            // so the total field of view is itsOverSample times larger than the
            // original field of view.
            /// Work space
            casa::Matrix<casa::DComplex> thisPlane(cnx, cny);
            thisPlane.set(0.0);

            // use either support determined for this particular plane or a generic one,
            // determined from the first plane (largest support as we have the largest w-term)
            const int support = (int(itsConvFunc[plane].nrow()) - 1) / 2;
            ASKAPDEBUGASSERT(itsConvFunc[plane].nrow() % 2 == 1);
            ASKAPDEBUGASSERT(itsConvFunc[plane].nrow() == itsConvFunc[plane].ncolumn());

            const std::pair<int, int> cfOffset = getConvFuncOffset(iz);

            for (int iy = -support; iy < +support; ++iy) {
                for (int ix = -support; ix < +support; ++ix) {
                    const int xPos = ix + ccenx + cfOffset.first;
                    const int yPos = iy + cceny + cfOffset.second;

                    if ((xPos < 0) || (yPos < 0) || (xPos >= int(thisPlane.nrow())) || (yPos >= int(thisPlane.ncolumn()))) {
                        continue;
                    }

                    thisPlane(xPos, yPos) = itsConvFunc[plane](ix + support, iy + support);
                }
            }

            scimath::fft2d(thisPlane, false);
            thisPlane *= casa::DComplex(cnx * cny);

            // Now we need to cut out only the part inside the field of view
            for (int chan = 0; chan < nChan; chan++) {
                for (int pol = 0; pol < nPol; pol++) {
                    casa::IPosition ip(4, 0, 0, pol, chan);
                    const double wt = sumOfWeights()(iz, pol, chan);
                    ASKAPCHECK(!std::isnan(wt), "sumOfWeights returns NaN for row=" << iz <<
                               " pol=" << pol << " chan=" << chan);

                    for (int ix = 0; ix < cnx; ix++) {
                        ip(0) = ix;

                        for (int iy = 0; iy < cny; iy++) {
                            ip(1) = iy;
                            const casa::DComplex val = thisPlane(ix, iy);
                            cOut(ip) += double(wt) * casa::real(val * conj(val));//*ccfx(ix)*ccfy(iy);
                        }
                    }
                }
            }
        }
    }

    //SphFuncVisGridder::correctConvolution(cOut);
    scimath::PaddingUtils::fftPad(cOut, out, paddingFactor());

    ASKAPLOG_DEBUG_STR(logger,
        "Finished finalising the weights, the sum over all convolution functions is "
        << totSumWt);
}


/// Correct for gridding convolution function
/// @param image image to be corrected
void AWProjectVisGridder::correctConvolution(casa::Array<double>& /*image*/)
{
    //SphFuncVisGridder::correctConvolution(image);

    // experiments with grid-correction, I didn't manage to bring this code
    // into working order so far, so it is commented out (and grid-correction is
    // temporary disabled during the generation of CFs

    /*

    // unlike for plain spheroidal function gridder or w-projection we multuply
    // by FT of the antialiasing filter rather than divide here. The reason is
    // because the weight also computed from acutal CFs with the antialiasing filter
    // applied (weight squared is in denominator). The alternative approach would be
    // to grid-correct the weight as well, but then searching for maximum becomes difficult

    ASKAPDEBUGASSERT(itsShape.nelements()>=2);
    ASKAPDEBUGASSERT(image.shape() == itsShape);
    const casa::Int xHalfSize = itsShape(0)/2;
    const casa::Int yHalfSize = itsShape(1)/2;
    const casa::Int nx = itsShape(0);
    const casa::Int ny = itsShape(1);
    casa::Vector<double> ccfx(itsShape(0));
    casa::Vector<double> ccfy(itsShape(1));
    ASKAPDEBUGASSERT(itsShape(0)>1);
    ASKAPDEBUGASSERT(itsShape(1)>1);

    // note grdsf(-1)=0.
    for (int ix=0; ix<nx; ++ix) {
         const double nux=std::abs(double(ix-xHalfSize))/double(xHalfSize);
         ccfx(ix) = grdsf(nux);
    }
    for (int iy=0; iy<ny; ++iy) {
         const double nuy=std::abs(double(iy-yHalfSize))/double(yHalfSize);
         ccfy(iy) = grdsf(nuy);
    }
    casa::Matrix<casa::DComplex> buffer(nx,ny);
    for (casa::Int ix = 0; ix < nx; ++ix) {
         for (casa::Int iy = 0; iy < ny; ++iy) {
              buffer(ix,iy) = ccfx(ix)*ccfy(iy);
         }
    }
    scimath::fft2d(buffer, true);
    buffer *= casa::DComplex(1./(double(nx)*double(ny)));
    for (casa::Int x = 0; x < nx; ++x) {
         for (casa::Int y = 0; y < ny; ++y) {
              buffer(x,y) *= conj(buffer(x,y));
         }
    }
    scimath::fft2d(buffer, false);
    buffer *= casa::DComplex(double(nx)*double(ny));


    for (casa::ArrayIterator<double> it(image, 2); !it.pastEnd(); it.next()) {
       casa::Matrix<double> mat(it.array());
       ASKAPDEBUGASSERT(int(mat.nrow()) == nx);
       ASKAPDEBUGASSERT(int(mat.ncolumn()) == ny);
       for (int ix=0; ix<nx; ++ix) {
            for (int iy=0; iy<ny; ++iy) {
                 const double val = ccfx(ix)*ccfy(iy);
                 if (casa::abs(val)<1e-20) {
                     mat(ix,iy) = 0.;
                 } else {
                    mat(ix, iy) *= real(buffer(ix,iy))/val;
                 }
            }
       }

       //casa::Array<float> img(mat.shape());
       //casa::convertArray<float, double>(img, mat);
       //SynthesisParamsHelper::saveAsCasaImage("dbg.img",img);
       //throw 1;
    }
    */
}

int AWProjectVisGridder::cIndex(int row, int pol, int chan)
{
    const int plane = itsCMap(row, pol, chan);
    if (nWPlanes() > 0) {
        notifyOfWPlaneUse(plane % nWPlanes());
    }

    return plane;
}

/// @brief assignment operator (not to be called)
/// @details It is made private, so we can't call it inadvertently
/// @param[in] other input object
AWProjectVisGridder& AWProjectVisGridder::operator=(const AWProjectVisGridder &)
{
    ASKAPTHROW(AskapError, "This method is not supposed to be called!");
    return *this;
}


/// @brief static method to create gridder
/// @details Each gridder should have a static factory method, which is
/// able to create a particular type of the gridder and initialise it with
/// the parameters taken form the given parset. It is assumed that the
/// method receives a subset of parameters where the gridder name is already
/// taken out.
/// @param[in] parset input parset file
/// @return a shared pointer to the gridder instance
IVisGridder::ShPtr AWProjectVisGridder::createGridder(const LOFAR::ParameterSet& parset)
{
    boost::shared_ptr<AWProjectVisGridder> gridder = createAProjectGridder<AWProjectVisGridder>(parset);
    gridder->configureGridder(parset);

    return gridder;
}


} // end namespace synthesis
} // end namespace askap
