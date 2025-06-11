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
#include <askap/askap_synthesis.h>

// ASKAPsoft includes
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.awprojectvisgridder");
#include <askap/gridding/AWProjectVisGridder.h>
#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <casacore/casa/Arrays/ArrayIter.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/casa/Quanta/MVDirection.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Quanta/MVTime.h>

// Local package includes
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/gridding/SupportSearcher.h>
#include <askap/profile/AskapProfiler.h>

// for debugging - to export intermediate images
#include <askap/scimath/utils/ImageUtils.h>

namespace askap {
namespace synthesis {

//std::vector<casa::Matrix<casa::Complex> > AWProjectVisGridder::theirCFCache;
//std::vector<std::pair<int,int> > AWProjectVisGridder::theirConvFuncOffsets;

/// @brief a helper method for a ref copy of casa arrays held in
/// stl vector
/// @param[in] in input array
/// @param[out] out output array (will be resized)
/// @return size of the cache in bytes (assuming Complex array elements)
//template<typename T>
//size_t deepRefCopyOfSTDVector(const std::vector<T> &in,
//                            std::vector<T> &out)
//{
//   out.resize(in.size());
//   size_t total = 0;
//   const typename std::vector<T>::const_iterator inEnd = in.end();
//   typename std::vector<T>::iterator outIt = out.begin();
//   for (typename std::vector<T>::const_iterator inIt = in.begin();
//       inIt != inEnd; ++inIt,++outIt) {
//       outIt->reference(*inIt);
//       total += outIt->nelements()*sizeof(casa::Complex)+sizeof(T);
//   }
//   return total;
//}

AWProjectVisGridder::AWProjectVisGridder(const boost::shared_ptr<IBasicIllumination const> &illum,
        const double wmax, const int nwplanes,
        const double cutoff, const int overSample,
        const int maxSupport, const int limitSupport,
        const int maxFeeds, const int maxFields, const double pointingTol,
        const double paTol, const double freqTol, const bool frequencyDependent,
        const bool spheroidalTaper, const double spheroidalWeightsCutoff, const std::string& name) :
        AProjectGridderBase(maxFeeds, maxFields, pointingTol, paTol, freqTol),
        WProjectVisGridder(wmax, nwplanes, cutoff, overSample, maxSupport, limitSupport, name),
        itsReferenceFrequency(0.0), itsIllumination(illum), itsFreqDep(frequencyDependent),
        itsSpheroidalTaper(spheroidalTaper), itsSpheroidalWeightsCutoff(spheroidalWeightsCutoff),
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
        itsSpheroidalTaper(other.itsSpheroidalTaper),
        itsSpheroidalWeightsCutoff(other.itsSpheroidalWeightsCutoff),
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

    const casacore::Vector<casacore::RigidVector<double, 3> > &rotatedUVW = acc.rotatedUVW(getTangentPoint());
    const casacore::Vector<casacore::Double> & chanFreq = acc.frequency();

    for (int i = 0; i < nSamples; ++i) {
        const int feed = acc.feed1()(i);
        ASKAPCHECK(feed < itsMaxFeeds, "Exceeded specified maximum number of feeds");
        ASKAPCHECK(feed > -1, "Illegal negative feed number");

        const double w = (rotatedUVW(i)(2)) / (casacore::C::c);

        for (int chan = 0; chan < nChan; ++chan) {
            const double freq = chanFreq[chan];
            const int iw = getWPlane(w * freq);

            for (int pol = 0; pol < nPol; ++pol) {
                /// Order is (iw, chan, feed)
                itsCMap(i, pol, chan) = iw;
                if (iw >= 0) {
                    if (itsFreqDep) {
                        itsCMap(i, pol, chan) += nWPlanes() * (chan + nChan * (feed + itsMaxFeeds * currentField()));
                        ASKAPCHECK(itsCMap(i, pol, chan) < nWPlanes()*itsMaxFeeds*itsMaxFields*nChan,
                                   "CMap index too large");
                        //ASKAPCHECK(itsCMap(i, pol, chan) > -1, "CMap index less than zero");
                    } else {
                        itsCMap(i, pol, chan) += nWPlanes() * (feed + itsMaxFeeds * currentField());
                        ASKAPCHECK(itsCMap(i, pol, chan) < nWPlanes()*itsMaxFeeds*itsMaxFields,
                                   "CMap index too large");
                        //ASKAPCHECK(itsCMap(i, pol, chan) > -1, "CMap index less than zero");
                    }
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
void AWProjectVisGridder::initialiseGrid(const scimath::Axes& axes,  const casacore::IPosition& shape,
        const bool dopsf, const bool dopcf)
{
    ASKAPTRACE("AWProjectVisGridder::initialiseGrid");
    WProjectVisGridder::initialiseGrid(axes, shape, dopsf, dopcf);

    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casacore::uInt nx = maxSupport();
    const casacore::uInt ny = maxSupport();

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
        const casacore::Array<imtype>& image)
{
    ASKAPTRACE("AWProjectVisGridder::initialiseDegrid");
    WProjectVisGridder::initialiseDegrid(axes, image);
    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casacore::uInt nx = maxSupport();
    const casacore::uInt ny = maxSupport();

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

    // sharecf assumes frequencies and feedpa stay constant throughout the data
    if (itsShareCF && itsSupport > 0 && !isPSFGridder()) return;

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

    // Note the Image and PSF CF differ (no phase slope for PSF), so we need to recalculate for PSF
    if (itsShareCF && theirCFCache.size()>0 && !isPSFGridder()) {
        itsSupport = 3;
        // we already have what we need for the image case, for PCF we can create it here
        if (isPCFGridder()) {
            size_t nplane = theirCFCache.size();
            ASKAPLOG_INFO_STR(logger,"Setting PCF gridder support from shared cache");
            for (size_t plane = 0; plane < nplane; plane++) {
                itsConvFunc[plane].resize(3,3);
                itsConvFunc[plane].set(0.0);
                int support = theirCFCache[plane].shape()(0);
                itsConvFunc[plane](1,1) = casacore::Complex(1.0,support);
            }
        } else {
            // residual or model gridder
            size_t size = deepRefCopyOfSTDVector(theirCFCache,itsConvFunc)/1024/1024;
            ASKAPLOG_INFO_STR(logger, "Using cached convolution functions ("<<size<<" MB)");
            if (isOffsetSupportAllowed()) {
                for (size_t i=0; i<theirConvFuncOffsets.size(); i++) {
                    setConvFuncOffset(i,theirConvFuncOffsets[i].first,theirConvFuncOffsets[i].second);
                }
            }
        }
        return;
    }

    casacore::MVDirection out = getImageCentre();
    const int nSamples = acc.nRow();

    ASKAPDEBUGASSERT(itsIllumination);
    // just to avoid a repeated call to a virtual function from inside the loop
    const bool hasSymmetricIllumination = itsIllumination->isSymmetric();

    // check whether the output pattern is image based. If so, inverse FFT is not needed
    const bool imageBasedPattern = itsIllumination->isImageBased();

    // check whether the pattern depends on the feed number
    const bool feedDependent = itsIllumination->isFeedDependent();

    validateCFCache(acc, hasSymmetricIllumination);

    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const int nx = maxSupport();
    const int ny = maxSupport();

    const casacore::uInt qnx = nx / itsOverSample;
    const casacore::uInt qny = ny / itsOverSample;

    // Find the actual cellsizes in x and y (radians)
    // corresponding to the limited support
    const double ccellx = 1.0 / (double(qnx) * itsUVCellSize(0));
    const double ccelly = 1.0 / (double(qny) * itsUVCellSize(1));

    casacore::Vector<double> ccfx;
    casacore::Vector<double> ccfy;
    if (itsSpheroidalTaper) {
        // Include spheroidal for anti-aliasing and general kernel robustness
        ccfx.resize(qnx);
        ccfy.resize(qny);
        for (casacore::uInt qix = 0; qix < qnx; ++qix) {
            const double nux = std::abs(double(qix) - double(qnx / 2)) / double(qnx / 2);
            ccfx(qix) = grdsf(nux);
        }
        for (casacore::uInt qiy = 0; qiy < qny; ++qiy) {
            const double nuy = std::abs(double(qiy) - double(qny / 2)) / double(qny / 2);
            ccfy(qiy) = grdsf(nuy);
        }
        if (doInterpolation()) {
          // The spheroidal is undefined and set to zero at nu=1, but that
          // is not the numerical limit. Estimate it from its neighbours.
          interpolateEdgeValues(ccfx);
          interpolateEdgeValues(ccfy);
        }
    }

    UVPattern &pattern = uvPattern();
    casacore::Matrix<imtypeComplex> thisPlane = getCFBuffer();
    ASKAPDEBUGASSERT(thisPlane.nrow() == nx);
    ASKAPDEBUGASSERT(thisPlane.ncolumn() == ny);
    const casacore::Vector<casacore::Double> & chanFreq = acc.frequency();

    int nDone = 0;
    // Limit number of fft threads to 8 (more is slower for our fft sizes)
    scimath::FFT2DWrapper<imtypeComplex> fft2d(true,8);
    for (int row = 0; row < nSamples; ++row) {
        const int feed = acc.feed1()(row);

        if (!isCFValid(feed, currentField())) {
            makeCFValid(feed, currentField());
            nDone++;
            casacore::MVDirection offset(acc.pointingDir1()(row).getAngle());
            const double parallacticAngle = hasSymmetricIllumination ? 0. : acc.feed1PA()(row);

            for (int chan = 0; chan < nChan; ++chan) {

                /// Extract illumination pattern for this channel
                itsIllumination->getPattern(chanFreq[chan], pattern, out, offset,
                                            parallacticAngle, isPSFGridder() || isPCFGridder(),
                                            feed);
                // If the pattern is in the Fourier domain, FFT to the image domain
                if( !imageBasedPattern ) {
                    fft2d(pattern.pattern(), false);
                }
                /// Calculate the total convolution function including
                /// the w term and the antenna convolution function

                for (int iw = 0; iw < nWPlanes(); ++iw) {
                    thisPlane.set(0.0);

                    // Loop over the central nx, ny region, setting it to the product
                    // of the phase screen and the a-projection function
                    double maxCF = 0.0;
                    const double w = 2.0f * casacore::C::pi * getWTerm(iw);
                    //std::cout<<"plane "<<iw<<" w="<<w<<std::endl;

                    for (int iy = 0; iy < int(ny); ++iy) {
                        const double y2 = casacore::square((double(iy) - double(ny) / 2) * ccelly);
                        const int qiy = iy + int(qny)/2 - int(ny)/2;

                        for (int ix = 0; ix < int(nx); ++ix) {
                            const double x2 = casacore::square((double(ix) - double(nx) / 2) * ccellx);
                            const int qix = ix + int(qnx)/2 - int(nx)/2;
                            const double r2 = x2 + y2;

                            // The spheroidal is only defined in the qnx,qny region
                            if ((r2 < 1.0) &&
                                (!itsSpheroidalTaper || ((qix>=0) && (qix<qnx) && (qiy>=0) && (qiy<qny)))) {
                                const double phase = w * (1.0 - sqrt(1.0 - r2));
                                const imtype taper = itsSpheroidalTaper ? ccfx(qix) * ccfy(qiy) : 1.0;
                                const imtypeComplex wt = pattern(ix, iy) * conj(pattern(ix, iy)) * taper;
                                // this ensures the oversampling is done
                                thisPlane(ix, iy) = wt * imtypeComplex(cos(phase), -sin(phase));
                                //thisPlane(ix, iy)=wt*imtypeComplex(cos(phase));
                                maxCF += casacore::abs(wt);
                            }
                        }
                    }

                    ASKAPCHECK(maxCF > 0.0, "Convolution function is empty");

                    // At this point, we have the phase screen multiplied by the illumination
                    // function, sampled on larger cellsize (itsOverSample larger) in image
                    // space. Only the inner qnx, qny pixels have a non-zero value

                    // Now we have to calculate the Fourier transform to get the
                    // convolution function in uv space
                    fft2d(thisPlane, true);

                    // Now correct for normalization of FFT
                    thisPlane *= imtypeComplex(1.0 / (double(nx) * double(ny)));
                    // use this norm later on during normalisation
                    const double thisPlaneNorm = sum(real(thisPlane));
                    ASKAPDEBUGASSERT(thisPlaneNorm > 0.);

                    const int zIndex = iw + nWPlanes() * (chan + nChan * (feed + itsMaxFeeds * currentField()));

                    // If the support is not yet set, find it and size the
                    // convolution function appropriately

                    // by default the common support without offset is used
                    CFSupport cfSupport(itsSupport);

                    if (isSupportPlaneDependent() || (itsSupport == 0)) {
                        cfSupport = extractSupport(thisPlane);
                        const int support = cfSupport.itsSize;

                        ASKAPCHECK((support+1)*itsOverSample < int(nx) / 2,
                                   "Overflowing convolution function - increase maxSupport or cutoff or decrease overSample. " <<
                                   "Current support size = " << support << " oversampling factor=" << itsOverSample <<
                                   " image size nx=" << nx)

                        cfSupport.itsSize = limitSupportIfNecessary(support);

                        if (itsSupport == 0) {
                            itsSupport = cfSupport.itsSize;
                            ASKAPLOG_DEBUG_STR(logger, "Number of planes in convolution function = " <<
                                               itsConvFunc.size() << " or " <<
                                               itsConvFunc.size() / itsOverSample / itsOverSample <<
                                               " before oversampling with factor " << itsOverSample);
                        }

                        if (isOffsetSupportAllowed()) {
                            setConvFuncOffset(zIndex, cfSupport.itsOffsetU, cfSupport.itsOffsetV);
                        }

                        // just for log output
                        const double cell = std::abs(itsUVCellSize(0) * (casacore::C::c
                                                     / acc.frequency()[chan]));
                        ASKAPLOG_DEBUG_STR(logger, "CF cache w-plane=" << iw << " feed=" << feed << " field=" <<
                                           currentField() << ": maximum extent = " << support*cell <<
                                           " (m) sampled at " << cell / itsOverSample << " (m)" << " offset (m): " <<
                                           cfSupport.itsOffsetU*cell << " " << cfSupport.itsOffsetV*cell);
                    }

                    // use either support determined for this particular plane or a generic one,
                    // determined from the first plane (largest support as we have the largest w-term)
                    const int support = isSupportPlaneDependent() ? cfSupport.itsSize : itsSupport;

                    // Since we are decimating, we need to rescale by the
                    // decimation factor
                    const double rescale = double(itsOverSample * itsOverSample);
                    const int cSize = 2 * support + 1;

                    // work out range of kx, ky and see if they will overflow the array
                    const int kxmin = (-support + cfSupport.itsOffsetU)*itsOverSample + nx/2;
                    const int kxmax = (support + cfSupport.itsOffsetU)*itsOverSample + itsOverSample-1 + nx/2;
                    const int kymin = (-support + cfSupport.itsOffsetV)*itsOverSample + ny/2;
                    const int kymax = (support + cfSupport.itsOffsetV)*itsOverSample + itsOverSample-1 + ny/2;
                    int overflow = 0;
                    if (kxmin<0) {
                        overflow = -kxmin;
                    }
                    if (kxmax>=nx) {
                        overflow = std::max(overflow, kxmax-(nx-1));
                    }
                    if (kymin<0) {
                        overflow = std::max(overflow, -kymin);
                    }
                    if (kymax>=ny) {
                        overflow = std::max(overflow, kymax-(ny-1));
                    }

                    ASKAPCHECK(overflow==0,"Convolution function overflowing - increase maxsupport or cutoff or decrease oversample, overflow="<<overflow);


                    for (int fracu = 0; fracu < itsOverSample; fracu++) {
                        for (int fracv = 0; fracv < itsOverSample; fracv++) {
                            const int plane = fracu + itsOverSample * (fracv + itsOverSample * zIndex);
                            ASKAPDEBUGASSERT(plane >= 0 && plane < int(itsConvFunc.size()));
                            if (isPCFGridder()) {
                                itsConvFunc[plane].resize(3,3);
                                itsConvFunc[plane].set(0.0);
                                itsConvFunc[plane](1,1) = casacore::Complex(1.0,support);
                                continue;
                            }
                            itsConvFunc[plane].resize(cSize, cSize);
                            itsConvFunc[plane].set(0.0);

                            // Now cut out the inner part of the convolution function and
                            // insert it into the convolution function
                            for (int iy = -support; iy <= support; iy++) {
                                const int ky = (iy + cfSupport.itsOffsetV)*itsOverSample + fracv + int(ny) / 2;
                                for (int ix = -support; ix <= support; ix++) {
                                    const int kx = (ix + cfSupport.itsOffsetU)*itsOverSample + fracu + int(nx) / 2;
                                    itsConvFunc[plane](ix + support, iy + support) = imtype(rescale) * thisPlane(kx, ky);
                                }
                            }

                            /*
                            // force normalization for all fractional offsets (or planes)
                            const double norm = sum(real(itsConvFunc[plane]));
                            //    ASKAPLOG_INFO_STR(logger, "Sum of convolution function = " << norm<<" for plane "<<plane<<
                            //       " full buffer has sum="<<thisPlaneNorm<<" ratio="<<norm/thisPlaneNorm);
                            ASKAPDEBUGASSERT(norm>0.);
                            if (norm>0.) {
                                //itsConvFunc[plane]*=casacore::Complex(norm/thisPlaneNorm);
                                itsConvFunc[plane] /= casacore::Complex(norm);
                            }
                            */

                        } // for fracv
                    } // for fracu
                } // w loop
            } // chan loop

            if (isPSFGridder() && !feedDependent && itsShareCF) {
                // All illumination patterns are identical across feeds for the PSF gridder,
                // so copy CFs across and make valid
                // Use sharecf=false to bypass all CF optimizations
                ASKAPLOG_INFO_STR(logger,"Using feed "<<feed<<" for all PSF gridding convolution functions");
                for (int ifeed = 0; ifeed< itsMaxFeeds; ifeed++) {
                    if (ifeed == feed) continue;
                    makeCFValid(ifeed, currentField());
                    for (int chan = 0; chan < nChan; chan++) {
                        for (int iw = 0; iw < nWPlanes(); ++iw) {
                            const int zIndex = iw + nWPlanes() * (chan + nChan * (ifeed +
                                itsMaxFeeds * currentField()));
                            const int refzIndex = iw + nWPlanes() * (chan + nChan * (feed +
                                itsMaxFeeds * currentField()));
                            for (int fracu = 0; fracu < itsOverSample; fracu++) {
                                for (int fracv = 0; fracv < itsOverSample; fracv++) {
                                    const int plane = fracu + itsOverSample * (fracv + itsOverSample * zIndex);
                                    const int refPlane = fracu + itsOverSample * (fracv + itsOverSample * refzIndex);
                                    itsConvFunc[plane].reference(itsConvFunc[refPlane]);
                                }
                            }
                        }
                    }
                }
            }
        } // isCFValid
    } // row of the accessor

    if (nDone == itsMaxFeeds*itsMaxFields) {
        if (isSupportPlaneDependent()) {
            ASKAPLOG_INFO_STR(logger, "Convolution function cache has " << itsConvFunc.size() << " planes");
            ASKAPLOG_INFO_STR(logger, "Maximum kernel size is " << itsConvFunc[0].shape());

            ASKAPLOG_DEBUG_STR(logger, "Variable support size is used:");
            const size_t step = casacore::max(itsConvFunc.size() / itsOverSample / itsOverSample / 10, 1);
            for (size_t plane = 0; plane < itsConvFunc.size(); plane += step * itsOverSample * itsOverSample) {
                ASKAPLOG_DEBUG_STR(logger, "CF cache plane " << plane << " (" << plane / itsOverSample / itsOverSample <<
                                   " prior to oversampling) shape is " << itsConvFunc[plane].shape());
            }
            size_t total = 0;
            for (auto matrix : itsConvFunc) {
                total += matrix.shape().product() * sizeof(imtypeComplex);
            }
            ASKAPLOG_INFO_STR(logger,"Convolution function cache size = "<< total/1024/1024 <<" MB");
        } else {
            ASKAPLOG_INFO_STR(logger, "Shape of convolution function = "
                                  << itsConvFunc[0].shape() << " by " << itsConvFunc.size() << " planes");
        }
    }

    ASKAPCHECK(itsSupport > 0, "Support not calculated correctly");
    updateStats(nDone);

    // Save the CF to the cache
    if (itsShareCF && !(isPSFGridder() || isPCFGridder())) {
        deepRefCopyOfSTDVector(itsConvFunc,theirCFCache);
        if (isOffsetSupportAllowed()) {
            theirConvFuncOffsets.resize(nWPlanes()*itsMaxFeeds*itsMaxFields*nChan);
            for (int nw=0; nw<theirConvFuncOffsets.size(); nw++) {
                theirConvFuncOffsets[nw]=getConvFuncOffset(nw);
            }
        }
    }

}


/// To finalize the transform of the weights, we use the following steps:
/// 1. For each plane of the convolution function, transform to image plane
/// and multiply by conjugate to get abs value squared.
/// 2. Sum all planes weighted by the weight for that convolution function.
void AWProjectVisGridder::finaliseWeights(casacore::Array<imtype>& out)
{
    ASKAPTRACE("AWProjectVisGridder::finaliseWeights");
    ASKAPLOG_DEBUG_STR(logger, "Calculating sum of weights image");
    ASKAPLOG_DEBUG_STR(logger, "Shape " << itsShape);
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

    casacore::Vector<double> ccfx;
    casacore::Vector<double> ccfy;
    if (itsSpheroidalTaper) {
        // the spheroidal ccf needs to be removed from the weights before they are stored
        // it will also be removed from the image in correctConvolution
        ccfx.resize(cnx);
        ccfy.resize(cny);
        for (int ix=0; ix<cnx; ++ix) {
             const double nux = std::abs(double(ix)-double(ccenx))/double(ccenx);
             const double val = grdsf(nux);
             //ccfx(ix) = val; //casacore::abs(val) > 1e-10 ? 1./val : 0.;
             ccfx(ix) = val > itsSpheroidalWeightsCutoff ? val : 0.;
        }
        for (int iy=0; iy<cny; ++iy) {
             const double nuy = std::abs(double(iy)-double(cceny))/double(cceny);
             const double val = grdsf(nuy);
             //ccfy(iy) = val; //casacore::abs(val) > 1e-10 ? 1./val : 0.;
             ccfy(iy) = val > itsSpheroidalWeightsCutoff ? val : 0.;
        }
        // this isn't really needed unless itsSpheroidalWeightsCutoff==0. But shouldn't hurt
        if (doInterpolation()) {
          // The spheroidal is undefined and set to zero at nu=1, but that
          // is not the numerical limit. Estimate it from its neighbours.
          interpolateEdgeValues(ccfx);
          interpolateEdgeValues(ccfy);
        }
    }

    /// This is the output array before sinc padding
    casacore::Array<imtype> cOut(casacore::IPosition(4, cnx, cny, nPol, nChan));
    cOut.set(0.0);
    // Limit number of fft threads to 8 (more is slower for our fft sizes)
    scimath::FFT2DWrapper<imtypeComplex> fft2d(true,8);

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
            casacore::Matrix<imtypeComplex> thisPlane(cnx, cny);
            thisPlane.set(0.0);

            // use either support determined for this particular plane or a generic one,
            // determined from the first plane (largest support as we have the largest w-term)
            const int support = (int(itsConvFunc[plane].nrow()) - 1) / 2;
            ASKAPDEBUGASSERT(itsConvFunc[plane].nrow() % 2 == 1);
            ASKAPDEBUGASSERT(itsConvFunc[plane].nrow() == itsConvFunc[plane].ncolumn());

            const std::pair<int, int> cfOffset = getConvFuncOffset(iz);

            for (int iy = -support; iy <= support; ++iy) {
                for (int ix = -support; ix <= support; ++ix) {
                    const int xPos = ix + ccenx + cfOffset.first;
                    const int yPos = iy + cceny + cfOffset.second;

                    if ((xPos < 0) || (yPos < 0) || (xPos >= int(thisPlane.nrow())) || (yPos >= int(thisPlane.ncolumn()))) {
                        continue;
                    }

                    thisPlane(xPos, yPos) = itsConvFunc[plane](ix + support, iy + support);
                }
            }
            fft2d(thisPlane, false);
            thisPlane *= imtypeComplex(cnx * cny);

            // Now we need to cut out only the part inside the field of view
            for (int chan = 0; chan < nChan; chan++) {
                for (int pol = 0; pol < nPol; pol++) {
                    casacore::IPosition ip(4, 0, 0, pol, chan);
                    const double wt = sumOfWeights()(iz, pol, chan);
                    ASKAPCHECK(!std::isnan(wt), "sumOfWeights returns NaN for row=" << iz <<
                               " pol=" << pol << " chan=" << chan);

                    if (itsSpheroidalTaper) {
                        // if itsConvFunc has an additional taper, remove it before updating the weights.
                        for (int ix = 0; ix < cnx; ix++) {
                            ip(0) = ix;
                            for (int iy = 0; iy < cny; iy++) {
                                ip(1) = iy;
                                imtype taper = ccfx(ix) * ccfy(iy);
                                const imtypeComplex val = taper > 0 ? thisPlane(ix, iy) / taper : 0.;
                                cOut(ip) += double(wt) * casacore::real(val * conj(val));
                            }
                        }
                    }
                    else {
                        // if itsConvFunc has an additional taper, remove it before updating the weights.
                        for (int ix = 0; ix < cnx; ix++) {
                            ip(0) = ix;
                            for (int iy = 0; iy < cny; iy++) {
                                ip(1) = iy;
                                const imtypeComplex val = thisPlane(ix, iy);
                                cOut(ip) += double(wt) * casacore::real(val * conj(val));
                            }
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
void AWProjectVisGridder::correctConvolution(casacore::Array<imtype>& image)
{

    // if there are any kernel factors that are not A or W terms, remove them now
    if (itsSpheroidalTaper) {
        SphFuncVisGridder::correctConvolution(image);
    }

    // experiments with grid-correction, I didn't manage to bring this code
    // into working order so far, so it is commented out (and grid-correction is
    // temporary disabled during the generation of CFs

    //
    // // unlike for plain spheroidal function gridder or w-projection we multuply
    // // by FT of the antialiasing filter rather than divide here. The reason is
    // // because the weight also computed from acutal CFs with the antialiasing filter
    // // applied (weight squared is in denominator). The alternative approach would be
    // // to grid-correct the weight as well, but then searching for maximum becomes difficult
    //
    // ASKAPDEBUGASSERT(itsShape.nelements()>=2);
    // ASKAPDEBUGASSERT(image.shape() == itsShape);
    // const casacore::Int xHalfSize = itsShape(0)/2;
    // const casacore::Int yHalfSize = itsShape(1)/2;
    // const casacore::Int nx = itsShape(0);
    // const casacore::Int ny = itsShape(1);
    // casacore::Vector<double> ccfx(itsShape(0));
    // casacore::Vector<double> ccfy(itsShape(1));
    // ASKAPDEBUGASSERT(itsShape(0)>1);
    // ASKAPDEBUGASSERT(itsShape(1)>1);
    //
    // // note grdsf(-1)=0.
    // for (int ix=0; ix<nx; ++ix) {
    //      const double nux=std::abs(double(ix-xHalfSize))/double(xHalfSize);
    //      ccfx(ix) = grdsf(nux);
    // }
    // for (int iy=0; iy<ny; ++iy) {
    //      const double nuy=std::abs(double(iy-yHalfSize))/double(yHalfSize);
    //      ccfy(iy) = grdsf(nuy);
    // }
    // casacore::Matrix<imtypeComplex> buffer(nx,ny);
    // for (casacore::Int ix = 0; ix < nx; ++ix) {
    //      for (casacore::Int iy = 0; iy < ny; ++iy) {
    //           buffer(ix,iy) = ccfx(ix)*ccfy(iy);
    //      }
    // }
    // //scimath::fft2d(buffer, true);
    // buffer *= imtypeComplex(1./(double(nx)*double(ny)));
    // for (casacore::Int x = 0; x < nx; ++x) {
    //      for (casacore::Int y = 0; y < ny; ++y) {
    //           buffer(x,y) *= conj(buffer(x,y));
    //      }
    // }
    // //scimath::fft2d(buffer, false);
    // buffer *= imtypeComplex(double(nx)*double(ny));
    //
    //
    // for (casacore::ArrayIterator<double> it(image, 2); !it.pastEnd(); it.next()) {
    //    casacore::Matrix<double> mat(it.array());
    //    ASKAPDEBUGASSERT(int(mat.nrow()) == nx);
    //    ASKAPDEBUGASSERT(int(mat.ncolumn()) == ny);
    //    for (int ix=0; ix<nx; ++ix) {
    //         for (int iy=0; iy<ny; ++iy) {
    //              const double val = ccfx(ix)*ccfy(iy);
    //              if (casacore::abs(val)<1e-20) {
    //                  mat(ix,iy) = 0.;
    //              } else {
    //                 mat(ix, iy) *= real(buffer(ix,iy))/val;
    //              }
    //         }
    //    }
    //
    //    //casacore::Array<float> img(mat.shape());
    //    //casacore::convertArray<float, double>(img, mat);
    //    //SynthesisParamsHelper::saveAsCasaImage("dbg.img",img);
    //    //throw 1;
    // }
    //
}

int AWProjectVisGridder::cIndex(int row, int pol, int chan)
{
    const int plane = itsCMap(row, pol, chan);
    if (nWPlanes() > 0 && plane >= 0) {
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
/// the parameters taken from the given parset. It is assumed that the
/// method receives a subset of parameters where the gridder name is already
/// taken out.
/// @param[in] parset input parset file
/// @return a shared pointer to the gridder instance
IVisGridder::ShPtr AWProjectVisGridder::createGridder(const LOFAR::ParameterSet& parset)
{
    boost::shared_ptr<AWProjectVisGridder> gridder = createAProjectGridder<AWProjectVisGridder>(parset);
    gridder->configureGridder(parset); // this calls WProjectVisGridder
    return gridder;
}


} // end namespace synthesis
} // end namespace askap
