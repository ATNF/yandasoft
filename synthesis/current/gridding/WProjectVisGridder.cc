/// @file WProjectVisGridder.cc
///
/// @copyright (c) 2007,2016 CSIRO
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

// System includes
#include <cmath>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <fft/FFTWrapper.h>
#include <profile/AskapProfiler.h>

// Local package includes
#include <gridding/WProjectVisGridder.h>
#include <gridding/SupportSearcher.h>

ASKAP_LOGGER(logger, ".gridding.wprojectvisgridder");

namespace askap {
namespace synthesis {

std::vector<casa::Matrix<casa::Complex> > WProjectVisGridder::theirCFCache;
std::vector<std::pair<int,int> > WProjectVisGridder::theirConvFuncOffsets;

/// @brief a helper method for a ref copy of casa arrays held in
/// stl vector
/// @param[in] in input array
/// @param[out] out output array (will be resized)
template<typename T>
void deepRefCopyOfSTDVector(const std::vector<T> &in,
                            std::vector<T> &out)
{
   out.resize(in.size());

   const typename std::vector<T>::const_iterator inEnd = in.end();
   typename std::vector<T>::iterator outIt = out.begin();
   for (typename std::vector<T>::const_iterator inIt = in.begin();
       inIt != inEnd; ++inIt,++outIt) {
       outIt->reference(*inIt);
   }
}

WProjectVisGridder::WProjectVisGridder(const double wmax,
                                       const int nwplanes,
                                       const double cutoff,
                                       const int overSample,
                                       const int maxSupport,
                                       const int limitSupport,
                                       const std::string& name,
                                       const float alpha,
                                       const bool useDouble,
                                       const bool shareCF) :
        WDependentGridderBase(wmax, nwplanes, alpha),
        itsMaxSupport(maxSupport), itsCutoff(cutoff), itsLimitSupport(limitSupport),
        itsPlaneDependentCFSupport(false), itsOffsetSupportAllowed(false), itsCutoffAbs(false),
        itsDoubleCF(useDouble), itsShareCF(shareCF)
{
    ASKAPCHECK(overSample > 0, "Oversampling must be greater than 0");
    ASKAPCHECK(maxSupport > 0, "Maximum support must be greater than 0")
    itsSupport = 0;
    itsOverSample = overSample;
    setTableName(name);
    itsConvFunc.resize(nWPlanes()*itsOverSample*itsOverSample);
}

WProjectVisGridder::~WProjectVisGridder()
{
}

/// @brief copy constructor
/// @details It is required to decouple internal arrays in the input
/// object and the copy.
/// @param[in] other input object
WProjectVisGridder::WProjectVisGridder(const WProjectVisGridder &other) :
        IVisGridder(other), WDependentGridderBase(other),
        itsCMap(other.itsCMap.copy()), itsMaxSupport(other.itsMaxSupport),
        itsCutoff(other.itsCutoff), itsLimitSupport(other.itsLimitSupport),
        itsPlaneDependentCFSupport(other.itsPlaneDependentCFSupport),
        itsOffsetSupportAllowed(other.itsOffsetSupportAllowed),
        itsCutoffAbs(other.itsCutoffAbs),itsDoubleCF(other.itsDoubleCF),
        itsShareCF(other.itsShareCF) {}


/// Clone a copy of this Gridder
IVisGridder::ShPtr WProjectVisGridder::clone()
{
    return IVisGridder::ShPtr(new WProjectVisGridder(*this));
}

/// @brief initialise sum of weights
/// @details We keep track the number of times each convolution function is used per
/// channel and polarisation (sum of weights). This method is made virtual to be able
/// to do gridder specific initialisation without overriding initialiseGrid.
/// This method accepts no parameters as itsShape, itsNWPlanes, etc should have already
/// been initialised by the time this method is called.
void WProjectVisGridder::initialiseSumOfWeights()
{
    resizeSumOfWeights(nWPlanes());
    zeroSumOfWeights();
}

/// Initialize the convolution function into the cube. If necessary this
/// could be optimized by using symmetries.
void WProjectVisGridder::initIndices(const accessors::IConstDataAccessor& acc)
{
    ASKAPTRACE("WProjectVisGridder::initIndices");
    /// We have to calculate the lookup function converting from
    /// row and channel to plane of the w-dependent convolution
    /// function
    const int nSamples = acc.nRow();
    const int nChan = acc.nChannel();
    const int nPol = acc.nPol();

    itsCMap.resize(nSamples, nPol, nChan);

#ifdef ASKAP_DEBUG
    // in the debug mode we check that all used indices are initialised.
    // negative value means an uninitialised index. In the production version we don't care
    // about uninitialised indices as long as they are not used.
    itsCMap.set(-1);
#endif

    const casa::Vector<casa::RigidVector<double, 3> > &rotatedUVW = acc.rotatedUVW(getTangentPoint());
    const casa::Vector<casa::Double> & chanFreq = acc.frequency();

    for (int i = 0; i < nSamples; ++i) {
        const double w = (rotatedUVW(i)(2)) / (casa::C::c);
        for (int chan = 0; chan < nChan; ++chan) {
            /// Calculate the index into the convolution functions
            const double freq = chanFreq[chan];
            const int wPlane = getWPlane(w * freq);
            for (int pol = 0; pol < nPol; ++pol) {
                itsCMap(i, pol, chan) = wPlane;
            }
        }
    }
}

/// Initialize the convolution function into the cube. If necessary this
/// could be optimized by using symmetries.
void WProjectVisGridder::initConvolutionFunction(const accessors::IConstDataAccessor&)
{
    ASKAPTRACE("WProjectVisGridder::initConvolutionFunction");
    /// We have to calculate the lookup function converting from
    /// row and channel to plane of the w-dependent convolution
    /// function

    if (itsSupport > 0) {
        return;
    }

    itsSupport = 0;

    if (isOffsetSupportAllowed()) {
        // this command is executed only once when itsSupport is not set.
        initConvFuncOffsets(nWPlanes());
    }

    if (isPCFGridder()) {

      // A simple grid kernel for use in setting up the preconditioner function.
      // Set up as a nearest neighbour gridder (based partly on the Box gridder).
      // Gridding weight is written to the real part, gridding support is written
      // to the imaginary part.
      // It is assumed that any non-zero cfSupport.itsOffsetU or
      // cfSupport.itsOffsetV will be taken care of in gridding/preconditioning
      // with an appropriate use of the wKernelPix data.

      itsSupport=1;
      const int cSize=2*itsSupport+1;
      const int cCenter=(cSize-1)/2;

      // Use the w & oversampling setup as the main grid kernels.
      for (int iw = 0; iw < nWPlanes(); ++iw) {

        // Store the main kernel size (in pixels) in the imag part.
        // What size should be stored? From some simple fits:
        //   1e-2 cutoff: support ~ sqrt( 7^2 + (w.theta^2)^2 )
        //   1e-3 cutoff: support ~ 6.02 + 1.14*w.theta^2
        // Preconditioning happens after "deconvolving" the anti-aliasing part
        // of the kernel. When w=0 this should reduce the support from the
        // nominal value of 7, however for larger w kernels it won't (in fact
        // it may result in the gridded data ringing out across the uv plane).
        // So do we want to limit the support when w.Theta ~ 0 but keep +/-3 for
        // larger kernels? I don't know that this is right, but it's a start.
        float wThetaPix = fabs(getWTerm(iw)) / (itsUVCellSize(0) * itsUVCellSize(0));
        float wKernelPix;
        if (wThetaPix < 1) {
          wKernelPix = 3;
        } else if (itsCutoff < 0.01) {
          wKernelPix = 6 + 1.14*wThetaPix;
        } else {
          wKernelPix = sqrt( 49 + wThetaPix*wThetaPix );
        }

        for (int fracu = 0; fracu < itsOverSample; ++fracu) {
            for (int fracv = 0; fracv < itsOverSample; ++fracv) {
                const int plane = fracu + itsOverSample * (fracv + itsOverSample * iw);
                ASKAPDEBUGASSERT(plane < int(itsConvFunc.size()));
                itsConvFunc[plane].resize(cSize, cSize);
                itsConvFunc[plane].set(0.0);
                // are fracu and fracv being correctly used here?
                // I think they should be -ve, since the offset in nux & nuy is +ve.
                const int ix = -float(fracu)/float(itsOverSample);
                const int iy = -float(fracv)/float(itsOverSample);
                itsConvFunc[plane](ix + cCenter, iy + cCenter) =  casa::Complex(1.0, wKernelPix);
            }
        }

      }

      return;
    }

    if (itsShareCF && theirCFCache.size()>0) {
        // we already have what we need
        itsSupport = 1;
        ASKAPLOG_INFO_STR(logger, "Using cached convolution functions");
        deepRefCopyOfSTDVector(theirCFCache,itsConvFunc);
        if (isOffsetSupportAllowed()) {
            for (size_t i=0; i<theirConvFuncOffsets.size(); i++) {
                setConvFuncOffset(i,theirConvFuncOffsets[i].first,theirConvFuncOffsets[i].second);
            }
        }
        return;
    }

    /// These are the actual cell sizes used
    const double cellx = 1.0 / (double(itsShape(0)) * itsUVCellSize(0));
    const double celly = 1.0 / (double(itsShape(1)) * itsUVCellSize(1));

    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    //      int nx=std::min(maxSupport(), itsShape(0));
    //      int ny=std::min(maxSupport(), itsShape(1));
    const int nx = maxSupport();
    const int ny = maxSupport();

    // initialise the buffer for full-sized CF
    ASKAPDEBUGASSERT((nx > 0) && (ny > 0));
    initCFBuffer(casa::uInt(nx), casa::uInt(ny));

    /// We want nx * ccellx = overSample * itsShape(0) * cellx

    const int qnx = nx / itsOverSample;
    const int qny = ny / itsOverSample;
    ASKAPDEBUGASSERT((qnx != 0) && (qny != 0));

    // Find the actual cellsizes in x and y (radians) after over
    // oversampling (in uv space)
    const double ccellx = double(itsShape(0)) * cellx / double(qnx);
    const double ccelly = double(itsShape(1)) * celly / double(qny);

    casa::Vector<float> ccfx(qnx);
    casa::Vector<float> ccfy(qny);

    for (int ix = 0; ix < qnx; ix++) {
        float nux = std::abs(float(ix - qnx / 2)) / float(qnx / 2);
        ccfx(ix) = grdsf(nux) / float(qnx);
    }

    for (int iy = 0; iy < qny; iy++) {
        float nuy = std::abs(float(iy - qny / 2)) / float(qny / 2);
        ccfy(iy) = grdsf(nuy) / float(qny);
    }

    if (itsInterp) {
      // The spheroidal is undefined and set to zero at nu=1, but that
      // is not the numerical limit. Estimate it from its neighbours.
      interpolateEdgeValues(ccfx);
      interpolateEdgeValues(ccfy);
    }

    // Now we step through the w planes, starting the furthest
    // out. We calculate the support for that plane and use it
    // for all the others.

    // We pad here to do sinc interpolation of the convolution
    // function in uv space
    casa::Matrix<casa::DComplex> thisPlane;
    if (itsDoubleCF) thisPlane.reference(getCFBuffer());
    casa::Matrix<casa::Complex> thisPlaneF;
    if (!itsDoubleCF) thisPlaneF.reference(getCFBufferF());
    ASKAPDEBUGASSERT(thisPlane.nrow() == casa::uInt(nx)||thisPlaneF.nrow() == casa::uInt(nx));
    ASKAPDEBUGASSERT(thisPlane.ncolumn() == casa::uInt(ny)||thisPlaneF.ncolumn() == casa::uInt(ny));

    for (int iw = 0; iw < nWPlanes(); ++iw) {
        if (itsDoubleCF) thisPlane.set(0.0);
        else thisPlaneF.set(0.0);

        //const double w = isPSFGridder() ? 0. : 2.0f*casa::C::pi*getWTerm(iw);
        const double w = 2.0f * casa::C::pi * getWTerm(iw);

        // Loop over the central nx, ny region, setting it to the product
        // of the phase screen and the spheroidal function
        for (int iy = 0; iy < qny; iy++) {
            double y2 = double(iy - qny / 2) * ccelly;
            y2 *= y2;

            for (int ix = 0; ix < qnx; ix++) {
                double x2 = double(ix - qnx / 2) * ccellx;
                x2 *= x2;
                const double r2 = x2 + y2;

                if (r2 < 1.0) {
                    const double phase = w * (1.0 - sqrt(1.0 - r2));
                    const float wt = ccfx(ix) * ccfy(iy);
                    ASKAPDEBUGASSERT(ix - qnx / 2 + nx / 2 < nx);
                    ASKAPDEBUGASSERT(iy - qny / 2 + ny / 2 < ny);
                    ASKAPDEBUGASSERT(ix + nx / 2 >= qnx / 2);
                    ASKAPDEBUGASSERT(iy + ny / 2 >= qny / 2);
                    if (itsDoubleCF) {
                        thisPlane(ix - qnx / 2 + nx / 2, iy - qny / 2 + ny / 2) =
                        casa::DComplex(wt * cos(phase), -wt * sin(phase));
                    } else {
                        thisPlaneF(ix - qnx / 2 + nx / 2, iy - qny / 2 + ny / 2) =
                        casa::Complex(wt * cos(phase), -wt * sin(phase));

                    }
                    //thisPlane(ix-qnx/2+nx/2, iy-qny/2+ny/2)=casa::DComplex(wt*cos(phase));
                }
            }
        }

        // At this point, we have the phase screen multiplied by the spheroidal
        // function, sampled on larger cellsize (itsOverSample larger) in image
        // space. Only the inner qnx, qny pixels have a non-zero value

        // Now we have to calculate the Fourier transform to get the
        // convolution function in uv space
        if (itsDoubleCF) scimath::fft2d(thisPlane, true);
        else scimath::fft2d(thisPlaneF, true);

        /*
            for (uint xx=0;xx<thisPlane.nrow();++xx) {
             for (uint yy=0;yy<thisPlane.ncolumn();++yy) {
                  ASKAPCHECK(!std::isinf(casa::abs(thisPlane(xx,yy))),
                      "Infinite value detected for plane="<<iw<<
                      " at "<<xx<<","<<yy<<" "<<thisPlane(xx,yy));
                 }
        }
        */

        // Now thisPlane is filled with convolution function
        // sampled on a finer grid in u,v
        //
        // If the support is not yet set, find it and size the
        // convolution function appropriately

        // by default the common support without offset is used
        CFSupport cfSupport(itsSupport);

        if (isSupportPlaneDependent() || (itsSupport == 0)) {
            cfSupport = (itsDoubleCF ? extractSupport(thisPlane) : extractSupport(thisPlaneF));
            const int support = cfSupport.itsSize;

            ASKAPCHECK(support*itsOverSample < nx / 2,
                       "Overflowing convolution function for w-plane " << iw <<
                       " - increase maxSupport or decrease overSample; support=" <<
                       support << " oversample=" << itsOverSample << " nx=" << nx);
            cfSupport.itsSize = limitSupportIfNecessary(support);

            if (itsSupport == 0) {
                itsSupport = cfSupport.itsSize;
            }

            if (isOffsetSupportAllowed()) {
                setConvFuncOffset(iw, cfSupport.itsOffsetU, cfSupport.itsOffsetV);
            }

        }

        ASKAPCHECK(itsConvFunc.size() > 0, "Convolution function not sized correctly");
        // use either support determined for this particular plane or a generic one,
        // determined from the first plane (largest support as we have the largest w-term)
        const int support = isSupportPlaneDependent() ? cfSupport.itsSize : itsSupport;

        const int cSize = 2 * support + 1;

        for (int fracu = 0; fracu < itsOverSample; ++fracu) {
            for (int fracv = 0; fracv < itsOverSample; ++fracv) {
                const int plane = fracu + itsOverSample * (fracv + itsOverSample * iw);
                ASKAPDEBUGASSERT(plane < int(itsConvFunc.size()));
                itsConvFunc[plane].resize(cSize, cSize);
                itsConvFunc[plane].set(0.0);

                // Now cut out the inner part of the convolution function and
                // insert it into the convolution function
                for (int iy = -support; iy <= support; ++iy) {
                    for (int ix = -support; ix <= support; ++ix) {
                        const int kx = (ix + cfSupport.itsOffsetU)*itsOverSample + fracu + nx / 2;
                        const int ky = (iy + cfSupport.itsOffsetV)*itsOverSample + fracv + ny / 2;
                        ASKAPDEBUGASSERT((ix + support >= 0) && (iy + support >= 0));
                        ASKAPDEBUGASSERT(ix + support < int(itsConvFunc[plane].nrow()));
                        ASKAPDEBUGASSERT(iy + support < int(itsConvFunc[plane].ncolumn()));
                        ASKAPDEBUGASSERT(kx >= 0);
                        ASKAPDEBUGASSERT(ky >= 0);
                        ASKAPDEBUGASSERT(kx < nx);
                        ASKAPDEBUGASSERT(ky < ny);
                        //if (w < 0) {
                        //    itsConvFunc[plane](ix + support, iy + support) =
                        //        conj(thisPlane((ix + cfSupport.itsOffsetU) * itsOverSample + fracu + nx / 2,
                        //              (iy + cfSupport.itsOffsetV) * itsOverSample + fracv + ny / 2));
                        //} else {
                        itsConvFunc[plane](ix + support, iy + support) =
                            (itsDoubleCF ? thisPlane(kx, ky): thisPlaneF(kx, ky));
                        //}
                    }
                }

            } // for fracv
        } // for fracu

    } // for iw

    // force normalization for all fractional offsets (or planes)
    for (size_t plane = 0; plane < itsConvFunc.size(); ++plane) {
        if (itsConvFunc[plane].nelements() == 0) {
            // this plane of the cache is unused
            continue;
        }

        const double norm = sum(casa::real(itsConvFunc[plane]));
        // ASKAPLOG_INFO_STR(logger, "Sum of convolution function = " << norm);
        ASKAPDEBUGASSERT(norm > 0.);

        if (norm > 0.) {
            const casa::Complex invNorm = casa::Complex(1.0/norm);
            itsConvFunc[plane] *= invNorm;
        }
    } // for plane

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

    ASKAPCHECK(itsSupport > 0, "Support not calculated correctly");
    // we can free up the memory because for WProject gridder this method is called only once!
    if (itsDoubleCF) itsCFBuffer.reset();
    else itsCFBufferF.reset();

    // Save the CF to the cache
    if (itsShareCF) {
        deepRefCopyOfSTDVector(itsConvFunc,theirCFCache);
        if (isOffsetSupportAllowed()) {
            theirConvFuncOffsets.resize(nWPlanes());
            for (int nw=0; nw<nWPlanes(); nw++) {
                theirConvFuncOffsets[nw]=getConvFuncOffset(nw);
            }
        }
    }
}

/// @brief search for support parameters
/// @details This method encapsulates support search operation, taking into account the
/// cutoff parameter and whether or not an offset is allowed.
/// @param[in] cfPlane const reference to 2D plane with the convolution function
/// @return an instance of CFSupport with support parameters
WProjectVisGridder::CFSupport WProjectVisGridder::extractSupport(const casa::Matrix<casa::DComplex> &cfPlane) const
{
    ASKAPDEBUGTRACE("WProjectVisGridder::extractSupport");
    CFSupport result(-1);
    SupportSearcher ss(itsCutoff);

    if (isCutoffAbsolute()) {
        ss.search(cfPlane, 1.);
    } else {
        ss.search(cfPlane);
    }

    if (isOffsetSupportAllowed()) {
        result.itsSize = ss.support();
        const casa::IPosition peakPos = ss.peakPos();
        ASKAPDEBUGASSERT(peakPos.nelements() == 2);
        result.itsOffsetU = (peakPos[0] - int(cfPlane.nrow()) / 2) / itsOverSample;
        result.itsOffsetV = (peakPos[1] - int(cfPlane.ncolumn()) / 2) / itsOverSample;
    } else {
        result.itsSize = ss.symmetricalSupport(cfPlane.shape());
        ASKAPCHECK(result.itsSize > 0, "Unable to determine support of convolution function");
    }

    result.itsSize /= 2 * itsOverSample;
    if (result.itsSize < 3) {
        result.itsSize = 3;
    }

    return result;
}

/// @brief search for support parameters
/// @details This method encapsulates support search operation, taking into account the
/// cutoff parameter and whether or not an offset is allowed.
/// @param[in] cfPlane const reference to 2D plane with the convolution function
/// @return an instance of CFSupport with support parameters
WProjectVisGridder::CFSupport WProjectVisGridder::extractSupport(const casa::Matrix<casa::Complex> &cfPlane) const
{
    ASKAPDEBUGTRACE("WProjectVisGridder::extractSupport");
    CFSupport result(-1);
    SupportSearcher ss(itsCutoff);

    if (isCutoffAbsolute()) {
        ss.search(cfPlane, 1.);
    } else {
        ss.search(cfPlane);
    }

    if (isOffsetSupportAllowed()) {
        result.itsSize = ss.support();
        const casa::IPosition peakPos = ss.peakPos();
        ASKAPDEBUGASSERT(peakPos.nelements() == 2);
        result.itsOffsetU = (peakPos[0] - int(cfPlane.nrow()) / 2) / itsOverSample;
        result.itsOffsetV = (peakPos[1] - int(cfPlane.ncolumn()) / 2) / itsOverSample;
    } else {
        result.itsSize = ss.symmetricalSupport(cfPlane.shape());
        ASKAPCHECK(result.itsSize > 0, "Unable to determine support of convolution function");
    }

    result.itsSize /= 2 * itsOverSample;
    if (result.itsSize < 3) {
        result.itsSize = 3;
    }

    return result;
}

/// @brief truncate support, if necessary
/// @details This method encapsulates all usage of itsLimitSupport. It truncates the support
/// if necessary and reports the new value back.
/// @param[in] support support size to truncate according to itsLimitSupport
/// @return support size to use (after possible truncation)
int WProjectVisGridder::limitSupportIfNecessary(int support) const
{
    if (itsLimitSupport > 0  &&  support > itsLimitSupport) {
        ASKAPLOG_INFO_STR(logger, "Convolution function support = "
                              << support << " pixels exceeds upper support limit; "
                              << "set to limit = " << itsLimitSupport << " pixels");
        support = itsLimitSupport;
    }

    const int cSize = 2 * support + 1;
    ASKAPLOG_DEBUG_STR(logger, "Convolution function support = "
                           << support << " pixels, convolution function size = "
                           << cSize << " pixels");
    return support;
}

int WProjectVisGridder::cIndex(int row, int pol, int chan)
{
    const int plane = itsCMap(row, pol, chan);
    //ASKAPDEBUGASSERT(plane >= 0);
    if (plane >=0) notifyOfWPlaneUse(plane);
    return plane;
}


/// @brief static method to create gridder
/// @details Each gridder should have a static factory method, which is
/// able to create a particular type of the gridder and initialise it with
/// the parameters taken form the given parset. It is assumed that the
/// method receives a subset of parameters where the gridder name is already
/// taken out.
/// @param[in] parset input parset file
/// @return a shared pointer to the gridder instance
IVisGridder::ShPtr WProjectVisGridder::createGridder(const LOFAR::ParameterSet& parset)
{
    const double wmax = parset.getDouble("wmax", 35000.0);
    const int nwplanes = parset.getInt32("nwplanes", 65);
    const double cutoff = parset.getDouble("cutoff", 1e-3);
    const int oversample = parset.getInt32("oversample", 8);
    const int maxSupport = parset.getInt32("maxsupport", 256);
    const int limitSupport = parset.getInt32("limitsupport", 0);
    const string tablename = parset.getString("tablename", "");
    const float alpha=parset.getFloat("alpha", 1.);
    const bool useDouble = parset.getBool("usedouble",false);
    const bool shareCF = parset.getBool("sharecf",false);

    ASKAPLOG_INFO_STR(logger, "Gridding using W projection with " << nwplanes << " w-planes");
    ASKAPLOG_INFO_STR(logger, "Using " << (useDouble ? "double":"single")<<
                      " precision to calculate convolution functions");
    boost::shared_ptr<WProjectVisGridder> gridder(new WProjectVisGridder(wmax, nwplanes,
            cutoff, oversample, maxSupport, limitSupport, tablename, alpha, useDouble, shareCF));
    gridder->configureGridder(parset);
    gridder->configureWSampling(parset);

    return gridder;
}

/// @brief additional operations to configure gridder
/// @details This method is supposed to be called from createGridder and could be
/// used in derived classes to avoid too much duplication of the code. For this
/// particular class it configures variable/offset support and cutoff behavior.
/// @param[in] parset input parset file
void WProjectVisGridder::configureGridder(const LOFAR::ParameterSet& parset)
{
    const bool planeDependentSupport = parset.getBool("variablesupport", false);

    if (planeDependentSupport) {
        ASKAPLOG_INFO_STR(logger, "Support size will be calculated separately for each w-plane");
    } else {
        ASKAPLOG_INFO_STR(logger, "Common support size will be used for all w-planes");
    }

    WProjectVisGridder::planeDependentSupport(planeDependentSupport);

    const bool offsetSupport = parset.getBool("offsetsupport", false);
    ASKAPCHECK((!offsetSupport && !planeDependentSupport) || planeDependentSupport,
               "offsetsupport option of the gridder should only be used together with variablesupport option");
    WProjectVisGridder::offsetSupport(offsetSupport);

    const bool absCutoff = parset.getBool("cutoff.absolute", false);

    if (absCutoff) {
        ASKAPLOG_INFO_STR(logger, "Cutoff value of " << itsCutoff <<
            " will be treated as an absolute threshold during CF generation");
    } else {
        ASKAPLOG_INFO_STR(logger, "Cutoff value of " << itsCutoff <<
            " will be treated as a threshold relative to the peak during CF generation");
        ASKAPCHECK(itsCutoff > 0.0, "Cutoff must be positive");
        ASKAPCHECK(itsCutoff < 1.0, "Cutoff must be less than 1.0");
    }

    setAbsCutoffFlag(absCutoff);
}


/// @brief obtain buffer used to create convolution functions
/// @return a reference to the buffer held as a shared pointer
casa::Matrix<casa::DComplex> WProjectVisGridder::getCFBuffer() const
{
    ASKAPDEBUGASSERT(itsCFBuffer);
    return *itsCFBuffer;
}

/// @brief obtain buffer used to create convolution functions
/// @return a reference to the buffer held as a shared pointer
casa::Matrix<casa::Complex> WProjectVisGridder::getCFBufferF() const
{
    ASKAPDEBUGASSERT(itsCFBufferF);
    return *itsCFBufferF;
}

/// @brief assignment operator
/// @details Defined as private, so it can't be called (to enforce usage of the
/// copy constructor
/// @param[in] other input object
/// @return reference to itself
WProjectVisGridder& WProjectVisGridder::operator=(const WProjectVisGridder &)
{
    ASKAPTHROW(AskapError, "This method is not supposed to be called!");
    return *this;
}


} // namespace askap
} // namespace synthesis
