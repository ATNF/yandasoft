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

// System includes
#include <cmath>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/BasicSL/Constants.h>
#include <fft/FFTWrapper.h>
#include <profile/AskapProfiler.h>

// Local package includes
#include <gridding/WProjectVisGridder.h>
#include <gridding/SupportSearcher.h>

ASKAP_LOGGER(logger, ".gridding.wprojectvisgridder");

namespace askap {
namespace synthesis {

WProjectVisGridder::WProjectVisGridder(const double wmax,
                                       const int nwplanes,
                                       const double cutoff,
                                       const int overSample,
                                       const int maxSupport,
                                       const int limitSupport,
                                       const std::string& name) :
        WDependentGridderBase(wmax, nwplanes),
        itsMaxSupport(maxSupport), itsCutoff(cutoff), itsLimitSupport(limitSupport),
        itsPlaneDependentCFSupport(false), itsOffsetSupportAllowed(false), itsCutoffAbs(false)
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
        itsCutoffAbs(other.itsCutoffAbs) {}


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

    for (int i = 0; i < nSamples; ++i) {
        const double w = (rotatedUVW(i)(2)) / (casa::C::c);
        for (int chan = 0; chan < nChan; ++chan) {
            for (int pol = 0; pol < nPol; ++pol) {
                const double freq = acc.frequency()[chan];
                /// Calculate the index into the convolution functions
                itsCMap(i, pol, chan) = getWPlane(w * freq);
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

    // Now we step through the w planes, starting the furthest
    // out. We calculate the support for that plane and use it
    // for all the others.

    // We pad here to do sinc interpolation of the convolution
    // function in uv space
    casa::Matrix<casa::DComplex> thisPlane = getCFBuffer();
    ASKAPDEBUGASSERT(thisPlane.nrow() == casa::uInt(nx));
    ASKAPDEBUGASSERT(thisPlane.ncolumn() == casa::uInt(ny));

    for (int iw = 0; iw < nWPlanes(); ++iw) {
        thisPlane.set(0.0);

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
                const float r2 = x2 + y2;

                if (r2 < 1.0) {
                    const double phase = w * (1.0 - sqrt(1.0 - r2));
                    const float wt = ccfx(ix) * ccfy(iy);
                    ASKAPDEBUGASSERT(ix - qnx / 2 + nx / 2 < nx);
                    ASKAPDEBUGASSERT(iy - qny / 2 + ny / 2 < ny);
                    ASKAPDEBUGASSERT(ix + nx / 2 >= qnx / 2);
                    ASKAPDEBUGASSERT(iy + ny / 2 >= qny / 2);
                    thisPlane(ix - qnx / 2 + nx / 2, iy - qny / 2 + ny / 2) = casa::DComplex(wt * cos(phase), -wt * sin(phase));
                    //thisPlane(ix-qnx/2+nx/2, iy-qny/2+ny/2)=casa::DComplex(wt*cos(phase));
                }
            }
        }

        // At this point, we have the phase screen multiplied by the spheroidal
        // function, sampled on larger cellsize (itsOverSample larger) in image
        // space. Only the inner qnx, qny pixels have a non-zero value

        // Now we have to calculate the Fourier transform to get the
        // convolution function in uv space
        //casa::Matrix<casa::DComplex> buffer(thisPlane.nrow(),thisPlane.ncolumn());
        //casa::convertArray<casa::DComplex,casa::Complex>(buffer,thisPlane);
        scimath::fft2d(thisPlane, true);
        //scimath::fft2d(buffer, true);
        //casa::convertArray<casa::Complex,casa::DComplex>(thisPlane,buffer);


        /*
            for (uint xx=0;xx<thisPlane.nrow();++xx) {
             for (uint yy=0;yy<thisPlane.ncolumn();++yy) {
                  ASKAPCHECK(!std::isinf(casa::abs(thisPlane(xx,yy))), "Infinite value detected for plane="<<iw<<" at "<<xx<<","<<yy<<" "<<thisPlane(xx,yy));
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
            cfSupport = extractSupport(thisPlane);
            const int support = cfSupport.itsSize;

            ASKAPCHECK(support*itsOverSample < nx / 2,
                       "Overflowing convolution function for w-plane " << iw <<
                       " - increase maxSupport or decrease overSample; support=" << support << " oversample=" << itsOverSample <<
                       " nx=" << nx);
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
                for (int iy = -support; iy < support; ++iy) {
                    for (int ix = -support; ix < support; ++ix) {
                        ASKAPDEBUGASSERT((ix + support >= 0) && (iy + support >= 0));
                        ASKAPDEBUGASSERT(ix + support < int(itsConvFunc[plane].nrow()));
                        ASKAPDEBUGASSERT(iy + support < int(itsConvFunc[plane].ncolumn()));
                        ASKAPDEBUGASSERT((ix + cfSupport.itsOffsetU)*itsOverSample + fracu + nx / 2 >= 0);
                        ASKAPDEBUGASSERT((iy + cfSupport.itsOffsetV)*itsOverSample + fracv + ny / 2 >= 0);
                        ASKAPDEBUGASSERT((ix + cfSupport.itsOffsetU)*itsOverSample + fracu + nx / 2 < int(thisPlane.nrow()));
                        ASKAPDEBUGASSERT((iy + cfSupport.itsOffsetV)*itsOverSample + fracv + ny / 2 < int(thisPlane.ncolumn()));
                        //if (w < 0) {
                        //    itsConvFunc[plane](ix + support, iy + support) =
                        //        conj(thisPlane((ix + cfSupport.itsOffsetU) * itsOverSample + fracu + nx / 2,
                        //              (iy + cfSupport.itsOffsetV) * itsOverSample + fracv + ny / 2));
                        //} else {
                        itsConvFunc[plane](ix + support, iy + support) =
                            thisPlane((ix + cfSupport.itsOffsetU) * itsOverSample + fracu + nx / 2,
                                   (iy + cfSupport.itsOffsetV) * itsOverSample + fracv + ny / 2);
                        //}
                    } // for ix
                } // for iy
            } // for fracv
        } // for fracu

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
                itsConvFunc[plane] /= casa::Complex(norm);
            }
        } // for plane
    } // for iw

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
    itsCFBuffer.reset();
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
    ASKAPDEBUGASSERT(plane >= 0);
    notifyOfWPlaneUse(plane);
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
    ASKAPLOG_INFO_STR(logger, "Gridding using W projection with " << nwplanes << " w-planes");
    boost::shared_ptr<WProjectVisGridder> gridder(new WProjectVisGridder(wmax, nwplanes, cutoff, oversample,
            maxSupport, limitSupport, tablename));
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
        ASKAPLOG_INFO_STR(logger, "Cutoff value of " << itsCutoff << " will be treated as an absolute threshold during CF generation");
    } else {
        ASKAPLOG_INFO_STR(logger, "Cutoff value of " << itsCutoff << " will be treated as a threshold relative to the peak during CF generation");
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
