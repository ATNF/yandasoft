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
#include <askap/askap_synthesis.h>

// System includes
#include <cmath>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <askap/scimath/fft/FFTWrapper.h>
#include <askap/profile/AskapProfiler.h>

// Local package includes
#include <askap/gridding/MPIWProjectVisGridder.h>
#include <askap/gridding/SupportSearcher.h>

ASKAP_LOGGER(logger, ".gridding.MPIwprojectvisgridder");

namespace askap {
namespace synthesis {


MPI_Aint MPIWProjectVisGridder::itsWindowSize;
int      MPIWProjectVisGridder::itsWindowDisp;
MPI_Win  MPIWProjectVisGridder::itsWindowTable;
MPI_Comm MPIWProjectVisGridder::itsNodeComms;
MPI_Comm MPIWProjectVisGridder::itsNonRankZeroComms;
MPI_Group MPIWProjectVisGridder::itsWorldGroup;
MPI_Group MPIWProjectVisGridder::itsGridderGroup;
int MPIWProjectVisGridder::itsNodeSize;
int MPIWProjectVisGridder::itsNodeRank;
imtypeComplex* MPIWProjectVisGridder::itsMpiSharedMemory = nullptr;
std::mutex MPIWProjectVisGridder::ObjCountMutex;
unsigned int MPIWProjectVisGridder::ObjCount = 0;

std::vector<casa::Matrix<casa::Complex> > MPIWProjectVisGridder::theirCFCache;
std::vector<std::pair<int,int> > MPIWProjectVisGridder::theirConvFuncOffsets;

/// @brief a helper method for a ref copy of casa arrays held in
/// stl vector
/// @param[in] in input array
/// @param[out] out output array (will be resized)
/// @return size of the cache in bytes (assuming Complex array elements)

template<typename T>
size_t deepRefCopyOfSTDVector(const std::vector<T> &in,
                            std::vector<T> &out)
{
   out.resize(in.size());
   size_t total = 0;
   
   const typename std::vector<T>::const_iterator inEnd = in.end();
   typename std::vector<T>::iterator outIt = out.begin();
   for (typename std::vector<T>::const_iterator inIt = in.begin();
       inIt != inEnd; ++inIt,++outIt) {
       outIt->reference(*inIt);
       total += outIt->nelements()*sizeof(casa::Complex)+sizeof(T);
   }
   
   return total;
}

MPIWProjectVisGridder::MPIWProjectVisGridder(const double wmax,
                                       const int nwplanes,
                                       const double cutoff,
                                       const int overSample,
                                       const int maxSupport,
                                       const int limitSupport,
                                       const std::string& name,
                                       const float alpha,
                                       const bool shareCF) :
        WDependentGridderBase(wmax, nwplanes, alpha),
        itsMaxSupport(maxSupport), itsCutoff(cutoff), itsLimitSupport(limitSupport),
        itsPlaneDependentCFSupport(false), itsOffsetSupportAllowed(false), itsCutoffAbs(false),
        itsShareCF(shareCF)
{
    ASKAPCHECK(overSample > 0, "Oversampling must be greater than 0");
    ASKAPCHECK(maxSupport > 0, "Maximum support must be greater than 0")
    itsSupport = 0;
    itsOverSample = overSample;
    setTableName(name);
    itsConvFunc.resize(nWPlanes()*itsOverSample*itsOverSample);

    std::lock_guard<std::mutex> lk(ObjCountMutex);
    ObjCount += 1;
}

MPIWProjectVisGridder::~MPIWProjectVisGridder()
{
    {
        std::lock_guard<std::mutex> lk(ObjCountMutex);
        ObjCount -= 1;
    }
    if ( ObjCount == 0 ) {
        MPI_Group_free(&itsWorldGroup);
        MPI_Group_free(&itsGridderGroup);
        MPI_Comm_free(&itsNodeComms);
        MPI_Comm_free(&itsNonRankZeroComms);
        MPI_Win_free(&itsWindowTable);
    }
}

/// @brief copy constructor
/// @details It is required to decouple internal arrays in the input
/// object and the copy.
/// @param[in] other input object
MPIWProjectVisGridder::MPIWProjectVisGridder(const MPIWProjectVisGridder &other) :
        IVisGridder(other), WDependentGridderBase(other),
        itsCMap(other.itsCMap.copy()), itsMaxSupport(other.itsMaxSupport),
        itsCutoff(other.itsCutoff), itsLimitSupport(other.itsLimitSupport),
        itsPlaneDependentCFSupport(other.itsPlaneDependentCFSupport),
        itsOffsetSupportAllowed(other.itsOffsetSupportAllowed),
        itsCutoffAbs(other.itsCutoffAbs),
        itsShareCF(other.itsShareCF) 
{
	std::lock_guard<std::mutex> lk(ObjCountMutex);
	ObjCount += 1;
}


/// Clone a copy of this Gridder
IVisGridder::ShPtr MPIWProjectVisGridder::clone()
{
    return IVisGridder::ShPtr(new MPIWProjectVisGridder(*this));
}

/// @brief initialise sum of weights
/// @details We keep track the number of times each convolution function is used per
/// channel and polarisation (sum of weights). This method is made virtual to be able
/// to do gridder specific initialisation without overriding initialiseGrid.
/// This method accepts no parameters as itsShape, itsNWPlanes, etc should have already
/// been initialised by the time this method is called.
void MPIWProjectVisGridder::initialiseSumOfWeights()
{
    resizeSumOfWeights(nWPlanes());
    zeroSumOfWeights();
}

/// Initialize the convolution function into the cube. If necessary this
/// could be optimized by using symmetries.
void MPIWProjectVisGridder::initIndices(const accessors::IConstDataAccessor& acc)
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

    const casacore::Vector<casacore::RigidVector<double, 3> > &rotatedUVW = acc.rotatedUVW(getTangentPoint());
    const casacore::Vector<casacore::Double> & chanFreq = acc.frequency();

    for (int i = 0; i < nSamples; ++i) {
        const double w = (rotatedUVW(i)(2)) / (casacore::C::c);
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
void MPIWProjectVisGridder::initConvolutionFunction(const accessors::IConstDataAccessor&)
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
                itsConvFunc[plane](ix + cCenter, iy + cCenter) =  casacore::Complex(1.0, wKernelPix);
            }
        }

      }

      return;
    }

    ASKAPLOG_INFO_STR(logger, "theirCFCache size: " << theirCFCache.size()
		    << ", itsShareCF: " << itsShareCF);
    if (itsShareCF && theirCFCache.size()>0) {
        // we already have what we need
        itsSupport = 1;
        size_t size = deepRefCopyOfSTDVector(theirCFCache,itsConvFunc)/1024/1024;
        ASKAPLOG_INFO_STR(logger, "Using cached convolution functions ("<<size<<" MB)");
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
    initCFBuffer(casacore::uInt(nx), casacore::uInt(ny));

    /// We want nx * ccellx = overSample * itsShape(0) * cellx

    const int qnx = nx / itsOverSample;
    const int qny = ny / itsOverSample;
    ASKAPDEBUGASSERT((qnx != 0) && (qny != 0));

    // Find the actual cellsizes in x and y (radians) after over
    // oversampling (in uv space)
    const double ccellx = double(itsShape(0)) * cellx / double(qnx);
    const double ccelly = double(itsShape(1)) * celly / double(qny);

    casacore::Vector<float> ccfx(qnx);
    casacore::Vector<float> ccfy(qny);

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
    casacore::Matrix<imtypeComplex> thisPlane(getCFBuffer());
    ASKAPDEBUGASSERT(thisPlane.nrow() == casacore::uInt(nx));
    ASKAPDEBUGASSERT(thisPlane.ncolumn() == casacore::uInt(ny));

    for (int iw = 0; iw < nWPlanes(); ++iw) {
        thisPlane.set(0.0);

        //const double w = isPSFGridder() ? 0. : 2.0f*casacore::C::pi*getWTerm(iw);
        const double w = 2.0f * casacore::C::pi * getWTerm(iw);

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
                    thisPlane(ix - qnx / 2 + nx / 2, iy - qny / 2 + ny / 2) =
                    imtypeComplex(wt * cos(phase), -wt * sin(phase));
                }
                    //thisPlane(ix-qnx/2+nx/2, iy-qny/2+ny/2)=imtypeComplex(wt*cos(phase));
            }
        }

        // At this point, we have the phase screen multiplied by the spheroidal
        // function, sampled on larger cellsize (itsOverSample larger) in image
        // space. Only the inner qnx, qny pixels have a non-zero value

        // Now we have to calculate the Fourier transform to get the
        // convolution function in uv space
        scimath::fft2d(thisPlane, true);

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
            // fail here if the cutoff level is on the edge of the image
            ASKAPCHECK((support+1)*itsOverSample < nx / 2,
                       "Overflowing convolution function for w-plane " << iw <<
                       " - increase maxSupport or cutoff or decrease overSample; support=" <<
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
                        //if (w < 0) {
                        //    itsConvFunc[plane](ix + support, iy + support) =
                        //        conj(thisPlane(kx, ky);
                        //} else {
                        itsConvFunc[plane](ix + support, iy + support) =thisPlane(kx, ky);
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

        const double norm = sum(casacore::real(itsConvFunc[plane]));
        // ASKAPLOG_INFO_STR(logger, "Sum of convolution function = " << norm);
        ASKAPDEBUGASSERT(norm > 0.);

        if (norm > 0.) {
            const casacore::Complex invNorm = casacore::Complex(1.0/norm);
            itsConvFunc[plane] *= invNorm;
        }
    } // for plane

    if (isSupportPlaneDependent()) {
        ASKAPLOG_DEBUG_STR(logger, "Convolution function cache has " << itsConvFunc.size() << " planes");
        ASKAPLOG_DEBUG_STR(logger, "Variable support size is used:");
        const size_t step = casacore::max(itsConvFunc.size() / itsOverSample / itsOverSample / 10, 1);

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

    // Save the CF to the cache
    ASKAPLOG_INFO_STR(logger,"itsShareCF: " << itsShareCF);
    if (itsShareCF) {
    	ASKAPLOG_INFO_STR(logger,"copy itsConvFunc vector to cache");
	
        if ( itsNodeRank == 0 ) {
            ASKAPLOG_INFO_STR(logger, "Copy to shared memory etc ...");
            ASKAPLOG_INFO_STR(logger, "number of elements in itsConvFunc: " << itsConvFunc.size());
        }
        // workout the size of itsConvFunc
        size_t total = 0; // in bytes
	    unsigned int nplane = itsConvFunc.size();
        for (auto it = itsConvFunc.begin();
            it != itsConvFunc.end(); ++it) {
            total += it->nelements() * sizeof(imtypeComplex);
		    // matrixSizePerPlaneVect.push_back(std::make_pair(it->nrow(),it->ncolumn()));
        }
        // create MPI shared memory
        setupMpiMemory(total);

        // copy itsConvFunc to shared memory
        imtypeComplex* shareMemPtr = itsMpiSharedMemory; // a contiguous chunk of shared memory
	    // only itsNodeRank 0 does the copy 
	    if ( itsNodeRank == 0 ) {
      	    ASKAPLOG_INFO_STR(logger, "copy itsConvFunc to shared memory");
            for ( auto it = itsConvFunc.begin();
            	it != itsConvFunc.end(); ++it) {
            	std::copy(it->data(),it->data() + it->nelements(),shareMemPtr);
            	shareMemPtr += it->nelements();
       	    }
	    }
	    MPI_Barrier(itsNodeComms);

    	// copy shared memory back to itsConvFunc
    	ASKAPLOG_INFO_STR(logger, "copy shared memory back to itsConvFunc");
     	shareMemPtr = itsMpiSharedMemory;
        for (int plane = 0; plane < nplane; plane++) {
            casacore::IPosition pos(2,itsConvFunc[plane].nrow(),itsConvFunc[plane].ncolumn());
            casacore::Matrix<imtypeComplex> m(pos,shareMemPtr,casacore::SHARE);
            itsConvFunc[plane].reference(m);
            shareMemPtr += m.nelements();
        }	

        // this is the original code from WProjectVisGridder
        total = deepRefCopyOfSTDVector(itsConvFunc,theirCFCache);
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
MPIWProjectVisGridder::CFSupport MPIWProjectVisGridder::extractSupport(const casacore::Matrix<casacore::DComplex> &cfPlane) const
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
        const casacore::IPosition peakPos = ss.peakPos();
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
MPIWProjectVisGridder::CFSupport MPIWProjectVisGridder::extractSupport(const casacore::Matrix<casacore::Complex> &cfPlane) const
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
        const casacore::IPosition peakPos = ss.peakPos();
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
int MPIWProjectVisGridder::limitSupportIfNecessary(int support) const
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

int MPIWProjectVisGridder::cIndex(int row, int pol, int chan)
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
IVisGridder::ShPtr MPIWProjectVisGridder::createGridder(const LOFAR::ParameterSet& parset)
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

    ASKAPLOG_INFO_STR(logger, "Gridding using W projection with " << nwplanes << " w-planes");
    ASKAPLOG_INFO_STR(logger, "Gridding using maxsupport: " << maxSupport );
    ASKAPLOG_INFO_STR(logger, "Using " << (useDouble ? "double":"single")<<
                      " precision to calculate convolution functions");
    boost::shared_ptr<MPIWProjectVisGridder> gridder(new MPIWProjectVisGridder(wmax, nwplanes,
            cutoff, oversample, maxSupport, limitSupport, tablename, alpha, useDouble));
    gridder->configureGridder(parset);
    gridder->configureWSampling(parset);

    return gridder;
}

/// @brief additional operations to configure gridder
/// @details This method is supposed to be called from createGridder and could be
/// used in derived classes to avoid too much duplication of the code. For this
/// particular class it configures variable/offset support and cutoff behavior.
/// @param[in] parset input parset file
void MPIWProjectVisGridder::configureGridder(const LOFAR::ParameterSet& parset)
{
    const bool planeDependentSupport = parset.getBool("variablesupport", false);

    MPIWProjectVisGridder::planeDependentSupport(planeDependentSupport);

    const bool offsetSupport = parset.getBool("offsetsupport", false);
    ASKAPCHECK((!offsetSupport && !planeDependentSupport) || planeDependentSupport,
               "offsetsupport option of the gridder should only be used together with variablesupport option");
    MPIWProjectVisGridder::offsetSupport(offsetSupport);

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

    //itsShareCF = parset.getBool("sharecf",false);
    itsShareCF = parset.getBool("sharecf",true);

    // Each rank should only run the code below once.
    static bool doOnce = false;
    if ( doOnce ) return;
    doOnce = true;

    ASKAPLOG_INFO_STR(logger, "Setup MPI subgroup communicator for MPI Shared Memory");

    // The code below setup a MPI communicator for each node where the 
    // rank 0 (in the COMM_WORLD) of the first node is not included.
    int r;
    r = MPI_Comm_group(MPI_COMM_WORLD, &itsWorldGroup);
    ASKAPASSERT(r == MPI_SUCCESS);
    const int exclude_ranks[1] = {0};
    r = MPI_Group_excl(itsWorldGroup, 1, exclude_ranks, &itsGridderGroup);
    ASKAPASSERT(r == MPI_SUCCESS);
    r = MPI_Comm_create_group(MPI_COMM_WORLD, itsGridderGroup, 0, &itsNonRankZeroComms);
    ASKAPASSERT(r == MPI_SUCCESS);
    r = MPI_Comm_split_type(itsNonRankZeroComms, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &itsNodeComms);
	ASKAPASSERT(r == MPI_SUCCESS);
    // number of ranks within the node
    MPI_Comm_size(itsNodeComms, &itsNodeSize);
    // rank within a node
    MPI_Comm_rank(itsNodeComms, &itsNodeRank);
    MPI_Barrier(itsNodeComms);
}


/// @brief obtain buffer used to create convolution functions
/// @return a reference to the buffer held as a shared pointer
casacore::Matrix<imtypeComplex> MPIWProjectVisGridder::getCFBuffer() const
{
    ASKAPDEBUGASSERT(itsCFBuffer);
    return *itsCFBuffer;
}

/// @brief assignment operator
/// @details Defined as private, so it can't be called (to enforce usage of the
/// copy constructor
/// @param[in] other input object
/// @return reference to itself
MPIWProjectVisGridder& MPIWProjectVisGridder::operator=(const MPIWProjectVisGridder &)
{
    ASKAPTHROW(AskapError, "This method is not supposed to be called!");
    return *this;
}

void  MPIWProjectVisGridder::setupMpiMemory(size_t bufferSize /* in bytes */)
{
    static bool runOnce = false;
    if ( runOnce  ) {
	    ASKAPLOG_INFO_STR(logger,"itsNodeRank: " << itsNodeRank << " - mpi shared memory already setup");
        return;
    }

    runOnce = true;

    ASKAPLOG_INFO_STR(logger,"itsNodeRank: " << itsNodeRank << ", bufferSize: " << bufferSize);

    if ( itsNodeComms != MPI_COMM_NULL ) {
        // the case where every rank has a gridder object
        size_t memSizeInBytes = 0;
        if ( itsNodeRank == 0 ) {
            memSizeInBytes = bufferSize;
        }
        int r = MPI_Win_allocate_shared(memSizeInBytes,sizeof(imtypeComplex),
                                MPI_INFO_NULL, itsNodeComms, &itsMpiSharedMemory,
                                &itsWindowTable);
        ASKAPCHECK(r == MPI_SUCCESS, "MPI_Win_allocate_shared failed.");
        // For itsNodeRanks != 0, get their itsMpiSharedMemory pointer variable to point the
        // start of the MPI shared memory allocated by itsNodeRank = 0
        if ( itsNodeRank != 0 ) {
            int r = MPI_Win_shared_query(itsWindowTable, 0, &itsWindowSize, &itsWindowDisp, &itsMpiSharedMemory);
            ASKAPCHECK(r == MPI_SUCCESS, "MPI_Win_shared_query failed.");
        }
        MPI_Barrier(itsNodeComms);
	    //unsigned long* val;
	    //int flag = 1;
	    //r = MPI_Win_get_attr(itsWindowTable,MPI_WIN_SIZE,&val,&flag);
    }
}

} // namespace askap
} // namespace synthesis
