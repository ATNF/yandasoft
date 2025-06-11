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

#include <askap/measurementequation/WienerPreconditioner.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.wienerpreconditioner");

#include <askap/askap/AskapError.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/profile/AskapProfiler.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Arrays/Vector.h>
#include <askap/scimath/fft/FFT2DWrapper.h>

// for debugging - to export intermediate images
#include <askap/scimath/utils/ImageUtils.h>

#include <Common/OpenMP.h>

using namespace casa;

#include <iostream>
#include <cmath>
using std::abs;

namespace askap
{
  namespace synthesis
  {

    bool WienerPreconditioner::itsUseCachedPcf = false;
    double WienerPreconditioner::itsAveWgtSum = 0.0;
    casacore::Matrix<float> WienerPreconditioner::itsPcf;

    WienerPreconditioner::WienerPreconditioner() :
        itsParameter(0.0), itsDoNormalise(false), itsUseRobustness(false)
    {
    }

    /// @brief constructor with explicitly defined noise power
    /// @param[in] noisepower parameter of the
    /// @param[in] normalise if true, PSF is normalised during filter construction
    WienerPreconditioner::WienerPreconditioner(float noisepower, bool normalise) :
            itsParameter(noisepower), itsDoNormalise(normalise),
            itsUseRobustness(false) {}

    /// @brief constructor with explicitly defined robustness
    /// @details In this version, the noise power is calculated from
    /// the robustness parameter
    /// @param[in] robustness robustness parameter (roughly matching Briggs' weighting)
    /// @note Normalisation of PSF is always used when noise power is defined via robustness
    WienerPreconditioner::WienerPreconditioner(float robustness) : itsParameter(robustness),
            itsDoNormalise(true), itsUseRobustness(true) {}


    /// @brief copy constructor
    /// @param[in] other an opject to copy from
    WienerPreconditioner::WienerPreconditioner(const WienerPreconditioner &other) :
          IImagePreconditioner(other),
          itsParameter(other.itsParameter), itsDoNormalise(other.itsDoNormalise),
          itsUseRobustness(other.itsUseRobustness)
    {
       if (other.itsTaperCache) {
           itsTaperCache.reset(new GaussianTaperCache(*(other.itsTaperCache)));
       }
    }

    IImagePreconditioner::ShPtr WienerPreconditioner::clone()
    {
        return IImagePreconditioner::ShPtr(new WienerPreconditioner(*this));
    }

    bool WienerPreconditioner::doPreconditioning(casacore::Array<float>& psf,
                                                 casacore::Array<float>& dirty,
                                                 casacore::Array<float>& pcfIn) const
    {
      // Note this routine uses generic arrays, but they are really always 2d Matrices, all with the same shape
      ASKAPTRACE("WienerPreconditioner::doPreconditioning");
      if (!itsUseRobustness && (itsParameter < 1e-6)) {
          return false;
      }

      ASKAPCHECK(psf.shape().conform(dirty.shape()),
          "Dirty image and PSF do not conform - shapes: " <<
          dirty.shape() << " & " << psf.shape());

      const casacore::IPosition shape(2,psf.shape()(0),psf.shape()(1));

      bool useCachedFilter = true;
      if (itsWienerfilter.shape() == 0 || itsWienerfilter.shape() != shape) {
        useCachedFilter = false;
        ASKAPLOG_DEBUG_STR(logger,"Can't reuse Wiener filter: filter shape="<<itsWienerfilter.shape()<<", psf shape="<<shape);
        itsWienerfilter.resize(shape);
      }

      if (itsPcf.shape() != shape) {
          ASKAPLOG_DEBUG_STR(logger,"Can't reuse Wiener filter: PCF shape="<<itsPcf.shape()<<", psf shape="<<shape);
          itsUseCachedPcf = false;
      }

      bool newFilter = true;
      if (itsUseCachedPcf) {
          // If the preconditioner function is cached, use that.
          if (useCachedFilter) {
              ASKAPLOG_INFO_STR(logger, "Applying cached Wiener filter");
          } else {
              ASKAPLOG_INFO_STR(logger, "Applying Wiener filter generated from cached PCF weights");
          }
      } else if (pcfIn.shape() == 0) {
          // If the preconditioner function is undefined, use the psf
          ASKAPLOG_INFO_STR(logger, "Applying old-style Wiener filter");
          itsPcf.reference(psf.nonDegenerate());
          // could also use the stored version, but leave for now.
          useCachedFilter = false;
          newFilter = false;
      } else {
          ASKAPLOG_INFO_STR(logger, "Generating and applying Wiener filter of size " << pcfIn.nonDegenerate().shape());
          // Usually only get here if itsPcf is empty. But with Nyquist gridding that may not be the case.
          // So check, and if it has already be set but is a different size, reset it.
          if (!itsPcf.empty() && !itsPcf.shape().isEqual(pcfIn.nonDegenerate().shape())) {
              itsPcf.resize(pcfIn.nonDegenerate().shape());
          }
          //itsPcf.reference(pcfIn);
          itsPcf = pcfIn.nonDegenerate();
      }
      ASKAPCHECK(itsPcf.shape().conform(shape),
          "PSF and preconditioner function do not conform - shapes: " <<
          shape << " & " << itsPcf.shape());

      if (itsUseRobustness) {
          ASKAPLOG_INFO_STR(logger,
              "Wiener filter noise power defined via robustness = " << itsParameter);
      } else {
          ASKAPLOG_INFO_STR(logger,
              "Wiener filter noise power = " << itsParameter);
      }

      float maxPSFBefore = casacore::max(psf);
      ASKAPCHECK(maxPSFBefore > 0.0, "Peak of PSF before Wiener filtering, " <<
          maxPSFBefore << ", is less than or equal to zero");
      ASKAPLOG_INFO_STR(logger, "Peak of PSF before Wiener filtering = " << maxPSFBefore);

      /*
       * Not using the PSF to generate the filter, so this isn't needed and complicates subsequent normalisation
       *
      if (itsDoNormalise) {
          ASKAPLOG_INFO_STR(logger, "The PSF will be normalised to 1 before filter construction");
          psf=psf/maxPSFBefore;
          // dirty=dirty/maxPSFBefore;
          maxPSFBefore=1.0;
      }
       */

      // MV: there is some technical debt here. First, using statics in this context is a no-no. Second, the code should've worked with Matrix
      // upfront. At least, nonDegenerate should have starting axis index passed to it - but leave as is for now.
      // also, this method is too long - more structure is necessary

      casacore::Matrix<casacore::Float> psf2D(psf.nonDegenerate());
      casacore::Matrix<casacore::Float> dirty2D(dirty.nonDegenerate());

      // Make the scratch array into which we will calculate the Wiener filter
      casacore::Matrix<casacore::Complex> scratch(shape);
      convertArray<casacore::Complex,casacore::Float>(scratch, itsPcf);

      // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
      // skip it here as don't care about particular harmonic placement in the Fourier space in most cases. However, tapering and filtering
      // would need to be altered to support different harmonic placement - this would add to code complexity although may provide some performance benefits
      // MV: also note we only work with single-precision complex here for some reason (set by the interface of the method, the class is not a template)
      // Limit number of fft threads to 8 (more is slower for our fft sizes)
      scimath::FFT2DWrapper<casacore::Complex> fft2d(true,8);

      if (!itsUseCachedPcf) {

        // the filter has been accumulated in the image domain, so transform to uv

        fft2d(scratch, true);

        // The filter is rescaling Fourier components of dirty and psf based on
        // the value of non-zero Fourier components of pcf. So set small
        // erroneous components to zero.
        const float scratchThreshold = pcfThreshold(scratch);
        if (scratchThreshold > 0.0) {
          // "storage order is the same as in FORTRAN, i.e. memory varies most rapidly with the first index"
          for (int y=0; y<shape(1); ++y) {
            for (int x=0; x<shape(0); ++x) {
              if (real(scratch(x,y)) < scratchThreshold) {
                scratch(x,y) = 0.0;
              }
            }
          }
        }

        if (newFilter) {

          // The PCF should be set up to have nearest-neighbour gridding of
          // SNR weights in the real part and nearest-neighbour gridding of
          // weights*kernel-size (in pixels) in the imaginary part. From this
          // the idea is to build a filter that alters the local density of
          // weighted data without breaking local projections.

          // get weighted sum of kernel size from imaginary part of the data
          casacore::Matrix<casacore::Float> kernelWidthMatrix(shape(0),shape(1),0.0);
          // "storage order is the same as in FORTRAN, i.e. memory varies most rapidly with the first index"
          for (int y=0; y<shape(1); ++y) {
            for (int x=0; x<shape(0); ++x) {
              casacore::Complex c = scratch(x,y);
              if (abs(real(c)) > scratchThreshold) {
                kernelWidthMatrix(x,y)=abs(imag(c)/real(c));
              }
            }
          }

          // for debugging - to export kernel widths
          //scimath::saveAsCasaImage("kernelwidths.uv",kernelWidthMatrix);
          //throw 1;

          // Build the filter as follows:
          // * set a large box (maximum kernel size) around each pixel
          // * find largest local kernel size
          // * make non-zero if there is data within this radius
          // * estimate weight based on the average of near by points,
          //   where "near by" can be different to the kernel width.

          const int maxKernelWidth = ceil(max(kernelWidthMatrix));
          ASKAPCHECK(maxKernelWidth > 0, "Maximum kernel width in the Wiener filter, " <<
              maxPSFBefore << ", is less than or equal to zero");
          const int boxWidth = maxKernelWidth;
          ASKAPDEBUGASSERT(boxWidth>0);

          // the factor for running mean (relative to the local kernel width)
          const int extra = 2;

          // Are boxes faster than simply searching scratch? Test. - No

          // As with the threshold search above, here we assume that
          // zero-padding is used to over sample the PSF, meaning that we
          // can be sure that visibilities are not gridded to the edge of
          // the uv plane. Should test for the presence of zero-padding...

          // set pcf to to contain the new weights, so they don't have to be regenerated.
          // note that it is going from the image domain to the uv domain.
          itsPcf = 0.0;
#ifdef HAVE_MPI
          double timerStart = MPI_Wtime();
#endif
          #pragma omp parallel
          {
          #pragma omp master
          {
              uint nthreads = LOFAR::OpenMP::numThreads();
              if (nthreads>1) ASKAPLOG_INFO_STR(logger, "Computing Wiener filter using "<<nthreads<< " threads");
          }
          #pragma omp for
          for (int y=extra*boxWidth/2; y<shape[1]-extra*boxWidth/2; ++y) {
              int boxStart1 = y - boxWidth/2;
              for (int x=extra*boxWidth/2; x<shape[0]-extra*boxWidth/2; ++x) {

                int boxStart0 = x - boxWidth/2;

                int localCount = 0;
                // double localSum = 0.0;
                int regionCount = 0;
                double regionSum = 0.0;

                //const int kernelWidth = ceil(max(
                //    kernelWidthMatrix(Slice(boxStart0,boxWidth),Slice(boxStart1,boxWidth))));

                //try writing out max
                casacore::Float kernelW = 0;
                for (int yb=boxStart1; yb < boxStart1+boxWidth; yb++) {
                    for (int xb=boxStart0; xb < boxStart0+boxWidth; xb++) {
                        kernelW = max(kernelW,kernelWidthMatrix(xb,yb));
                    }
                }
                // ceil was experiencing some roundoff error, so changing to round.
                //const int kernelWidth = ceil(kernelW);
                const int kernelWidth = round(kernelW);

                if (kernelWidth>0) {

                  // reset box to the kernelWidth
                  const int regionWidth = 1 + extra*(kernelWidth-1);
                  ASKAPDEBUGASSERT(regionWidth>=kernelWidth);
                  boxStart0 = x - regionWidth/2;
                  boxStart1 = y - regionWidth/2;

                    const float localRadiusSq = 0.25 * kernelWidth*kernelWidth;
                    const float regionRadiusSq = 0.25 * regionWidth*regionWidth;

                    for (int yb=boxStart1; yb<boxStart1+regionWidth; ++yb) {
                      const int dy = yb - boxStart1 - regionWidth/2;
                      const int dy2 = dy * dy;
                      for (int xb=boxStart0; xb<boxStart0+regionWidth; ++xb) {
                        const int dx = xb - boxStart0 - regionWidth/2;
                        const float val = real(scratch(xb,yb));
                        if ( val > scratchThreshold) {
                          const float rsq = dx*dx + dy2;
                          if (rsq<=regionRadiusSq) {
                            regionCount += 1;
                            regionSum += val;
                            if (rsq<=localRadiusSq) {
                              localCount += 1;
                              //localSum += localBox(xb,yb); // Unused, why?
                            }
                          }
                        }
                      }
                  }
                }


                if (localCount > 0) {
                  itsPcf(x,y) = regionSum/double(regionCount);
                }

              } // x
          } // y
          } // parallel
#ifdef HAVE_MPI
          double timerStop = MPI_Wtime();
          ASKAPLOG_INFO_STR(logger,"Wiener filter calculation took "<< timerStop-timerStart<<"s");
#else
          ASKAPLOG_INFO_STR(logger,"Cant compute how long Wiener filter took because MPI library is not installed (i.e cant use MPI_Wtime())");
#endif

          // Calc ave SNR weight-sum *over visibilities* (not pixels).
          // const casacore::Array<float> wgts(real(scratch.asArray()));
          casacore::Array<double> wgts(shape);
          casacore::convertArray<double,float>(wgts, real(scratch));
          // The following approx assumes that vis have equal SNR weights.
          // * see D. Briggs thesis
          itsAveWgtSum = sum(wgts*wgts) / sum(wgts);

          // update the scratch array
          convertArray<casacore::Complex,casacore::Float>(scratch, itsPcf);

          itsUseCachedPcf = true;

        } // if (newFilter)

      } // if (!itsUseCachedPcf)

      if (!useCachedFilter) {

        double aveWgtSum = itsAveWgtSum;

        if (itsTaperCache) {
          // Apply Gaussian taper to the preconditioner function (in the image domain)
          // The resulting uv convolution has various effects:
          // 1. it averages over neighbouring pixels, and as such has similarities with super-uniform weighting,
          // 2. it smears the PCF sampling fn out into unsampled regions, ending will small values,
          // 3. it can result in extra low-level rumble across the uv plane.
          // Effect #1 is useful, the others are not, particularly when robustness is set towards uniform
          // weighting. The latter two effects can be removed by reapplying the original thresholding.
          // So, save the original thresholding
          casacore::Matrix<casacore::Bool> scratch_mask(shape, casacore::True);
          for (int y=0; y<shape[1]; ++y) {
            for (int x=0; x<shape[0]; ++x) {
              if (real(scratch(x,y)) == 0.0) {
                scratch_mask(x,y) = casacore::False;
              }
            }
          }
          ASKAPLOG_INFO_STR(logger, "Applying Gaussian taper to the preconditioner function in the image domain");

          fft2d(scratch, false);

          casacore::Matrix<float> taperArray(itsTaperCache->taper(shape));
          scratch *= taperArray;

          fft2d(scratch, true);

          // Redo thresholding
          for (int y=0; y<shape[1]; ++y) {
            for (int x=0; x<shape[0]; ++x) {
              if (!scratch_mask(x,y)) {
                scratch(x,y)=0.0;
              }
            }
          }
          // Recalculate ave SNR weight-sum
          // const casacore::Array<float> wgts(real(scratch.asArray()));
          casacore::Array<double> wgts(shape);
          casacore::convertArray<double,float>(wgts, real(scratch));
          aveWgtSum = sum(wgts*wgts) / sum(wgts);
        }

        // Calculate the Wiener filter
        if (newFilter) {
          const float noisePower = itsUseRobustness ?
              (1.0/aveWgtSum)*25.0*std::pow(10., -2.0*itsParameter) : itsParameter;
          ASKAPLOG_INFO_STR(logger, "Effective noise power of the Wiener filter = " << noisePower);
          // The filter should turn itself off as wgt->0, however we don't trust those
          // regions and will explicitly set them to zero.
          itsWienerfilter.set(0.0);
          for (int y=0; y<shape[1]; ++y) {
            for (int x=0; x<shape[0]; ++x) {
              if (real(scratch(x,y)) != 0) {
                itsWienerfilter(x,y)=(1.0 / (noisePower*real(scratch(x,y)) + 1.0));
              }
            }
          }
        } else {
          const float normFactor = itsDoNormalise ? maxPSFBefore : 1.;
          const float noisePower = (itsUseRobustness ?
              std::pow(10., 4.*itsParameter) : itsParameter)*normFactor*normFactor;
          ASKAPLOG_INFO_STR(logger, "Effective noise power of the Wiener filter = " << noisePower);
          itsWienerfilter = conj(scratch)/(real(scratch*conj(scratch)) + noisePower);
          if (itsDoNormalise) itsWienerfilter *= normFactor;
        }

      }

      // for debugging - to export Wiener filter
      //scimath::saveAsCasaImage("wienerfilter.uv",real(itsWienerfilter.asArray()));
      //throw 1;

      // Apply the Wiener filter to the xfr and transform to the filtered PSF
      convertArray(scratch, psf2D);

      fft2d(scratch, true);

      // Estimate the normalised point source sensitivity (AKA the normalised thermal RMS)
      // * see D. Briggs thesis
      if (!useCachedFilter) {
          casacore::Matrix<Float> realWienerArr = real(itsWienerfilter);
          casacore::Matrix<Float> realScrArr = real(scratch);
          casacore::Float sumScr = sum(realScrArr);
          realScrArr *= realWienerArr;
          ASKAPLOG_INFO_STR(logger,
              "Normalisation factor for the point source sensitivity " <<
              "(theorectical image RMS / naturally weighted RMS): " <<
              sqrt( sum( realScrArr * realWienerArr) * sumScr ) / sum(realScrArr) );
      }

      scratch *= itsWienerfilter;

      fft2d(scratch, false);

      psf2D = real(scratch);
      const float maxPSFAfter=casacore::max(psf2D);
      ASKAPLOG_INFO_STR(logger,
          "Peak of PSF after Wiener filtering = " << maxPSFAfter <<
          " = " << 100.0*maxPSFAfter/maxPSFBefore << "%");
      psf2D *= maxPSFBefore/maxPSFAfter;
      ASKAPLOG_INFO_STR(logger, "Normalized to unit peak");

      // Apply the filter to the dirty image
      convertArray<casacore::Complex,casacore::Float>(scratch, dirty2D);

      fft2d(scratch, true);

      scratch *= itsWienerfilter;

      fft2d(scratch, false);

      dirty2D = real(scratch);
      dirty2D *= maxPSFBefore/maxPSFAfter;

      return true;

    }

    /// @brief calculate a threshold to use to clean up the PCF
    /// @details The preconditioner function can have a rumble of low-level
    /// Fourier components due to various issues (sharp edges in the image,
    /// tapering, de-warping). These need to be dealt with before forming the
    /// Wiener filter. One approach to this is thresholding out small values.
    /// @param[in] input array from which to determine threshold
    /// @return threshold value
    float WienerPreconditioner::pcfThreshold(casacore::Matrix<casacore::Complex>& pcf) const
    {

      float xfrMax = max(abs(real(pcf)));

      // The default threshold. Will be overwritten if useAutoThreshold==true.
      float scratchThreshold = 1e-6 * xfrMax;

      // Auto thresholding in this way will not work with Nyquist gridding.
      // Upcoming changes to the PCF should remove the necessity of the PCF.
      ASKAPLOG_INFO_STR(logger, "Bypassing statistical thresholding of the uv sampling function and instead using " <<
                        scratchThreshold << " (" << 100.0*scratchThreshold/xfrMax << "% of max) instead");
      return scratchThreshold;

      // set a test value that we expect the noise to be below
      float thresholdTestVal = 1e-2 * xfrMax;

      // Take slices across the zero-padding-region of the image
      size_t nx = pcf.shape()[0];
      size_t ny = pcf.shape()[1];
      // horizontal slice:
      float thresh1 = max(abs(real(pcf(Slice(nx/4,nx/2),Slice(0,ny/8)))));
      // vertical slice:
      float thresh2 = max(abs(real(pcf(Slice(0,nx/8),Slice(ny/4,ny/2)))));

      float autoThreshold = 3.0 * max(thresh1,thresh2);

      // test that the threshold is at a low level. If this fails, then
      // it is likely that visibilities have been gridded to the edge of
      // the uv plane. In this case, just set to something sensible.
      // It would be better to test the maximum baseline length either here
      // or earlier (and set a flag).
      if (autoThreshold > thresholdTestVal) {
        ASKAPLOG_INFO_STR(logger, "Auto-threshold seems too high. Using default.");
      } else {
        scratchThreshold = autoThreshold;
      }

      ASKAPLOG_INFO_STR(logger,
          "Thresholding the uv sampling function at " << scratchThreshold
          << " (" << 100.0*scratchThreshold/xfrMax << "% of max)");

      return scratchThreshold;

    }

    /// @brief static factory method to create preconditioner from a parset
    /// @details
    /// @param[in] parset subset of parset file (with preconditioner.Wiener. removed)
    /// @return shared pointer
    boost::shared_ptr<WienerPreconditioner>
        WienerPreconditioner::createPreconditioner(const LOFAR::ParameterSet &parset,
                                                   const bool useCachedPcf)
    {
      ASKAPCHECK(parset.isDefined("noisepower") != parset.isDefined("robustness"),
           "Exactly one parameter, either noisepower or robustness parameter must be given. You gave either none or both of them.");

      boost::shared_ptr<WienerPreconditioner> result;
      if (parset.isDefined("noisepower")) {
          const float noisepower = parset.getFloat("noisepower");
          const bool normalise = parset.getBool("normalise",false);
          result.reset(new WienerPreconditioner(noisepower,normalise));
      } else {

          ASKAPDEBUGASSERT(parset.isDefined("robustness"));

          const float robustness = parset.getFloat("robustness");
          ASKAPCHECK((robustness >= -2.00001) && (robustness <= 2.0001),
                     "Robustness parameter is supposed to be between -2 and 2, you have = "<<robustness);
          ASKAPCHECK(!parset.isDefined("normalise"),
                     "Normalise option of the Wiener preconditioner is not compatible with the "
                     "preconditioner definition via robustness (as normalisation of PSF is always done in this case)");
          result.reset(new WienerPreconditioner(robustness));
      }
      ASKAPASSERT(result);
      // configure tapering
      if (parset.isDefined("taper")) {
          const double fwhm = parset.getDouble("taper");
          result->enableTapering(fwhm);
      }
      //

      if (itsPcf.shape() == 0) {
          // there is no cache, so don't change anything.
          itsUseCachedPcf = false;
      } else if (!itsUseCachedPcf) {
          // weren't using the cache before, so don't start now.
          itsUseCachedPcf = false;
      } else {
          itsUseCachedPcf = useCachedPcf;
      }
      return result;
    }

    /// @brief assignment operator, to ensure it is not called
    WienerPreconditioner& WienerPreconditioner::operator=(const WienerPreconditioner &)
    {
      ASKAPTHROW(AskapError, "Assignment operator is not supposed to be used");
      return *this;
    }

    /// @brief enable Filter tapering
    /// @details Wiener filter can optionally be tapered in the image domain, so it is not extended over
    /// the whole field of view.
    /// @param[in] fwhm full width at half maximum of the taper given in image cells
    void WienerPreconditioner::enableTapering(double fwhm)
    {
      ASKAPLOG_INFO_STR(logger, "Wiener filter will be tapered by a circular Gaussian with FWHM="<<fwhm<<
                        " pixels in the image plane");
      itsTaperCache.reset(new GaussianTaperCache(fwhm));
    }

  }
}
