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

#include <measurementequation/WienerPreconditioner.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.wienerpreconditioner");

#include <askap/AskapError.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <utils/PaddingUtils.h>
#include <profile/AskapProfiler.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/lattices/Lattices/SubLattice.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#include <casacore/lattices/LEL/LatticeExpr.h>

// for debugging - to export intermediate images
#include <utils/ImageUtils.h>

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
    casa::Array<float> WienerPreconditioner::itsPcf;

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
    
    bool WienerPreconditioner::doPreconditioning(casa::Array<float>& psf,
                                                 casa::Array<float>& dirty,
                                                 casa::Array<float>& pcfIn) const
    {
      ASKAPTRACE("WienerPreconditioner::doPreconditioning");
      if (!itsUseRobustness && (itsParameter < 1e-6)) {
          return false;
      }

      ASKAPCHECK(psf.shape().conform(dirty.shape()),
          "Dirty image and PSF do not conform - shapes: " <<
          dirty.shape() << " & " << psf.shape());

      const casa::IPosition shape = psf.shape();

      bool useCachedFilter = true;
      if (itsWienerfilter.shape() == 0) {
        useCachedFilter = false;
        itsWienerfilter = casa::ArrayLattice<casa::Complex>(shape);
      }

      bool newFilter = true;
      //casa::Array<float> pcf;
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
          itsPcf.reference(psf);
          // could also use the stored version, but leave for now.
          useCachedFilter = false;
          newFilter = false;
      } else {
          ASKAPLOG_INFO_STR(logger, "Applying Wiener filter");
          //itsPcf.reference(pcfIn);
          itsPcf = pcfIn.copy();
      }
      ASKAPCHECK(itsPcf.shape().conform(psf.shape()),
          "PSF and preconditioner function do not conform - shapes: " <<
          psf.shape() << " & " << itsPcf.shape());

      if (itsUseRobustness) {
          ASKAPLOG_INFO_STR(logger,
              "Wiener filter noise power defined via robustness = " << itsParameter);
      } else {
          ASKAPLOG_INFO_STR(logger,
              "Wiener filter noise power = " << itsParameter);
      }

      IPosition pos(shape.nelements());
      for (uInt k=0; k<pos.nelements(); ++k) {
        pos(k) = 0;
      }

      float maxPSFBefore = casa::max(psf);
      ASKAPLOG_INFO_STR(logger, "Peak of PSF before Wiener filtering = " << maxPSFBefore);

      if (itsDoNormalise) {
          ASKAPLOG_INFO_STR(logger, "The PSF will be normalised to 1 before filter construction");
          psf=psf/maxPSFBefore;
          // dirty=dirty/maxPSFBefore;
          maxPSFBefore=1.0;
      }

      casa::ArrayLattice<float> lpsf(psf);
      casa::ArrayLattice<float> ldirty(dirty);
      casa::ArrayLattice<float> lpcf(itsPcf);

      // Make the scratch array into which we will calculate the Wiener filter
      casa::ArrayLattice<casa::Complex> scratch(shape);
      scratch.copyData(casa::LatticeExpr<casa::Complex>(toComplex(lpcf)));

      if (!itsUseCachedPcf) {

        // the filter has been accumulated in the image domain, so transform to uv
        casa::LatticeFFT::cfft2d(scratch, casa::True);

        // The filter is rescaling Fourier components of dirty and psf based on
        // the value of non-zero Fourier components of pcf. So set small
        // erroneous components to zero.
        float scratchThreshold = 1e-6;
    
        // Reset the threshold based on data in the zero-padding region of the uv
        // plane. Should test for the presence of zero-padding...
        bool autoThreshold = true;
        if (autoThreshold) {
    
          // Take slices across the zero-padding-region of the image
          IPosition sliceStart(scratch.shape().nelements());
          IPosition sliceShape(scratch.shape().nelements());
          for (uInt k=0; k<sliceStart.nelements(); ++k) {
            sliceStart(k) = 0;
            sliceShape(k) = 1;
          }
          // horizontal slice:
          sliceStart(0) = scratch.shape()[0] / 4;
          sliceShape(0) = scratch.shape()[0] / 2;
          sliceStart(1) = 0;
          sliceShape(1) = scratch.shape()[1] / 8;
          scratchThreshold = max(scratchThreshold,
              3.0 * max(abs(real(scratch.getSlice(sliceStart,sliceShape)))));
          // vertical slice:
          sliceStart(0) = 0;
          sliceShape(0) = scratch.shape()[0] / 8;
          sliceStart(1) = scratch.shape()[1] / 4;
          sliceShape(1) = scratch.shape()[1] / 2;
          scratchThreshold = max(scratchThreshold,
              3.0 * max(abs(real(scratch.getSlice(sliceStart,sliceShape)))));
    
          ASKAPLOG_INFO_STR(logger,
              "Thresholding the input uv sampling function at " <<
              scratchThreshold << " (" <<
              100.0*scratchThreshold/max(abs(real(scratch.asArray()))) <<
              "% of max)");
    
          for (int x=0; x<scratch.shape()[0]; ++x) {
            for (int y=0; y<scratch.shape()[1]; ++y) {
              pos(0) = x;
              pos(1) = y;
              if (abs(real(scratch.getAt(pos))) < scratchThreshold) {
                scratch.putAt(0.0, pos);
              }
            }
          }

        } // if (autoThreshold)
 
        if (newFilter) {
       
          // The PCF should be set up to have nearest-neighbour gridding of
          // SNR weights in the real part and nearest-neighbour gridding of
          // weights*kernel-size (in pixels) in the imaginary part. From this
          // the idea is to build a filter that alters the local density of
          // weighted data without breaking local projections.
       
          // get weighted sum of kernel size from imaginary part of the data
          casa::ArrayLattice<casa::Float> kernelWidthArray(shape);
          kernelWidthArray.set(0.0);
          for (int x=0; x<scratch.shape()[0]; ++x) {
            for (int y=0; y<scratch.shape()[1]; ++y) {
              pos(0) = x;
              pos(1) = y;
              if (real(scratch.getAt(pos)) > scratchThreshold) {
                kernelWidthArray.putAt(abs(imag(scratch.getAt(pos))/real(scratch.getAt(pos))), pos);
              }
            }
          }
       
          // for debugging - to export kernel widths
          //scimath::saveAsCasaImage("kernelwidths.uv",kernelWidthArray.asArray());
          //throw 1;
       
          // Build the filter as follows:
          // * set a large box (maximum kernel size) around each pixel
          // * find largest local kernel size
          // * make non-zero if there is data within this radius
          // * estimate weight based on the average of near by points,
          //   where "near by" can be different to the kernel width.
       
          const int maxKernelWidth = ceil(max(kernelWidthArray.asArray()));
          ASKAPDEBUGASSERT(maxKernelWidth>0);
          const int boxWidth = maxKernelWidth;
          ASKAPDEBUGASSERT(boxWidth>0);
       
          // the factor for running mean (relative to the local kernel width)
          const int extra = 2;
       
          // Initialise a box.
          // Are boxes faster than simply searching scratch? Test.
          IPosition boxStart(scratch.shape().nelements());
          IPosition boxShape(scratch.shape().nelements());
          for (uInt k=0; k<boxStart.nelements(); ++k) {
            boxStart(k) = 0;
            boxShape(k) = 1;
          }
       
          // As with the threshold search above, here we assume that 
          // zero-padding is used to over sample the PSF, meaning that we
          // can be sure that visibilities are not gridded to the edge of
          // the uv plane. Should test for the presence of zero-padding...

          // set pcf to to contain the new weights, so they don't have to be regenerated.
          // note that it is going from the image domain to the uv domain.
          itsPcf = 0.0;
       
          for (int x=extra*boxWidth/2; x<scratch.shape()[0]-extra*boxWidth/2; ++x) {
            for (int y=extra*boxWidth/2; y<scratch.shape()[1]-extra*boxWidth/2; ++y) {
       
              boxStart(0) = x - boxWidth/2;
              boxStart(1) = y - boxWidth/2;
              boxShape(0) = boxWidth;
              boxShape(1) = boxWidth;
       
              int localCount = 0;
              double localSum = 0.0;
              int regionCount = 0;
              double regionSum = 0.0;
       
              const int kernelWidth = ceil(max(kernelWidthArray.getSlice(boxStart,boxShape)));
       
              if (kernelWidth>0) {
       
                // reset box to the kernelWidth
                const int regionWidth = 1 + extra*(kernelWidth-1);
                ASKAPDEBUGASSERT(regionWidth>=kernelWidth);
                boxStart(0) = x - regionWidth/2;
                boxStart(1) = y - regionWidth/2;
                boxShape(0) = regionWidth;
                boxShape(1) = regionWidth;
                const int x0 = x - boxStart(0);
                const int y0 = y - boxStart(1);
                const casa::Array<casa::Float> localBox = real(scratch.getSlice(boxStart,boxShape));
       
                const float localRadiusSq = 0.25 * kernelWidth*kernelWidth;
                const float regionRadiusSq = 0.25 * regionWidth*regionWidth;
       
                // check max(localBox) > 0 before continuing?
                // there are efficient ways of doing a running mean, but with a varying box size?
                for (int xb=0; xb<boxShape(0); ++xb) {
                  pos(0) = xb;
                  const int dx = xb - x0;
                  for (int yb=0; yb<boxShape(1); ++yb) {
                    pos(1) = yb;
                    const int dy = yb - y0;
                    if (localBox(pos) > scratchThreshold) {
                      const float rsq = dx*dx + dy*dy;
                      if (rsq<=regionRadiusSq) {
                        regionCount += 1;
                        regionSum += localBox(pos);
                        if (rsq<=localRadiusSq) {
                          localCount += 1;
                          localSum += localBox(pos);
                        }
                      }
                    }
                  }
                }
       
              }
       
              pos(0) = x;
              pos(1) = y;
       
              if (localCount > 0) {
                itsPcf(pos) = regionSum/double(regionCount);
              }
       
            } // y
          } // x
      
          // Calc ave SNR weight-sum *over visibilities* (not pixels).
          // * do this before the taper?
          // * using the pcf or the psf? The psf may be normalised, so pcf.
          const casa::Array<float> wgts(real(scratch.asArray()));
          // The following approx assumes that vis have equal SNR weights.
          // * see D. Briggs thesis
          itsAveWgtSum = sum(wgts*wgts) / sum(wgts);
 
          // update the scratch array
          scratch.copyData(casa::LatticeExpr<casa::Complex>(toComplex(lpcf)));

          itsUseCachedPcf = true;
 
        } // if (newFilter)

      } // if (!itsUseCachedPcf)

      if (!useCachedFilter) {

        // Apply Gaussian taper to the preconditioner function (in the image domain)
        // this needs to happen after the weights are set up
        if (itsTaperCache) {
          ASKAPLOG_INFO_STR(logger, "Applying Gaussian taper to the preconditioner function in the image domain");
          casa::LatticeFFT::cfft2d(scratch, casa::False);
          casa::Array<casa::Complex> taperArray(itsTaperCache->taper(shape));
          casa::ArrayLattice<casa::Complex> taperLattice(taperArray);
          scratch.copyData(casa::LatticeExpr<casa::Complex>(scratch * taperLattice));
          casa::LatticeFFT::cfft2d(scratch, casa::True);
        }

        // Calculate the Wiener filter
        if (newFilter) {
          const float noisePower = itsUseRobustness ?
              (1.0/itsAveWgtSum)*25.0*std::pow(10., -2.0*itsParameter) : itsParameter;
          ASKAPLOG_INFO_STR(logger, "Effective noise power of the Wiener filter = " << noisePower);
          // the filter should turn itself off as wgt->0, however we don't trust those
          // regions and will explicitly set them to zero.
          itsWienerfilter.set(0.0);
          for (int x=0; x<scratch.shape()[0]; ++x) {
            for (int y=0; y<scratch.shape()[1]; ++y) {
              pos(0) = x;
              pos(1) = y;
              if (real(scratch.getAt(pos)) != 0) {
                itsWienerfilter.putAt(1.0 / (noisePower*real(scratch.getAt(pos)) + 1.0), pos);
              }
            }
          }
        } else {
          const float normFactor = itsDoNormalise ? maxPSFBefore : 1.;
          const float noisePower = (itsUseRobustness ?
              std::pow(10., 4.*itsParameter) : itsParameter)*normFactor*normFactor;
          ASKAPLOG_INFO_STR(logger, "Effective noise power of the Wiener filter = " << noisePower);     
          itsWienerfilter.copyData(casa::LatticeExpr<casa::Complex>(
              normFactor*conj(scratch)/(real(scratch*conj(scratch)) + noisePower)));
        }

      }

      // for debugging - to export Wiener filter
      //scimath::saveAsCasaImage("wienerfilter.uv",real(itsWienerfilter.asArray()));
      //throw 1;

      // Apply the Wiener filter to the xfr and transform to the filtered PSF
      scratch.copyData(casa::LatticeExpr<casa::Complex>(toComplex(lpsf)));
      casa::LatticeFFT::cfft2d(scratch, casa::True);
      scratch.copyData(casa::LatticeExpr<casa::Complex> (itsWienerfilter * scratch));
      casa::LatticeFFT::cfft2d(scratch, casa::False);       
      lpsf.copyData(casa::LatticeExpr<float>(real(scratch)));
      const float maxPSFAfter=casa::max(psf);
      ASKAPLOG_INFO_STR(logger,
          "Peak of PSF after Wiener filtering = " << maxPSFAfter <<
          " = " << 100.0*maxPSFAfter/maxPSFBefore << "%");
      psf *= maxPSFBefore/maxPSFAfter;
      ASKAPLOG_INFO_STR(logger, "Normalized to unit peak");

      // Apply the filter to the dirty image
      scratch.copyData(casa::LatticeExpr<casa::Complex>(toComplex(ldirty)));       
      casa::LatticeFFT::cfft2d(scratch, casa::True);
      scratch.copyData(casa::LatticeExpr<casa::Complex> (itsWienerfilter * scratch));
      casa::LatticeFFT::cfft2d(scratch, casa::False);
      ldirty.copyData(casa::LatticeExpr<float>(real(scratch)));
      dirty *= maxPSFBefore/maxPSFAfter;

      return true;

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


