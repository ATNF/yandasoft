/// @file
///
/// @brief apply gaussian taper
/// @details This preconditioner applies gaussian taper in the uv-space to the normal
/// equations.
///
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

#include <askap/measurementequation/GaussianTaperPreconditioner.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.gaussiantaperpreconditioner");

#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>

#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>

namespace askap {

namespace synthesis {


/// @brief set up the preconditioner
/// @details This constructor just sets the taper size. The size is full width at
/// half maximum expressed in the units of uv-cell.
/// @param[in] majFWHM full width at half maximum of the major axis in the uv-plane
/// (given as a fraction of the uv-cell size).
/// @param[in] minFWHM full width at half maximum of the minor axis in the uv-plane
/// (given as a fraction of the uv-cell size).
/// @param[in] pa position angle in radians
/// @param[in] cutoff, the cutoff used to determine the support area for beam fitting
/// @param[in] maxsupport Max support size of beam above cutoff level
/// @param[in] tolerance, the fractional tolerance in fitted beam size when using isPsfSize=true
/// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
/// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is
/// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
/// cellsize and FWHM are image-plane cell size and FWHM in angular units.
GaussianTaperPreconditioner::GaussianTaperPreconditioner(double majFWHM,
    double minFWHM, double pa, bool isPsfSize, double cutoff, int maxsupport,
    double tolerance) :
    GaussianTaperCache(majFWHM, minFWHM, pa),itsFitBeam(isPsfSize),itsCutoff(cutoff),
    itsMaxSupport(maxsupport), itsTolerance(tolerance)
{}

/// @brief set up the preconditioner for the circularly symmetric taper
/// @details This constructor just sets the taper size, same for both axis.
/// The size is full width at half maximum expressed in the units of uv-cell.
/// @param[in] fwhm full width at half maximum of the taper in the uv-plane
/// (given as a fraction of the uv-cell size).
/// @param[in] cutoff, the cutoff used to determine the support area for beam fitting
/// @param[in] maxsupport Max support size of beam above cutoff level
/// @param[in] tolerance, the fractional tolerance in fitted beam size when using isPsfSize=true
/// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
/// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is
/// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
/// cellsize and FWHM are image-plane cell size and FWHM in angular units.
GaussianTaperPreconditioner::GaussianTaperPreconditioner(double fwhm, bool isPsfSize,
    double cutoff, int maxsupport, double tolerance) :
    GaussianTaperCache(fwhm),itsFitBeam(isPsfSize),itsCutoff(cutoff),itsMaxSupport(maxsupport),
    itsTolerance(tolerance)
{}

/// @brief Clone this object
/// @return shared pointer to a cloned copy
IImagePreconditioner::ShPtr GaussianTaperPreconditioner::clone()
{
  return IImagePreconditioner::ShPtr(new GaussianTaperPreconditioner(*this));
}

/// @brief Apply preconditioning to Image Arrays
/// @details This is the actual method, which does preconditioning.
/// It is applied to the PSF as well as the current residual image.
/// @param[in] psf array with PSF
/// @param[in] dirty array with dirty image
/// @return true if psf and dirty have been altered
bool GaussianTaperPreconditioner::doPreconditioning(casacore::Array<float>& psf,
                                                    casacore::Array<float>& dirty,
                                                    casacore::Array<float>& pcf) const
{
  ASKAPTRACE("GaussianTaperPreconditioner::doPreconditioning");

  const float maxPSFBefore=casacore::max(psf);
  ASKAPLOG_INFO_STR(logger, "Peak of PSF before Gaussian taper = " << maxPSFBefore);

  ASKAPLOG_INFO_STR(logger, "Applying Gaussian taper "<<majorAxis()*sqrt(8.*log(2.))<<" x "<<
                    minorAxis()*sqrt(8.*log(2.))<<" uv cells at the position angle of "<<posAngle()/M_PI*180.<<" degrees");
  ASKAPDEBUGASSERT(psf.shape().isEqual(dirty.shape()));

  if (!itsFitBeam) {
      applyTaper(psf);
  } else {
      ASKAPLOG_DEBUG_STR(logger, "Tuning taper to achieve final resolution");
      casacore::Array<float> origPsf;
      // just to initialise
      taper(psf.shape());
      origPsf = psf;
      GaussianTaperCache::TuningStatus status = GaussianTaperCache::NOTCONVERGED;
      int count = 0;

      while (status == NOTCONVERGED && count++ < 20) {
          ASKAPLOG_DEBUG_STR(logger, "Taper tuning iteration: "<<count);
          casacore::Vector<double> beam = fitPsf(psf);
          float tolerance = itsTolerance;
          status = tuneTaper(beam, tolerance, count);
          if (status == GaussianTaperCache::NOTCONVERGED) {
              psf = origPsf;
              applyTaper(psf);
          }
      }
      if (status != GaussianTaperCache::CONVERGED) {
          ASKAPLOG_WARN_STR(logger,"Failed to achieve requested beam size");
      }
      itsFitBeam = false; // Do we need this? Yes - avoids refitting in restore phase
  }

  const float maxPSFAfter = casacore::max(psf);
  ASKAPLOG_INFO_STR(logger, "Peak of PSF after Gaussian taper  = " << maxPSFAfter);
  ASKAPCHECK(maxPSFAfter>0., "Peak of the PSF is supposed to be a positive number");

  psf*=maxPSFBefore/maxPSFAfter;
  ASKAPLOG_INFO_STR(logger, "Normalized to unit peak");

  applyTaper(dirty);
  dirty*=maxPSFBefore/maxPSFAfter;

  return true;
}

casacore::Vector<double> GaussianTaperPreconditioner::fitPsf(casacore::Array<float>& psfArray) const {
    ASKAPTRACE("GaussianTaperPreconditioner::fitPsf");
    #ifdef ASKAP_FLOAT_IMAGE_PARAMS
    return SynthesisParamsHelper::fitBeam(psfArray,itsCutoff,itsMaxSupport);
    #else
    casa::Array<double> psfDArray(psfArray.shape());
    casa::convertArray<double, float>(psfDArray, psfArray);
    return SynthesisParamsHelper::fitBeam(psfDArray,itsCutoff,itsMaxSupport);
    #endif
}
/// @brief a helper method to apply the taper to one given array
/// @details We need exactly the same operation for psf and dirty image. This method
/// encapsulates the code which is actually doing the job. It is called twice from
/// doPreconditioning.
/// @param[in] image an image to apply the taper to
void GaussianTaperPreconditioner::applyTaper(casacore::Array<float> &image) const
{
  // This will throw exception if there are !=2 non degenerate axes
  casacore::Matrix<float> matImage(image.nonDegenerate());
  // Limit number of fft threads to 8 (more is slower for our fft sizes)
  scimath::FFT2DWrapper<casacore::Complex> fft2d(true,8);
  casacore::Matrix<casacore::Complex> grid(matImage.shape(),casacore::Complex(0.f));

  // copy image into a complex scratch space
  casacore::setReal(grid, matImage);

  // fft to uv
  fft2d(grid, true);

  // multiply by taper
  grid *= taper(matImage.shape());

  // ifft back to image and return the real part
  fft2d(grid, false);
  casacore::real(matImage, grid);
}

} // namespace synthesis

} // namespace askap
