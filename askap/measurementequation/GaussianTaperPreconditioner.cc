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

#include <askap/gridding/SupportSearcher.h>
#include <askap/scimath/utils/PaddingUtils.h>

#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#include <casacore/lattices/LEL/LatticeExpr.h>
#include <casacore/lattices/LatticeMath/Fit2D.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/scimath/Mathematics/SquareMatrix.h>
#include <casacore/scimath/Mathematics/RigidVector.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>

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
/// @param[in] tolerance, the fractional tolerance in fitted beam size when using isPsfSize=true
/// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
/// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is
/// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
/// cellsize and FWHM are image-plane cell size and FWHM in angular units.
GaussianTaperPreconditioner::GaussianTaperPreconditioner(double majFWHM,
    double minFWHM, double pa, bool isPsfSize, double cutoff, double tolerance) :
     GaussianTaperCache(majFWHM, minFWHM, pa),itsFitBeam(isPsfSize),itsCutoff(cutoff),
     itsTolerance(tolerance) {}

/// @brief set up the preconditioner for the circularly symmetric taper
/// @details This constructor just sets the taper size, same for both axis.
/// The size is full width at half maximum expressed in the units of uv-cell.
/// @param[in] fwhm full width at half maximum of the taper in the uv-plane
/// (given as a fraction of the uv-cell size).
/// @param[in] cutoff, the cutoff used to determine the support area for beam fitting
/// @param[in] tolerance, the fractional tolerance in fitted beam size when using isPsfSize=true
/// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
/// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is
/// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
/// cellsize and FWHM are image-plane cell size and FWHM in angular units.
GaussianTaperPreconditioner::GaussianTaperPreconditioner(double fwhm, bool isPsfSize,
    double cutoff, double tolerance) :
     GaussianTaperCache(fwhm),itsFitBeam(isPsfSize),itsCutoff(cutoff),itsTolerance(tolerance) { }

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

  casacore::Array<float> origPsf;
  if (itsFitBeam) {
       origPsf = psf;
       ASKAPLOG_DEBUG_STR(logger, "Tuning taper to achieve final resolution");
  }
  applyTaper(psf);

  bool converged = !itsFitBeam;
  int count = 20;
  while (!converged && count-- >0) {
      ASKAPLOG_DEBUG_STR(logger, "Taper tuning iteration: "<<20-count);
      casacore::Vector<double> beam = fitPsf(psf);
      float tolerance = itsTolerance;
      converged = tuneTaper(beam,tolerance);
      if (!converged) {
          psf = origPsf;
          applyTaper(psf);
      }
  }
  if (itsFitBeam && !converged) {
      ASKAPLOG_WARN_STR(logger,"Failed to achieve requested beam size");
  }
  itsFitBeam = false; // Do we need this?

  const float maxPSFAfter = casacore::max(psf);
  ASKAPLOG_INFO_STR(logger, "Peak of PSF after Gaussian taper  = " << maxPSFAfter);
  ASKAPCHECK(maxPSFAfter>0., "Peak of the PSF is supposed to be a positive number");

  psf*=maxPSFBefore/maxPSFAfter;
  ASKAPLOG_INFO_STR(logger, "Normalized to unit peak");

  applyTaper(dirty);
  dirty*=maxPSFBefore/maxPSFAfter;

  return true;
}

/// @brief fit the beam using similar code used to get the beam fit in the imager (SynthesisParamsHelper::fitBeam)
casacore::Vector<double> GaussianTaperPreconditioner::fitPsf(casacore::Array<float>& psfArray) const {
    ASKAPTRACE("GaussianTaperPreconditioner::fitPsf");
    ASKAPCHECK(psfArray.shape().nelements()>=2,"PSF image is supposed to be at least 2-dimensional, shape="
    <<psfArray.shape());
    casacore::Matrix<float> psf = imagemath::MultiDimArrayPlaneIter::getFirstPlane(psfArray).nonDegenerate();
    SupportSearcher ss(itsCutoff);
    ss.search(psf);
    // search only looks in x and y direction, now extend the support to diagonals
    ss.extendedSupport(psf);
    casacore::uInt support = ss.symmetricalSupport(psf.shape());
    support = 2 * (support/2) + 1; // if even, move up to next odd.

    ASKAPLOG_DEBUG_STR(logger, "Extracting support of "<<support
      <<" pixels for 2D gaussian fitting");
    const casacore::IPosition newShape(2,support,support);
    for (int dim=0; dim<2; ++dim) {
         ASKAPCHECK(psf.shape()[dim] >= int(support), "Support is greater than the original size, shape="<<
                    psf.shape());
    }
    // avoid making a reference copy since we rescale in next step
    casacore::Matrix<float> psfSub;
    psfSub = scimath::PaddingUtils::centeredSubArray(psf,newShape);

    // normalise to 1
    const float maxPSF = casacore::max(psfSub);
    if (fabs(maxPSF-1.)>1e-6) {
        psfSub /= maxPSF;
    }

    // actual fitting
    casacore::LogIO os;
    casacore::Fit2D fitter(os);
    casacore::Vector<casacore::Double> initialEstimate =  fitter.estimate(casacore::Fit2D::GAUSSIAN,psfSub);
    initialEstimate[0]=1.; // PSF peak is always 1
    initialEstimate[1]=psfSub.shape()[0]/2; // centre
    initialEstimate[2]=psfSub.shape()[1]/2; // centre
    casacore::Vector<casacore::Bool> parameterMask(6,casacore::False);
    parameterMask[3] = casacore::True; // fit maj
    parameterMask[4] = casacore::True; // fit min
    parameterMask[5] = casacore::True; // fit pa

    ASKAPLOG_DEBUG_STR(logger,"initial estimate: "<<initialEstimate );
    fitter.addModel(casacore::Fit2D::GAUSSIAN,initialEstimate,parameterMask);
    casacore::Array<casacore::Float> sigma(psfSub.shape(),1.);
    const casacore::Fit2D::ErrorTypes fitError = fitter.fit(psfSub,sigma);
    ASKAPCHECK(fitError == casacore::Fit2D::OK, "Error fitting the beam. fitError="<<fitError<<
               " message: "<<fitter.errorMessage());
    casacore::Vector<casacore::Double> result = fitter.availableSolution();
    ASKAPLOG_DEBUG_STR(logger, "Got fit result (in pixels) "<<result<<" and uncertainties "<<fitter.availableErrors());
    ASKAPCHECK(result.nelements() == 6, "Expect 6 parameters for 2D gaussian, result vector has "<<result.nelements());
    casacore::Vector<double> beam(3);
    beam[0] = result[3];
    beam[1] = result[4];
    // position angle in radians
    //double pa = increments[0]<0 ? result[5] - casacore::C::pi/2 : casacore::C::pi/2 - result[5];
    // Note: Code below assumes square RA/Dec grid, as we have no access to the axes increments here
    double pa = result[5] - casacore::C::pi/2;
    if (pa < -casacore::C::pi/2) {
        pa += casacore::C::pi;
    }
    beam[2] = pa ; //casacore::Quantum<double>(pa,"rad");
    return beam;
}
/// @brief a helper method to apply the taper to one given array
/// @details We need exactly the same operation for psf and dirty image. This method
/// encapsulates the code which is actually doing the job. It is called twice from
/// doPreconditioning.
/// @param[in] image an image to apply the taper to
void GaussianTaperPreconditioner::applyTaper(casacore::Array<float> &image) const
{
  casacore::ArrayLattice<float> lattice(image);

  /*
  casacore::IPosition paddedShape = image.shape();
  ASKAPDEBUGASSERT(paddedShape.nelements()>=2);
  paddedShape[0] *= 2;
  paddedShape[1] *= 2;
  */

  // Setup work arrays.
  const casacore::IPosition shape = lattice.shape();
  //const casacore::IPosition shape = paddedShape;
  casacore::ArrayLattice<casacore::Complex> scratch(shape);

  // fft to transform the image into uv-domain
  scratch.copyData(casacore::LatticeExpr<casacore::Complex>(toComplex(lattice)));
  casacore::LatticeFFT::cfft2d(scratch, true);

  // apply the taper
  casacore::ArrayLattice<float> taperLattice(taper(shape));

  scratch.copyData(casacore::LatticeExpr<casacore::Complex> (taperLattice * scratch));

  // transform back to the image domain
  casacore::LatticeFFT::cfft2d(scratch, false);
  lattice.copyData(casacore::LatticeExpr<float> ( real(scratch) ));
}

} // namespace synthesis

} // namespace askap
