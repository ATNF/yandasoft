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

#include <measurementequation/GaussianTaperPreconditioner.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.gaussiantaperpreconditioner");

#include <askap/AskapError.h>
#include <profile/AskapProfiler.h>

#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/Lattices/LatticeFFT.h>
#include <lattices/Lattices/LatticeExpr.h>
#include <casa/BasicSL/Constants.h>
#include <casa/Arrays/ArrayMath.h>
#include <scimath/Mathematics/SquareMatrix.h>
#include <scimath/Mathematics/RigidVector.h>

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
/// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
/// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is 
/// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
/// cellsize and FWHM are image-plane cell size and FWHM in angular units.
GaussianTaperPreconditioner::GaussianTaperPreconditioner(double majFWHM, double minFWHM, double pa) :
     GaussianTaperCache(majFWHM, minFWHM, pa) {}
   
/// @brief set up the preconditioner for the circularly symmetric taper 
/// @details This constructor just sets the taper size, same for both axis.
/// The size is full width at half maximum expressed in the units of uv-cell.
/// @param[in] fwhm full width at half maximum of the taper in the uv-plane
/// (given as a fraction of the uv-cell size).
/// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
/// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is 
/// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
/// cellsize and FWHM are image-plane cell size and FWHM in angular units.
GaussianTaperPreconditioner::GaussianTaperPreconditioner(double fwhm) : 
     GaussianTaperCache(fwhm) { }
   
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
bool GaussianTaperPreconditioner::doPreconditioning(casa::Array<float>& psf, casa::Array<float>& dirty) const
{
  ASKAPTRACE("GaussianTaperPreconditioner::doPreconditioning");

  const float maxPSFBefore=casa::max(psf);
  ASKAPLOG_INFO_STR(logger, "Peak of PSF before Gaussian taper = " << maxPSFBefore);
  
  ASKAPLOG_INFO_STR(logger, "Applying Gaussian taper "<<majorAxis()*sqrt(8.*log(2.))<<" x "<<
                    minorAxis()*sqrt(8.*log(2.))<<" uv cells at the position angle of "<<posAngle()/M_PI*180.<<" degrees");
  ASKAPDEBUGASSERT(psf.shape().isEqual(dirty.shape()));
  
  applyTaper(psf);

  const float maxPSFAfter=casa::max(psf);

  ASKAPLOG_INFO_STR(logger, "Peak of PSF after Gaussian taper  = " << maxPSFAfter);
  ASKAPCHECK(maxPSFAfter>0., "Peak of the PSF is supposed to be a positive number");
 
  psf*=maxPSFBefore/maxPSFAfter;
  ASKAPLOG_INFO_STR(logger, "Normalized to unit peak");

  applyTaper(dirty);
  dirty*=maxPSFBefore/maxPSFAfter;

  return true;
}

/// @brief a helper method to apply the taper to one given array
/// @details We need exactly the same operation for psf and dirty image. This method
/// encapsulates the code which is actually doing the job. It is called twice from
/// doPreconditioning.
/// @param[in] image an image to apply the taper to
void GaussianTaperPreconditioner::applyTaper(casa::Array<float> &image) const
{
  casa::ArrayLattice<float> lattice(image);
  
  /*
  casa::IPosition paddedShape = image.shape();
  ASKAPDEBUGASSERT(paddedShape.nelements()>=2);
  paddedShape[0] *= 2;
  paddedShape[1] *= 2;
  */
  
  // Setup work arrays.
  const casa::IPosition shape = lattice.shape();
  //const casa::IPosition shape = paddedShape;
  casa::ArrayLattice<casa::Complex> scratch(shape);
      
  // fft to transform the image into uv-domain
  scratch.copyData(casa::LatticeExpr<casa::Complex>(toComplex(lattice)));
  casa::LatticeFFT::cfft2d(scratch, true);
  
  // apply the taper
  casa::Array<casa::Complex> taperCache = taper(shape); 
  casa::ArrayLattice<casa::Complex> taperLattice(taperCache);
  
  scratch.copyData(casa::LatticeExpr<casa::Complex> (taperLattice * scratch));
  
  // transform back to the image domain
  casa::LatticeFFT::cfft2d(scratch, false);
  lattice.copyData(casa::LatticeExpr<float> ( real(scratch) ));
}

} // namespace synthesis

} // namespace askap

