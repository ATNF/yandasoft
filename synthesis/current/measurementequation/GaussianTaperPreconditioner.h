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

#ifndef GAUSSIAN_TAPER_PRECONDITIONER_H
#define GAUSSIAN_TAPER_PRECONDITIONER_H

#include <measurementequation/IImagePreconditioner.h>
#include <measurementequation/GaussianTaperCache.h>

namespace askap {

namespace synthesis {

/// @brief apply gaussian taper
/// @details This preconditioner applies gaussian taper in the uv-space to the normal
/// equations.
/// @ingroup measurementequation
class GaussianTaperPreconditioner : public IImagePreconditioner, private GaussianTaperCache {
public:
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
   GaussianTaperPreconditioner(double majFWHM, double minFWHM, double pa);
   
   /// @brief set up the preconditioner for the circularly symmetric taper 
   /// @details This constructor just sets the taper size, same for both axis.
   /// The size is full width at half maximum expressed in the units of uv-cell.
   /// @param[in] fwhm full width at half maximum of the taper in the uv-plane
   /// (given as a fraction of the uv-cell size).
   /// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
   /// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is 
   /// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
   /// cellsize and FWHM are image-plane cell size and FWHM in angular units.
   GaussianTaperPreconditioner(double fwhm);
   
   /// @brief Clone this object
   /// @return shared pointer to a cloned copy
   virtual IImagePreconditioner::ShPtr clone();
        
   /// @brief Apply preconditioning to Image Arrays
   /// @details This is the actual method, which does preconditioning.
   /// It is applied to the PSF as well as the current residual image.
   /// @param[in] psf array with PSF
   /// @param[in] dirty array with dirty image
   /// @return true if psf and dirty have been altered
   virtual bool doPreconditioning(casa::Array<float>& psf, casa::Array<float>& dirty) const;

   /// @brief a helper method to apply the taper to one given array
   /// @details We need exactly the same operation for psf and dirty image. This method
   /// encapsulates the code which is actually doing the job. It is called twice from
   /// doPreconditioning. In addition, this method is used in Wiener preconditioner to 
   /// smooth PSF before filter construction 
   /// @param[in] image an image to apply the taper to
   void applyTaper(casa::Array<float> &image) const;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef GAUSSIAN_TAPER_PRECONDITIONER_H

