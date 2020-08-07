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

#include <askap/measurementequation/IImagePreconditioner.h>
#include <askap/measurementequation/GaussianTaperCache.h>

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
   /// @param[in] isPsfSize the specified taper size is the required fitted beam size of the output
   /// @param[in] cutoff, the cutoff used to determine the support area for beam fitting
   /// @param[in] tolerance, the fractional tolerance in fitted beam size when using isPsfSize=true
   /// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
   /// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is
   /// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
   /// cellsize and FWHM are image-plane cell size and FWHM in angular units.
   GaussianTaperPreconditioner(double majFWHM, double minFWHM, double pa, bool isPsfSize = false,
       double cutoff = 0.5, double tolerance = 0.005);

   /// @brief set up the preconditioner for the circularly symmetric taper
   /// @details This constructor just sets the taper size, same for both axis.
   /// The size is full width at half maximum expressed in the units of uv-cell.
   /// @param[in] fwhm full width at half maximum of the taper in the uv-plane
   /// (given as a fraction of the uv-cell size).
   /// @param[in] isPsfSize the specified taper size is the required fitted beam size of the output
   /// @param[in] cutoff, the cutoff used to determine the support area for beam fitting
   /// @param[in] tolerance, the fractional tolerance in fitted beam size when using isPsfSize=true
   /// @note Gaussian taper is set up in the uv-space. Constructors accept sizes given as FWHM expressed
   /// as fractions of uv-cell size. The relation between FWHMs in fourier and image plane is
   /// uvFWHM = (Npix*cellsize / FWHM) * (4*log(2)/pi), where Npix is the number of pixels
   /// cellsize and FWHM are image-plane cell size and FWHM in angular units.
   GaussianTaperPreconditioner(double fwhm, bool isPsfSize = false, double cutoff = 0.5,
        double tolerance = 0.005);

   /// @brief Clone this object
   /// @return shared pointer to a cloned copy
   virtual IImagePreconditioner::ShPtr clone();

   /// @brief Apply preconditioning to Image Arrays
   /// @details This is the actual method, which does preconditioning.
   /// It is applied to the PSF as well as the current residual image.
   /// @param[in] psf array with PSF
   /// @param[in] dirty array with dirty image
   /// @return true if psf and dirty have been altered
   virtual bool doPreconditioning(casacore::Array<float>& psf,
                                  casacore::Array<float>& dirty,
                                  casacore::Array<float>& pcf) const;

   /// @brief a helper method to apply the taper to one given array
   /// @details We need exactly the same operation for psf and dirty image. This method
   /// encapsulates the code which is actually doing the job. It is called twice from
   /// doPreconditioning. In addition, this method is used in Wiener preconditioner to
   /// smooth PSF before filter construction
   /// @param[in] image an image to apply the taper to
   void applyTaper(casacore::Array<float> &image) const;

   /// @ brief fit the psf with a gaussian and return size and pa
   casacore::Vector<double> fitPsf(casacore::Array<float>& psf) const;

private:
    mutable bool itsFitBeam;
    double itsCutoff;
    double itsTolerance;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef GAUSSIAN_TAPER_PRECONDITIONER_H
