/// @file
///
/// NormWienerPreconditioner: Precondition the normal equations
///                       by applying a Wiener filter
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
/// @author Urvashi Rau <rurvashi@aoc.nrao.edu>
///
#ifndef SYN_NORMWEINER_PRECONDITIONER_H_
#define SYN_NORMWEINER_PRECONDITIONER_H_

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <fitting/Axes.h>

#include <boost/shared_ptr.hpp>

#include <measurementequation/IImagePreconditioner.h>

namespace askap
{
  namespace synthesis
  {
    /// @brief Precondition the normal equations via a Wiener filter
    ///
    /// @details It constructs a Wiener filter from the PSF 
    /// and applies it to the PSF and current Residual image
    /// The robust parameter scaling is used to aid 
    /// use and comparison with robust weighting
    /// @ingroup measurementequation
    class NormWienerPreconditioner : public IImagePreconditioner
    {
      public:

        /// @brief default constructor with zero robust parameter
        /// @note do we really need it?
        NormWienerPreconditioner();
        
        /// @brief constructor with explicitly defined robust 
        /// @param[in] robust parameter of the fileter
        NormWienerPreconditioner(const float& robust);
        
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

      private:
  	    /// @brief Noise Power Spectrum
	    float itsRobust;
   };

  }
}
#endif
