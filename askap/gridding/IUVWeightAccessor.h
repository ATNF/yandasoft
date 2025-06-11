/// @file
/// @brief Interface class to access UVWeights
/// @details This interface is intended to be used inside the gridder to structure access 
/// to uv weights when they are applied (in traditional weighting). In the simplest form
/// it can be simply the weight collection. But other options may involve index translation
/// (e.g. if we want to apply the weight obtained for one beam to another or interpolate,
/// etc).
///
/// @copyright (c) 2023 CSIRO
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_ACCESSOR_H
#define ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_ACCESSOR_H

// casa includes
#include <casacore/casa/aipstype.h>

// own includes
#include <askap/gridding/UVWeight.h>

namespace askap {

namespace synthesis {

/// @brief Interface class to access UVWeights
/// @details This interface is intended to be used inside the gridder to structure access 
/// to uv weights when they are applied (in traditional weighting). In the simplest form
/// it can be simply the weight collection. But other options may involve index translation
/// (e.g. if we want to apply the weight obtained for one beam to another or interpolate,
/// etc).
/// @ingroup gridding
struct IUVWeightAccessor { 

   /// @brief virtual destructor to keep the compiler happy
   ~IUVWeightAccessor();

   /// @brief obtain weight grid for a given index
   /// @details index interpretation is left for the implementation which can be different
   /// for different gridders. A derived interface provides index translation to standardise this operation. 
   /// Trying to make it more general, we pass beam, field and source indices to this method,
   /// as these indices are currently used in various gridders. This can change in the future, if necessary.
   /// @note Strictly speaking, detection of field changes will be subject to some assumptions. And the 
   /// numbering will be per gridder (figured out from the accessor supplied), other gridders may have
   /// different indexing. The field index is expected to match currentField() output in AProjectGridderBase.
   /// @param[in] beam beam index (from accessor for the given row, it is assumed that we don't have cross-beam correlations)
   /// @param[in] field field index if the gridder is a mosaicing one, zero otherwise
   /// @param[in] source source index used to form direction-dependent offset index (not sure if it is needed, 
   /// add it here just to keep things general as it is used in the gridder code)
   /// @return UVWeight object with the selected grid of weights
   virtual UVWeight getWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) const = 0;
};

} // namespace synthesis

} // namespace askap

#endif // #define ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_ACCESSOR_H

