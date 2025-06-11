/// @file
/// @brief Interface class to translate indices for uv weight access
/// @details Implementation of traditional weighting works with flat indices which may cover
/// different beams, fields, facets, etc. Moreover, it is worth not to design out the possibility
/// to apply different index translation for the case of building weights and applying them. 
/// This interface class encapsulate such a translation. It can be enabled via the uv weight accessor
/// or builder interfaces.
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_INDEX_TRANSLATOR_H
#define ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_INDEX_TRANSLATOR_H

// casa includes
#include <casacore/casa/aipstype.h>

namespace askap {

namespace synthesis {

/// @brief Interface class to translate indices for uv weight access
/// @details Implementation of traditional weighting works with flat indices which may cover
/// different beams, fields, facets, etc. Moreover, it is worth not to design out the possibility
/// to apply different index translation for the case of building weights and applying them. 
/// This interface class encapsulate such a translation. It can be enabled via the uv weight accessor
/// or builder interfaces.
/// @ingroup gridding
struct IUVWeightIndexTranslator { 

   /// @brief virtual destructor to keep the compiler happy
   ~IUVWeightIndexTranslator();


   /// @brief obtain uv weight index for given gridder indices
   /// @details index interpretation is left for the implementation which can be different
   /// for different gridders. An implementation of this interface defines a possible translation
   /// from the indices passed to IUVWeightAccessor or IUVWeightBuilder interfaces to the flat index used in
   /// collection of weight grids.  Trying to make the framework more general, we pass beam, field and source 
   /// indices to these interfaces as they are currently used in various griddersi to select convolution functions. 
   /// This can be changed in the future, if necessary.
   /// @note Strictly speaking, detection of field changes will be subject to some assumptions. And the 
   /// numbering will be per gridder (figured out from the accessor supplied), other gridders may have
   /// different indexing. The field index is expected to match currentField() output in AProjectGridderBase.
   /// @param[in] beam beam index (from accessor for the given row, it is assumed that we don't have cross-beam correlations)
   /// @param[in] field field index if the gridder is a mosaicing one, zero otherwise
   /// @param[in] source source index used to form direction-dependent offset index (not sure if it is needed, 
   /// add it here just to keep things general as it is used in the gridder code)
   /// @return flat index into the weight collection
   virtual casacore::uInt indexOf(casacore::uInt beam, casacore::uInt field, casacore::uInt source) const = 0;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_INDEX_TRANSLATOR_H

