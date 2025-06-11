/// @file
/// @brief Simple linear translation of indices
/// @details It converts the indices used inside the gridder (through the interface to access uv weights)
/// into a flat index used by weight collection class by applying linear factors and summing across (i.e. like a dot product).
/// This is sufficient for ignoring some or all indices (via multiplying them by zero) or translate
/// into a flat index with the desired order and stride. This class is a basic implementation of IUVWeightIndexTranslator
/// interface.
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

// own includes
#include <askap/gridding/GenericUVWeightIndexTranslator.h>

namespace askap {

namespace synthesis {

/// @brief construct the translator
/// @details The translation of indices is set up by this method. By default, all coefficients are zero
/// which means that only index 0 will be used from the collection regardless of the indices supplied by the gridder. 
/// Otherwise, the resulting index is coeffBeam * beam + coeffField * field + coeffSource * source
/// @note At this stage, it is not envisaged that we'd change the indices during the lifetime of this object, so the
/// data members are defined as const to help with the optimisation (i.e. the values are fixed at construction)
/// @param[in] coeffBeam beam index coefficient
/// @param[in] coeffField field index coefficient
/// @param[in] coeffSource source index coefficient
GenericUVWeightIndexTranslator::GenericUVWeightIndexTranslator(casacore::uInt coeffBeam, casacore::uInt coeffField, 
         casacore::uInt coeffSource) : itsBeamCoeff(coeffBeam), itsFieldCoeff(coeffField), itsSourceCoeff(coeffSource) {}

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
casacore::uInt GenericUVWeightIndexTranslator::indexOf(casacore::uInt beam, casacore::uInt field, casacore::uInt source) const 
{
   return itsBeamCoeff * beam + itsFieldCoeff * field + itsSourceCoeff * source;
}

} // namespace synthesis

} // namespace askap

