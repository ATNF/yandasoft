/// @file
/// @brief Simple implementation of IUVWeightAccessor interface
/// @details It holds the collection of weights with flat index and converts the indices
/// required by the interace by applying some custom index translator (a linear translation is
/// setup by default via GenericUVWeightIndexTranslator). 
/// This is sufficient for ignoring some or all indices (via multiplying them by zero) or translate
/// into a flat index with the desired order.
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_GENERIC_UV_WEIGHT_ACCESSOR_H
#define ASKAP_SYNTHESIS_GRIDDING_GENERIC_UV_WEIGHT_ACCESSOR_H

// casa includes
#include <casacore/casa/aipstype.h>

// own includes
#include <askap/gridding/IUVWeightAccessor.h>
#include <askap/gridding/UVWeightIndexTranslationHelper.h>
#include <askap/gridding/GenericUVWeightIndexTranslator.h>
#include <askap/gridding/UVWeightCollection.h>

// boost includes
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief Simple implementation of IUVWeightAccessor interface
/// @details It holds the collection of weights with flat index and converts the indices
/// required by the interace by applying some custom index translator (a linear translation is
/// setup by default via GenericUVWeightIndexTranslator). 
/// This is sufficient for ignoring some or all indices (via multiplying them by zero) or translate
/// into a flat index with the desired order.
/// @ingroup gridding
struct GenericUVWeightAccessor : virtual public UVWeightIndexTranslationHelper<IUVWeightAccessor> { 

   /// @brief construct weight accessor from the given collection
   /// @details The default translation of indices is set up by this method. Also, by default, all 
   /// coefficients are zero which means that only index 0 will be used from the collection. 
   /// A different index translation class can be assigned via the appropriate setter method of UVWeightIndexTranslationHelper.
   /// The constructor is not "explicit" by intention. This way we can cast it from the collection type directly.
   /// Otherwise, the resulting index is coeffBeam * beam + coeffField * field + coeffSource * source
   /// @param[in] wts weight collection
   /// @param[in] coeffBeam beam index coefficient
   /// @param[in] coeffField field index coefficient
   /// @param[in] coeffSource source index coefficient
   GenericUVWeightAccessor(const UVWeightCollection &wts, casacore::uInt coeffBeam = 0,
                           casacore::uInt coeffField = 0, casacore::uInt coeffSource = 0);

   /// @brief construct weight accessor via shared pointer interfaces
   /// @details This constructor should be used when weight accessor is expected to absorb the ownership of the collection, i.e.
   /// it is not used as an adapter. This is handy if say a collection is built based on the model parameters. Care must be taken to
   /// manage the life-cycle, i.e. the weight data will be preserved as long as all copies of the accessor object are alive.
   /// @param[in] wts shared pointer to weight collection
   /// @param[in] translator shared pointer to the index translator object 
   GenericUVWeightAccessor(const boost::shared_ptr<UVWeightCollection const> &wts, 
                           const boost::shared_ptr<IUVWeightIndexTranslator const> &translator);

   /// @brief obtain weight grid for a given index
   /// @details index interpretation is left for the implementation which can be different
   /// for different gridders. Trying to make it more general, we pass beam, field and source indices
   /// as they are currently used in various gridders. This can change in the future, if necessary.
   /// @note Strictly speaking, detection of field changes will be subject to some assumptions. And the 
   /// numbering will be per gridder (figured out from the accessor supplied), other gridders may have
   /// different indexing. The field index is expected to match currentField() output in AProjectGridderBase.
   /// @param[in] beam beam index (from accessor for the given row, it is assumed that we don't have cross-beam correlations)
   /// @param[in] field field index if the gridder is a mosaicing one, zero otherwise
   /// @param[in] source source index used to form direction-dependent offset index (not sure if it is needed, 
   /// add it here just to keep things general as it is used in the gridder code)
   /// @return UVWeight object with the selected grid of weights
   UVWeight getWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) const final;

private:

   /// @brief shared pointer to weight collection
   /// @note reference semantics is used by casacore arrays, so the weight collection itself is not very large (although grows with the
   /// number of elements) and we could've held it by value (but it would need to be made copyable). Having a reference is handy when this
   /// object is used as a wrapper. Having a shared pointer allows us to absorb ownership (although note that in this case one must be careful
   /// to avoid having extra copy hanging out). The reference can be stored in a shared pointer with NullDeleter essentially enabling both interfaces.
   const boost::shared_ptr<UVWeightCollection const> itsWeights;
};

} // namespace synthesis

} // namespace askap

#endif // #define ASKAP_SYNTHESIS_GRIDDING_GENERIC_UV_WEIGHT_ACCESSOR_H

