/// @file
/// @brief A helper template to manage uv weight index translator   
/// @details Implementation of traditional weighting works with flat indices which may cover
/// different beams, fields, facets, etc. Moreover, it is worth not to design out the possibility
/// to apply different index translation for the case of building weights and applying them. 
/// The IUVWeightIndexTranslator interface class encapsulate such a translation. This helper template
/// represents an adapter allowing us to avoid the need to inherit from the index translator interface
/// and have the aggregation relationship instead.
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_INDEX_TRANSLATION_HELPER_H
#define ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_INDEX_TRANSLATION_HELPER_H

// own includes
#include <askap/gridding/IUVWeightIndexTranslator.h>
#include <askap/askap/AskapError.h>

// boost icnludes
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief A helper template to manage uv weight index translator   
/// @details Implementation of traditional weighting works with flat indices which may cover
/// different beams, fields, facets, etc. Moreover, it is worth not to design out the possibility
/// to apply different index translation for the case of building weights and applying them. 
/// The IUVWeightIndexTranslator interface class encapsulate such a translation. This helper template
/// represents an adapter allowing us to avoid the need to inherit from the index translator interface
/// and have the aggregation relationship instead.
/// @note No copy constructor/assignment operator is defined, so if copy is performed, reference semantics
/// will be used for the translator class as it will be copied by copying shared pointer
/// @ingroup gridding
template<typename Base>
struct UVWeightIndexTranslationHelper : virtual public Base { 
   /// @brief default constructor - does nothing 
   /// @note It is assumed that the Base class is default constructable
   UVWeightIndexTranslationHelper() {}

   /// @brief constructor setting translator class up front
   /// @param[in] translator shared pointer to the translator class to be used with this adapter
   /// @note It is assumed that the Base class is default constructable
   explicit UVWeightIndexTranslationHelper(const boost::shared_ptr<IUVWeightIndexTranslator const> &translator) :
           itsTranslator(translator) {}

   /// @brief set the translator class
   /// @param[in] translator shared pointer to the translator class to be used with this adapter
   void setTranslator(const boost::shared_ptr<IUVWeightIndexTranslator const> &translator) {
        itsTranslator = translator;
   }

   /// @brief obtain flat index for the given set of indices used inside the gridder
   /// @details index interpretation is left for the implementation which can be different
   /// for different gridders. This template class wraps around IUVWeightIndexTranslator and provides 
   /// the translation service via this method
   /// @note Strictly speaking, detection of field changes will be subject to some assumptions. And the 
   /// numbering will be per gridder (figured out from the accessor supplied), other gridders may have
   /// different indexing. The field index is expected to match currentField() output in AProjectGridderBase.
   /// Also, if the translation is not setup an exception is thrown.
   /// @param[in] beam beam index (from accessor for the given row, it is assumed that we don't have cross-beam correlations)
   /// @param[in] field field index if the gridder is a mosaicing one, zero otherwise
   /// @param[in] source source index used to form direction-dependent offset index (not sure if it is needed, 
   /// add it here just to keep things general as it is used in the gridder code)
   /// @return flat index into the weight collection
   casacore::uInt indexOf(casacore::uInt beam, casacore::uInt field, casacore::uInt source) const {
      ASKAPCHECK(itsTranslator, "Index translation has not been set up in UVWeightIndexTranslatonHelper");
      return itsTranslator->indexOf(beam, field, source);
   };
private:
   /// @brief shared pointer for the actual index translation class
   boost::shared_ptr<IUVWeightIndexTranslator const> itsTranslator;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_INDEX_TRANSLATION_HELPER_H

