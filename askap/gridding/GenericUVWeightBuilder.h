/// @file
/// @brief Simple implementation of IUVWeightBuilder interface
/// @details It holds the collection of weights with flat index and converts the indices
/// required by the interace by applying some custom index translator (a linear translation is
/// setup by default via GenericUVWeightIndexTranslator) in a way similar to GenericUVWeightAccessor.
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_GENERIC_UV_WEIGHT_BUILDER_H
#define ASKAP_SYNTHESIS_GRIDDING_GENERIC_UV_WEIGHT_BUILDER_H

// casa includes
#include <casacore/casa/aipstype.h>

// own includes
#include <askap/gridding/IUVWeightBuilder.h>
#include <askap/gridding/UVWeightIndexTranslationHelper.h>
#include <askap/gridding/GenericUVWeightIndexTranslator.h>
#include <askap/gridding/UVWeightCollection.h>

// other 3rd party
#include <Blob/BlobIStream.h>
#include <Blob/BlobOStream.h>


namespace askap {

namespace synthesis {

/// @brief Simple implementation of IUVWeightBuilder interface
/// @details It holds the collection of weights with flat index and converts the indices
/// required by the interace by applying some custom index translator (a linear translation is
/// setup by default via GenericUVWeightIndexTranslator) in a way similar to GenericUVWeightAccessor.
/// This is sufficient for ignoring some or all indices (via multiplying them by zero) or translate
/// into a flat index with the desired order. 
/// @ingroup gridding
struct GenericUVWeightBuilder : virtual public UVWeightIndexTranslationHelper<IUVWeightBuilder>,
                                public boost::noncopyable { 

   // need to think whether we should allow passing collection from outside or it just should be managed inside the class 
 
   /// @brief construct weight builder from the given collection
   /// @details The default translation of indices is set up by this method. Also, by default, all 
   /// coefficients are zero which means that only index 0 will be used from the collection. 
   /// A different index translation class can be assigned via the appropriate setter method of UVWeightIndexTranslationHelper.
   /// Unlike with weight accessor class, the builder owes the collection. It is setup on demand using parameters passed in 
   /// the initialised method. The translated index is coeffBeam * beam + coeffField * field + coeffSource * source
   /// @param[in] coeffBeam beam index coefficient
   /// @param[in] coeffField field index coefficient
   /// @param[in] coeffSource source index coefficient
   explicit GenericUVWeightBuilder(casacore::uInt coeffBeam = 0, casacore::uInt coeffField = 0, casacore::uInt coeffSource = 0);

   /// @brief reset the object to a pristine state like after the default constructor
   /// @details This method is a bit of the technical debt. It is not necessary for the uv-weight framework itself, but
   /// needed to be able to reuse normal equations reduction code via the EstimatorAdapter. As defined in the normal equations
   /// class hierarchy, it is supposed to restore the state of object as it would be after the default constructor. Note, in
   /// the case of this class, the default constructor actually implies using zero for all coefficients (coeffBeam, etc). So if the
   /// original object was setup with non-trivial mapping it would change. It is, however, expected that we won't add new data
   /// to the builder following either reset or merge call (in our use case, reset could happen as part of the merge and only after
   /// the appropriate portion of data has already been accumulated).
   void reset();

   // a copy of the accessor method. May need changes, at least to the doco as it is not expected to be used
   // it may be possible to have the same behaviour as addWeight  

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

   /// @brief setup the shape for the uv weight grid
   /// @details This method can be called from a gridder when the grid shape is setup to
   /// ensure uv weights are collected on the grid with the same dimensions. At this stage we assume
   /// that the gridder has the same grid size for all potential weight grids (i.e. for different beams).
   /// This is the only case we support anyway inside the gridder.
   /// @param[in] uSize size in the direction of u-coordinate
   /// @param[in] vSize size in the direction of v-coordinate
   /// @param[in] nPlanes number of planes for the given weight (3rd dimension)
   /// @note This class is agnostic about the physical pixel sizes (e.g. uv cell size) -
   /// we assume it is always the same as in the setup of the gridder
   void initialise(casacore::uInt uSize, casacore::uInt vSize, casacore::uInt nPlanes) final;

   /// @brief obtain weight grid for writing for the given metadata (indices)
   /// @details index interpretation is left for the implementation which can be different
   /// for different gridders. Trying to make it more general, we pass beam, field and source indices
   /// as they are currently used in various gridders. This can change in the future, if necessary.
   /// @note Strictly speaking, detection of field changes will be subject to some assumptions. And the 
   /// numbering will be per gridder (figured out from the accessor supplied), other gridders may have
   /// different indexing. The field index is expected to match currentField() output in AProjectGridderBase.
   /// In this class, the index translator (initially setup to the generic version with parameters passed in the constructor)
   /// is used to convert generic indices passed to this method to the flat index of the collection.
   /// @param[in] beam beam index (from accessor for the given row, it is assumed that we don't have cross-beam correlations)
   /// @param[in] field field index if the gridder is a mosaicing one, zero otherwise
   /// @param[in] source source index used to form direction-dependent offset index (not sure if it is needed, 
   /// add it here just to keep things general as it is used in the gridder code)
   /// @return UVWeight object with the selected grid of weights
   /// @note This method can be called multiple times for the same indices, if necessary. The new grid will be created on demand.
   UVWeight addWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) final;

   /// @brief merge weight information from the other builder
   /// @details This method is expected to be used in conjunction with the EstimatorAdapter and normal equation tree reduction
   /// routines. It enables merging multiple weight grids (corresponding to the same indices) together in distributed data
   /// processing.
   /// @param[in] src other builder to merge from
   /// @note The source should be the object which can be cast to the same type
   void merge(const IUVWeightBuilder &src) final;

   /// @brief finalise uv weight computation by applying some post-processing to accumulated weight grids
   /// @details This method is expected to call the given object function for each weight in the collection.
   /// Examples are calculation of robust weights, ensure Hermitian symmetry, uv-taper or a combination of those.
   /// It is envisaged that we'll work with a chain of algorithms, but keep this complexity outside of this interface.
   /// @param[in] calc processing algorithm/object function 
   /// @return a reference to the uv-weight collection which is the final result (may be useful if the weight collection is
   /// owed by implementation of this interface instead of being supplied by reference).
   /// MV: should it be the const reference? The class itself is non-copyable by design.
   UVWeightCollection& finalise(const IUVWeightCalculator &calc) final;

   // do we need a getter method to obtain a const reference to weight collection outside of the finalise call?
   // if so, it needs to be added in the interface

   // serialisation / deserialisation

   /// @brief write the object to a blob stream
   /// @param[in] os the output stream
   void writeToBlob(LOFAR::BlobOStream& os) const;

   /// @brief read the object from a blob stream
   /// @param[in] is the input stream
   /// @note Not sure whether the parameter should be made const or not
   void readFromBlob(LOFAR::BlobIStream& is);
 
private:

   /// @brief weight collection
   /// @note Unlike for the accessor, the builder ows weight collection. 
   UVWeightCollection itsWeights;

   // default shape of the weight grid, set via initialise method

   /// @brief default size in the direction of u-coordinate
   casacore::uInt itsUSize;
   
   /// @brief default size in the direction of v-coordinate
   casacore::uInt itsVSize;

   /// @brief default number of planes in the weight grid (3rd dimension)
   casacore::uInt itsNPlanes;

   /// @brief payload version for the serialisation
   static constexpr int theirPayloadVersion = 1;
};

} // namespace synthesis

} // namespace askap

#endif // #define ASKAP_SYNTHESIS_GRIDDING_GENERIC_UV_WEIGHT_BUILDER_H

