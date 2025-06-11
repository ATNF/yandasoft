/// @file
/// @brief Interface class to create/write UVWeights
/// @details This interface is to be used inside the gridder to structure access 
/// to uv weights when they are created (in traditional weighting). An alternative gridder-less
/// implementation is also possible (essentially replicating the functionality of the box gridder
/// without the rest of the overhead of a gridder class). The weight access itself 
/// is similar to IUVWeightAccessor interface. The only difference in behaviour is that a new element
/// is supposed to be created if not already known instead of generating an exception. Also, the derived
/// classes are expected to owe weight collection (as opposed to refer to it by reference in the case of
/// the accessor class, although this is not the requirement and is hidden behind the interface) to simplify
/// operations in the distributed data environment.
/// In addition, the interface provides other methods essential for the full job of creating weights.
/// For the time being, this interface is derived from IUVWeightAccessor, but it can be changed.
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_BUILDER_H
#define ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_BUILDER_H

// casa includes
#include <casacore/casa/aipstype.h>

// own includes
#include <askap/gridding/UVWeight.h>
#include <askap/gridding/IUVWeightAccessor.h>
#include <askap/gridding/IUVWeightCalculator.h>
#include <askap/gridding/UVWeightCollection.h>

namespace askap {

namespace synthesis {

/// @brief Interface class to create/write UVWeights
/// @details This interface is to be used inside the gridder to structure access 
/// to uv weights when they are created (in traditional weighting). An alternative gridder-less
/// implementation is also possible (essentially replicating the functionality of the box gridder
/// without the rest of the overhead of a gridder class). The weight access itself 
/// is similar to IUVWeightAccessor interface. The only difference in behaviour is that a new element
/// is supposed to be created if not already known instead of generating an exception. Also, the derived
/// classes are expected to owe weight collection (as opposed to refer to it by reference in the case of
/// the accessor class, although this is not the requirement and is hidden behind the interface) to simplify
/// operations in the distributed data environment.
/// In addition, the interface provides other methods essential for the full job of creating weights.
/// For the time being, this interface is derived from IUVWeightAccessor, but it can be changed.
/// @ingroup gridding
struct IUVWeightBuilder : virtual public IUVWeightAccessor { 

   /// @brief setup shape for the uv weight grid
   /// @details This method can be called from a gridder when the grid shape is setup to
   /// ensure uv weights are collected on the grid with the same dimensions. At this stage we assume
   /// that the gridder has the same grid size for all potential weight grids (i.e. for different beams).
   /// This is the only case we support anyway inside the gridder.
   /// @param[in] uSize size in the direction of u-coordinate
   /// @param[in] vSize size in the direction of v-coordinate
   /// @param[in] nPlanes number of planes for the given weight (3rd dimension)
   /// @note This class is agnostic about the physical pixel sizes (e.g. uv cell size) -
   /// we assume it is always the same as in the setup of the gridder
   virtual void initialise(casacore::uInt uSize, casacore::uInt vSize, casacore::uInt nPlanes) = 0;

   /// @brief obtain weight grid for writing for the given metadata (indices)
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
   /// @note Due to reference semantics of casacore arrays, we can access data for writing. It may be possible to
   /// reuse IUVWeightAccessor::getWeight method for this purpose. But for the sake of clarity, make a separate method for now.
   /// It can be called multiple times for the same indices, if necessary. The new grid will be created on demand.
   virtual UVWeight addWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) = 0;

   /// @brief merge with other builder
   /// @details This method is expected to be used in conjunction with the EstimatorAdapter and normal equation tree reduction
   /// routines. It enables merging multiple weight grids (corresponding to the same indices) together in distributed data
   /// processing.
   /// @param[in] src other builder to merge from
   virtual void merge(const IUVWeightBuilder &src) = 0;

   /// @brief finalise uv weight computation by applying some post-processing to accumulated weight grids
   /// @details This method is expected to call the given object function for each weight in the collection.
   /// Examples are calculation of robust weights, ensure Hermitian symmetry, uv-taper or a combination of those.
   /// It is envisaged that we'll work with a chain of algorithms, but keep this complexity outside of this interface.
   /// @param[in] calc processing algorithm/object function 
   /// @return a reference to the uv-weight collection which is the final result (may be useful if the weight collection is
   /// owed by implementation of this interface instead of being supplied by reference).
   virtual UVWeightCollection& finalise(const IUVWeightCalculator &calc) = 0;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_BUILDER_H

