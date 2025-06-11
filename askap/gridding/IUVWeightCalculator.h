/// @file
/// @brief Interface class to an object function working on UVWeights
/// @details The derived classes of this interface are the ones which actually calculate
/// robust or other weights from the accumulated grid of weights. All operations are expected
/// to be in situ. The object function is called in the finalise method of the builder class.
/// At this stage we assume that all frequency planes and all indices can be processed independently
/// (and potentially in parallel). Therefore, this method works with the matrix interface rather than
/// the cube (low-level access is needed for performance here).
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_CALCULATOR_H
#define ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_CALCULATOR_H

// casa includes
#include <casacore/casa/aipstype.h>
#include <casacore/casa/Arrays/Matrix.h>

namespace askap {

namespace synthesis {

/// @brief Interface class to an object function working on UVWeights
/// @details The derived classes of this interface are the ones which actually calculate
/// robust or other weights from the accumulated grid of weights. All operations are expected
/// to be in situ. The object function is called in the finalise method of the builder class.
/// At this stage we assume that all frequency planes and all indices can be processed independently
/// (and potentially in parallel). Therefore, this method works with the matrix interface rather than
/// the cube (low-level access is needed for performance here).
/// @ingroup gridding
struct IUVWeightCalculator { 

   /// @brief virtual destructor to keep the compiler happy
   ~IUVWeightCalculator();

   /// @brief perform processing for the given weight (single grid slice along the 3rd axis)
   /// @details For performance reasons, slices along the 3rd axis are taken inside finalise method
   /// of the builder (this can be changed if we ever had any effect where frequency dependence matter).
   /// At this stage, we can guarantee that supplied matrix has contiguous storage.
   /// @param[in] wt weight to work with (it is modified in situ).
   /// @note The shape is supposed to stay intact.
   virtual void process(casacore::Matrix<float> &wt) const = 0;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_I_UV_WEIGHT_CALCULATOR_H
