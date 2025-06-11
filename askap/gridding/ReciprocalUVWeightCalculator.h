/// @file
/// @brief Class calculating inverse of the weights for all elements of the grid.
/// @details This calculator simply fills the grid with reciprocal for each value.
/// The optional parameter is a threshold. All weights below threshold are assumed to be 0.
/// 
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_RECIPROCAL_UV_WEIGHT_CALCULATOR_H
#define ASKAP_SYNTHESIS_GRIDDING_RECIPROCAL_UV_WEIGHT_CALCULATOR_H

// own includes
#include <askap/gridding/IUVWeightCalculator.h>

namespace askap {

namespace synthesis {

/// @brief Class calculating inverse of the weights for all elements of the grid.
/// @details This calculator simply fills the grid with reciprocal for each value.
/// The optional parameter is a threshold. All weights below threshold are assumed to be 0.
/// @ingroup gridding
struct ReciprocalUVWeightCalculator : virtual public IUVWeightCalculator {

   /// @brief initialise the calculator, set desired robustness
   /// @details
   /// @param[in] threshold threshold to clip small weight pixels
   ReciprocalUVWeightCalculator(float threshold) : itsThreshold(threshold) {}

   /// @brief perform processing for the given weight (single grid slice along the 3rd axis)
   /// @details For performance reasons, slices along the 3rd axis are taken inside finalise method
   /// of the builder (this can be changed if we ever had any effect where frequency dependence matter).
   /// At this stage, we can guarantee that supplied matrix has contiguous storage.
   /// @param[in] wt weight to work with (it is modified in situ).
   /// @note The shape is supposed to stay intact.
   void process(casacore::Matrix<float> &wt) const final;

private:
   /// @brief threshold parameter, all weights below this value are replaced by zero
   /// @note it is made const because we don't need to change it after construction (can be changed if needed)
   const float itsThreshold;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_RECIPROCAL_UV_WEIGHT_CALCULATOR_H
