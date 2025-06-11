/// @file
/// @brief Class calculating robust weights from the grid of weights
/// @details This is an implementation of UVWeight calculator interface making robust weights.
/// The math follows the original code by Daniel Mitchell, who presumably checked it against casa.
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_ROBUST_UV_WEIGHT_CALCULATOR_H
#define ASKAP_SYNTHESIS_GRIDDING_ROBUST_UV_WEIGHT_CALCULATOR_H

// own includes
#include <askap/gridding/IUVWeightCalculator.h>

namespace askap {

namespace synthesis {

/// @brief Class calculating robust weights from the grid of weights
/// @details This is an implementation of UVWeight calculator interface making robust weights.
/// The math follows the original code by Daniel Mitchell, who presumably checked it against casa.
/// @ingroup gridding
struct RobustUVWeightCalculator : virtual public IUVWeightCalculator {

   /// @brief initialise the calculator, set desired robustness
   /// @details
   /// @param[in] robustness robustness parameter
   RobustUVWeightCalculator(float robustness) : itsRobustness(robustness) {}

   /// @brief perform processing for the given weight (single grid slice along the 3rd axis)
   /// @details For performance reasons, slices along the 3rd axis are taken inside finalise method
   /// of the builder (this can be changed if we ever had any effect where frequency dependence matter).
   /// At this stage, we can guarantee that supplied matrix has contiguous storage.
   /// @param[in] wt weight to work with (it is modified in situ).
   /// @note The shape is supposed to stay intact.
   void process(casacore::Matrix<float> &wt) const final;
private:
   /// @brief robustness parameter
   /// @note it is made const because we don't need to change it after construction (can be changed if needed)
   const float itsRobustness;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_ROBUST_UV_WEIGHT_CALCULATOR_H
