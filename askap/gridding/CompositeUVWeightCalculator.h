/// @file
/// @brief weight calculator which calls others in a sequence
/// @details This is an implementation of UVWeight calculator interface which calls a list of other
/// algorithms in a sequence (e.g. adding conjugates, robust weighting, tapering).
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_COMPOSITE_UV_WEIGHT_CALCULATOR_H
#define ASKAP_SYNTHESIS_GRIDDING_COMPOSITE_UV_WEIGHT_CALCULATOR_H

// own includes
#include <askap/gridding/IUVWeightCalculator.h>

// std includes
#include <vector>

// boost includes
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief weight calculator which calls others in a sequence
/// @details This is an implementation of UVWeight calculator interface which calls a list of other
/// algorithms in a sequence (e.g. adding conjugates, robust weighting, tapering).
/// @ingroup gridding
struct CompositeUVWeightCalculator : virtual public IUVWeightCalculator {

   /// @brief default constructor - empty task list
   CompositeUVWeightCalculator() {}

   /// @brief generic constructor from begin/end iterators
   /// @details This constructor builds task list straight away using begin and end iterators. The
   /// value type must be convertible to shared pointer to IUVWeightCalculator type.
   /// @param[in] begin start iterator
   /// @param[in] end end iterator (i.e. pointing to an element after last)
   template<typename Iter>
   CompositeUVWeightCalculator(const Iter &begin, const Iter &end) : itsTasks(begin, end) {}

   /// @brief add a new procedure to the task list 
   /// @details The given shared pointer is appended to the task list. So this method should
   /// be called in the order the algorithms are supposed to be executed.
   /// @param[in] task new procedure to add
   void add(const boost::shared_ptr<IUVWeightCalculator> &task);

   /// @brief perform processing for the given weight (single grid slice along the 3rd axis)
   /// @details For performance reasons, slices along the 3rd axis are taken inside finalise method
   /// of the builder (this can be changed if we ever had any effect where frequency dependence matter).
   /// At this stage, we can guarantee that supplied matrix has contiguous storage.
   /// @param[in] wt weight to work with (it is modified in situ).
   /// @note The shape is supposed to stay intact.
   void process(casacore::Matrix<float> &wt) const final;

private:
   
   /// @brief individual procedures we want to apply in the order we want them
   std::vector<boost::shared_ptr<IUVWeightCalculator> > itsTasks;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_COMPOSITE_UV_WEIGHT_CALCULATOR_H

