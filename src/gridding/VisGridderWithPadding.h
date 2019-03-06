/// @file
/// @brief Base class for gridders supporting padding option
/// @details The gridder factory tries to cast the type of a concrete gridder to this type
/// to determine whether padding is supported. It is also handy take some padding-related
/// functionality from TableVisGridder into this class, so we have access to it from the snap-shot imaging
/// adapter.
///
/// @copyright (c) 2007 CSIRO
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
///

#ifndef VIS_GRIDDER_WITH_PADDING_H
#define VIS_GRIDDER_WITH_PADDING_H

#include <gridding/IVisGridder.h>

namespace askap {

namespace synthesis {


/// @brief Base class for gridders supporting padding option
/// @details The gridder factory tries to cast the type of a concrete gridder to this type
/// to determine whether padding is supported. It is also handy take some padding-related
/// functionality from TableVisGridder into this class, so we have access to it from the snap-shot imaging
/// adapter.
/// @ingroup gridding
struct VisGridderWithPadding : virtual public IVisGridder {
  /// @brief constructor
  /// @details Optionally sets padding factor, defaults to no padding
  /// @param[in] padding padding factor
  explicit VisGridderWithPadding(const float padding = 1.) : itsPaddingFactor(padding) {}
  
  /// @brief return padding factor
  /// @return current padding factor
  float inline paddingFactor() const { return itsPaddingFactor;}

  /// @brief set padding factor
  /// @param[in] padding new padding factor
  void inline setPaddingFactor(const float padding) { itsPaddingFactor = padding;}
private:
  
  /// @brief internal padding factor, 1 by default
  float itsPaddingFactor;    
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef VIS_GRIDDER_WITH_PADDING_H


