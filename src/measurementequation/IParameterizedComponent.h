/// @file
///
/// @brief An abstract interface for a component depending on a number of 
///        parameters
/// @details
///     This interface does not implement any method. It is a structural type
///     allowing to hold a number of derived objects, which potentially
///     depend on a different number of parameters, in a container.
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

#ifndef I_PARAMETERIZED_COMPONENT_H
#define I_PARAMETERIZED_COMPONENT_H

// own includes
#include <measurementequation/IComponent.h>

// std includes
#include <string>

namespace askap {

namespace synthesis {

/// @brief An abstract interface for a component depending on a number of 
///        parameters
/// @details
///     This interface does not implement any method. It is a structural type
///     allowing to hold a number of derived objects, which potentially
///     depend on a different number of parameters, in a container.
/// @ingroup measurementequation  
struct IParameterizedComponent : virtual public IComponent {
  /// @brief get number of parameters
  /// @return a number of parameters  this component depends upon. 
  virtual size_t nParameters() const throw() = 0;
  
  /// @brief get the name of the given parameter
  /// @details All parameters are handled in the synthesis code using their
  /// string name, which allows to fix or free any of them easily. This method
  /// allows to obtain this string name using a integer index
  /// @param[in] index an integer index of the parameter (should be less than
  /// nParameters).
  /// @return a const reference to the string name of the parameter 
  virtual const std::string& parameterName(size_t index) const = 0;
};

} // namespace synthesis

} // namespace askap


#endif // #ifndef I_PARAMETERIZED_COMPONENT_H
