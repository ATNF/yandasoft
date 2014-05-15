/// @file
///
/// @brief An unpolarized component (stokes Q,U and V always give 0)
/// @details
///     This is a derived template from ParameterizedComponent, which
///     represents an unpolarized component. It implements calculate
///     methods via new calculate methods, which don't have pol parameter
///     in their interface and always return stokes I. Having a separate
///     type allows to avoid unnecessary loops in polarization in
///     ComponentEquation, by testing the type with dynamic_cast. 
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

#ifndef UNPOLARIZED_COMPONENT_H
#define UNPOLARIZED_COMPONENT_H

#include <measurementequation/ParameterizedComponent.h>
#include <measurementequation/IUnpolarizedComponent.h>

namespace askap {

namespace synthesis {

/// @brief An unpolarized component (stokes Q,U and V always give 0)
/// @details
///     This is a derived template from ParameterizedComponent, which
///     represents an unpolarized component. It implements calculate
///     methods via new calculate methods, which don't have pol parameter
///     in their interface and always return stokes I. Having a separate
///     type allows to avoid unnecessary loops in polarization in
///     ComponentEquation, by testing the type with dynamic_cast. 
/// @ingroup measurementequation  
template<size_t NComp>
struct UnpolarizedComponent : public ParameterizedComponent<NComp>,
                              virtual public IUnpolarizedComponent {
  
  /// @brief construct the object with a given parameters
  /// @details
  /// @param[in] param parameters of the component (meaning is defined in the
  /// derived classes)
  UnpolarizedComponent(const casa::RigidVector<double, NComp> &param) :
            ParameterizedComponent<NComp>(param) {}
  
};

} // namespace synthesis

} // namespace askap


#endif // #ifndef UNPOLARIZED_COMPONENT_H
