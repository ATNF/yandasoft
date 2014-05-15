/// @file
/// 
/// @brief Interface for any unary object function acting on the visibility cube
/// @details Unlike calibration interfaces, this one doesn't pass any metadata and
/// therefore can work with a casa cube instead of the accessor. The main motivation
/// was the optional interrank communication for parallel measurement equation.
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

#ifndef I_VIS_CUBE_UPDATE_H
#define I_VIS_CUBE_UPDATE_H

// casa includes
#include <casa/Arrays/Cube.h>
#include <casa/BasicSL/Complex.h>

namespace askap {

namespace synthesis {

/// @brief Interface for any unary object function acting on the visibility cube
/// @details Unlike calibration interfaces, this one doesn't pass any metadata and
/// therefore can work with a casa cube instead of the accessor. The main motivation
/// was the optional interrank communication for parallel measurement equation.
/// @ingroup measurementequation
struct IVisCubeUpdate {

  /// @brief virtual destructor (does nothing in this class)
  virtual ~IVisCubeUpdate();
  
  /// @brief update visibility cube
  /// @param[in,out] cube reference to visiblity cube to update 
  virtual void update(casa::Cube<casa::Complex> &cube) const = 0;
  
  /// @brief aggregate flag with the logical or operation
  /// @param[in,out] flag flag to reduce
  virtual void aggregateFlag(bool &flag) const = 0;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef I_VIS_CUBE_UPDATE_H


