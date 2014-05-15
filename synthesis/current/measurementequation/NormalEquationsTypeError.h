/// @file
/// 
/// @brief An error class passed in the exception during an attempt to use 
/// a wrong type of the normal equations class.
/// @details Previously, an instance of AskapError was thrown if 
/// ImagingMultiChunkEquation or GenericMultiChunkEquation encountered a
/// wrong type of normal equations. However, in composite equations
/// this particular error must be ignored, while other occurences also 
/// producing AskapError must not. Here is a special exception class for
/// type error related to NormalEquations which allows a separate handling
/// of this particular type of exception

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

#ifndef NORMAL_EQUATIONS_TYPE_ERROR_H
#define NORMAL_EQUATIONS_TYPE_ERROR_H

#include <askap/AskapError.h>
#include <string>

namespace askap {

namespace synthesis {

/// @brief An error class passed in the exception during an attempt to use 
/// a wrong type of the normal equations class.
/// @details Previously, an instance of AskapError was thrown if 
/// ImagingMultiChunkEquation or GenericMultiChunkEquation encountered a
/// wrong type of normal equations. However, in composite equations
/// this particular error must be ignored, while other occurences also 
/// producing AskapError must not. This is a special exception class for
/// type error related to NormalEquations which allows a separate handling
/// of this particular type of exception
struct NormalEquationsTypeError: public askap::AskapError
{
  /// Constructor taking a message
  /// @param[in] message Message string
  explicit NormalEquationsTypeError(const std::string& message);
};


} // namespace synthesis

} // namespace askap

#endif // #ifndef NORMAL_EQUATIONS_TYPE_ERROR_H
