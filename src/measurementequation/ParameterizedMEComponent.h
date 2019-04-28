/// @file
/// 
/// @brief A calibration component (i.e. individual effect).
/// @details The easiest way of creating individual components is
/// by deriving from this class. This class is mainly a structural
/// unit, but it holds the reference to parameters, which are passed
/// around to all components of the measurement equation.
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

#ifndef PARAMETERIZED_ME_COMPONENT_H
#define PARAMETERIZED_ME_COMPONENT_H

#include <string>

// own includes
#include <fitting/Params.h>
#include <measurementequation/MEComponent.h>
#include <fitting/ComplexDiff.h>
#include <askap/AskapError.h>


namespace askap {

namespace synthesis {

/// @brief A Calibration component (i.e. individual effect).
/// @details The easiest way of creating individual components is
/// by deriving from this class. This class is mainly a structural
/// unit, but it holds a reference to parameters, which are passed
/// around to all components of the measurement equation.
/// @ingroup measurementequation
template<bool FDP>
struct ParameterizedMEComponent : public MEComponent<FDP> {
   
   /// @brief constructor, store reference to paramters
   /// @param[in] par shared pointer to parameters
   inline explicit ParameterizedMEComponent(const scimath::Params::ShPtr &par) : itsParameters(par) {}
   
protected:
   /// @return reference to parameters
   inline scimath::Params::ShPtr parameters() const { return itsParameters; }
   
   /// @brief obtain a value of the parameter
   /// @details This helper method returns the parameter wrapped 
   /// around in a ComplexDiff class. An exception is thrown if the parameter
   /// is not defined.
   /// @param[in] name parameter name
   /// @return value of the parameter wrapped in a complex diff object
   inline scimath::ComplexDiff getParameter(const std::string &name) const;
   
private:
   /// @brief shared pointer to paramters
   scimath::Params::ShPtr itsParameters;
};

} // namespace synthesis

} // namespace askap

#include <measurementequation/ParameterizedMEComponent.tcc>

#endif // #define PARAMETERIZED_ME_COMPONENT_H
