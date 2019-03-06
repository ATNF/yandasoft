/// @file
/// 
/// @brief Calibration effect: Identity Mueller matrix.
/// @details This is a simple effect which doesn't change anything.
/// it is used mainly for debugging. 
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

#ifndef IDENTITY_COMPONENT_H
#define IDENTITY_COMPONENT_H

// own includes
#include <fitting/ComplexDiffMatrix.h>
#include <fitting/ComplexDiff.h>
#include <fitting/Params.h>
#include <dataaccess/IConstDataAccessor.h>
#include <askap/AskapError.h>
#include <measurementequation/MEComponent.h>


// std includes
#include <string>

namespace askap {

namespace synthesis {

/// @brief Calibration effect: antenna gains without cross-pol
/// @details This is a simple effect which can be used in conjunction
/// with the CalibrationME template (as its template argument)
/// @ingroup measurementequation
struct IdentityComponent : public MEComponent<false> {
   
   /// @brief constructor, parameters are actually ignored
   /// @param[in] par const reference to parameters
   inline explicit IdentityComponent(const scimath::Params::ShPtr &par)  {}
   
   /// @brief main method returning Mueller matrix and derivatives
   /// @details This method has to be overloaded (in the template sense) for
   /// all classes representing various calibration effects. CalibrationME
   /// template will call it when necessary. 
   /// @param[in] chunk accessor to work with
   /// @param[in] row row of the chunk to work with
   /// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
   /// this effect
   inline scimath::ComplexDiffMatrix get(const accessors::IConstDataAccessor &chunk, 
                                casa::uInt row) const;
};

/// @brief main method returning Mueller matrix and derivatives
/// @details This method has to be overloaded (in the template sense) for
/// all classes representing various calibration effects. CalibrationME
/// template will call it when necessary. 
/// @param[in] chunk accessor to work with
/// @param[in] row row of the chunk to work with
/// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
/// this effect
inline scimath::ComplexDiffMatrix IdentityComponent::get(const accessors::IConstDataAccessor &chunk, 
                                      casa::uInt row) const
{
   
   scimath::ComplexDiffMatrix calFactor(chunk.nPol(), chunk.nPol(), 0.);

   for (casa::uInt pol=0; pol<chunk.nPol(); ++pol) {            
        calFactor(pol,pol) = scimath::ComplexDiff(1.);            
   }
   return calFactor;
}

} // namespace synthesis

} // namespace askap



#endif // #ifndef IDENTITY_COMPONENT_H
