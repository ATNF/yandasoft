/// @file
/// 
/// @brief Calibration effect: polarisation leakage
/// @details This is a simple effect which can be used in conjunction
/// with the CalibrationME template (as its template argument)
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

#ifndef LEAKAGE_TERM_H
#define LEAKAGE_TERM_H

// own includes
#include <fitting/ComplexDiffMatrix.h>
#include <fitting/ComplexDiff.h>
#include <fitting/Params.h>
#include <dataaccess/IConstDataAccessor.h>
#include <askap/AskapError.h>
#include <measurementequation/ParameterizedMEComponent.h>
#include <utils/PolConverter.h>

// std includes
#include <string>
#include <utility>

namespace askap {

namespace synthesis {

/// @brief Calibration effect: polarisation leakage
/// @details This is a simple effect which can be used in conjunction
/// with the CalibrationME template (as its template argument)
/// @note The relation between leakage terms and Jones matrices matches 
/// the definition of Hamaker, Bregman & Sault. See their equation 
/// (14) for details. Our parameters d12 (corresponding to Stokes:XY) and
/// d21 (corresponding to Stokes::YX) correspond to d_{Ap} and d_{Aq} from
/// Hamaker, Bregman & Sault, respectively.
/// @ingroup measurementequation
struct LeakageTerm : public ParameterizedMEComponent<false> {
   
   /// @brief constructor, store reference to paramters
   /// @param[in] par shared pointer to parameters
   inline explicit LeakageTerm(const scimath::Params::ShPtr &par) : 
                              ParameterizedMEComponent<false>(par) {}
   
   /// @brief main method returning Mueller matrix and derivatives
   /// @details This method has to be overloaded (in the template sense) for
   /// all classes representing various calibration effects. CalibrationME
   /// template will call it when necessary. It returns 
   /// @param[in] chunk accessor to work with
   /// @param[in] row row of the chunk to work with
   /// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
   /// this effect
   inline scimath::ComplexDiffMatrix get(const accessors::IConstDataAccessor &chunk, 
                                casa::uInt row) const;   

};

} // namespace synthesis

} // namespace askap

#include <measurementequation/LeakageTerm.tcc>

#endif // #ifndef LEAKAGE_TERM_H
