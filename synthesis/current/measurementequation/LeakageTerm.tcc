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

#ifndef LEAKAGE_TERM_TCC
#define LEAKAGE_TERM_TCC

// own includes
#include <utils/PolConverter.h>
#include <calibaccess/CalParamNameHelper.h>

// casa includes
#include <measures/Measures/Stokes.h>

// std includes
#include <string>
#include <utility>

namespace askap {

namespace synthesis {


/// @brief main method returning Mueller matrix and derivatives
/// @details This method has to be overloaded (in the template sense) for
/// all classes representing various calibration effects. CalibrationME
/// template will call it when necessary. It returns 
/// @param[in] chunk accessor to work with
/// @param[in] row row of the chunk to work with
/// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
/// this effect
inline scimath::ComplexDiffMatrix LeakageTerm::get(const accessors::IConstDataAccessor &chunk, 
                                      casa::uInt row) const
{
   const casa::uInt nPol = chunk.nPol();
   ASKAPDEBUGASSERT(nPol != 0);   
   const casa::Vector<casa::Stokes::StokesTypes> stokes = chunk.stokes();   
   ASKAPDEBUGASSERT(stokes.nelements() == nPol);
   ASKAPDEBUGASSERT(!scimath::PolConverter::isStokes(stokes));
   
   
   const casa::uInt ant1 = chunk.antenna1()[row];
   const casa::uInt ant2 = chunk.antenna2()[row];
   
   const casa::uInt beam1 = chunk.feed1()[row];
   const casa::uInt beam2 = chunk.feed2()[row];
  
   // main diagonal is always 1.
   scimath::ComplexDiffMatrix calFactor(4, 4, 1.);
   
   // flag showing that the polarisation products are present
   // in the canonic form (e.g. XX,XY,YX,YY for linears)
   bool canonicPolOrder = (nPol == 4);
   
   calFactor(3, 1) = -1.*getParameter(accessors::CalParamNameHelper::paramName(ant1, beam1, casa::Stokes::YX));
   calFactor(1, 3) = getParameter(accessors::CalParamNameHelper::paramName(ant1, beam1, casa::Stokes::XY));
   
   calFactor(3, 2) = -1.*conj(getParameter(accessors::CalParamNameHelper::paramName(ant2, beam2, casa::Stokes::YX)));
   calFactor(2, 3) = conj(getParameter(accessors::CalParamNameHelper::paramName(ant2, beam2, casa::Stokes::XY)));
   
   for (casa::uInt pol=0; pol<4; ++pol) {
        
        if (pol<nPol) {
            const casa::uInt polIndex = scimath::PolConverter::getIndex(stokes[pol]);
            ASKAPDEBUGASSERT(polIndex<4);
            // polIndex is index in the polarisation frame, i.e.
            // XX is 0, XY is 1, YX is 2 and YY is 3
            // we need an index into matrix 
            if (polIndex != pol) {
                canonicPolOrder = false;
            }
        } else {
              canonicPolOrder = false;
        }

        // cross-diagonal terms                   
        calFactor(pol, 3 - pol) = (pol % 3 == 0 ? 1. : -1.)*getParameter(accessors::CalParamNameHelper::paramName(ant1, 
                            beam1, pol < 2 ? casa::Stokes::XY : casa::Stokes::YX))*
                            conj(getParameter(accessors::CalParamNameHelper::paramName(ant2, beam2, 
                            pol % 2 == 0 ? casa::Stokes::XY : casa::Stokes::YX)));
        if (pol % 3 != 0) {
            // middle rows and columns of the 4x4 matrix (index is 0-based)
            // exploit the symmetries
            calFactor(0, pol) = calFactor(3 - pol, 3);
            calFactor(pol,0) = calFactor(3, 3 - pol);            
        }
   }
   ASKAPCHECK(canonicPolOrder, "Only canonic order of polarisation products (e.g. XX,XY,YX,YY) is currently supported");
   return calFactor;
}

} // namespace synthesis

} // namespace askap


#endif // #ifndef LEAKAGE_TERM_TCC
