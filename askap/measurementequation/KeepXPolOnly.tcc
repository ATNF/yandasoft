/// @file
/// 
/// @brief Calibration effect: zeros parallel-hand products and keeps cross-pols
/// @details This effect does not itroduce any parameters, it simply zeros parallel hand products.
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

#ifndef KEEP_X_POL_ONLY_TCC
#define KEEP_X_POL_ONLY_TCC

// own includes
#include <askap/scimath/utils/PolConverter.h>
#include <askap/calibaccess/CalParamNameHelper.h>

// casa includes
#include <casacore/measures/Measures/Stokes.h>

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
inline scimath::ComplexDiffMatrix KeepXPolOnly::get(const accessors::IConstDataAccessor &chunk, 
                                      casacore::uInt row) const
{
   const casacore::uInt nPol = chunk.nPol();
   ASKAPDEBUGASSERT(nPol != 0);   
   const casacore::Vector<casacore::Stokes::StokesTypes> stokes = chunk.stokes();   
   ASKAPDEBUGASSERT(stokes.nelements() == nPol);
   ASKAPDEBUGASSERT(!scimath::PolConverter::isStokes(stokes));
   
   scimath::ComplexDiffMatrix calFactor(nPol, nPol, 0.);

   for (casacore::uInt pol=0; pol<nPol; ++pol) {
        
        const casacore::uInt polIndex = scimath::PolConverter::getIndex(stokes[pol]);
        // polIndex is index in the polarisation frame, i.e.
        // XX is 0, XY is 1, YX is 2 and YY is 3
        // pol is the index into matrix 

        // keep only cross-pols
        if ((polIndex == 1) || (polIndex == 2)) {
            calFactor(pol, pol) = 1.;
        }        
   }
   return calFactor;
}

} // namespace synthesis

} // namespace askap


#endif // #ifndef KEEP_X_POL_ONLY_TCC
