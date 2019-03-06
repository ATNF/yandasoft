/// @file
/// 
/// @brief Calibration component (or individual effect).
/// @details This class is mainly a structural unit. A derived class
/// ParameterizedMEComponent holds a reference to parameters.
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

#ifndef ME_COMPONENT_H
#define ME_COMPONENT_H

// own includes
#include <fitting/Params.h>

namespace askap {

namespace synthesis {

/// @brief Calibration component (or individual effect).
/// @details This class is mainly a structural unit. A derived class
/// ParameterizedMEComponent holds a reference to parameters.
/// @ingroup measurementequation
template <bool FDP>
struct MEComponent {

   /// @brief check whether the effect is frequency-dependent
   /// @details For frequency-dependent effects the "get" method returns block matrix with 
   /// one block corresponding to every channel (i.e. the size is nPol x nPol*nChannel 
   /// instead of simply nPol x nPol). We need this flag to unroll the matrix multiplication
   /// correctly.
   /// @return true, if the effect is frequency-dependent
   static inline bool isFrequencyDependent() { return FDP; }
      
   /// @brief frequency dependency flag
   static const bool theirFDPFlag = FDP;   
};

} // namespace synthesis

} // namespace askap


#endif // #define ME_COMPONENT_H
