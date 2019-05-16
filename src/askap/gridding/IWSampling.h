/// @file
///
/// @brief Interface for generic w-samplign helper 
/// @details
/// W-dependent gridders support non-linear sampling in the w-space (through WDependentGridderBase).
/// A class derived from this interface is responsible for mapping [-1,1] domain into [-1,1] range and
/// back to facilitate plane to w-term conversion. The WDependentGridderBase does the appropriate scaling
/// (planes from 0 to nplanes-1 are mapped to w-terms from -wmax to wmax). For example, a simple identity
/// transform (y=x) is equivalent to the linear sampling.
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

#ifndef I_W_SAMPLING_H
#define I_W_SAMPLING_H

namespace askap {

namespace synthesis {

/// @brief Interface for generic w-samplign helper 
/// @details
/// W-dependent gridders support non-linear sampling in the w-space (through WDependentGridderBase).
/// A class derived from this interface is responsible for mapping [-1,1] domain into [-1,1] range and
/// back to facilitate plane to w-term conversion. The WDependentGridderBase does the appropriate scaling
/// (planes from 0 to nplanes-1 are mapped to w-terms from -wmax to wmax). For example, a simple identity
/// transform (y=x) is equivalent to the linear sampling.
/// @ingroup gridding
struct IWSampling {
  
  /// @brief just to keep the compiler happy
  inline virtual ~IWSampling() {};
  
  /// @brief plane to w-term conversion (mapping)
  /// @details This is a forward abstract method mapping scaled w-plane to scaled w-term.
  /// @param[in] plane plane number scaled down to interval [-1;1]
  /// @return w-term scaled down to interval [-1;1]Œ
  /// @note The result is unpredictable, if the plane is outside [-1;1] interval
  virtual double map(double plane) const = 0;

  /// @brief w-term to plane conversion (indexing)
  /// @details This is a reverse abstract method indexing dimensionless w-tern to get dimensionless w-plane.
  /// @param[in] wterm scaled down to interval [-1;1]
  /// @return w-term scaled down to interval [-1;1]
  /// @note The result is unpredictable, if the input w-term is outside [-1;1] interval
  virtual double index(double wterm) const = 0;  
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef I_W_SAMPLING_H

