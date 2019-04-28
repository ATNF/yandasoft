/// @file
///
/// @brief Power-law w-sampling 
/// @details
/// W-dependent gridders support non-linear sampling in the w-space (through WDependentGridderBase).
/// This class implements IWSampling interface and provides power-law sampling in the w-space. The
/// class is parameterised with a single parameter being the index of the power function. An index 
/// value of 1. is equivalent to the default linear sampling.
///
///     y = sign(x)*pow(abs(x),alpha)
///
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

#ifndef POWER_W_SAMPLING_H
#define POWER_W_SAMPLING_H

#include <gridding/IWSampling.h>

namespace askap {

namespace synthesis {

/// @brief Power-law w-sampling 
/// @details
/// W-dependent gridders support non-linear sampling in the w-space (through WDependentGridderBase).
/// This class implements IWSampling interface and provides power-law sampling in the w-space. The
/// class is parameterised with a single parameter being the index of the power function. An index 
/// value of 1. is equivalent to the default linear sampling.
///
///     y = sign(x)*pow(abs(x),alpha)
///
/// @ingroup gridding
struct PowerWSampling : public IWSampling {

  /// @brief initialise the class
  /// @details
  /// @param[in] exponent exponent of the desired power law. The largest w-plane number always corresponds to wmax,
  /// the centre always correspond to 0.
  explicit PowerWSampling(const double exponent);
    
  /// @brief plane to w-term conversion (mapping)
  /// @details This is a forward method mapping scaled w-plane to scaled w-term.
  /// @param[in] plane plane number scaled down to interval [-1;1]
  /// @return w-term scaled down to interval [-1;1]
  /// @note The result is unpredictable, if the plane is outside [-1;1] interval
  virtual double map(double plane) const;

  /// @brief w-term to plane conversion (indexing)
  /// @details This is a reverse method indexing dimensionless w-tern to get dimensionless w-plane.
  /// @param[in] wterm scaled down to interval [-1;1]
  /// @return w-term scaled down to interval [-1;1]
  /// @note The result is unpredictable, if the input w-term is outside [-1;1] interval
  virtual double index(double wterm) const;  
private:
  /// @brief exponent of the power law
  double itsExp;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef POWER_W_SAMPLING_H

