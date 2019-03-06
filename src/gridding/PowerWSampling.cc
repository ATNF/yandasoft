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

#include <gridding/PowerWSampling.h>
#include <askap/AskapError.h>
#include <cmath>

using namespace askap;
using namespace askap::synthesis;

/// @brief initialise the class
/// @details
/// @param[in] exponent exponent of the desired power law. The largest w-plane number always corresponds to wmax,
/// the centre always correspond to 0.
PowerWSampling::PowerWSampling(const double exponent) : itsExp(exponent) 
{
  ASKAPCHECK(exponent != 0., "Zero exponent in power law doesn't make sense. Use a single w-plane instead");
}


/// @brief plane to w-term conversion (mapping)
/// @details This is a forward method mapping scaled w-plane to scaled w-term.
/// @param[in] plane plane number scaled down to interval [-1;1]
/// @return w-term scaled down to interval [-1;1]
/// @note The result is unpredictable, if the plane is outside [-1;1] interval
double PowerWSampling::map(double plane) const
{
  if (plane > 0) {
      return std::pow(plane,itsExp);
  } else if (plane<0) {
      return -std::pow(-plane,itsExp);
  } 
  return 0.;
}

/// @brief w-term to plane conversion (indexing)
/// @details This is a reverse  method indexing dimensionless w-tern to get dimensionless w-plane.
/// @param[in] wterm scaled down to interval [-1;1]
/// @return w-term scaled down to interval [-1;1]
/// @note The result is unpredictable, if the input w-term is outside [-1;1] interval
double PowerWSampling::index(double wterm) const
{
  if (wterm > 0) {
      return std::pow(wterm, 1. / itsExp);
  } else if (wterm < 0) {
      return -std::pow(-wterm, 1. / itsExp);
  } 
  return 0.;  
}

