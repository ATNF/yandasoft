/// @file
///
/// @brief Generic calibration information for one antenna
/// @details This is essentially a data holding structure which is 
/// intended to be used with OpCalImpl and related solvers. The main
/// OpCalImpl class obtains delays and complex phases for each interval 
/// identified bases on some criteria. These results are aggregated in a
/// single structure per antenna for each interval and then are passed to
/// a high-level "solver" for analysis. The main reason behind these classes
/// is experimentation with calibration which otherwise would require too much
/// scripting and/or hand work if ccalibrator is used directly.
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>


#include <opcal/GenericCalInfo.h>
#include <askap/AskapError.h>

namespace askap {

namespace synthesis {

/// @brief default constructor - all parameters are undefined
GenericCalInfo::GenericCalInfo() : itsGainDefined(false), itsDelayDefined(false) {}
  
/// @brief constructor setting parameters explicitly
/// @details
/// @param[in] gain complex gain to set
/// @param[in] delay delay to set  
GenericCalInfo::GenericCalInfo(const casa::Complex &gain, double delay) : itsGain(gain),
   itsGainDefined(true), itsDelay(delay), itsDelayDefined(true) {}
  
/// @brief obtain delay
/// @return delay (in seconds)
/// @note an exception is thrown if the value is undefined
double GenericCalInfo::delay() const 
{
  ASKAPCHECK(itsDelayDefined, "Attempting to access delay which has not been defined");
  return itsDelay;
}
  
/// @brief obtain complex gain
/// @return complex gain
/// @note an exception is thrown if the value is undefined
const casa::Complex& GenericCalInfo::gain() const
{
  ASKAPCHECK(itsGainDefined, "Attempting to access complex gain which has not been defined");
  return itsGain;
}
  
  
/// @brief invalidate the values
void GenericCalInfo::invalidate()
{
  itsDelayDefined = false;
  itsGainDefined = false;
}
  
/// @brief set complex gain
/// @param[in] gain complex gain to set
void GenericCalInfo::setGain(const casa::Complex &gain)
{
  itsGain = gain;
  itsGainDefined = true;
}
  
/// @brief set delay
/// @param[in] delay delay in seconds
void GenericCalInfo::setDelay(double delay)
{
  itsDelay = delay;
  itsDelayDefined = true;  
}

} // namespace synthesis

} // namespace askap

