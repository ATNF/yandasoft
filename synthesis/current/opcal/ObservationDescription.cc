/// @file
///
/// @brief A class representing an observation
/// @details It is intended for use in observation-specific utilities where the
/// measurement set may be quite inhomogeneous. For real astronomy data our
/// approach is to keep flexibility to a minimum, so building this kind of 
/// information is not necessary.
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
/// @author Max Voronkov <maxim.voronkov@csiro.au

#include <opcal/ObservationDescription.h>
#include <askap/AskapError.h>

namespace askap {

namespace synthesis {

/// @brief default constructor
ObservationDescription::ObservationDescription() : itsBeam(-1) {}
  
/// @brief construct as a single cycle observation
/// @details The start cycle and time are set to be equal to end cycle and time, respectively.
/// Data fields are initialised from the parameters of the constructor.
/// @param[in] name file or other key name defined by the user
/// @param[in] cycle current cycle (i.e. iteration number)
/// @param[in] time current time (in accessor units)
/// @param[in] beam beam ID
/// @param[in] dir current direction (in accessor frame)
/// @param[in] freq centre frequency (in accessor units)
ObservationDescription::ObservationDescription(const std::string &name, casa::uInt cycle, double time, 
       casa::uInt beam, const casa::MVDirection &dir, double freq) : itsName(name), itsStartCycle(cycle), itsEndCycle(cycle),
       itsStartTime(time), itsEndTime(time), itsBeam(static_cast<int>(beam)), itsDirection(dir), itsFreq(freq) {}
          
  
/// @brief check whether the structure has been initialised
/// @return true, if the structure has data, false otherwise
/// @note the convention is that non-negative beam ID is a signature of initialised structure
bool ObservationDescription::isValid() const
{
  return itsBeam >= 0;
}

/// @brief update observation parameters for a single cycle observation
/// @param[in] name file or other key name defined by the user
/// @param[in] cycle current cycle (i.e. iteration number)
/// @param[in] time current time (in accessor units)
/// @param[in] beam beam ID
/// @param[in] dir current direction (in accessor frame)
/// @param[in] freq centre frequency (in accessor units)
void ObservationDescription::set(const std::string &name, casa::uInt cycle, double time, casa::uInt beam, 
                                 const casa::MVDirection &dir, double freq)
{
  itsName = name;
  itsStartCycle = itsEndCycle = cycle;
  itsStartTime = itsEndTime = time;
  itsBeam = static_cast<int>(beam);
  itsDirection = dir;
  itsFreq = freq;
}
  
// getter methods

/// @brief obtain beam
/// @return beam number
/// @note an exception is thrown if the structure is uninitialised
casa::uInt ObservationDescription::beam() const
{
  ASKAPCHECK(itsBeam >= 0, "An attempt to get beam for an undefined observation structure");
  return static_cast<casa::uInt>(itsBeam);
}
  
/// @brief obtain name
/// @return file name or other user-defined string key identifying this dataset
/// @note an exception is thrown if the structure is uninitialised
const std::string& ObservationDescription::name() const
{
  ASKAPCHECK(isValid(), "An attempt to get name for an undefined observation structure");
  return itsName;
}
  
/// @brief obtain start cycle 
/// @return iteration number of the start of the observation
/// @note an exception is thrown if the structure is uninitialised
casa::uInt ObservationDescription::startCycle() const
{
  ASKAPCHECK(isValid(), "An attempt to get start cycle for an undefined observation structure");
  return itsStartCycle;
}
  
/// @brief obtain end cycle 
/// @return iteration number of the end of the observation
/// @note an exception is thrown if the structure is uninitialised
casa::uInt ObservationDescription::endCycle() const
{
  ASKAPCHECK(isValid(), "An attempt to get end cycle for an undefined observation structure");
  return itsEndCycle;
}
  
/// @brief obtain start time
/// @return time of the start of this observation (in units set up in accessor)
/// @note an exception is thrown if the structure is uninitialised
double ObservationDescription::startTime() const
{
  ASKAPCHECK(isValid(), "An attempt to get start time for an undefined observation structure");
  return itsStartTime;
}
  
/// @brief obtain end time
/// @return time of the end of this observation (in units set up in accessor)
/// @note an exception is thrown if the structure is uninitialised
double ObservationDescription::endTime() const
{
  ASKAPCHECK(isValid(), "An attempt to get end time for an undefined observation structure");
  return itsEndTime;
}  

/// @brief extend existing observation for additional cycle(s)  
/// @param[in] cycle current cycle (i.e. iteration number)
/// @param[in] time current time (in accessor units)
/// @note an exception is thrown if the structure is uninitialised
void ObservationDescription::update(casa::uInt cycle, double time)
{
  ASKAPCHECK(isValid(), "An attempt to extend an undefined observation");
  ASKAPCHECK(cycle > itsEndCycle, "New cycle number is supposed to be greater than the previous value");
  ASKAPCHECK(time > itsEndTime, "New time stamp is supposed to be later than the previous one");
  itsEndCycle = cycle;
  itsEndTime = time;
}

  
/// @brief observed direction
/// @return phase centre of the observation (same frame as used in the accessor)
/// @note an exception is thrown if the structure is uninitialised 
const casa::MVDirection& ObservationDescription::direction() const
{
  ASKAPCHECK(isValid(), "An attempt to get direction for an undefined observation structure");
  return itsDirection;
}


/// @brief centre frequency
/// @return frequency corresponding to the centre of the band (in units used in the accessor)
/// @note an exception is thrown if the structure is uninitialised 
double ObservationDescription::frequency() const
{
  ASKAPCHECK(isValid(), "An attempt to get frequency for an undefined observation structure");
  return itsFreq;
}

} // namespace synthesis

} // namespace askap

