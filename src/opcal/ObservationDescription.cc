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
ObservationDescription::ObservationDescription(const std::string &name, casacore::uInt cycle, double time, 
       casacore::uInt beam, const casacore::MVDirection &dir, double freq) : itsName(name), itsStartCycle(cycle), itsEndCycle(cycle),
       itsStartTime(time), itsEndTime(time), itsBeam(static_cast<int>(beam)), 
       itsScanID(0u), itsFieldID(0u), itsDirection(dir), itsFreq(freq) {}
          
  
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
void ObservationDescription::set(const std::string &name, casacore::uInt cycle, double time, casacore::uInt beam, 
                                 const casacore::MVDirection &dir, double freq)
{
  itsName = name;
  itsStartCycle = itsEndCycle = cycle;
  itsStartTime = itsEndTime = time;
  itsBeam = static_cast<int>(beam);
  itsDirection = dir;
  itsFreq = freq;
}

/// @brief set scan and field IDs
/// @details This information can only be obtained for table-based datasets. Therefore, we deal with these misc fields
/// in a separate method.
/// @param[in] scanID scan ID to set
/// @param[in] fieldID field ID to set
void ObservationDescription::setScanAndFieldIDs(casacore::uInt scanID, casacore::uInt fieldID)
{
  itsScanID = scanID;
  itsFieldID = fieldID;
}

// getter methods

/// @brief obtain beam
/// @return beam number
/// @note an exception is thrown if the structure is uninitialised
casacore::uInt ObservationDescription::beam() const
{
  ASKAPCHECK(itsBeam >= 0, "An attempt to get beam for an undefined observation structure");
  return static_cast<casacore::uInt>(itsBeam);
}
  
/// @brief obtain name
/// @return file name or other user-defined string key identifying this dataset
/// @note an exception is thrown if the structure is uninitialised
const std::string& ObservationDescription::name() const
{
  ASKAPCHECK(isValid(), "An attempt to get name for an undefined observation structure");
  return itsName;
}

/// @brief obtain scan ID
/// @return scan ID in the appropriate dataset. Zero if the input iterator is not table-based
/// @note an exception is thrown if the structure is uninitialised
casacore::uInt ObservationDescription::scanID() const
{
  ASKAPCHECK(isValid(), "An attempt to get scan ID for an undefined observation structure");
  return itsScanID;
}

/// @brief obtain field ID
/// @return field ID in the appropriate dataset. Zero if the input iterator is not table-based
/// @note an exception is thrown if the structure is uninitialised
casacore::uInt ObservationDescription::fieldID() const
{
  ASKAPCHECK(isValid(), "An attempt to get field ID for an undefined observation structure");
  return itsFieldID;
}

  
/// @brief obtain start cycle 
/// @return iteration number of the start of the observation
/// @note an exception is thrown if the structure is uninitialised
casacore::uInt ObservationDescription::startCycle() const
{
  ASKAPCHECK(isValid(), "An attempt to get start cycle for an undefined observation structure");
  return itsStartCycle;
}
  
/// @brief obtain end cycle 
/// @return iteration number of the end of the observation
/// @note an exception is thrown if the structure is uninitialised
casacore::uInt ObservationDescription::endCycle() const
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
void ObservationDescription::update(casacore::uInt cycle, double time)
{
  ASKAPCHECK(isValid(), "An attempt to extend an undefined observation");
  ASKAPCHECK(cycle > itsEndCycle, "New cycle number is supposed to be greater than the previous value");
  ASKAPCHECK(time > itsEndTime, "New time stamp is supposed to be later than the previous one");
  itsEndCycle = cycle;
  itsEndTime = time;
}

/// @brief extend exsting observation by merging in another structure
/// @details Unlike update, it also processes flagging informaton.
/// @param[in] other structure to merge in
/// @note an exception is thrown if added structure does not follow in time or cycle
void ObservationDescription::merge(const ObservationDescription &other)
{
  ASKAPCHECK(isValid(), "An attempt to extend an undefined observation");
  ASKAPCHECK(other.isValid(), "An attempt merge in an undefined observation");
  ASKAPCHECK(other.itsStartCycle > itsEndCycle, "Merged in structure should start with tje cycle number which is greater than the previous value");
  ASKAPCHECK(other.itsStartTime > itsEndTime, "New chunk is supposed to be later in time than the previous one");
  ASKAPCHECK(other.itsStokes.nelements() == itsStokes.nelements(), "Merged observations are supposed to have matching polarisations");
  ASKAPDEBUGASSERT(other.itsEndCycle > itsEndCycle);
  ASKAPDEBUGASSERT(other.itsEndTime > itsEndTime);
  // could've also double check directions, fequencies, etc 
  itsEndCycle = other.itsEndCycle;
  itsEndTime = other.itsEndTime;
  // now aggregate flags
  for (size_t pol=0; pol<itsAntennasWithValidData.size(); ++pol) {
       ASKAPDEBUGASSERT(pol < itsStokes.nelements());
       ASKAPCHECK(itsStokes[pol] == other.itsStokes[pol], "Merged observations are supposed to have matcing polarisations, product "<<pol+1<<" is different");
       itsAntennasWithValidData[pol].insert(other.itsAntennasWithValidData[pol].begin(),other.itsAntennasWithValidData[pol].end());
  }
}

  
/// @brief observed direction
/// @return phase centre of the observation (same frame as used in the accessor)
/// @note an exception is thrown if the structure is uninitialised 
const casacore::MVDirection& ObservationDescription::direction() const
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

/// @brief stokes vector
/// @return vector with polarisation descriptors for each recorded product
/// @note an exception is thrown if the structure is uninitialised
const casacore::Vector<casacore::Stokes::StokesTypes>& ObservationDescription::stokes() const
{
   ASKAPCHECK(isValid(), "An attempt to get stokes vector for an undefined observation structure");
   return itsStokes;
}

/// @brief initialise the stokes vector
/// @param[in] stokes vector with polarisation descriptors for each recorded product
void ObservationDescription::setStokes(const casacore::Vector<casacore::Stokes::StokesTypes> &stokes)
{
   itsStokes.reference(stokes.copy());
   ASKAPCHECK(itsAntennasWithValidData.size() == 0, 
        "Flagging information can only be set after initialisation of the stokes vector");
   itsAntennasWithValidData.resize(itsStokes.nelements());
}

/// @brief update antennas with valid data
/// @details This method processes a flag vector for the given baseline and updates
///          the list of antennas with valid data.
/// @param[in] ant1 first antenna of the baseline
/// @param[in] ant2 second antenna of the baseline
/// @param[in] flags per-polarisation flags for the given baseline (should match the 
///                  length of stokes vector)
/// @note This method is not supposed to be called by the reader, it is called when
///       the structure is populated
void ObservationDescription::processBaselineFlags(casacore::uInt ant1, casacore::uInt ant2, 
                              const casacore::Vector<casacore::Bool> &flags)
{
   // note, a better performance implementation is possible if we iterate over row for
   // each polarisation. But it is not a bottleneck at the moment
   ASKAPCHECK(flags.nelements() == itsStokes.nelements(), "Flag vector should have the same dimension as the Stokes vector. Most likely attempting to process baseline flags before initialisation of Stokes vector");
   ASKAPDEBUGASSERT(flags.nelements() == itsAntennasWithValidData.size());
   for (size_t i=0; i<itsAntennasWithValidData.size(); ++i) {
        if (!flags[i]) {
            // good data for this product
            itsAntennasWithValidData[i].insert(ant1);
            itsAntennasWithValidData[i].insert(ant2);
        }
   }
}

// access to the valid antenna information - we can add other methods if necessary
  
/// @brief obtain a set of flagged antennas
/// @details This method returns a set of antenna indices corresponding to antennas completely
/// flagged in the data 'scan' described by this structure for the given polarisation. It is 
/// handy to have bad antennas listed rather than good ones because by the nature of this tool,
/// little input data should be flagged. The total number of antennas is a parameter (antennas with
/// higher indices may be present but completely flagged). The returned set may contain indices up to
/// the total number of antennas minus 1.
/// @param[in] stokes polarisation product of interest
/// @param[in] nAnt total number of antennas, indices probed go from 0 to nAnt-1
/// @return set flagged antennas represented by their indices
std::set<casacore::uInt> ObservationDescription::flaggedAntennas(casacore::Stokes::StokesTypes stokes, casacore::uInt nAnt) const
{
   ASKAPCHECK(isValid(), "An attempt to get stokes vector for an undefined observation structure");
   ASKAPCHECK(itsAntennasWithValidData.size() > 0, "Perhaps, the Stokes vector is uninitialised");
   ASKAPDEBUGASSERT(itsStokes.nelements() == itsAntennasWithValidData.size());
   size_t pol = 0;
   for (; pol < itsStokes.nelements(); ++pol) {
        if (itsStokes[pol] == stokes) {
            break;
        }
   }
   ASKAPCHECK(pol < itsStokes.nelements(), "Requested stokes type has not been observed");
   
   // now check all antennas from 0 to nAnt-1
   std::set<casacore::uInt> result;
   const std::set<casacore::uInt>& flags = itsAntennasWithValidData[pol];
   for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
        if (flags.find(ant) == flags.end()) {
            result.insert(ant);
        }       
   }
   return result;
}

/// @brief copy constructor
/// @details casa arrays use reference semantics, so need a copy constructor
/// @param[in] src input object
ObservationDescription::ObservationDescription(const ObservationDescription &src) :
    itsName(src.itsName), itsStartCycle(src.itsStartCycle), itsEndCycle(src.itsEndCycle),
    itsStartTime(src.itsStartTime), itsEndTime(src.itsEndTime), itsBeam(src.itsBeam),
    itsScanID(src.itsScanID), itsFieldID(src.itsFieldID), 
    itsDirection(casacore::MVDirection(src.itsDirection)), itsFreq(src.itsFreq),
    itsStokes(src.itsStokes.copy()), itsAntennasWithValidData(src.itsAntennasWithValidData)
{
}

/// @brief assignment operator
/// @details casa arrays use reference semantics, so need a copy constructor
/// @param[in] src input object
ObservationDescription& ObservationDescription::operator=(const ObservationDescription &src)
{
   if (this != &src) {
    itsName=src.itsName;
    itsStartCycle=src.itsStartCycle;
    itsEndCycle=src.itsEndCycle;
    itsStartTime=src.itsStartTime;
    itsEndTime=src.itsEndTime; 
    itsBeam=src.itsBeam;
    itsScanID = src.itsScanID;
    itsFieldID = src.itsFieldID; 
    itsDirection = casacore::MVDirection(src.itsDirection);
    itsFreq = src.itsFreq;
    itsStokes.reference(src.itsStokes.copy());
    itsAntennasWithValidData = src.itsAntennasWithValidData;
   }
   return *this;
}

  

} // namespace synthesis

} // namespace askap

