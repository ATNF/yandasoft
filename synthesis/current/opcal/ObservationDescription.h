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

#ifndef SYNTHESIS_OBSERVATION_DESCRIPTION_H
#define SYNTHESIS_OBSERVATION_DESCRIPTION_H

// casa includes
#include <casacore/measures/Measures/MDirection.h>

// std includes
#include <string>

namespace askap {

namespace synthesis {

/// @brief description of a single observation
struct ObservationDescription {

  /// @brief default constructor
  ObservationDescription();
  
  /// @brief construct as a single cycle observation
  /// @details The start cycle and time are set to be equal to end cycle and time, respectively.
  /// Data fields are initialised from the parameters of the constructor.
  /// @param[in] name file or other key name defined by the user
  /// @param[in] cycle current cycle (i.e. iteration number)
  /// @param[in] time current time (in accessor units)
  /// @param[in] beam beam ID
  /// @param[in] dir current direction (in accessor frame)
  /// @param[in] freq centre frequency (in accessor units)
  ObservationDescription(const std::string &name, casa::uInt cycle, double time, casa::uInt beam, 
                         const casa::MVDirection &dir, double freq);
  
  /// @brief check whether the structure has been initialised
  /// @return true, if the structure has data, false otherwise
  /// @note the convention is that non-negative beam ID is a signature of initialised structure
  bool isValid() const;

  /// @brief update observation parameters for a single cycle observation
  /// @param[in] name file or other key name defined by the user
  /// @param[in] cycle current cycle (i.e. iteration number)
  /// @param[in] time current time (in accessor units)
  /// @param[in] beam beam ID
  /// @param[in] dir current direction (in accessor frame)
  /// @param[in] freq centre frequency (in accessor units)
  void set(const std::string &name, casa::uInt cycle, double time, casa::uInt beam, const casa::MVDirection &dir, double freq);
  
  /// @brief extend existing observation for additional cycle(s)  
  /// @param[in] cycle current cycle (i.e. iteration number)
  /// @param[in] time current time (in accessor units)
  /// @note an exception is thrown if the structure is uninitialised
  void update(casa::uInt cycle, double time);
  
  /// @brief set scan and field IDs
  /// @details This information can only be obtained for table-based datasets. Therefore, we deal with these misc fields
  /// in a separate method.
  /// @param[in] scanID scan ID to set
  /// @param[in] fieldID field ID to set
  void setScanAndFieldIDs(casa::uInt scanID, casa::uInt fieldID);
  
  // getter methods

  /// @brief obtain beam
  /// @return beam number
  /// @note an exception is thrown if the structure is uninitialised
  casa::uInt beam() const;
  
  /// @brief obtain name
  /// @return file name or other user-defined string key identifying this dataset
  /// @note an exception is thrown if the structure is uninitialised
  const std::string& name() const;
  
  /// @brief obtain start cycle 
  /// @return iteration number of the start of the observation
  /// @note an exception is thrown if the structure is uninitialised
  casa::uInt startCycle() const;
  
  /// @brief obtain end cycle 
  /// @return iteration number of the end of the observation
  /// @note an exception is thrown if the structure is uninitialised
  casa::uInt endCycle() const;
  
  /// @brief obtain start time
  /// @return time of the start of this observation (in units set up in accessor)
  /// @note an exception is thrown if the structure is uninitialised
  double startTime() const;
  
  /// @brief obtain end time
  /// @return time of the end of this observation (in units set up in accessor)
  /// @note an exception is thrown if the structure is uninitialised
  double endTime() const;

  /// @brief obtain scan ID
  /// @return scan ID in the appropriate dataset. Zero if the input iterator is not table-based
  /// @note an exception is thrown if the structure is uninitialised
  casa::uInt scanID() const;

  /// @brief obtain field ID
  /// @return field ID in the appropriate dataset. Zero if the input iterator is not table-based
  /// @note an exception is thrown if the structure is uninitialised
  casa::uInt fieldID() const;
    
  /// @brief observed direction
  /// @return phase centre of the observation (same frame as used in the accessor)
  /// @note an exception is thrown if the structure is uninitialised 
  const casa::MVDirection& direction() const;
  
  /// @brief centre frequency
  /// @return frequency corresponding to the centre of the band (in units used in the accessor)
  /// @note an exception is thrown if the structure is uninitialised 
  double frequency() const;
  
private:  
  /// @brief file name or other string key defined by the user
  std::string itsName;
  
  /// @brief start cycle (the cycle is just the sequence number of iteration)
  casa::uInt itsStartCycle;
  
  /// @brief end cycle (the cycle is just the sequence number of iteration)
  casa::uInt itsEndCycle;
  
  /// @brief start time (same units/frame as used in the accessor to get the data)
  double itsStartTime;
  
  /// @brief end time (same units/frame as used in the accessor to get the data)
  double itsEndTime;
  
  /// @brief beam ID (we treat here different beams as different observations done in parallel)
  /// @note negative value means the class is uninitialised
  int itsBeam;
  
  /// @brief scan ID
  casa::uInt itsScanID;
  
  /// @brief field ID
  casa::uInt itsFieldID;
    
  /// @brief phase centre direction (same coordinate system as the accessor is set up with)
  casa::MVDirection itsDirection;
  
  /// @brief centre frequency (in accessor units). 
  /// @note At this stage we assume that all observations are taken in a standard mode, so no
  /// catering for different bandwidth/resolution, etc. It can be changed in the future if there is
  /// a different use case on the horizon
  double itsFreq;  
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_OBSERVATION_DESCRIPTION_H


