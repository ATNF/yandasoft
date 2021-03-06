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
#include <casacore/measures/Measures/Stokes.h>
#include <casacore/casa/Arrays/Vector.h>

// std includes
#include <string>
#include <set>
#include <vector>

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
  ObservationDescription(const std::string &name, casacore::uInt cycle, double time, casacore::uInt beam, 
                         const casacore::MVDirection &dir, double freq);
  
  
  /// @brief copy constructor
  /// @details casa arrays use reference semantics, so need a copy constructor
  /// @param[in] src input object
  ObservationDescription(const ObservationDescription &src);

  /// @brief assignment operator
  /// @details casa arrays use reference semantics, so need a copy constructor
  /// @param[in] src input object
  ObservationDescription& operator=(const ObservationDescription &src);

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
  void set(const std::string &name, casacore::uInt cycle, double time, casacore::uInt beam, const casacore::MVDirection &dir, double freq);
  
  /// @brief extend existing observation for additional cycle(s)  
  /// @param[in] cycle current cycle (i.e. iteration number)
  /// @param[in] time current time (in accessor units)
  /// @note an exception is thrown if the structure is uninitialised. The flags are not updated!
  void update(casacore::uInt cycle, double time);

  /// @brief extend exsting observation by merging in another structure
  /// @details Unlike update, it also processes flagging informaton.
  /// @param[in] other structure to merge in
  /// @note an exception is thrown if added structure does not follow in time or cycle
  void merge(const ObservationDescription &other);
  
  /// @brief set scan and field IDs
  /// @details This information can only be obtained for table-based datasets. Therefore, we deal with these misc fields
  /// in a separate method.
  /// @param[in] scanID scan ID to set
  /// @param[in] fieldID field ID to set
  void setScanAndFieldIDs(casacore::uInt scanID, casacore::uInt fieldID);
  
  // getter methods

  /// @brief obtain beam
  /// @return beam number
  /// @note an exception is thrown if the structure is uninitialised
  casacore::uInt beam() const;
  
  /// @brief obtain name
  /// @return file name or other user-defined string key identifying this dataset
  /// @note an exception is thrown if the structure is uninitialised
  const std::string& name() const;
  
  /// @brief obtain start cycle 
  /// @return iteration number of the start of the observation
  /// @note an exception is thrown if the structure is uninitialised
  casacore::uInt startCycle() const;
  
  /// @brief obtain end cycle 
  /// @return iteration number of the end of the observation
  /// @note an exception is thrown if the structure is uninitialised
  casacore::uInt endCycle() const;
  
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
  casacore::uInt scanID() const;

  /// @brief obtain field ID
  /// @return field ID in the appropriate dataset. Zero if the input iterator is not table-based
  /// @note an exception is thrown if the structure is uninitialised
  casacore::uInt fieldID() const;
    
  /// @brief observed direction
  /// @return phase centre of the observation (same frame as used in the accessor)
  /// @note an exception is thrown if the structure is uninitialised 
  const casacore::MVDirection& direction() const;
  
  /// @brief centre frequency
  /// @return frequency corresponding to the centre of the band (in units used in the accessor)
  /// @note an exception is thrown if the structure is uninitialised 
  double frequency() const;


  /// @brief stokes vector
  /// @return vector with polarisation descriptors for each recorded product
  /// @note an exception is thrown if the structure is uninitialised
  const casacore::Vector<casacore::Stokes::StokesTypes>& stokes() const;

  /// @brief initialise the stokes vector
  /// @param[in] stokes vector with polarisation descriptors for each recorded product
  void setStokes(const casacore::Vector<casacore::Stokes::StokesTypes> &stokes);

  /// @brief update antennas with valid data
  /// @details This method processes a flag vector for the given baseline and updates
  ///          the list of antennas with valid data.
  /// @param[in] ant1 first antenna of the baseline
  /// @param[in] ant2 second antenna of the baseline
  /// @param[in] flags per-polarisation flags for the given baseline (should match the 
  ///                  length of stokes vector)
  /// @note This method is not supposed to be called by the reader, it is called when
  ///       the structure is populated
  void processBaselineFlags(casacore::uInt ant1, casacore::uInt ant2, const casacore::Vector<casacore::Bool> &flags);

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
  std::set<casacore::uInt> flaggedAntennas(casacore::Stokes::StokesTypes stokes, casacore::uInt nAnt) const;
  
private:  
  /// @brief file name or other string key defined by the user
  std::string itsName;
  
  /// @brief start cycle (the cycle is just the sequence number of iteration)
  casacore::uInt itsStartCycle;
  
  /// @brief end cycle (the cycle is just the sequence number of iteration)
  casacore::uInt itsEndCycle;
  
  /// @brief start time (same units/frame as used in the accessor to get the data)
  double itsStartTime;
  
  /// @brief end time (same units/frame as used in the accessor to get the data)
  double itsEndTime;
  
  /// @brief beam ID (we treat here different beams as different observations done in parallel)
  /// @note negative value means the class is uninitialised
  int itsBeam;
  
  /// @brief scan ID
  casacore::uInt itsScanID;
  
  /// @brief field ID
  casacore::uInt itsFieldID;
    
  /// @brief phase centre direction (same coordinate system as the accessor is set up with)
  casacore::MVDirection itsDirection;
  
  /// @brief centre frequency (in accessor units). 
  /// @note At this stage we assume that all observations are taken in a standard mode, so no
  /// catering for different bandwidth/resolution, etc. It can be changed in the future if there is
  /// a different use case on the horizon
  double itsFreq;  

  /// @brief measured polarisations
  /// @details The size of the vector is the number of polarisations, content describes frame and which
  /// product is where. Largely added for an extra consistency check - we don't expect this code to be
  /// used for polarisation work
  casacore::Vector<casacore::Stokes::StokesTypes> itsStokes;

  /// @brief set of antennas with valid data for each polarisation
  /// @details If a given antenna has some unflagged data for the given polarisation it will be
  /// added to the apppropriate set. Antennas are referred to by their indices.
  std::vector<std::set<casacore::uInt> > itsAntennasWithValidData;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_OBSERVATION_DESCRIPTION_H


