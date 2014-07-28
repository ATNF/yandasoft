/// @file
///
/// @brief A class to build summary of what is in the dataset
/// @details For operations-specific calibration it is often handy to observe
/// a number of fields in the same measurement set in semi-random order. As a result,
/// to simplify processing, we want to be able to extract this information from the 
/// measurement set itself. The idea is similar to cadvise and VisMetaDataStats, but
/// we don't need a parallel version and don't aggregate statistics in any way. This
/// particular code will not be ran in "online" software and more intended for 
/// offline analysis and operations-specific tools. Note, the "scans" here do not
/// necessarily map to the SCAN_ID in the measurement set. The definition is more like
/// a unique observation (pointing and/or frequency)
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

// own includes
#include <measurementequation/ScanStats.h>
#include <dataaccess/IConstDataIterator.h>

// std includes
#include <set>
#include <map>
#include <algorithm>

namespace askap {

namespace synthesis {

// helper class to check the beam in ObservationDescription
struct BeamCheckHelper {
  explicit BeamCheckHelper(const casa::uInt beam) : itsBeam(beam) {}
  bool operator()(const ObservationDescription& obs) const { return obs.beam() == itsBeam; }
private:
  casa::uInt itsBeam;  
};

/// @brief constructor
/// @details The user can optionally limit the maximum duration of the scan
/// both in cycles and time (the most restrictive limit applies)
/// @param[in] timeLimit maximum duration of observation in time units, negative value means no limit (accessor units)
/// @param[in] cycleLimit maximum duration of observation in cycles, negative value means no limit
ScanStats::ScanStats(double timeLimit, int cycleLimit) : itsTimeLimit(timeLimit), itsCycleLimit(cycleLimit) {} 

/// @brief inspect data referred by a given iterator
/// @details
/// @param[in] name string key identifying the data (e.g. file name)
/// @param[in] iter iterator to work with (note, it is changed in process)
void ScanStats::inspect(const std::string &name, const accessors::IConstDataSharedIter &iter)
{
  ASKAPASSERT(iter);
  const double dirTolerance = 1e-7;
  const double freqTolerance = 1.;
  casa::uInt cycle = 0;
  for (accessors::IConstDataSharedIter it = iter; it!=it.end(); ++it,++cycle) {
       const casa::Vector<casa::uInt>& beam1 = it->feed1();
       const casa::Vector<casa::uInt>& beam2 = it->feed2();
       
       // buffer for observations of the current integration, one per beam
       std::map<casa::uInt, ObservationDescription> currentObs;
       
       ASKAPDEBUGASSERT(it->nRow() == beam1.nelements());
       ASKAPDEBUGASSERT(it->nRow() == beam2.nelements());
       const casa::uInt centreChan = it->nChannel() / 2;
       ASKAPDEBUGASSERT(it->frequency().nelements() > centreChan);
       const double freq = it->frequency()[centreChan];
       // go through all rows and aggregate baselines
       for (casa::uInt row=0; row<it->nRow(); ++row) {
            const casa::uInt beam = beam1[row];
            ASKAPCHECK(beam == beam2[row], "Cross-correlations of different beams are not supported");
            
            // either existing or a brand new element in the map
            ObservationDescription& obs = currentObs[beam];
            
            ASKAPDEBUGASSERT(it->pointingDir1().nelements() > row);
            ASKAPDEBUGASSERT(it->pointingDir2().nelements() > row);
            
            if (obs.isValid()) {
                // existing element - just check for consistency
                ASKAPCHECK(obs.direction().separation(it->pointingDir1()[row]) < dirTolerance, 
                          "Pointing direction is different for different baselines of the same integration, unsupported scenario");
                ASKAPCHECK(obs.direction().separation(it->pointingDir2()[row]) < dirTolerance, 
                          "Pointing direction is different for two antennas, unsupported scenario");                
            } else {
                // brand new element in the map
                obs.set(name, cycle, it->time(), beam, it->pointingDir1()[row],freq);
                ASKAPCHECK(obs.direction().separation(it->pointingDir2()[row]) < dirTolerance, 
                          "Pointing direction is different for two antennas, unsupported scenario");
            }        
       } // loop over row
       // aggregate the map into the final buffer
       for (std::map<casa::uInt, ObservationDescription>::const_iterator ci = currentObs.begin(); ci != currentObs.end(); ++ci) {
            std::vector<ObservationDescription>::reverse_iterator lastObsThisBeam = 
                  std::find_if(itsObs.rbegin(), itsObs.rend(), BeamCheckHelper(ci->first));
            bool sameScan = (lastObsThisBeam != itsObs.rend());
            if (sameScan) {            
                sameScan = (lastObsThisBeam->endCycle() + 1 == ci->second.startCycle());
                sameScan &= (lastObsThisBeam->name() == name);
                if (sameScan && (itsCycleLimit >= 0)) {
                    sameScan = (ci->second.endCycle() - lastObsThisBeam->startCycle() < static_cast<casa::uInt>(itsCycleLimit));
                } 
                if (sameScan && (itsTimeLimit >= 0.)) {
                    sameScan = (ci->second.endTime() - lastObsThisBeam->startTime() < itsTimeLimit);
                } 
                if (sameScan) {
                    sameScan = (fabs(ci->second.frequency() - freq) < freqTolerance);
                } 
                sameScan &= (lastObsThisBeam->direction().separation(ci->second.direction()) < dirTolerance);
            }
            if (sameScan) {
                // continuing the same scan
                ASKAPDEBUGASSERT(lastObsThisBeam != itsObs.rend());
                lastObsThisBeam->update(ci->second.endCycle(), ci->second.endTime());
            } else {
                // start a new scan
                itsObs.push_back(ci->second);
            }
       } // loop over existing scans
       
  }
}


} // namespace synthesis

} // namespace askap

