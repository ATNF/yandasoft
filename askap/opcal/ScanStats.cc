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
#include <askap/opcal/ScanStats.h>
#include <askap/dataaccess/IConstDataIterator.h>
#include <askap/dataaccess/TableConstDataIterator.h>

#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>

ASKAP_LOGGER(logger, ".ScanStats");


// std includes
#include <set>
#include <map>
#include <algorithm>

namespace askap {

namespace synthesis {

// helper class to check the beam in ObservationDescription
struct BeamCheckHelper {
  explicit BeamCheckHelper(const casacore::uInt beam) : itsBeam(beam) {}
  bool operator()(const ObservationDescription& obs) const { return obs.beam() == itsBeam; }
private:
  casacore::uInt itsBeam;  
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
  // correctly typed iterator - for specific table-based functions
  const boost::shared_ptr<accessors::TableConstDataIterator> tableIt = iter.dynamicCast<accessors::TableConstDataIterator>();
  if (!tableIt) {
       ASKAPLOG_WARN_STR(logger, "Dataset "<<name<<" appears not to be table-based. Unable to access scan and field information");
  }
  
  const double dirTolerance = 1e-7;
  const double freqTolerance = 1.;
  casacore::uInt cycle = 0;
  for (accessors::IConstDataSharedIter it = iter; it!=it.end(); ++it,++cycle) {
       const casacore::Vector<casacore::uInt>& beam1 = it->feed1();
       const casacore::Vector<casacore::uInt>& beam2 = it->feed2();
       const casacore::Vector<casacore::uInt>& antenna1 = it->antenna1();
       const casacore::Vector<casacore::uInt>& antenna2 = it->antenna2();
       
       // buffer for observations of the current integration, one per beam
       std::map<casacore::uInt, ObservationDescription> currentObs;
       
       ASKAPDEBUGASSERT(it->nRow() == beam1.nelements());
       ASKAPDEBUGASSERT(it->nRow() == beam2.nelements());
       const casacore::uInt centreChan = it->nChannel() / 2;
       ASKAPDEBUGASSERT(it->frequency().nelements() > centreChan);
       const double freq = it->frequency()[centreChan];
       const casacore::Matrix<casacore::Bool> flags = allChannelsFlagged(it->flag());
       const casacore::Vector<casacore::Stokes::StokesTypes> stokes = it->stokes();
       ASKAPDEBUGASSERT(flags.ncolumn() == stokes.nelements());
       ASKAPDEBUGASSERT(flags.nrow() == it->nRow());
       // go through all rows and aggregate baselines
       for (casacore::uInt row=0; row<it->nRow(); ++row) {
            const casacore::uInt beam = beam1[row];
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

                const casacore::Vector<casacore::Stokes::StokesTypes>& stokes = obs.stokes();
                ASKAPCHECK(stokes.nelements() == it->nPol(), "Number of polarisations appears to have changed");
                for (casacore::uInt pol = 0; pol < stokes.nelements(); ++pol) {
                     ASKAPCHECK(stokes[pol] == it->stokes()[pol], "Polarisation product "<<pol+1<<" appears to have changed type");
                }
            } else {
                // brand new element in the map
                obs.set(name, cycle, it->time(), beam, it->pointingDir1()[row],freq);
                ASKAPCHECK(obs.direction().separation(it->pointingDir2()[row]) < dirTolerance, 
                          "Pointing direction is different for two antennas, unsupported scenario");
                if (tableIt) {
                    // fill implementation-specific fields (scan and field IDs)
                    obs.setScanAndFieldIDs(tableIt->currentScanID(), tableIt->currentFieldID());
                }          
                obs.setStokes(it->stokes());
            }        
            // now process flagging information for the given row
            obs.processBaselineFlags(antenna1[row], antenna2[row], flags.row(row));
       } // loop over row
       // aggregate the map into the final buffer
       for (std::map<casacore::uInt, ObservationDescription>::const_iterator ci = currentObs.begin(); ci != currentObs.end(); ++ci) {
            std::vector<ObservationDescription>::reverse_iterator lastObsThisBeam = 
                  std::find_if(itsObs.rbegin(), itsObs.rend(), BeamCheckHelper(ci->first));
            bool sameScan = (lastObsThisBeam != itsObs.rend());
            if (sameScan) {            
                sameScan = (lastObsThisBeam->endCycle() + 1 == ci->second.startCycle());
                sameScan &= (lastObsThisBeam->name() == name);
                if (sameScan && (itsCycleLimit >= 0)) {
                    sameScan = (ci->second.endCycle() - lastObsThisBeam->startCycle() < static_cast<casacore::uInt>(itsCycleLimit));
                } 
                if (sameScan && (itsTimeLimit >= 0.)) {
                    sameScan = (ci->second.endTime() - lastObsThisBeam->startTime() < itsTimeLimit);
                } 
                if (sameScan) {
                    sameScan = (fabs(ci->second.frequency() - freq) < freqTolerance);
                }
                sameScan &= (lastObsThisBeam->scanID() == ci->second.scanID());
                sameScan &= (lastObsThisBeam->fieldID() == ci->second.fieldID());                
                sameScan &= (lastObsThisBeam->direction().separation(ci->second.direction()) < dirTolerance);
            }
            if (sameScan) {
                // continuing the same scan
                ASKAPDEBUGASSERT(lastObsThisBeam != itsObs.rend());
                //lastObsThisBeam->update(ci->second.endCycle(), ci->second.endTime());
                lastObsThisBeam->merge(ci->second);
            } else {
                // start a new scan
                itsObs.push_back(ci->second);
            }
       } // loop over existing scans
       
  }
}

/// @brief helper method to compact flags across frequency axis
/// @details Each row each polarisation is considered flagged if all corresponding frequency channels
/// are flagged
/// @param[in] flags nPol x nChan x nRow cube as provided by the accessor
/// @return nPol x nRow matrix with aggregated flags 
casacore::Matrix<casacore::Bool> ScanStats::allChannelsFlagged(const casacore::Cube<casacore::Bool> &flags)
{
   ASKAPDEBUGASSERT(flags.nelements() > 0);
   casacore::Matrix<casacore::Bool> result(flags.nplane(), flags.nrow(), true);
   for (casacore::uInt row = 0; row < flags.nplane(); ++row) {
        for (casacore::uInt pol = 0; pol < flags.nrow(); ++pol) {
             for (casacore::uInt chan = 0; chan < flags.ncolumn(); ++chan) {
                  if (!flags(pol,chan,row)) {
                      result(row,pol) = false;
                      break;
                  }
             }
        }
   }
   return result;
}

/// @brief access to the selected scan
/// @param[in] index scan index (from 0 to size()-1)
/// @return description of the selected scan
const ObservationDescription& ScanStats::operator[](size_t index) const
{
   ASKAPCHECK(index < size(), "Index "<<index<<" is outside the range of available scans, size="<<size());
   return itsObs[index];
}


} // namespace synthesis

} // namespace askap

