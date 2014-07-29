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

#ifndef SYNTHESIS_SCAN_STATS_H
#define SYNTHESIS_SCAN_STATS_H

#include <dataaccess/SharedIter.h>
#include <opcal/ObservationDescription.h>

#include <string>
#include <vector>

namespace askap {

namespace synthesis {


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
/// @ingroup measurementequation
class ScanStats {
public:

   /// @brief constructor
   /// @details The user can optionally limit the maximum duration of the scan
   /// both in cycles and time (the most restrictive limit applies)
   /// @param[in] timeLimit maximum duration of observation in time units, negative value means no limit (accessor units)
   /// @param[in] cycleLimit maximum duration of observation in cycles, negative value means no limit
   explicit ScanStats(double timeLimit = -1, int cycleLimit = -1);

   // iterator type
   typedef std::vector<ObservationDescription>::const_iterator const_iterator;
   
   // access iterators
   
   /// @brief start iterator
   /// @return iterator pointing to the first element
   inline const_iterator begin() const {return itsObs.begin();}
   
   /// @brief end iterator
   /// @return iterator pointing to the element after last
   inline const_iterator end() const  {return itsObs.end();}
   
   /// @brief number of observations
   /// @return total number of observations
   /// @note each beam is present as a separate observation as beams can have different phase centres
   inline size_t size() const {return itsObs.size();}
   
   /// @brief inspect data referred by a given iterator
   /// @details
   /// @param[in] name string key identifying the data (e.g. file name)
   /// @param[in] iter iterator to work with (note, it is changed in process)
   void inspect(const std::string &name, const accessors::IConstDataSharedIter &iter);
   
private:
   /// collection of individual observations
   std::vector<ObservationDescription> itsObs;    

   /// @brief time limit on the scan duration, a negative value means no limit is imposed (accessor units)
   double itsTimeLimit;
   
   /// @brief limit on the scan duration in cycles, a negative value means no limit is imposed (accessor units)
   int itsCycleLimit;
   
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_SCAN_STATS_H


