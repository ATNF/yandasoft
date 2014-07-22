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

#include <dataaccess/IConstDataAccessor.h>

#include <string>

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
   /// @details builds an empty map
   //ScanStats();

   
private:
       
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_SCAN_STATS_H


