/// @file MSFlaggingSummary.h
///
/// @copyright (c) 2013 CSIRO
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>

#ifndef ASKAP_SYNTHESIS_MSFLAGGINGSUMMARY_H
#define ASKAP_SYNTHESIS_MSFLAGGINGSUMMARY_H

// ASKAPsoft includes
#include "casacore/ms/MeasurementSets/MSColumns.h"

namespace askap {
namespace synthesis {

/// @brief Writes a Measurement Set flagging summary to the log.
class MSFlaggingSummary {
    public:

        /// @brief Prints a summary of the measurement set to the log
        /// @param[in] msc  Measurement set columns
        static void printToLog(const casacore::MSColumns& msc);

        /// @brief Print summary for a single chunk of data.
        /// This is a utility function used by printToLog()
        /// A chunk is defined as a contiguous series of rows with the same
        /// observation id, scan id, field id, and data description id.
        ///
        /// @param[in] msc  Measurement set columns
        /// @param[in] row  The row to start reading from
        /// @param[in] chunkId The chunk ID (just used for printing)
        ///
        /// @return the row index for "1" row past the end of the current chunk.
        ///         (i.e. the next chunk starts at this row, or this row is
        ///         past the end of the table in the case the last chunk
        ///         was being processed).
        static casacore::uInt summariseChunk(const casacore::MSColumns& msc,
                                         casacore::uInt row, casacore::uInt chunkId);
};

}
}

#endif
