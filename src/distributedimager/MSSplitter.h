/// @file MSSplitter.h
///
/// @copyright (c) 2012 CSIRO
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

#ifndef ASKAP_CP_MSSPLITTER_H
#define ASKAP_CP_MSSPLITTER_H


// System includes
#include <string>
#include <set>
#include <utility>
#include <stdint.h>

// ASKAPsoft includes
#include "askap/AskapError.h"
#include "askap/AskapLogging.h"
#include "askap/Application.h"
#include "askap/AskapUtil.h"
#include "boost/shared_ptr.hpp"
#include "boost/optional.hpp"
#include "Common/ParameterSet.h"
#include "casacore/casa/aips.h"
#include "casacore/ms/MeasurementSets/MeasurementSet.h"

namespace askap {
namespace cp {


class MSSplitter {
    public:
        /// Constructor
        MSSplitter(LOFAR::ParameterSet& parset);
        /// Entry point method
        int split(const std::string& invis, const std::string& outvis,
              const uint32_t startChan,
              const uint32_t endChan,
              const uint32_t width,
              const LOFAR::ParameterSet& parset);

    

    private:

        static boost::shared_ptr<casa::MeasurementSet> create(
            const std::string& filename, const casa::Bool addSigmaSpec,
            casa::uInt bucketSize, casa::uInt tileNcorr, casa::uInt tileNchan);

        static void copyAntenna(const casa::MeasurementSet& source, casa::MeasurementSet& dest);

        static void copyDataDescription(const casa::MeasurementSet& source, casa::MeasurementSet& dest);

        static void copyFeed(const casa::MeasurementSet& source, casa::MeasurementSet& dest);

        static void copyField(const casa::MeasurementSet& source, casa::MeasurementSet& dest);

        static void copyObservation(const casa::MeasurementSet& source, casa::MeasurementSet& dest);

        static void copyPointing(const casa::MeasurementSet& source, casa::MeasurementSet& dest);

        static void copyPolarization(const casa::MeasurementSet& source, casa::MeasurementSet& dest);

        /// @brief add non-standard column to POINTING table
        /// @details We use 3 non-standard columns to capture
        /// actual pointing on all three axes. This method creates one such
        /// column.
        /// @throws AskapError if column already exist in destPointing
        /// @param[in] name column name
        /// @param[in] srcPointing source MS POINTING table
        /// @param[in] destPointing destination MS POINTING table
        static void addNonStandardPointingColumn(const std::string &name,
                                                 const casa::MSPointing &srcPointing,
                                                 casa::MSPointing &destPointing);

        /// @throws AskapError  if all rows in the main table don't refer to the
        ///                     same spectral window
        /// @return the spectral window id refered to by all rows in the main table,
        ///         or -1 if the main table how no rows;
        static casa::Int findSpectralWindowId(const casa::MeasurementSet& ms);

        /// Writes a new row to the spectral window table of the destination measurement
        /// set which the correct information describing the output spectral window.
        static void splitSpectralWindow(const casa::MeasurementSet& source,
                                 casa::MeasurementSet& dest,
                                 const uint32_t startChan,
                                 const uint32_t endChan,
                                 const uint32_t width,
                                 const casa::Int spwId);

        void splitMainTable(const casa::MeasurementSet& source,
                            casa::MeasurementSet& dest,
                            const uint32_t startChan,
                            const uint32_t endChan,
                            const uint32_t width);

    
        // Returns true if row filtering is enabled, otherwise false.
        bool rowFiltersExist() const;

        // Returns true if the the row should be filtered (i.e excluded), otherwise
        // true.
        bool rowIsFiltered(uint32_t scanid, uint32_t fieldid, uint32_t feed1,
                           uint32_t feed2, double time) const;

        // Helper method for the configuration of the time range filters.
        // Parses the parset value associated with "key" (using MVTime::read()),
        // sets "var" to MVTime::second(), and logs a message "msg".
        // @throws AskapError is thrown if the time string cannot be parsed by
        // MVTime::read()
        void configureTimeFilter(const std::string& key, const std::string& msg,
                                 double& var);

        // Helper method for the configuration of the field name filters.
        // Parses the parset value associated with "fieldnames" and
        // returns all fields in invis with one of these names.
        // @throws AskapError is thrown if none of the field names are present.
        std::vector<uint32_t>
            configureFieldNameFilter(const std::vector<std::string>& names,
                                     const std::string invis);

        /// Set of beam IDs to include in the new measurement set, or empty
        /// if all beams are to be included
        std::set<uint32_t> itsBeams;

        /// Set of scan IDs to include in the new measurement set, or empty
        /// if all scans are to be included
        std::set<uint32_t> itsScans;

        /// Set of fields to include in the new measurement set, or empty
        /// if all scans are to be included
        std::set<uint32_t> itsFieldIds;

        // Optional begin time filter. Rows with TIME < this value will be
        // excluded
        double itsTimeBegin;

        // Optional end time filter. Rows with TIME > this value will be
        // excluded
        double itsTimeEnd;
    
        //Parset - this class came from the MsSplitApp - and needs a config
    LOFAR::ParameterSet& itsParset;
};



}
}
#endif
