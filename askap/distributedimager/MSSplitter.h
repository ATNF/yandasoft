/// @file MSSplitter.h
/// MV: Note, this file is not currently used (and excluded from cmake). It is left checked in just in case we need it in the future.
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
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/Application.h>
#include <askap/askap/AskapUtil.h>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <Common/ParameterSet.h>
#include <casacore/casa/aips.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

namespace askap {
namespace cp {


class MSSplitter {
    public:
        /// @brief Constructor - setup tiling parameters and bucket size of the output MS
        /// @param[in] bucketSize bucket size for the output measurement set
        /// @param[in] tileNCorr number of correlations in the tile for the output measurement set
        /// @param[in] tileNChan number of channels in the tile for the output measurement set
        explicit MSSplitter(casacore::uInt bucketSize = 65536u, casacore::uInt tileNCorr = 4u, 
                            casacore::uInt tileNChan = 1u);

        /// @brief Entry point method
        /// @param[in] invis input measurement set name
        /// @param[in] outvis the name of the measurement set to write
        /// @param[in] startChan start channel to extract from the input MS
        /// @param[in] endChan end channel to extract from the input MS
        /// @param[in] width number of consecutive channels to average
        void split(const std::string& invis, const std::string& outvis,
              const uint32_t startChan,
              const uint32_t endChan,
              const uint32_t width);
        
        /// @brief set beam selection
        /// @details By default, all beams are passed. This method can be used to narrow down the selection.
        /// @param[in] beams set of beams to select
        void chooseBeams(const std::set<uint32_t>& beams);

        /// @brief set scan selection
        /// @details By default, all scans are passed. This method can be used to narrow down the selection.
        /// @param[in] scans set of scans to select
        void chooseScans(const std::set<uint32_t>& scans);

        /// @brief configure class from the parset
        /// @deails This is the legacy interface configuring the class from the parset. It is not expected to be
        /// used long term. All usage of the parset is confined to this method.
        /// @param[in] parset configuration parameters
        void configure(const LOFAR::ParameterSet& parset);
    

    private:

        boost::shared_ptr<casacore::MeasurementSet> create(
            const std::string& filename, const casacore::Bool addSigmaSpec);

        static void copyAntenna(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest);

        static void copyDataDescription(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest);

        static void copyFeed(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest);

        static void copyField(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest);

        static void copyObservation(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest);

        static void copyPointing(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest);

        static void copyPolarization(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest);

        /// @brief add non-standard column to POINTING table
        /// @details We use 3 non-standard columns to capture
        /// actual pointing on all three axes. This method creates one such
        /// column.
        /// @throws AskapError if column already exist in destPointing
        /// @param[in] name column name
        /// @param[in] srcPointing source MS POINTING table
        /// @param[in] destPointing destination MS POINTING table
        static void addNonStandardPointingColumn(const std::string &name,
                                                 const casacore::MSPointing &srcPointing,
                                                 casacore::MSPointing &destPointing);

        /// @throws AskapError  if all rows in the main table don't refer to the
        ///                     same spectral window
        /// @return the spectral window id refered to by all rows in the main table,
        ///         or -1 if the main table how no rows;
        static casacore::Int findSpectralWindowId(const casacore::MeasurementSet& ms);

        /// Writes a new row to the spectral window table of the destination measurement
        /// set which the correct information describing the output spectral window.
        static void splitSpectralWindow(const casacore::MeasurementSet& source,
                                 casacore::MeasurementSet& dest,
                                 const uint32_t startChan,
                                 const uint32_t endChan,
                                 const uint32_t width,
                                 const casacore::Int spwId);

        void splitMainTable(const casacore::MeasurementSet& source,
                            casacore::MeasurementSet& dest,
                            const uint32_t startChan,
                            const uint32_t endChan,
                            const uint32_t width);

    
        // Returns true if row filtering is enabled, otherwise false.
        bool rowFiltersExist() const;

        // Returns true if the the row should be filtered (i.e excluded), otherwise
        // true.
        bool rowIsFiltered(uint32_t scanid, uint32_t fieldid, uint32_t feed1,
                           uint32_t feed2, double time) const;

        /*
        // MV - unused method. If we need something like this it should probably goes elsewhere. Otherwise,
        // forces the dependency of this class on the parset (which is a technical debt causing problems elsewhere).
        

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
        */

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
    
        /// @brief bucket size for the output measurement set
        casacore::uInt itsBucketSize;

        /// @brief number of correlations in the tile for the output measurement set
        casacore::uInt itsTileNCorr;
 
        /// @brief number of channels in the tile for the output measurement set
        casacore::uInt itsTileNChan;
};



}
}
#endif
