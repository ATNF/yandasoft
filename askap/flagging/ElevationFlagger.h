/// @file ElevationFlagger.h
///
/// @copyright (c) 2013,2014 CSIRO
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

#ifndef ASKAP_SYNTHESIS_ELEVATIONFLAGGER_H
#define ASKAP_SYNTHESIS_ELEVATIONFLAGGER_H

// System includes
#include <vector>

// ASKAPsoft includes
#include "Common/ParameterSet.h"
//#include "casacore/casa/Quanta/Quantum.h"

// Local package includes
#include "askap/flagging/IFlagger.h"
#include "askap/flagging/FlaggingStats.h"

namespace askap {
namespace synthesis {

/// @brief Applies flagging based on elevation. This flagger will flag any visibilities
/// where one or both of the antennas have an elevation either lower than the lower threshold
/// or higher than the upper threshold.
class ElevationFlagger : public IFlagger {
    public:

        /// @brief Constructs zero or more instances of the ElevationFlagger.
        /// The flagger is responsible for reading the "parset" and constructing
        /// zero or more instances of itself, depending on the configuration.
        static std::vector<std::shared_ptr<IFlagger> > build(
                const LOFAR::ParameterSet& parset);

        /// @brief Constructor
        ElevationFlagger(const LOFAR::ParameterSet& parset);

        /// @see IFlagger::processRows()
        virtual void processRows(const accessors::IDataSharedIter& di,
                                 const casacore::Vector<bool>& rowFlag,
                                 const casacore::uInt pass, const bool dryRun) override;

        /// @see IFlagger::stats()
        virtual FlaggingStats stats(void) const override;

        /// @see IFlagger::stats()
        virtual casacore::Bool processingRequired(const casacore::uInt pass) const override;

    private:

        // Elevations are cached in "itsAntennaElevations" for a given timestamp
        // (itsTimeElevCalculated). This method updates the elevations.
        void updateElevations(const accessors::IDataSharedIter& di);

        // Utility method to flag the current row.
        void flagRow(casacore::Cube<casacore::Bool>& flag, const casacore::uInt row, const bool dryRun);

        // Flagging statistics
        FlaggingStats itsStats;

        // Flagging threshold. If the elevation of an antenna is
        // larger than this then the row will be flagged.
        casacore::Quantity itsHighLimit;

        // Flagging threshold. If the elevation of an antenna is
        // less than this then the row will be flagged.
        casacore::Quantity itsLowLimit;

        // Timestamp that the antenna elevations std::vector was updated
        casacore::Double itsTimeElevCalculated;

        // Antenna elevations, as calculated at time "itsTimeElevCalculated"
        casacore::Vector<casacore::Quantity> itsAntennaElevations;
};

}
}

#endif
