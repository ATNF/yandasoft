/// @file SelectionFlagger.h
///
/// @copyright (c) 2012-2014 CSIRO
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

#ifndef ASKAP_SYNTHESIS_SELECTIONFLAGGER_H
#define ASKAP_SYNTHESIS_SELECTIONFLAGGER_H

// ASKAPsoft includes
#include "Common/ParameterSet.h"
#include "askap/dataaccess/TableDataAccessor.h"
#include "casacore/tables/TaQL/ExprNode.h"

// Local package includes
#include "askap/flagging/IFlagger.h"
#include "askap/flagging/FlaggingStats.h"

namespace askap {
namespace synthesis {

/// @brief A selection based flagging flagger. This allows flagging based on:
/// - Baseline (i.e. an antenna or a pair of antennas)
/// - Field index number
/// - Time range
/// - Scan index number
/// - Feed/beam index number
/// - UVRange
/// - Autocorrelations only
/// - Spectral (e.g. channel index number or frequency)
class SelectionFlagger : public IFlagger {
    public:

        /// @brief Constructs zero or more instances of the SelectionFlagger.
        /// The flagger is responsible for reading the "parset" and constructing
        /// zero or more instances of itself, depending on the configuration.
        ///
        /// @throw AskapError   If no selection criteria is specified
        ///                     in the parset for a listed rule.
        static std::vector<std::shared_ptr<IFlagger> > build(
                const LOFAR::ParameterSet& parset, const std::string &ms);

        /// @brief Constructor
        /// @throw AskapError   If no selection criteria is specified
        ///                     in the parset for a listed rule.
        SelectionFlagger(const LOFAR::ParameterSet& parset,
                         const std::string &ms);

        /// @see IFlagger::processRows()
        virtual void processRows(const accessors::IDataSharedIter& di,
                                 const casacore::Vector<bool>& rowFlag,
                                 const casacore::uInt pass, const bool dryRun) override;

        /// @see IFlagger::stats()
        virtual FlaggingStats stats(void) const override;

        /// @see IFlagger::stats()
        virtual casacore::Bool processingRequired(const casacore::uInt pass) const override;

    private:
        enum SelectionCriteria {
            BASELINE,
            FIELD,
            TIMERANGE,
            SCAN,
            FEED,
            AUTOCORR,
            TEN
        };

        bool checkBaseline(const accessors::IDataSharedIter& di, const casacore::uInt row);
        bool checkField(const casacore::uInt fieldId);
        bool checkTimerange(const casacore::Double time);
        bool checkScan(const casacore::Int scan);
        bool checkFeed(const accessors::IDataSharedIter& di, const casacore::uInt row) const;
        bool checkAutocorr(const accessors::IDataSharedIter& di, const casacore::uInt row) const;

        bool dispatch(const std::vector<SelectionCriteria>& v,
                      const accessors::IDataSharedIter& di, const casacore::uInt row);

        void checkDetailed(const accessors::IDataSharedIter& di, casacore::Cube<casacore::Bool>& flag, const casacore::uInt row,
                           const bool dryRun);

        // Sets the row flag to true, and also sets the flag true for each visibility
        void flagRow(casacore::Cube<casacore::Bool>& flag, const casacore::uInt row, const bool dryRun);

        // Flagging statistics
        FlaggingStats itsStats;

        // True if auto-correlations should be flagged.
        bool itsFlagAutoCorr;

        // Set to true if per channel or per polarisation product flagging
        // criteria is specified
        bool itsDetailedCriteriaExists;

        // A list indicating which of the row based selection criteria has
        // been specified. The criteria which are more granular than whole row
        // are indicated via the itsDetailedCriteriaExists attribute.
        std::vector<SelectionCriteria> itsRowCriteria;

        // The selections to apply to each row
        casacore::Matrix<casacore::Int> itsBaselines;
        casacore::Vector<casacore::Int> itsFields;
        casacore::Matrix<casacore::Double> itsTimeList;
        casacore::Vector<casacore::Int> itsScans;
        casacore::Matrix<casacore::Int> itsChanList;

        // A set containing the feeds that should be flagged.
        std::set<uint32_t> itsFeedsFlagged;

        // Table row selection expression
        casacore::TableExprNode itsTEN;
};

}
}

#endif
