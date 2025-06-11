/// @file SelectionFlagger.cc
///
/// @copyright (c) 2011-2014 CSIRO
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

// Include own header file first
#include "SelectionFlagger.h"

// Include package level header file
#include "askap/askap_synthesis.h"

// ASKAPsoft includes
#include "askap/askap/AskapLogging.h"

#include "casacore/ms/MSSel/MSSelection.h"

ASKAP_LOGGER(logger, ".SelectionFlagger");

using namespace askap;
using namespace casacore;
using namespace askap::synthesis;
using namespace askap::accessors;

std::vector<std::shared_ptr<IFlagger> > SelectionFlagger::build(
        const LOFAR::ParameterSet& parset,
        const std::string &ms)
{
    vector<std::shared_ptr<IFlagger> > flaggers;
    const string key = "selection_flagger.rules";
    if (parset.isDefined(key)) {
        const vector<string> rules = parset.getStringVector(key);
        vector<string>::const_iterator it;
        for (it = rules.begin(); it != rules.end(); ++it) {
            ASKAPLOG_DEBUG_STR(logger, "Processing rule: " << *it);
            string s = "selection_flagger." + *it + ".";
            const LOFAR::ParameterSet subset = parset.makeSubset(s);
            flaggers.push_back(std::shared_ptr<IFlagger>(new SelectionFlagger(subset, ms)));
        }
    }
    return flaggers;
}

SelectionFlagger:: SelectionFlagger(const LOFAR::ParameterSet& parset,
                                    const std::string &ms)
        : itsStats("SelectionFlagger"), itsFlagAutoCorr(false),
        itsDetailedCriteriaExists(false)
{
    casacore::MeasurementSet myMS;
    if (ms.size()>0) {
        myMS = casacore::MeasurementSet(ms);
    }
    casacore::MSSelection mySelection;
    mySelection.resetMS(myMS);
    if (parset.isDefined("field")) {
        mySelection.setFieldExpr(parset.getString("field"));
        itsRowCriteria.push_back(FIELD);
    }

    if (parset.isDefined("spw")) {
        mySelection.setSpwExpr(parset.getString("spw"));
        itsDetailedCriteriaExists = true;
    }

    if (parset.isDefined("antenna")) {
        mySelection.setAntennaExpr(parset.getString("antenna"));
        itsRowCriteria.push_back(BASELINE);
    }

    if (parset.isDefined("timerange")) {
        mySelection.setTimeExpr(parset.getString("timerange"));
        itsRowCriteria.push_back(TIMERANGE);
    }

    if (parset.isDefined("correlation")) {
        mySelection.setPolnExpr(parset.getString("correlation"));
        ASKAPTHROW(AskapError, "Correlation selection not yet implemented");
        itsDetailedCriteriaExists = true;
    }

    if (parset.isDefined("scan")) {
        mySelection.setScanExpr(parset.getString("scan"));
        itsRowCriteria.push_back(SCAN);
    }

    if (parset.isDefined("feed")) {
        const std::vector<uint32_t> v = parset.getUint32Vector("feed");
        itsFeedsFlagged.insert(v.begin(), v.end());
        itsRowCriteria.push_back(FEED);
    }

    if (parset.isDefined("uvrange")) {
        mySelection.setUvDistExpr(parset.getString("uvrange"));

        // Specifying a uvrange results in row selection
        itsRowCriteria.push_back(TEN);

        // Create a table expression over a MS representing the selection
        itsTEN = mySelection.toTableExprNode(&myMS);
    }

    if (parset.isDefined("autocorr")) {
        itsFlagAutoCorr = parset.getBool("autocorr");
        if (itsFlagAutoCorr) {
            itsRowCriteria.push_back(AUTOCORR);
        }
    }

    if (itsRowCriteria.empty() && !itsDetailedCriteriaExists) {
        ASKAPTHROW(AskapError, "No selection criteria for rule specified");
    }
    // Get all the selection info we need
    if (!myMS.isNull()) {
        itsBaselines = mySelection.getBaselineList();
        itsFields = mySelection.getFieldList();
        itsTimeList = mySelection.getTimeList();
        itsScans = mySelection.getScanList();
        itsChanList = mySelection.getChanList();
    }
}

FlaggingStats SelectionFlagger::stats(void) const
{
    return itsStats;
}

casacore::Bool SelectionFlagger::processingRequired(const casacore::uInt pass) const
{
    return (pass==0);
}

void SelectionFlagger::processRows(const IDataSharedIter& di,
    const casacore::Vector<bool>& rowFlag,
    const casacore::uInt pass, const bool dryRun)
{
    IFlagDataAccessor &fda= dynamic_cast<IFlagDataAccessor&>(*di);
    casacore::Cube<casacore::Bool> flag(dryRun ? fda.flag() : fda.rwFlag());
    for (uint row = 0; row < rowFlag.size(); row++) {
        if (!rowFlag(row)) {
            const bool rowCriteriaMatches = dispatch(itsRowCriteria, di, row);

            // 1: Handle the case where all row criteria match and no detailed criteria
            // exists
            if (rowCriteriaMatches && !itsDetailedCriteriaExists) {
                flagRow(flag, row, dryRun);
            }

            // 2: Handle the case where there is no row criteria, but there is detailed
            // criteria. Or, where the row criteria exists and match.
            if (itsDetailedCriteriaExists &&
                (itsRowCriteria.empty() || rowCriteriaMatches) ) {
                checkDetailed(di, flag, row, dryRun);
            }
        }
    }
}

bool SelectionFlagger::checkBaseline(const IDataSharedIter& di, const casacore::uInt row)
{
    const casacore::Matrix<casacore::Int>& m = itsBaselines;
    if (m.empty()) {
        return false;
    }
    ASKAPCHECK(m.ncolumn() == 2, "Expected two columns");

    const casacore::Int ant1 = di->antenna1()(row);
    const casacore::Int ant2 = di->antenna2()(row);
    for (size_t i = 0; i < m.nrow(); ++i) {
        if ((m(i, 0) == ant1 && m(i, 1) == ant2)
                || (m(i, 0) == ant2 && m(i, 1) == ant1)) {
            return true;
        }
    }

    return false;
}

bool SelectionFlagger::checkField(const casacore::uInt fieldId)
{
    const casacore::Vector<casacore::Int> & v = itsFields;
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] == fieldId) {
            return true;
        }
    }
    return false;
}

bool SelectionFlagger::checkTimerange(const casacore::Double time)
{
    const casacore::Matrix<casacore::Double>& timeList = itsTimeList;
    if (timeList.empty()) {
        ASKAPLOG_DEBUG_STR(logger, "Time list is EMPTY");
        return false;
    }
    ASKAPCHECK(timeList.nrow() == 2, "Expected two rows");
    ASKAPCHECK(timeList.ncolumn() == 1,
            "Only a single time range specification is supported");
    if (time > timeList(0, 0) && time < timeList(1, 0)) {
        return true;
    } else {
        return false;
    }
}

bool SelectionFlagger::checkScan(casacore::Int scan)
{
    const casacore::Vector<casacore::Int>& v = itsScans;
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] == scan) {
            return true;
        }
    }
    return false;
}

bool SelectionFlagger::checkFeed(const IDataSharedIter& di, const casacore::uInt row) const
{
    const casacore::Int feed1 = di->feed1()(row);
    const casacore::Int feed2 = di->feed2()(row);

    if ((itsFeedsFlagged.find(feed1) != itsFeedsFlagged.end())
            || (itsFeedsFlagged.find(feed2) != itsFeedsFlagged.end())) {
        return true;
    } else {
        return false;
    }
}

bool SelectionFlagger::checkAutocorr(const IDataSharedIter& di, const casacore::uInt row) const
{
    ASKAPDEBUGASSERT(itsFlagAutoCorr);

    const casacore::Int ant1 = di->antenna1()(row);
    const casacore::Int ant2 = di->antenna2()(row);
    return (ant1 == ant2);
}

bool SelectionFlagger::dispatch(const std::vector<SelectionCriteria>& v,
                                const IDataSharedIter& di, const casacore::uInt row)
{
    std::vector<SelectionCriteria>::const_iterator it;
    auto tdi = di.dynamicCast<TableConstDataIterator>();

    for (it = v.begin(); it != v.end(); ++it) {
        switch (*it) {
            case SelectionFlagger::BASELINE:
                if (!checkBaseline(di, row)) return false;
                break;
            case SelectionFlagger::FIELD:
                if (!checkField(tdi->currentFieldID())) return false;
                break;
            case SelectionFlagger::TIMERANGE:
                if (!checkTimerange(tdi->getTime())) return false;
                break;
            case SelectionFlagger::SCAN:
                if (!checkScan(tdi->currentScanID())) return false;
                break;
            case SelectionFlagger::FEED:
                if (!checkFeed(di, row)) return false;
                break;
            case SelectionFlagger::AUTOCORR:
                if (!checkAutocorr(di, row)) return false;
                break;
            case SelectionFlagger::TEN:
                if (!tdi->isSelected(itsTEN,row)) return false;
                break;
            default:
                break;
        }
    }
    return true;
}

void SelectionFlagger::checkDetailed(const IDataSharedIter& di,
    casacore::Cube<casacore::Bool>& flag, const casacore::uInt row, const bool dryRun)
{
    const casacore::Matrix<casacore::Int>& chanList = itsChanList;
    if (chanList.empty()) {
        ASKAPLOG_DEBUG_STR(logger, "Channel flagging list is EMPTY");
        return;
    }
    ASKAPCHECK(chanList.ncolumn() == 4, "Expected four columns");
    auto tdi = di.dynamicCast<TableConstDataIterator>();
    const casacore::Int curSpwId = tdi->currentSpWindowID();
    //ASKAPLOG_DEBUG_STR(logger, "Channel flagging list size: " << chanList.nrow());
    for (size_t i = 0; i < chanList.nrow(); ++i) {
        const casacore::Int spwID = chanList(i, 0);
        const casacore::Int startCh = chanList(i, 1);
        const casacore::Int stopCh = chanList(i, 2);
        const casacore::Int step = chanList(i, 3);
        //ASKAPLOG_DEBUG_STR(logger, "spwID: " << spwID
        //                       << ", startCh: " << startCh
        //                       << ", stopCh: " << stopCh
        //                       << ", step: " << step);
        ASKAPCHECK(step > 0, "Step must be greater than zero to avoid infinite loop");
        if (curSpwId != spwID) {
            continue;
        }

        if (!dryRun) {
            flag(Slice(), Slice(startCh,stopCh,1,false), row) = true;
        }
        itsStats.visFlagged += (stopCh-startCh+1)*flag.nrow();
    }
}

void SelectionFlagger::flagRow(casacore::Cube<casacore::Bool>& flag, const casacore::uInt row, const bool dryRun)
{

    itsStats.visFlagged += flag.shape()(1) * flag.shape()(0);
    itsStats.rowsFlagged++;

    if (!dryRun) {
        flag(casacore::Slice(),casacore::Slice(),casacore::Slice(row)) = true;
    }
}
