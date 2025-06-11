/// @file
///
/// FlagParallel: Support for parallel flagging
///
/// @copyright (c) 2023 CSIRO
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
/// @author Mark Wieringa
///

#include <askap/parallel/FlagParallel.h>
#include <askap/askap_synthesis.h>

#include <askap/scimath/fitting/NormalEquationsStub.h>
#include <askap/dataaccess/IFlagDataAccessor.h>
#include <askap/dataaccess/TableDataIterator.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/flagging/FlaggerFactory.h>
#include <askap/flagging/MSFlaggingSummary.h>
#include <askap/utils/TilingUtils.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

// logging stuff
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".FlagParallel");

#include <iomanip>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;
using namespace casacore;

/// @brief Constructor from ParameterSet
/// @details The parset is used to construct the internal state. We could
/// also support construction from a python dictionary (for example).
/// The command line inputs are needed solely for MPI - currently no
/// application specific information is passed on the command line.
/// @param comms communication object
/// @param parset ParameterSet for inputs
FlagParallel::FlagParallel(askap::askapparallel::AskapParallel& comms,
      const LOFAR::ParameterSet& parset) : MEParallelApp(comms,parset,true)
{
  // the stub allows to reuse MEParallelApp code although we're not solving
  // for the normal equations here
  itsNe.reset(new scimath::NormalEquationsStub);
}

void FlagParallel::flagOne(const std::string &ms, bool distributeByTile)
{
    // Print a summary if needed
    if (itsComms.isMaster() && parset().getBool("summary", true)) {
        MeasurementSet mset(ms);
        MSFlaggingSummary::printToLog(mset);
    }
    // don't start flagging before summary is done
    itsComms.barrier();

    // Is this a dry run?
    const bool dryRun = parset().getBool("dryrun", false);
    if (dryRun) {
        ASKAPLOG_INFO_STR(logger, "!!!!! DRY RUN ONLY - MeasurementSet will not be updated !!!!!");
    }
    // Create a vector of all the flagging strategies specified in the parset
    itsFlaggers = FlaggerFactory::build(parset(), ms);
    ASKAPCHECK(!itsFlaggers.empty(), "No flaggers configured - Aborting");


    // Open readonly, accessor will reopen table r/w when needed
    TableDataSource ds(ms, TableDataSource::MEMORY_BUFFERS | TableDataSource::WRITE_DATA_ONLY,
         dataColumn());
    ds.configureUVWMachineCache(uvwMachineCacheSize(),uvwMachineCacheTolerance());
    IDataSelectorPtr sel=ds.createSelector();
    if (distributeByTile) {
        utils::distributeByTile(sel, dataColumn(),nWorkers(),workerRank());
    }
    sel << parset();
    IDataConverterPtr conv=ds.createConverter();
    conv->setEpochFrame();
    conv->setFrequencyFrame(getFreqRefFrame(), "Hz");
    conv->setDirectionFrame(MDirection::Ref(MDirection::J2000));
    IDataSharedIter dataIt = ds.createIterator(sel, conv);

    rownr_t rowsAlreadyFlagged = 0;
    rownr_t nRows = 0;
    Bool passRequired = True;
    uInt pass = 0;
    while (passRequired) {
        for (dataIt.init(); dataIt.hasMore(); dataIt.next()) {
            const Cube<Bool>& flag = dataIt->flag();
            rownr_t nRow = flag.nplane();

            // Count flagged rows and keep a list
            // (accessor reads "FLAG_ROW" column and applies it to flag)
            uInt flagged = 0;
            Vector<Bool> rowFlag(nRow, False);
            for (rownr_t j = 0; j < nRow ; j++) {
                rowFlag(j) = allEQ(flag.xyPlane(j),True);
                if (rowFlag(j)) {
                    flagged++;
                }
            }
            if (pass == 0) {
                nRows += nRow;
                rowsAlreadyFlagged += flagged;
            }
            // If there are unflagged rows, do more flagging
            if (flagged < nRow) {
                // Invoke each flagger
                for (auto it : itsFlaggers) {
                    if (it->processingRequired(pass)) {
                        it->processRows(dataIt, rowFlag, pass, dryRun);
                    }
                }
            }
        }
        pass++;
        passRequired = False;
        for (auto it : itsFlaggers) {
            if (it->processingRequired(pass)) {
                passRequired = True;
            }
        }
    }

    // Write out flagging statistics
    ASKAPLOG_INFO_STR(logger, "Summary:");
    float rowPercent = static_cast<float>(rowsAlreadyFlagged) / nRows * 100.0;
    ASKAPLOG_INFO_STR(logger, "  Rows already flagged: " << rowsAlreadyFlagged
            << " (" << std::setprecision(3) << rowPercent << "%)");
    for (auto it : itsFlaggers) {
        const FlaggingStats stats = it->stats();
        rowPercent = static_cast<float>(stats.rowsFlagged) / nRows * 100.0;
        ASKAPDEBUGASSERT(rowPercent <= 100.0);
        ASKAPLOG_INFO_STR(logger, "  " << stats.name
                              << " - Entire rows flagged: " << stats.rowsFlagged
                              << " (" << std::setprecision(3) << rowPercent << "%)"
                              << ", Visibilities flagged: " << stats.visFlagged);
    }

    //stats.logSummary();
    //RODataManAccessor(ms, "TiledData", False).showCacheStatistics (cout);
}

/// @brief perform the subtraction
/// @details This method iterates over one or more datasets, flagging visibilities according to
/// the parset. In parallel mode with a single dataset we can optionally distribute work over tiles.
void FlagParallel::doFlag()
{
    if (itsComms.isParallel()) {
        if (doWork()) {
            const uint rank = workerRank();
            ASKAPLOG_INFO_STR(logger, "Worker "<<rank<< " is processing "<<measurementSets()[rank]);
            // do automatic distribution over tiles if requested and
            //   if all measurementset names are the same
            const std::vector<std::string>& v = measurementSets();
            const bool allEqual = (std::adjacent_find(v.begin(), v.end(),
                std::not_equal_to<std::string>()) == v.end());
            const bool distributeByTile = parset().isDefined("Tiles") &&
                parset().getString("Tiles")=="auto" && allEqual;
            if (allEqual && !parset().isDefined("Tiles")) {
                ASKAPLOG_WARN_STR(logger, "Multiple ranks specified with single dataset, but no distribution by Tiles - ranks may be doing the same work");
            }
            flagOne(measurementSets()[rank], distributeByTile);
        }
    } else {
        for (size_t iMs=0; iMs<measurementSets().size(); ++iMs) {
            flagOne(measurementSets()[iMs]);
        }
    }
}
