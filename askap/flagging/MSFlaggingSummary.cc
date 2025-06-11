/// @file MSFlaggingSummary.cc
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

// Include own header file first
#include "askap/flagging/MSFlaggingSummary.h"

// Include package level header file
#include "askap/askap_synthesis.h"

// ASKAPsoft includes
#include "askap/askap/AskapLogging.h"
#include "casacore/ms/MeasurementSets/MSColumns.h"
#include "casacore/measures/Measures/Stokes.h"

#include <iomanip>

// Using
using namespace std;
using namespace askap;
using namespace askap::synthesis;
using namespace casacore;

ASKAP_LOGGER(logger, ".MSFlaggingSummary");

casacore::uInt MSFlaggingSummary::summariseChunk(const casacore::MSColumns& msc, casacore::uInt start, casacore::uInt chunkId)
{
    // These attributes are used to delineate chunks (i.e. if they change
    // it is the end of this chunk)
    const casacore::Int scanId = msc.scanNumber()(start);
    const casacore::Int obsId = msc.observationId()(start);
    const casacore::Int dataDescId = msc.dataDescId()(start);
    const casacore::Int fieldId = msc.fieldId()(start);

    //
    const casacore::ROMSDataDescColumns& ddc = msc.dataDescription();
    const casacore::ROMSPolarizationColumns& polc = msc.polarization();
    const casacore::ROMSFieldColumns& fieldc = msc.field();
    const casacore::ROMSSpWindowColumns& spwc = msc.spectralWindow();

    const casacore::uInt descPolId = ddc.polarizationId()(dataDescId);
    const casacore::uInt descSpwId = ddc.spectralWindowId()(dataDescId);

    // Stats to capture
    uint64_t nRowsFlagged = 0;
    uint64_t nVisFlagged = 0;
    uint64_t nVis= 0;
    set<casacore::Double> timeset;
    set< pair<casacore::Int, casacore::Int> > baselineset;
    set<casacore::Int> feedset;

    // Process rows until a new chunk is found or the table ends
    casacore::uInt row = start;
    while (row < msc.nrow()
            && scanId == msc.scanNumber()(row)
            && obsId == msc.observationId()(row)
            && dataDescId == msc.dataDescId()(row)
            && fieldId == msc.fieldId()(row)) {

        // Count times
        timeset.insert(msc.time()(row));

        // Count baselines seen
        baselineset.insert(make_pair(msc.antenna1()(row), msc.antenna2()(row)));

        // Count feeds seen
        feedset.insert(msc.feed1()(row));
        feedset.insert(msc.feed2()(row));

        // Row flagging
        if (msc.flagRow()(row)) nRowsFlagged++;

        // Visibility flagging
        const Matrix<casacore::Bool> flags = msc.flag()(row);
        nVis += flags.size();
        for (casacore::uInt m = 0; m < flags.nrow(); ++m) {
            for (casacore::uInt n = 0; n < flags.ncolumn(); ++n) {
                if (flags(m, n)) nVisFlagged++;
            }
        }

        ++row;
    }


    // Print chunk details
    ASKAPLOG_INFO_STR(logger, "Chunk " << chunkId
            << " (ObsID: " << obsId
            << ", ScanID: " << scanId
            << ", FieldID: " << fieldId
            << ", Field Name: " << fieldc.name()(fieldId)
            << ", Spectral Window: " << descSpwId
            << ")");

    // Build a string array from stokes types
    const casacore::Vector<casacore::Int> stokesTypesInt = polc.corrType()(descPolId);
    string stokesList = "[";
    for (size_t i = 0; i < stokesTypesInt.size(); ++i) {
        if (i != 0) stokesList += ", ";
        stokesList += Stokes::name(Stokes::type(stokesTypesInt(i)));
    }
    stokesList += "]";

    const unsigned int nChan = spwc.numChan()(descSpwId);

    ASKAPLOG_INFO_STR(logger, stokesList << ", " << nChan << " channels, "
            << timeset.size() << " times, "
            << baselineset.size() << " baselines, "
            << feedset.size() << " beams, "
            << row - start << " rows");

    const float rowFlagPercent = static_cast<float>(nRowsFlagged) / (row - start) * 100.0;
    ASKAPLOG_INFO_STR(logger, nRowsFlagged << " out of " << row - start << " ("
            << setprecision(3) << rowFlagPercent << "%) rows are flagged");
    const float visFlagPercent = static_cast<float>(nVisFlagged) / nVis * 100.0;
    ASKAPLOG_INFO_STR(logger, nVisFlagged << " out of " << nVis << " ("
            << setprecision(3) << visFlagPercent << "%) visibilities are flagged");

    return row;
}

void MSFlaggingSummary::printToLog(const casacore::MSColumns& msc)
{
        ASKAPLOG_INFO_STR(logger, "Pre-flagging Measurement Set Summary:");

        const casacore::uInt nrow = msc.nrow();
        if (nrow == 0) {
            ASKAPLOG_INFO_STR(logger, "No rows");
            return;
        }

        // Print a summary for each chunk of data
        ASKAPDEBUGASSERT(nrow > 0);
        casacore::uInt row = 0;
        casacore::uInt chunkId = 1;
        const string line = "-------------------------------------------------------------------------------";
            ASKAPLOG_INFO_STR(logger, line);
        do {
            row = summariseChunk(msc, row, chunkId);
            ASKAPLOG_INFO_STR(logger, line);
        } while (row < nrow);
}
