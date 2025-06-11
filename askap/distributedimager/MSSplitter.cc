/// @file MSSplitter.cc
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


// Include own header file
#include "MSSplitter.h"

// System includes
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <stdint.h>

// ASKAPsoft includes
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/Application.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askap/StatReporter.h>
#include <askap/askap/Log4cxxLogSink.h>
#include <boost/shared_ptr.hpp>
#include <Common/ParameterSet.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/DataMan/IncrementalStMan.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/tables/DataMan/TiledShapeStMan.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>


ASKAP_LOGGER(logger, ".mssplitter");

using namespace askap;
using namespace casa;
using namespace std;
using namespace cp;

/// @brief Constructor - setup tiling parameters and bucket size of the output MS
/// @param[in] bucketSize bucket size for the output measurement set
/// @param[in] tileNCorr number of correlations in the tile for the output measurement set
/// @param[in] tileNChan number of channels in the tile for the output measurement set
MSSplitter::MSSplitter(casacore::uInt bucketSize, casacore::uInt tileNCorr, casacore::uInt tileNChan)
    : itsTimeBegin(std::numeric_limits<double>::min()),
    itsTimeEnd(std::numeric_limits<double>::max()),
    itsBucketSize(bucketSize), itsTileNCorr(tileNCorr), itsTileNChan(tileNChan) {}

/// @brief set beam selection
/// @details By default, all beams are passed. This method can be used to narrow down the selection.
/// @param[in] beams set of beams to select
void MSSplitter::chooseBeams(const std::set<uint32_t>& beams)
{
   itsBeams = beams;
}

/// @brief set scan selection
/// @details By default, all scans are passed. This method can be used to narrow down the selection.
/// @param[in] scans set of scans to select
void MSSplitter::chooseScans(const std::set<uint32_t>& scans)
{
   itsScans = scans;
}

/// @brief configure class from the parset
/// @deails This is the legacy interface configuring the class from the parset. It is not expected to be
/// used long term. All usage of the parset is confined to this method.
/// @param[in] parset configuration parameters
void MSSplitter::configure(const LOFAR::ParameterSet& parset)
{
   // tiling parameters for the output MS
   itsBucketSize = parset.getUint32("stman.bucketsize", 64 * 1024);
   itsTileNCorr = parset.getUint32("stman.tilencorr", 4);
   itsTileNChan = parset.getUint32("stman.tilenchan", 1);

   // Read beam selection parameters
   if (parset.isDefined("beams")) {
       const std::vector<uint32_t> v = parset.getUint32Vector("beams", true);
       itsBeams.clear();
       itsBeams.insert(v.begin(), v.end());
       ASKAPLOG_DEBUG_STR(logger, "Including ONLY beams: " << v);
   }

   // Read scan id selection parameters
   if (parset.isDefined("scans")) {
       const vector<uint32_t> v = parset.getUint32Vector("scans", true);
       itsScans.clear();
       itsScans.insert(v.begin(), v.end());
       ASKAPLOG_DEBUG_STR(logger, "Including ONLY scan numbers: " << v);
   }

   // MV - I am not sure why this was commented out. But leave as is.
   // Read field name selection parameters
   //if (parset.isDefined("fieldnames")) {
   //    const vector<string> names = parset.getStringVector("fieldnames", true);
   //    ASKAPLOG_DEBUG_STR(logger, "Including ONLY fields with names: " << names);
   //    const vector<uint32_t> v = configureFieldNameFilter(names,invis);
   //    itsFieldIds.clear();
   //    itsFieldIds.insert(v.begin(), v.end());
   //    ASKAPLOG_DEBUG_STR(logger, "  fields: " << v);
   //}

   // MV - we can add time selection here, although the long term plan is not to rely on this method
   // and carry parameters explicitly.
}

boost::shared_ptr<casacore::MeasurementSet> MSSplitter::create(
    const std::string& filename, const casacore::Bool addSigmaSpec)
{
    const casacore::uInt bucketSize = itsBucketSize < 8192u ? 8192u: itsBucketSize;

    const casacore::uInt tileNcorr = itsTileNCorr < 1u ? 1u : itsTileNCorr;

    const casacore::uInt tileNchan = itsTileNChan < 1u ? 1u : itsTileNChan;

    ASKAPLOG_DEBUG_STR(logger, "Creating dataset " << filename);

    // Make MS with standard columns
    TableDesc msDesc(MS::requiredTableDesc());

    // Add the DATA column.
    MS::addColumnToDesc(msDesc, MS::DATA, 2);

    // Add the SIGMA_SPECTRUM column?
    if (addSigmaSpec) {
        MS::addColumnToDesc(msDesc, MS::SIGMA_SPECTRUM, 2);
    }

    SetupNewTable newMS(filename, msDesc, Table::New);

    // Set the default Storage Manager to be the Incr one
    {
        IncrementalStMan incrStMan("ismdata", bucketSize);
        newMS.bindAll(incrStMan, True);
    }

    // Bind ANTENNA1, and ANTENNA2 to the StandardStMan
    // as they may change sufficiently frequently to make the
    // incremental storage manager inefficient for these columns.
    {
        // NOTE: The addition of the FEED columns here is a bit unusual.
        // While the FEED columns are perfect candidates for the incremental
        // storage manager, for some reason doing so results in a huge
        // increase in I/O to the file (see ticket: 4094 for details).
        StandardStMan ssm("ssmdata", bucketSize);
        newMS.bindColumn(MS::columnName(MS::ANTENNA1), ssm);
        newMS.bindColumn(MS::columnName(MS::ANTENNA2), ssm);
        newMS.bindColumn(MS::columnName(MS::FEED1), ssm);
        newMS.bindColumn(MS::columnName(MS::FEED2), ssm);
        newMS.bindColumn(MS::columnName(MS::UVW), ssm);
    }

    // These columns contain the bulk of the data so save them in a tiled way
    {
        // Get nr of rows in a tile.
        const int bytesPerRow = sizeof(std::complex<float>) * tileNcorr * tileNchan;
        const int nrowTile = std::max(1u, bucketSize / bytesPerRow);
        TiledShapeStMan dataMan("TiledData",
                                IPosition(3, tileNcorr, tileNchan, nrowTile));
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::DATA),
                         dataMan);
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::FLAG),
                         dataMan);
        if (addSigmaSpec) {
            newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::SIGMA_SPECTRUM),
                             dataMan);
        }
    }
    {
        const int bytesPerRow = sizeof(float) * tileNcorr;
        const int nrowTile = std::max(1u, bucketSize / bytesPerRow);
        TiledShapeStMan dataMan("TiledWeight",
                                IPosition(2, tileNcorr, nrowTile));
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::SIGMA),
                         dataMan);
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::WEIGHT),
                         dataMan);
    }

    // Now we can create the MeasurementSet and add the (empty) subtables
    boost::shared_ptr<casacore::MeasurementSet> ms(new MeasurementSet(newMS, 0));
    ms->createDefaultSubtables(Table::New);
    ms->flush();

    // Set the TableInfo
    {
        TableInfo& info(ms->tableInfo());
        info.setType(TableInfo::type(TableInfo::MEASUREMENTSET));
        info.setSubType(String(""));
        info.readmeAddLine("This is a MeasurementSet Table holding simulated astronomical observations");
    }
    return ms;
}

void MSSplitter::copyAntenna(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSAntennaColumns& sc = srcMsc.antenna();

    MSColumns destMsc(dest);
    MSAntennaColumns& dc = destMsc.antenna();

    // Add new rows to the destination and copy the data
    dest.antenna().addRow(sc.nrow());

    dc.name().putColumn(sc.name());
    dc.station().putColumn(sc.station());
    dc.type().putColumn(sc.type());
    dc.mount().putColumn(sc.mount());
    dc.position().putColumn(sc.position());
    dc.dishDiameter().putColumn(sc.dishDiameter());
    dc.flagRow().putColumn(sc.flagRow());
}

void MSSplitter::copyDataDescription(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSDataDescColumns& sc = srcMsc.dataDescription();

    MSColumns destMsc(dest);
    MSDataDescColumns& dc = destMsc.dataDescription();

    // Add new rows to the destination and copy the data
    dest.dataDescription().addRow(sc.nrow());

    dc.flagRow().putColumn(sc.flagRow());
    dc.spectralWindowId().putColumn(sc.spectralWindowId());
    dc.polarizationId().putColumn(sc.polarizationId());
}

void MSSplitter::copyFeed(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSFeedColumns& sc = srcMsc.feed();

    MSColumns destMsc(dest);
    MSFeedColumns& dc = destMsc.feed();

    // Add new rows to the destination and copy the data
    dest.feed().addRow(sc.nrow());

    dc.antennaId().putColumn(sc.antennaId());
    dc.feedId().putColumn(sc.feedId());
    dc.spectralWindowId().putColumn(sc.spectralWindowId());
    dc.beamId().putColumn(sc.beamId());
    dc.numReceptors().putColumn(sc.numReceptors());
    dc.position().putColumn(sc.position());
    dc.beamOffset().putColumn(sc.beamOffset());
    dc.polarizationType().putColumn(sc.polarizationType());
    dc.polResponse().putColumn(sc.polResponse());
    dc.receptorAngle().putColumn(sc.receptorAngle());
    dc.time().putColumn(sc.time());
    dc.interval().putColumn(sc.interval());
}

void MSSplitter::copyField(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSFieldColumns& sc = srcMsc.field();

    MSColumns destMsc(dest);
    MSFieldColumns& dc = destMsc.field();

    // Add new rows to the destination and copy the data
    dest.field().addRow(sc.nrow());

    dc.name().putColumn(sc.name());
    dc.code().putColumn(sc.code());
    dc.time().putColumn(sc.time());
    dc.numPoly().putColumn(sc.numPoly());
    dc.sourceId().putColumn(sc.sourceId());
    dc.delayDir().putColumn(sc.delayDir());
    dc.phaseDir().putColumn(sc.phaseDir());
    dc.referenceDir().putColumn(sc.referenceDir());
}

void MSSplitter::copyObservation(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSObservationColumns& sc = srcMsc.observation();

    MSColumns destMsc(dest);
    MSObservationColumns& dc = destMsc.observation();

    // Add new rows to the destination and copy the data
    dest.observation().addRow(sc.nrow());

    dc.timeRange().putColumn(sc.timeRange());
    //dc.log().putColumn(sc.log());
    //dc.schedule().putColumn(sc.schedule());
    dc.flagRow().putColumn(sc.flagRow());
    dc.observer().putColumn(sc.observer());
    dc.telescopeName().putColumn(sc.telescopeName());
    dc.project().putColumn(sc.project());
    dc.releaseDate().putColumn(sc.releaseDate());
    dc.scheduleType().putColumn(sc.scheduleType());
}

void MSSplitter::copyPointing(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSPointingColumns& sc = srcMsc.pointing();

    MSColumns destMsc(dest);
    MSPointingColumns& dc = destMsc.pointing();

    // Add new rows to the destination and copy the data
    dest.pointing().addRow(sc.nrow());

    // Create and copy non-standard columns, if they exist.
    // dest row order is different to src when copes come after required columns, so do them first.
    if (source.pointing().actualTableDesc().isColumn("AZIMUTH")) {
        addNonStandardPointingColumn("AZIMUTH", source.pointing(), dest.pointing());
    }
    if (source.pointing().actualTableDesc().isColumn("ELEVATION")) {
        addNonStandardPointingColumn("ELEVATION", source.pointing(), dest.pointing());
    }
    if (source.pointing().actualTableDesc().isColumn("POLANGLE")) {
        addNonStandardPointingColumn("POLANGLE", source.pointing(), dest.pointing());
    }

    // Copy required columns

    // Improved way of copying the target and direction arrays.
    // Need to copy the entire column (not just the first row), but
    // doing it with a single putColumn was found to be too slow (not
    // hanging, just taking a long time with the second call).
    ASKAPCHECK(sc.direction().nrow()==sc.target().nrow(),
               "Different numbers of rows for POINTING table's DIRECTION & TARGET columns. Exiting.");
    for(unsigned int i=0;i<sc.direction().nrow();i++){
        dc.direction().put(i,sc.direction().get(i));
        dc.target().put(i,sc.target().get(i));
    }

    dc.antennaId().putColumn(sc.antennaId());
    dc.interval().putColumn(sc.interval());
    dc.name().putColumn(sc.name());
    dc.numPoly().putColumn(sc.numPoly());
    dc.time().putColumn(sc.time());
    dc.timeOrigin().putColumn(sc.timeOrigin());
    dc.tracking().putColumn(sc.tracking());

}

void MSSplitter::copyPolarization(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSPolarizationColumns& sc = srcMsc.polarization();

    MSColumns destMsc(dest);
    MSPolarizationColumns& dc = destMsc.polarization();

    // Add new rows to the destination and copy the data
    dest.polarization().addRow(sc.nrow());

    dc.flagRow().putColumn(sc.flagRow());
    dc.numCorr().putColumn(sc.numCorr());
    dc.corrType().putColumn(sc.corrType());
    dc.corrProduct().putColumn(sc.corrProduct());
}

void MSSplitter::addNonStandardPointingColumn(const std::string &name,
                                              const MSPointing &srcPointing,
                                              MSPointing &destPointing)
{
    ASKAPDEBUGASSERT(!destPointing.actualTableDesc().isColumn(name));
    destPointing.addColumn(srcPointing.actualTableDesc().columnDesc(name));
    casacore::ScalarColumn<casacore::Float> destCol(destPointing, name);
    destCol.putColumn(ROScalarColumn<casacore::Float>(srcPointing, name));
}

casacore::Int MSSplitter::findSpectralWindowId(const casacore::MeasurementSet& ms)
{
    const ROMSColumns msc(ms);
    const casacore::uInt nrows = msc.nrow();
    ASKAPCHECK(nrows > 0, "No rows in main table");
    const casacore::ROMSDataDescColumns& ddc = msc.dataDescription();

    casacore::Int r0 = -1; // Row zero SpWindow id

    for (casacore::uInt row = 0; row < nrows; ++row) {
        const casacore::Int dataDescId = msc.dataDescId()(row);
        const casacore::Int spwId = ddc.spectralWindowId()(dataDescId);

        if (row == 0) {
            r0 = spwId;
        } else {
            ASKAPCHECK(spwId == r0, "All rows must be of the same spectral window");
        }
    }

    return r0;
}

void MSSplitter::splitSpectralWindow(const casacore::MeasurementSet& source,
        casacore::MeasurementSet& dest,
        const uint32_t startChan,
        const uint32_t endChan,
        const uint32_t width,
        const casacore::Int spwId)
{
    MSColumns destCols(dest);
    const ROMSColumns srcCols(source);

    MSSpWindowColumns& dc = destCols.spectralWindow();
    const ROMSSpWindowColumns& sc = srcCols.spectralWindow();
    const casacore::Int srow = spwId;
    const casacore::Int drow = dc.nrow();

    dest.spectralWindow().addRow();

    // 1: Copy over the simple cells (i.e. those not needing splitting/averaging)
    dc.measFreqRef().put(drow, sc.measFreqRef()(srow));
    dc.refFrequency().put(drow, sc.refFrequency()(srow));
    dc.flagRow().put(drow, sc.flagRow()(srow));
    dc.freqGroup().put(drow, sc.freqGroup()(srow));
    dc.freqGroupName().put(drow, sc.freqGroupName()(srow));
    dc.ifConvChain().put(drow, sc.ifConvChain()(srow));
    dc.name().put(drow, sc.name()(srow));
    dc.netSideband().put(drow, sc.netSideband()(srow));

    // 2: Now process each source measurement set, building up the arrays
    const uInt nChanIn = endChan - startChan + 1;
    const uInt nChanOut = nChanIn / width;
    vector<double> chanFreq;
    vector<double> chanWidth;
    vector<double> effectiveBW;
    vector<double> resolution;
    chanFreq.resize(nChanOut);
    chanWidth.resize(nChanOut);
    effectiveBW.resize(nChanOut);
    resolution.resize(nChanOut);
    double totalBandwidth = 0.0;

    for (uInt destChan = 0; destChan < nChanOut; ++destChan) {
        chanFreq[destChan] = 0.0;
        chanWidth[destChan] = 0.0;
        effectiveBW[destChan] = 0.0;
        resolution[destChan] = 0.0;

        // The offset for the first input channel for this destination channel
        const uInt chanOffset = startChan - 1 + (destChan * width);

        for (uInt i = chanOffset; i < chanOffset + width; ++i) {
            chanFreq[destChan] += sc.chanFreq()(srow)(casacore::IPosition(1, i));
            chanWidth[destChan] += sc.chanWidth()(srow)(casacore::IPosition(1, i));
            effectiveBW[destChan] += sc.effectiveBW()(srow)(casacore::IPosition(1, i));
            resolution[destChan] += sc.resolution()(srow)(casacore::IPosition(1, i));
            totalBandwidth += sc.chanWidth()(srow)(casacore::IPosition(1, i));
        }

        // Finally average chanFreq
        chanFreq[destChan] = chanFreq[destChan] / width;
    }

    // 3: Add those splitting/averaging cells
    dc.numChan().put(drow, nChanOut);
    dc.chanFreq().put(drow, casacore::Vector<double>(chanFreq));
    dc.chanWidth().put(drow, casacore::Vector<double>(chanWidth));
    dc.effectiveBW().put(drow, casacore::Vector<double>(effectiveBW));
    dc.resolution().put(drow, casacore::Vector<double>(resolution));
    dc.totalBandwidth().put(drow, totalBandwidth);
}

bool MSSplitter::rowFiltersExist() const
{
    return !itsBeams.empty() || !itsScans.empty() || !itsFieldIds.empty()
        || itsTimeBegin > std::numeric_limits<double>::min()
        || itsTimeEnd < std::numeric_limits<double>::max();
}

bool MSSplitter::rowIsFiltered(uint32_t scanid, uint32_t fieldid,
                               uint32_t feed1, uint32_t feed2,
                               double time) const
{
    // Include all rows if no filters exist
    if (!rowFiltersExist()) return false;

    if (time < itsTimeBegin || time > itsTimeEnd) return true;

    if (!itsScans.empty() && itsScans.find(scanid) == itsScans.end()) return true;

    if (!itsFieldIds.empty() && itsFieldIds.find(fieldid) == itsFieldIds.end()) return true;

    if (!itsBeams.empty()) {
        if (itsBeams.find(feed1) == itsBeams.end() ||
            itsBeams.find(feed2) == itsBeams.end()) {
            // beam not found in any element of row
            return true;
        }

        else {
            return false;
        }
    }

    return false;
}

void MSSplitter::splitMainTable(const casacore::MeasurementSet& source,
                                casacore::MeasurementSet& dest,
                                const uint32_t startChan,
                                const uint32_t endChan,
                                const uint32_t width)
{
    // Pre-conditions
    ASKAPDEBUGASSERT(endChan >= startChan);
    ASKAPDEBUGASSERT((endChan - startChan + 1) % width == 0);

    const ROMSColumns sc(source);
    MSColumns dc(dest);

    // Add rows upfront if no row based filters exist
    const casacore::uInt nRows = sc.nrow();
    if (!rowFiltersExist()) dest.addRow(nRows);

    // Work out how many channels are to be actual input and which output
    // and how many polarisations are involved.
    const uInt nChanIn = endChan - startChan + 1;
    const uInt nChanOut = nChanIn / width;
    const uInt nPol = sc.data().shape(0)(0);
    ASKAPDEBUGASSERT(nPol > 0);

    // Test to see whether SIGMA_SPECTRUM has been added
    casacore::Bool haveInSigmaSpec = source.isColumn(MS::SIGMA_SPECTRUM);
    casacore::Bool haveOutSigmaSpec = dest.isColumn(MS::SIGMA_SPECTRUM);
    if (haveInSigmaSpec) {
        ASKAPLOG_DEBUG_STR(logger, "Reading and using the spectra of sigma values");
    }
    if (haveOutSigmaSpec) {
        ASKAPLOG_DEBUG_STR(logger, "Calculating and storing spectra of sigma values");
    }

    // Decide how many rows to process simultaneously. This needs to fit within
    // a reasonable amount of memory, because all visibilities will be read
    // in for possible averaging. Assumes 128MB working space.
    std::size_t inDataSize = sizeof(casacore::Complex) + sizeof(casacore::Bool);
    std::size_t outDataSize = inDataSize;
    if (haveInSigmaSpec) {
        inDataSize += sizeof(casacore::Float);
    }
    if (haveOutSigmaSpec) {
        outDataSize += sizeof(casacore::Float);
    }
    uInt maxSimultaneousRows =  (128 * 1024 * 1024) / nPol / (nChanIn * inDataSize) / (nChanOut * outDataSize);
    if (maxSimultaneousRows<1) maxSimultaneousRows = 1;

    // However, if there is row-based filtering only one row can be copied
    // at a time.
    if (rowFiltersExist()) maxSimultaneousRows = 1;

    // Set a 64MB maximum cache size for the large columns
    const casacore::uInt cacheSize = 64 * 1024 * 1024;

    sc.data().setMaximumCacheSize(cacheSize);
    dc.data().setMaximumCacheSize(cacheSize);
    sc.flag().setMaximumCacheSize(cacheSize);
    dc.flag().setMaximumCacheSize(cacheSize);
    if (haveInSigmaSpec) {
        sc.sigmaSpectrum().setMaximumCacheSize(cacheSize);
    }
    if (haveOutSigmaSpec) {
        dc.sigmaSpectrum().setMaximumCacheSize(cacheSize);
    }

    uInt progressCounter = 0; // Used for progress reporting
    const uInt PROGRESS_INTERVAL_IN_ROWS = nRows / 100;

    // Row in destination table may differ from source table if row based
    // filtering is used
    uInt dstRow = 0;
    uInt row = 0;
    while (row < nRows) {
        // Number of rows to process for this iteration of the loop; either
        // maxSimultaneousRows or the remaining rows.
        const uInt nRowsThisIteration = min(maxSimultaneousRows, nRows - row);
        const Slicer srcrowslicer(IPosition(1, row), IPosition(1, nRowsThisIteration),
                Slicer::endIsLength);
        Slicer dstrowslicer = srcrowslicer;

        // Report progress at intervals and on completion
        progressCounter += nRowsThisIteration;
        if (progressCounter >= PROGRESS_INTERVAL_IN_ROWS ||
                (row >= nRows - 1)) {
            ASKAPLOG_DEBUG_STR(logger,  "Processed row " << row + 1 << " of " << nRows);
            progressCounter = 0;
        }

        // Debugging for chunk copying only
        if (nRowsThisIteration > 1) {
            ASKAPLOG_DEBUG_STR(logger,  "Processing " << nRowsThisIteration
                    << " rows this iteration");
        }

        // Skip this row if it is filtered out
        if (rowIsFiltered(sc.scanNumber()(row),
                    sc.fieldId()(row),
                    sc.feed1()(row),
                    sc.feed2()(row),
                    sc.time()(row))) {
            row += nRowsThisIteration;
            continue;
        }

        // Rows have been pre-added if no row based filtering is done
        if (rowFiltersExist()) {
            dest.addRow();
            dstrowslicer = Slicer(IPosition(1, dstRow), IPosition(1, nRowsThisIteration),
                    Slicer::endIsLength);
        }

        // Copy over the simple cells (i.e. those not needing averaging/merging)
        dc.scanNumber().putColumnRange(dstrowslicer, sc.scanNumber().getColumnRange(srcrowslicer));
        dc.fieldId().putColumnRange(dstrowslicer, sc.fieldId().getColumnRange(srcrowslicer));
        dc.dataDescId().putColumnRange(dstrowslicer, sc.dataDescId().getColumnRange(srcrowslicer));
        dc.time().putColumnRange(dstrowslicer, sc.time().getColumnRange(srcrowslicer));
        dc.timeCentroid().putColumnRange(dstrowslicer, sc.timeCentroid().getColumnRange(srcrowslicer));
        dc.arrayId().putColumnRange(dstrowslicer, sc.arrayId().getColumnRange(srcrowslicer));
        dc.processorId().putColumnRange(dstrowslicer, sc.processorId().getColumnRange(srcrowslicer));
        dc.exposure().putColumnRange(dstrowslicer, sc.exposure().getColumnRange(srcrowslicer));
        dc.interval().putColumnRange(dstrowslicer, sc.interval().getColumnRange(srcrowslicer));
        dc.observationId().putColumnRange(dstrowslicer, sc.observationId().getColumnRange(srcrowslicer));
        dc.antenna1().putColumnRange(dstrowslicer, sc.antenna1().getColumnRange(srcrowslicer));
        dc.antenna2().putColumnRange(dstrowslicer, sc.antenna2().getColumnRange(srcrowslicer));
        dc.feed1().putColumnRange(dstrowslicer, sc.feed1().getColumnRange(srcrowslicer));
        dc.feed2().putColumnRange(dstrowslicer, sc.feed2().getColumnRange(srcrowslicer));
        dc.uvw().putColumnRange(dstrowslicer, sc.uvw().getColumnRange(srcrowslicer));
        dc.flagRow().putColumnRange(dstrowslicer, sc.flagRow().getColumnRange(srcrowslicer));
        dc.weight().putColumnRange(dstrowslicer, sc.weight().getColumnRange(srcrowslicer));
        dc.sigma().putColumnRange(dstrowslicer, sc.sigma().getColumnRange(srcrowslicer)/sqrt(width));

        // Set the shape of the destination arrays
        for (uInt i = dstRow; i < dstRow + nRowsThisIteration; ++i) {
            dc.data().setShape(i, IPosition(2, nPol, nChanOut));
            dc.flag().setShape(i, IPosition(2, nPol, nChanOut));
            if (haveOutSigmaSpec) {
                dc.sigmaSpectrum().setShape(i, IPosition(2, nPol, nChanOut));
            }
        }

        //  Average (if applicable) then write data into the output MS
        const Slicer srcarrslicer(IPosition(2, 0, startChan - 1),
                                  IPosition(2, nPol, nChanIn), Slicer::endIsLength);
        const Slicer destarrslicer(IPosition(2, 0, 0),
                                   IPosition(2, nPol, nChanOut), Slicer::endIsLength);

        if (width == 1) {
            dc.data().putColumnRange(dstrowslicer, destarrslicer,
                sc.data().getColumnRange(srcrowslicer, srcarrslicer));
            dc.flag().putColumnRange(dstrowslicer, destarrslicer,
                sc.flag().getColumnRange(srcrowslicer, srcarrslicer));
            if (haveInSigmaSpec && haveOutSigmaSpec) {
                dc.sigmaSpectrum().putColumnRange(dstrowslicer, destarrslicer,
                    sc.sigmaSpectrum().getColumnRange(srcrowslicer, srcarrslicer));
            }
        } else {
            // Get (read) the input data/flag/sigma
            const casacore::Cube<casacore::Complex> indata = sc.data().getColumnRange(srcrowslicer, srcarrslicer);
            const casacore::Cube<casacore::Bool> inflag = sc.flag().getColumnRange(srcrowslicer, srcarrslicer);
            // This is only needed if generating sigmaSpectra, but that should be the
            // case with width>1, and this avoids testing in the tight loops below
            casacore::Cube<casacore::Float> insigma;
            if (haveInSigmaSpec) {
                insigma = sc.sigmaSpectrum().getColumnRange(srcrowslicer, srcarrslicer);
            } else {
                // There's only 1 sigma per pol & row, so spread over channels
                insigma = casacore::Cube<casacore::Float>(indata.shape());
                casacore::IPosition arrayShape(3, nPol,1,nRowsThisIteration);
                casacore::Array<casacore::Float> sigmaArray =
                    sc.sigma().getColumnRange(srcrowslicer).reform(arrayShape);
                for (uInt i = 0; i < nChanIn; ++i) {
                    const Slicer blockSlicer(IPosition(3, 0,i,0),
                                             arrayShape, Slicer::endIsLength);
                    insigma(blockSlicer) = sigmaArray;
                }
            }

            // Create the output data/flag/sigma
            casacore::Cube<casacore::Complex> outdata(nPol, nChanOut, nRowsThisIteration);
            casacore::Cube<casacore::Bool> outflag(nPol, nChanOut, nRowsThisIteration);
            // This is only needed if generating sigmaSpectra, but that should be the
            // case with width>1, and this avoids testing in the tight loops below
            casacore::Cube<casacore::Float> outsigma(nPol, nChanOut, nRowsThisIteration);

            // Average data and combine flag information
            for (uInt pol = 0; pol < nPol; ++pol) {
                for (uInt destChan = 0; destChan < nChanOut; ++destChan) {
                    for (uInt r = 0; r < nRowsThisIteration; ++r) {
                        casacore::Complex sum(0.0, 0.0);
                        casacore::Float varsum = 0.0;
                        casacore::uInt sumcount = 0;

                        // Starting at the appropriate offset into the source data, average "width"
                        // channels together
                        for (uInt i = (destChan * width); i < (destChan * width) + width; ++i) {
                            ASKAPDEBUGASSERT(i < nChanIn);
                            if (inflag(pol, i, r)) continue;
                            sum += indata(pol, i, r);
                            varsum += insigma(pol, i, r) * insigma(pol, i, r);
                            sumcount++;
                        }

                        // Now the input channels have been averaged, write the data to
                        // the output cubes
                        if (sumcount > 0) {
                            outdata(pol, destChan, r) = casacore::Complex(sum.real() / sumcount,
                                                                      sum.imag() / sumcount);
                            outflag(pol, destChan, r) = false;
                            outsigma(pol, destChan, r) = sqrt(varsum) / sumcount;
                        } else {
                            outflag(pol, destChan, r) = true;
                        }

                    }
                }
            }

            // Put (write) the output data/flag
            dc.data().putColumnRange(dstrowslicer, destarrslicer, outdata);
            dc.flag().putColumnRange(dstrowslicer, destarrslicer, outflag);
            if (haveOutSigmaSpec) {
                dc.sigmaSpectrum().putColumnRange(dstrowslicer, destarrslicer, outsigma);
            }
        }

        row += nRowsThisIteration;
        dstRow += nRowsThisIteration;
    }
}

/// @brief Entry point method
/// @param[in] invis input measurement set name
/// @param[in] outvis the name of the measurement set to write
/// @param[in] startChan start channel to extract from the input MS
/// @param[in] endChan end channel to extract from the input MS
/// @param[in] width number of consecutive channels to average
void MSSplitter::split(const std::string& invis, const std::string& outvis,
                      const uint32_t startChan,
                      const uint32_t endChan,
                      const uint32_t width)
{
    ASKAPLOG_DEBUG_STR(logger,  "Splitting out channel range " << startChan << " to "
                          << endChan << " (inclusive)");

    if (width > 1) {
        ASKAPLOG_DEBUG_STR(logger,  "Averaging " << width << " channels to form 1");
    } else {
        ASKAPLOG_DEBUG_STR(logger,  "No averaging");
    }

    // Verify split parameters
    const uInt nChanIn = endChan - startChan + 1;

    if ((width < 1) || (nChanIn % width != 0)) {
        ASKAPTHROW(AskapError, "Width must equally divide the channel range");
    }

    // Open the input measurement set
    const casacore::MeasurementSet in(invis);

    // Verify split parameters that require input MS info
    const casacore::uInt totChanIn = ROScalarColumn<casacore::Int>(in.spectralWindow(),"NUM_CHAN")(0);
    if ((startChan<1) || (endChan > totChanIn)) {
        ASKAPTHROW(AskapError, "Input channel range is inconsistent with input spectra: ["<<
            startChan<<","<<endChan<<"] is outside [1,"<<totChanIn<<"]");
    }

    // Create the output measurement set
    if (casacore::File(outvis).exists()) {
        ASKAPTHROW(AskapError, "File or table " << outvis << " already exists!");
    }

    // Add a sigma spectrum to the output measurement set?
    casacore::Bool addSigmaSpec = false;
    if ((width > 1) || in.isColumn(MS::SIGMA_SPECTRUM)) {
        addSigmaSpec = true;
    }

    boost::shared_ptr<casacore::MeasurementSet> out(create(outvis, addSigmaSpec));

    // Copy ANTENNA
    ASKAPLOG_DEBUG_STR(logger,  "Copying ANTENNA table");
    copyAntenna(in, *out);

    // Copy DATA_DESCRIPTION
    ASKAPLOG_DEBUG_STR(logger,  "Copying DATA_DESCRIPTION table");
    copyDataDescription(in, *out);

    // Copy FEED
    ASKAPLOG_DEBUG_STR(logger,  "Copying FEED table");
    copyFeed(in, *out);

    // Copy FIELD
    ASKAPLOG_DEBUG_STR(logger,  "Copying FIELD table");
    copyField(in, *out);

    // Copy OBSERVATION
    ASKAPLOG_DEBUG_STR(logger,  "Copying OBSERVATION table");
    copyObservation(in, *out);

    // Copy POINTING
    ASKAPLOG_DEBUG_STR(logger,  "Copying POINTING table");
    copyPointing(in, *out);

    // Copy POLARIZATION
    ASKAPLOG_DEBUG_STR(logger,  "Copying POLARIZATION table");
    copyPolarization(in, *out);

    // Get the spectral window id (must be common for all main table rows)
    const casacore::Int spwId = findSpectralWindowId(in);

    // Split SPECTRAL_WINDOW
    ASKAPLOG_DEBUG_STR(logger,  "Splitting SPECTRAL_WINDOW table");
    splitSpectralWindow(in, *out, startChan, endChan, width, spwId);

    // Split main table
    ASKAPLOG_DEBUG_STR(logger,  "Splitting main table");
    splitMainTable(in, *out, startChan, endChan, width);
}

/*
// MV - unused

void MSSplitter::configureTimeFilter(const std::string& key, const std::string& msg,
                                 double& var)
{

    const string ts = itsParset.getString(key);
    casacore::Quantity tq;
    if(!casacore::MVTime::read(tq, ts)) {
        ASKAPTHROW(AskapError, "Unable to convert " << ts << " to MVTime");
    }

    const casacore::MVTime t(tq);
    var = t.second();
    ASKAPLOG_DEBUG_STR(logger, msg << ts << " (" << var << " sec)");

}



std::vector<uint32_t> MSSplitter::configureFieldNameFilter(
                          const std::vector<std::string>& names,
                          const std::string invis)
{
    std::vector<uint32_t> fieldIds;
    if (!names.empty()) {
        const casacore::MeasurementSet in(invis);
        const ROMSColumns srcMsc(in);
        const ROMSFieldColumns& sc = srcMsc.field();
        const casacore::Vector<casacore::String> fieldNames = sc.name().getColumn();
        // Step through each field and find IDs for the filter.
        // Could set fieldIds in the following loop, but this seems easier.
        for (uInt i = 0; i < sc.nrow(); ++i) {
            if (find(names.begin(), names.end(), (std::string)fieldNames[i]) !=
                     names.end()) {
                fieldIds.push_back(i);
            }
        }
        // print a warning for any missing fields
        for (uInt i = 0; i < names.size(); ++i) {
            if (find(fieldNames.begin(), fieldNames.end(),
                       (casacore::String)names[i]) == fieldNames.end()) {
                ASKAPLOG_WARN_STR(logger, "  cannot find field name " <<
                    names[i] << " in ms "<< invis);
            }
        }
    }
    if (fieldIds.empty()) {
        ASKAPTHROW(AskapError, "Cannot find any of the field names " <<
            names << " in ms "<< invis);
    }
    return fieldIds;
}

*/

