/// @file msmerge.cc
///
/// @brief
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

// Package level header file
#include "askap/askap_synthesis.h"

// System includes
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

// ASKAPsoft includes
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/StatReporter.h>
#include <askap/askap/Log4cxxLogSink.h>
#include <boost/shared_ptr.hpp>
#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Quanta.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/DataMan/IncrementalStMan.h>
#include <casacore/tables/DataMan/StandardStMan.h>
#include <casacore/tables/DataMan/TiledShapeStMan.h>
//#include <casacore/tables/DataMan/DataManAccessor.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>
#include <casacore/casa/Quanta/MVTime.h>

ASKAP_LOGGER(logger, ".msmerge");

using namespace casa;

namespace askap {

boost::shared_ptr<casa::MeasurementSet> create(const std::string& filename, const casacore::StorageOption & option,
    casa::uInt tileNcorr = 4, casa::uInt tileNchan = 1, casa::uInt tileNrow = 0)
{
    casa::uInt bucketSize = 1024*1024;

    if (tileNcorr < 1) {
        tileNcorr = 1;
    }
    if (tileNchan < 1) {
        tileNchan = 1;
    }
    if (tileNrow < 1) {
        tileNrow = 1;
    }

    ASKAPLOG_DEBUG_STR(logger, "Creating dataset " << filename);

    // Make MS with standard columns
    TableDesc msDesc(MS::requiredTableDesc());

    // Add the DATA column.
    MS::addColumnToDesc(msDesc, MS::DATA, 2);

    SetupNewTable newMS(filename, msDesc, Table::New, option);

    // Set the default Storage Manager to be the Incr one
    {
        IncrementalStMan incrStMan("ismdata", bucketSize);
        newMS.bindAll(incrStMan, True);
    }

    // Bind ANTENNA1, and ANTENNA2 to the standardStMan
    // as they may change sufficiently frequently to make the
    // incremental storage manager inefficient for these columns.

    {
        StandardStMan ssm("ssmdata", bucketSize);
        newMS.bindColumn(MS::columnName(MS::ANTENNA1), ssm);
        newMS.bindColumn(MS::columnName(MS::ANTENNA2), ssm);
        newMS.bindColumn(MS::columnName(MS::UVW), ssm);
    }

    // These columns contain the bulk of the data so save them in a tiled way
    {
        TiledShapeStMan dataMan("TiledData",
                IPosition(3, tileNcorr, tileNchan, tileNrow));
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::DATA),
                dataMan);

        // Avoid storing less than 8 flag bits per row, it is slow
        casa::uInt tileNchanFlag = tileNchan;
        if (tileNcorr*tileNchanFlag<8) tileNchanFlag = 8 / tileNcorr;
        TiledShapeStMan dataManF("TiledFlag",
                IPosition(3, tileNcorr, tileNchanFlag, tileNrow));
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::FLAG),
                dataManF);
    }
    {
        const int nrowTile = std::max(1u, bucketSize / (4*8));
        TiledShapeStMan dataMan("TiledWeight",
                IPosition(2, 4, nrowTile));
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::SIGMA),
                dataMan);
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::WEIGHT),
                dataMan);
    }

    // Now we can create the MeasurementSet and add the (empty) subtables
    boost::shared_ptr<casa::MeasurementSet> ms(new MeasurementSet(newMS, 0));
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

void copyAntenna(const casa::MeasurementSet& source, casa::MeasurementSet& dest)
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

void mergeDataDescription(const casa::MeasurementSet& source, casa::MeasurementSet& dest)
{
    MSColumns destMsc(dest);
    MSDataDescColumns& dc = destMsc.dataDescription();
    dest.dataDescription().addRow();


    dc.flagRow().put(0, false);
    dc.spectralWindowId().put(0, 0);
    dc.polarizationId().put(0, 0);
}

void copyFeed(const casa::MeasurementSet& source, casa::MeasurementSet& dest)
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

void copyField(const casa::MeasurementSet& source, casa::MeasurementSet& dest)
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

void copyObservation(const casa::MeasurementSet& source, casa::MeasurementSet& dest)
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

void addNonStandardPointingColumn(const std::string &name,
                                    const MSPointing &srcPointing,
                                    MSPointing &destPointing,
                                    ROScalarColumn<casacore::Float> &src,
                                    ScalarColumn<casacore::Float> &dest)
{
    ASKAPDEBUGASSERT(!destPointing.actualTableDesc().isColumn(name));
    destPointing.addColumn(srcPointing.actualTableDesc().columnDesc(name));
    dest = casacore::ScalarColumn<casacore::Float>(destPointing, name);
    src = casacore::ROScalarColumn<casacore::Float>(srcPointing, name);
}

void copyPointing(const casacore::MeasurementSet& source, casacore::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSPointingColumns& sc = srcMsc.pointing();

    MSColumns destMsc(dest);
    MSPointingColumns& dc = destMsc.pointing();

    // Create and copy non-standard columns, if they exist.
    // dest row order is different to src when copes come after required columns, so do them first.
    const Bool doAz = source.pointing().actualTableDesc().isColumn("AZIMUTH");
    const Bool doEl = source.pointing().actualTableDesc().isColumn("ELEVATION");
    const Bool doPolAng = source.pointing().actualTableDesc().isColumn("POLANGLE");
    ScalarColumn<casacore::Float> dAz, dEl, dPolAng;
    ROScalarColumn<casacore::Float> sAz, sEl, sPolAng;
    if (doAz) {
        addNonStandardPointingColumn("AZIMUTH", source.pointing(), dest.pointing(), sAz, dAz);
    }
    if (doEl) {
        addNonStandardPointingColumn("ELEVATION", source.pointing(), dest.pointing(), sEl, dEl);
    }
    if (doPolAng) {
        addNonStandardPointingColumn("POLANGLE", source.pointing(), dest.pointing(), sPolAng, dPolAng);
    }

    // Copy required columns

    ASKAPCHECK(sc.direction().nrow()==sc.target().nrow(),
               "Different numbers of rows for POINTING table's DIRECTION & TARGET columns. Exiting.");

    // Pointing tables can get large. Doing whole column operations seems very slow (maybe due to
    // storage manager used). So doing row by row copy instead (finished in 6s vs >1h)
    const casacore::rownr_t nRow = sc.nrow();
    for(casacore::rownr_t i=0;i<nRow;i++){
        dest.pointing().addRow();
        dc.direction().put(i,sc.direction()(i));
        dc.target().put(i,sc.target()(i));
        dc.antennaId().put(i,sc.antennaId()(i));
        dc.interval().put(i,sc.interval()(i));
        dc.name().put(i,sc.name()(i));
        dc.numPoly().put(i,sc.numPoly()(i));
        dc.time().put(i,sc.time()(i));
        dc.timeOrigin().put(i,sc.timeOrigin()(i));
        dc.tracking().put(i,sc.tracking()(i));
        if (doAz) dAz.put(i,sAz(i));
        if (doEl) dEl.put(i,sEl(i));
        if (doPolAng) dPolAng.put(i,sPolAng(i));
    }
}

void copyPolarization(const casa::MeasurementSet& source, casa::MeasurementSet& dest)
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

void appendToVector(const casa::Vector<double>& src, std::vector<double>& dest)
{
    std::copy(src.begin(), src.end(), std::back_inserter(dest));
}

casa::Int findSpectralWindowId(const ROMSColumns& msc)
{
    const casa::uInt nrows = msc.nrow();
    ASKAPCHECK(nrows > 0, "No rows in main table");
    const casa::ROMSDataDescColumns& ddc = msc.dataDescription();

    casa::Int r0 = -1; // Row zero SpWindow id

    for (casa::uInt row = 0; row < nrows; ++row) {
        const casa::Int dataDescId = msc.dataDescId()(row);
        const casa::Int spwId = ddc.spectralWindowId()(dataDescId);

        if (row == 0) {
            r0 = spwId;
        } else {
            ASKAPCHECK(spwId == r0, "All rows must be of the same spectral window");
        }
    }

    return r0;
}

// Creates a single spectral window in the "dest" measurement set, which
// is the concatenation of the spectral windows used in the source measurement
// set. Each source measurement set is read through once to determine which
// spectral window id is used. All rows in a given source measurement set
// must refer to the same spectral window id.
void mergeSpectralWindow(const std::vector< boost::shared_ptr<const ROMSColumns> >& srcMscs,
                         casa::MeasurementSet& dest)
{
    ASKAPCHECK(!srcMscs.empty(), "Vector of source measurement sets is empty");

    MSColumns destMsc(dest);
    MSSpWindowColumns& dc = destMsc.spectralWindow();
    const casa::Int DEST_ROW = 0;

    // 1: Create a single spectral window in the destination measurement set.
    // Populate it with the simple cells (i.e. those not needing merging)
    const casa::Int spwIdForFirstMs = findSpectralWindowId(*srcMscs[0]);
    const ROMSSpWindowColumns& sc = srcMscs[0]->spectralWindow();
    dest.spectralWindow().addRow();
    dc.measFreqRef().put(DEST_ROW, sc.measFreqRef()(spwIdForFirstMs));
    dc.refFrequency().put(DEST_ROW, sc.refFrequency()(spwIdForFirstMs));
    dc.flagRow().put(DEST_ROW, false);
    dc.freqGroup().put(DEST_ROW, sc.freqGroup()(spwIdForFirstMs));
    dc.freqGroupName().put(DEST_ROW, sc.freqGroupName()(spwIdForFirstMs));
    dc.ifConvChain().put(DEST_ROW, sc.ifConvChain()(spwIdForFirstMs));
    dc.name().put(DEST_ROW, "Merged Window");
    dc.netSideband().put(DEST_ROW, sc.netSideband()(spwIdForFirstMs));

    // 2: Now process each source measurement sets, building up the arrays
    std::vector<double> chanFreq;
    std::vector<double> chanWidth;
    std::vector<double> effectiveBW;
    std::vector<double> resolution;
    uInt nChan = 0;
    double totalBandwidth = 0.0;

    for (uInt i = 0; i < srcMscs.size(); ++i) {
        // The below function obtains the spwid refered to by all rows in the
        // main table. This throws an exception if they are not all identical.
        const casa::Int srcSpwId = findSpectralWindowId(*srcMscs[i]);

        const ROMSSpWindowColumns& spwc = srcMscs[i]->spectralWindow();
        nChan += spwc.numChan()(srcSpwId);
        totalBandwidth += spwc.totalBandwidth()(srcSpwId);

        appendToVector(spwc.chanFreq()(srcSpwId), chanFreq);
        appendToVector(spwc.chanWidth()(srcSpwId), chanWidth);
        appendToVector(spwc.effectiveBW()(srcSpwId), effectiveBW);
        appendToVector(spwc.resolution()(srcSpwId), resolution);
    }

    // 3: Add those merged cells
    dc.numChan().put(DEST_ROW, nChan);
    dc.chanFreq().put(DEST_ROW, casa::Vector<double>(chanFreq));
    dc.chanWidth().put(DEST_ROW, casa::Vector<double>(chanWidth));
    dc.effectiveBW().put(DEST_ROW, casa::Vector<double>(effectiveBW));
    dc.resolution().put(DEST_ROW, casa::Vector<double>(resolution));
    dc.totalBandwidth().put(DEST_ROW, totalBandwidth);
}

Vector<double> collectTimes(const std::vector< boost::shared_ptr<const ROMSColumns> >& srcMscs) {
    const int n = srcMscs.size();
    casa::Vector<casa::Vector<double> > times(n);
    // work out first and last time present
    double first = 0, last=0;
    for (uInt i = 0; i < n; i++) {
        times(i) = srcMscs[i]->time().getColumn();
        if (i == 0) {
            first = min(times(i));
            last = max(times(i));
        } else {
            first = min(first,min(times(i)));
            last = max(last,max(times(i)));
        }
    }
    // get integration time
    double interval = 1e10;
    for (uInt i=1; i < times(0).nelements(); i++) {
        double diff = times(0)(i)-times(0)(i-1);
        if (diff > 0) interval = min (interval, diff);
    }

    // get max number of integrations present
    uint maxIntegrations = (last - first) / interval + 1;
    casa::Vector<double> mergedTimes(maxIntegrations);

    // get list of all times present
    uint count = 0;
    mergedTimes(0) = first;
    casa::Vector<uInt> ind(n,0);
    while (mergedTimes(count) < last) {
        double nextTime = last;
        for (uInt i = 0; i < n; i++) {
            uInt ntimes = times(i).nelements();
            while(ind(i) < ntimes && times(i)(ind(i)) <= mergedTimes(count)) ind(i)++;
            if (ind(i) < ntimes) {
                nextTime = min(times(i)(ind(i)), nextTime);
            }
        }
        mergedTimes(++count) = nextTime;
    }
    count++;
    mergedTimes.resize(count,True);

    return mergedTimes;
}

void mergeMainTable(const std::vector< boost::shared_ptr<const ROMSColumns> >& srcMscs,
                    boost::shared_ptr<casa::MeasurementSet>& dest, const IPosition& tileShape,
                    uInt nBaselines, bool dryRun, bool detectMissingIntegrations)
{
    uInt nMS = srcMscs.size();
    Vector<double> mergedTimes;
    uInt nRows = srcMscs[0]->nrow();
    // nTimes is really the number of tiles in the row direction
    // unless we are detecting missing integrations
    uInt nTimes = nRows / tileShape(2);
    if (nRows % tileShape(2) != 0) nTimes++;
    // baselines per tile
    uInt nBase = tileShape(2);
    // Check for missing integrations
    // This assumes ASKAP data with fixed number of baselines per integration throughout
    if (detectMissingIntegrations) {
        mergedTimes = collectTimes(srcMscs);
        nTimes = mergedTimes.nelements();
        nRows = nBaselines * nTimes;
        // baselines per integration
        nBase = nBaselines;
    }
    boost::shared_ptr<MSColumns> dc;
    if (!dryRun) {
        dc = boost::shared_ptr<MSColumns>(new MSColumns(*dest));
        // Add rows upfront
        dest->addRow(nRows);
    }
    ASKAPLOG_INFO_STR(logger,  "Number of rows in output:"<< nRows<<", tileshape = "<< tileShape);
    // 2: Size the matrix for data and flag
    const uInt nPol = srcMscs[0]->data().shape(0)(0);
    uInt nChanTotal = 0;
    casa::Vector<uInt> nChan(nMS);
    for (uInt i = 0; i < nMS; i++) {
        nChan(i) = srcMscs[i]->data().shape(0)(1);
        nChanTotal += nChan(i);
    }
    casa::Cube<casa::Complex> data(nPol, nChanTotal,tileShape(2));
    casa::Cube<casa::Bool> flag(nPol, nChanTotal,tileShape(2));

    int lastPerc = 0;
    uInt nIntperTile = tileShape(2)/nBase;
    uInt lastTileStart = nTimes - 1;
    uInt lastTileSize = nRows - lastTileStart * nBase;
    if (detectMissingIntegrations) {
        lastTileStart = (nTimes / nIntperTile) * nIntperTile;
        lastTileSize = (nTimes - lastTileStart) * nBase;
    }
    //ASKAPLOG_INFO_STR(logger, "last tile start "<< lastTileStart << ", last tile size = "<<lastTileSize);
    Vector<uInt> row(nMS,0);
    uInt dataRow = 0;
    const double tol = 1.0; // tolerance for comparing timestamps is 1s
    // For each integration
    for (uInt j = 0, outRow=0; j < nTimes; j++, outRow+=nBase) {
        const Slicer destRowSlicer(IPosition(1, outRow), IPosition(1, nBase), Slicer::endIsLength);
        if (10*outRow/nRows > lastPerc/10) {
            if (detectMissingIntegrations) {
                ASKAPLOG_INFO_STR(logger,  "Merging row " <<  outRow<< " of " << nRows <<", Integration "<<j<<"/"<<nTimes);
            } else {
                ASKAPLOG_INFO_STR(logger,  "Merging row " <<  outRow<< " of " << nRows);
            }
            lastPerc = 100*outRow/nRows;
        }

        // Find first MS that has this timeslot
        casa::Bool found = False;
        uInt k = 0;
        if (detectMissingIntegrations) {
            for (k = 0; k< srcMscs.size(); k++) {
                if (abs(srcMscs[k]->time()(row(k)) - mergedTimes(j)) < tol)  break;
            }
        }
        ASKAPCHECK(k<srcMscs.size(),"Logic error in mergeMainTable");
        const ROMSColumns& sc(*srcMscs[k]);
        const Slicer srcRowSlicer(IPosition(1,row(k)),IPosition(1,nBase), Slicer::endIsLength);
        //ASKAPLOG_INFO_STR(logger, "Integration # "<<j << ", outRow = "<<outRow << ", row("<<k<<")="<<row(k));
        // 1: Copy over the simple cells (i.e. those not needing merging)
        if (!dryRun) {
            dc->scanNumber().putColumnRange(destRowSlicer, sc.scanNumber().getColumnRange(srcRowSlicer));
            dc->fieldId().putColumnRange(destRowSlicer, sc.fieldId().getColumnRange(srcRowSlicer));
            dc->dataDescId().putColumnRange(destRowSlicer, sc.dataDescId().getColumnRange(srcRowSlicer));
            dc->time().putColumnRange(destRowSlicer, sc.time().getColumnRange(srcRowSlicer));
            dc->timeCentroid().putColumnRange(destRowSlicer, sc.timeCentroid().getColumnRange(srcRowSlicer));
            dc->arrayId().putColumnRange(destRowSlicer, sc.arrayId().getColumnRange(srcRowSlicer));
            dc->processorId().putColumnRange(destRowSlicer, sc.processorId().getColumnRange(srcRowSlicer));
            dc->exposure().putColumnRange(destRowSlicer, sc.exposure().getColumnRange(srcRowSlicer));
            dc->interval().putColumnRange(destRowSlicer, sc.interval().getColumnRange(srcRowSlicer));
            dc->observationId().putColumnRange(destRowSlicer, sc.observationId().getColumnRange(srcRowSlicer));
            dc->antenna1().putColumnRange(destRowSlicer, sc.antenna1().getColumnRange(srcRowSlicer));
            dc->antenna2().putColumnRange(destRowSlicer, sc.antenna2().getColumnRange(srcRowSlicer));
            dc->feed1().putColumnRange(destRowSlicer, sc.feed1().getColumnRange(srcRowSlicer));
            dc->feed2().putColumnRange(destRowSlicer, sc.feed2().getColumnRange(srcRowSlicer));
            dc->uvw().putColumnRange(destRowSlicer, sc.uvw().getColumnRange(srcRowSlicer));
            dc->flagRow().putColumnRange(destRowSlicer, sc.flagRow().getColumnRange(srcRowSlicer));
            dc->weight().putColumnRange(destRowSlicer, sc.weight().getColumnRange(srcRowSlicer));
            dc->sigma().putColumnRange(destRowSlicer, sc.sigma().getColumnRange(srcRowSlicer));
        }
        if (j == lastTileStart) {
            // for the last set of rows we need to resize
            data.resize(nPol,nChanTotal,lastTileSize);
            flag.resize(nPol,nChanTotal,lastTileSize);
        }

        // 3: Copy the data from each input into the output matrix
        uInt destChan = 0;
        uInt destRow = (j % nIntperTile) * nBase;
        for (uInt i = 0; i < nMS; ++i) {
            // check if data is available
            if (!detectMissingIntegrations || abs(srcMscs[i]->time()(row(i)) - mergedTimes(j)) < tol) {
                if (!dryRun) {
                    const Slicer srcDataSlicer(IPosition(1,row(i)),IPosition(1,nBase), Slicer::endIsLength);
                    const casa::Cube<casa::Complex> srcData = srcMscs[i]->data().getColumnRange(srcDataSlicer);
                    const casa::Cube<casa::Bool> srcFlag = srcMscs[i]->flag().getColumnRange(srcDataSlicer);
                    for (uInt irow = 0, our; irow < nBase; irow++) {
                        data(Slice(), Slice(destChan,nChan(i)), destRow + irow) = srcData(Slice(), Slice(), irow);
                        flag(Slice(), Slice(destChan,nChan(i)), destRow + irow) = srcFlag(Slice(), Slice(), irow);
                    }
                }
                row(i) += nBase;
            } else {
                // no data for this integration for this MS, set flags
                if (!dryRun) {
                    for (uInt irow = 0; irow < nBase; irow++) {
                        data(Slice(), Slice(destChan,nChan(i)), destRow + irow) = 0;
                        flag(Slice(), Slice(destChan,nChan(i)), destRow + irow) = True;
                    }
                }
                ASKAPLOG_INFO_STR(logger,  "Missing integration for input file " <<  srcMscs[i]->time().table().tableName()
                    << " at row " << row(i) << ", "<<MVTime(mergedTimes(j)/C::day).string(MVTime::YMD));
            }
            destChan += nChan(i);
        }
        // 4: Add those merged cells when we've filled a tile or reached the end
        if ((j+1) % nIntperTile == 0 || j == nTimes - 1) {
            //ASKAPLOG_INFO_STR(logger, "Writing tile at int# "<<j << "/"<<nTimes<<" data row ="<<
            //dataRow<< ", data shape = "<<data.shape());
            if (!dryRun) {
                const Slicer dataRowSlicer(IPosition(1, dataRow), IPosition(1, data.shape()(2)),
                        Slicer::endIsLength);
                dc->data().putColumnRange(dataRowSlicer, data);
                dc->flag().putColumnRange(dataRowSlicer, flag);
            }
            dataRow += nIntperTile * nBase;
        }
    }
}

void merge(const std::vector<std::string>& inFiles, const std::string& outFile,
    casacore::StorageOption& option, casa::uInt tileNcorr,
    casa::uInt tileNchan, casa::uInt tileNrow, bool dryRun, casa::uInt maxBuf)
{
    // Open the input measurement sets
    std::vector< boost::shared_ptr<const casa::MeasurementSet> > in;
    std::vector< boost::shared_ptr<const ROMSColumns> > inColumns;
    std::vector<std::string>::const_iterator it;
    uInt nChanOut = 0;
    uInt nBaselines = 0;
    uInt nRow = 0;
    bool sameSize = true;
    casa::String telescope;

    ASKAPASSERT(inFiles.size()>0);
    for (it = inFiles.begin(); it != inFiles.end(); ++it) {
        const boost::shared_ptr<const casa::MeasurementSet> p(new casa::MeasurementSet(*it));
        ASKAPCHECK(p->nrow(),"One of the input files has zero rows: "<<*it);
        in.push_back(p);
        inColumns.push_back(boost::shared_ptr<const ROMSColumns>(new ROMSColumns(*p)));
        nChanOut += inColumns.back()->data().shape(0)(1);
        if (nBaselines ==0) {
            uInt nAnt = inColumns.back()->antenna().nrow();
            // number of baselines including autocorrelations
            nBaselines = nAnt * (nAnt + 1) / 2;
        }
        // Check all inputs MSs have same number of rows
        if (nRow > 0) {
            sameSize = sameSize && (nRow == p->nrow());
        } else {
            nRow = p->nrow();
            telescope = inColumns.back()->observation().telescopeName()(0);
        }
    }
    // We can only deal with missing rows for ASKAP data
    ASKAPCHECK(sameSize || telescope == "ASKAP",
        "All input MeasurementSets should have the same number of rows");

    if (tileNcorr < 1) tileNcorr = 1;
    if (tileNchan < 1) tileNchan = 1;
    // Set tileNrow large, but not so large that caching uses more than maxBuf
    bool detectMissingIntegrations = false;
    if (tileNrow==0) {
        const casa::uInt nTilesPerRow = (nChanOut-1)/tileNchan+1;
        const casa::uInt bucketSize = std::max(8192u,maxBuf/nTilesPerRow);
        tileNrow = std::min(nRow,std::max(1u,bucketSize / (8 * tileNcorr * tileNchan)));
        // make it a multiple of the number of rows per integration for ASKAP
        if (nBaselines > 0 && telescope == "ASKAP") {
            tileNrow = std::max(1u, (tileNrow / nBaselines)) * nBaselines;
            detectMissingIntegrations = true;
        }
        ASKAPLOG_INFO_STR(logger, "Setting tileNrow to " << tileNrow);
        // Don't allow tiny size buckets
    } else if (tileNcorr * tileNchan * tileNrow * 8 < 8192u) {
        tileNrow = 1024u / (tileNcorr * tileNchan);
        ASKAPLOG_INFO_STR(logger, "Setting tileNrow to " << tileNrow);
    }
    // Refuse to proceed if it would produce corrupted data
    ASKAPCHECK(detectMissingIntegrations || sameSize,
        "All input MeasurementSets should have the same number of rows, \n" <<
        "for ASKAP data leave tileNrow unset to enable missing integration detection");
    if (detectMissingIntegrations) {
        ASKAPLOG_INFO_STR(logger,"Missing integration detection active - missing parts of the spectrum will be flagged");
    }
    // Create the output measurement set
    ASKAPCHECK(!casa::File(outFile).exists(), "File or table "
            << outFile << " already exists!");
    boost::shared_ptr<casa::MeasurementSet> out;
    if (!dryRun) {
        out = create(outFile,option,tileNcorr,tileNchan,tileNrow);

        ASKAPDEBUGASSERT(in.begin() != in.end());
        ASKAPDEBUGASSERT(out);

        ASKAPLOG_INFO_STR(logger,  "First copy " << inFiles[0]<< " into " << outFile);

        // Copy ANTENNA
        ASKAPLOG_INFO_STR(logger,  "Copying ANTENNA table");
        copyAntenna(**(in.begin()), *out);

        // Copy FEED
        ASKAPLOG_INFO_STR(logger,  "Copying FEED table");
        copyFeed(**(in.begin()), *out);

        // Copy FIELD
        ASKAPLOG_INFO_STR(logger,  "Copying FIELD table");
        copyField(**(in.begin()), *out);

        // Copy OBSERVATION
        ASKAPLOG_INFO_STR(logger,  "Copying OBSERVATION table");
        copyObservation(**(in.begin()), *out);

        // Copy POINTING
        ASKAPLOG_INFO_STR(logger,  "Copying POINTING table");
        copyPointing(**(in.begin()), *out);

        // Copy POLARIZATION
        ASKAPLOG_INFO_STR(logger,  "Copying POLARIZATION table");
        copyPolarization(**(in.begin()), *out);

        // Merge DATA_DESCRIPTION
        // This actually just creates a single row, for the single spectral window
        // that will be created in mergeSpectralWindow()
        ASKAPLOG_INFO_STR(logger,  "Merging DATA_DESCRIPTION table");
        mergeDataDescription(**(in.begin()), *out);

        // Merge SPECTRAL_WINDOW
        ASKAPLOG_INFO_STR(logger,  "Merging SPECTRAL_WINDOW table");
        mergeSpectralWindow(inColumns, *out);
    }
    // Merge main table
    ASKAPLOG_INFO_STR(logger,  "Merging main table");
    mergeMainTable(inColumns,out,IPosition(3,tileNcorr,tileNchan,tileNrow), nBaselines, dryRun,
                    detectMissingIntegrations);
    // Uncomment this to check if the caching is working
    //RODataManAccessor(**(in.begin()), "TiledData", False).showCacheStatistics (cout);
    //if (!dryRun) RODataManAccessor(*out, "TiledData", False).showCacheStatistics (cout);

}

// Main function
int main(int argc, const char** argv)
{
    // Now we have to initialize the logger before we use it
    // If a log configuration exists in the current directory then
    // use it, otherwise try to use the programs default one
    std::ifstream config("askap.log_cfg", std::ifstream::in);
    if (config) {
        ASKAPLOG_INIT("askap.log_cfg");
    } else {
        std::ostringstream ss;
        ss << argv[0] << ".log_cfg";
        ASKAPLOG_INIT(ss.str().c_str());
    }

    // Ensure that CASA log messages are captured
    casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
    casa::LogSink::globalSink(globalSink);

    try {
        StatReporter stats;
        int tileNcorr;
        int tileNchan;
        int tileNrow;
        int blockSize;
        uint maxBuf;
        std::string fileFormat;
        bool useODirect;
        bool dryRun;
        std::string outName;
        std::string paramFile;
        std::vector<std::string> inNamesVec;
        std::vector<std::string> inNames;
        namespace po = boost::program_options;
        // Declare the supported options.
        po::options_description command("Command line only options");
        command.add_options()
        ("help,h", "produce help message")
        ("parameter-file,p",po::value<std::string>(&paramFile),"Parameter file");

        po::options_description generic("Generic options");
        generic.add_options()
        ("tileNcorr,x", po::value<int>(&tileNcorr)->default_value(4), "Number of correlations per tile")
        ("tileNchan,c", po::value<int>(&tileNchan)->default_value(1), "Number of channels per tile")
        ("tileNrow,r", po::value<int>(&tileNrow)->default_value(0), "Number of rows per tile")
        ("blocksize,b", po::value<int>(&blockSize)->default_value(4*1024*1024), "Blocksize")
        ("bufferMB,B", po::value<uint>(&maxBuf)->default_value(2000u), "Max. buffer size in MB")
        ("format,f", po::value<string>(&fileFormat)->default_value("separate"),"File storage format (separate, combined, hdf5)")
        ("odirect", po::value<bool>(&useODirect)->default_value(false),"Use direct (unbuffered) IO if possible")
        ("dryrun,d", po::value<bool>(&dryRun)->default_value(false),"Don't produce any output")
        ("input-file,i", po::value< vector<string> >(&inNames), "Input file(s) - you can also just list them after the other options")
        ("output-file,o",po::value<std::string>(&outName)->default_value("out.ms"),"Output filename");

        po::positional_options_description p;
        p.add("input-file", -1);

        po::options_description desc("Program options");
        desc.add(command).add(generic);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

        po::notify(vm);
        if (vm.count("parameter-file")) {
            std::ifstream config_file(paramFile.c_str());
            if (config_file) {
                ASKAPLOG_INFO_STR(logger,"Reading parameter file: "<< paramFile);
                po::store(po::parse_config_file(config_file, generic), vm);
            }
        }

        po::notify(vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }


        if (dryRun) {
            ASKAPLOG_INFO_STR(logger,"Dry run only - not producing output");
        }

        if (vm.count("input-file"))
        {
            ASKAPLOG_INFO_STR(logger,"Input files are: "
            << vm["input-file"].as< vector<string> >());
        }

        ASKAPLOG_INFO_STR(logger,
                "This program merges given measurement sets and writes the output into `"
                << outName << "`" << "with a tiling: " << tileNcorr << "," << tileNchan);
        if (tileNrow > 0) {
            ASKAPLOG_INFO_STR(logger,"Number of rows per tile is set to "<< tileNrow);
        }
        // convert to bytes
        maxBuf *= 1024*1024;

        casacore::StorageOption option;
        if (fileFormat == "separate") {
            ASKAPLOG_INFO_STR(logger, "Setting up storage managers to write separate files");
            ASKAPCHECK(!useODirect, "ODirect option is not supported for storage managers using separate files");
            option = casacore::StorageOption(casacore::StorageOption::SepFile, static_cast<casacore::Int>(blockSize), 0);
        } else if (fileFormat == "combined") {
            ASKAPLOG_INFO_STR(logger, "Setting up storage managers to write combined file with "<<float(blockSize)/1024./1024.<<" Mb blocks and ODirect = "<<useODirect);
            option = casacore::StorageOption(casacore::StorageOption::MultiFile, static_cast<casacore::Int>(blockSize), useODirect ? 1 : 0);
        } else if (fileFormat == "hdf5") {
            ASKAPLOG_INFO_STR(logger, "Setting up storage managers to write combined HDF5 file with "<<float(blockSize)/1024./1024.<<" Mb blocks");
            ASKAPCHECK(!useODirect, "ODirect option is not supported for storage managers set up to write into an HDF5 file");
            option = casacore::StorageOption(casacore::StorageOption::MultiHDF5, static_cast<casacore::Int>(blockSize), 0);
        }

        // Turns inNames into vector<string>
        std::vector<std::string>::iterator it;
        for (it = inNames.begin(); it < inNames.end(); ++it) {
            // Check for file existence first - this array may be longer than required

            if (access((*it).c_str(),F_OK) != -1)
                inNamesVec.push_back(*it);
        }

        if (inNamesVec.size() > 0) {
            merge(inNamesVec, outName, option, tileNcorr, tileNchan, tileNrow, dryRun, maxBuf);
        } else {
            ASKAPLOG_WARN_STR(logger,"No input files could be found");
        }
        stats.logSummary();
        ///==============================================================================
    } catch (boost::exception &ex) {
        ASKAPLOG_FATAL_STR(logger, "Command line parser error, " << boost::diagnostic_information(ex));
        return 1;
    } catch (const askap::AskapError& x) {
        ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
        return 1;
    } catch (const std::exception& x) {
        ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
        return 1;
    }

    return 0;
}

} // end namespace askap

int main(int argc, const char** argv)
{
    return askap::main(argc, argv);
}
