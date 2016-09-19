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
#include "askap_synthesis.h"

// System includes
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

// ASKAPsoft includes
#include "askap/AskapError.h"
#include "askap/AskapLogging.h"
#include "askap/StatReporter.h"
#include "askap/Log4cxxLogSink.h"
#include "boost/shared_ptr.hpp"
#include "boost/exception/all.hpp"
#include "boost/program_options.hpp"
#include "casacore/casa/OS/File.h"
#include "casacore/casa/aips.h"
#include "casacore/casa/Quanta.h"
#include "casacore/casa/Arrays/Vector.h"
#include "casacore/casa/Arrays/MatrixMath.h"
#include "casacore/tables/Tables/TableDesc.h"
#include "casacore/tables/Tables/SetupNewTab.h"
#include "casacore/tables/DataMan/IncrementalStMan.h"
#include "casacore/tables/DataMan/StandardStMan.h"
#include "casacore/tables/DataMan/TiledShapeStMan.h"
#include "casacore/ms/MeasurementSets/MeasurementSet.h"
#include "casacore/ms/MeasurementSets/MSColumns.h"

ASKAP_LOGGER(logger, ".msmerge2");

using namespace askap;
using namespace casa;

boost::shared_ptr<casa::MeasurementSet> create(const std::string& filename, casa::uInt tileNcorr = 4, casa::uInt tileNchan = 1)
{
    // Get configuration first to ensure all parameters are present
    casa::uInt bucketSize =  128 * 1024;

    if (bucketSize < 8192) {
        bucketSize = 8192;
    }
    if (tileNcorr < 1) {
        tileNcorr = 1;
    }
    if (tileNchan < 1) {
        tileNchan = 1;
    }

    ASKAPLOG_DEBUG_STR(logger, "Creating dataset " << filename);

    // Make MS with standard columns
    TableDesc msDesc(MS::requiredTableDesc());

    // Add the DATA column.
    MS::addColumnToDesc(msDesc, MS::DATA, 2);

    SetupNewTable newMS(filename, msDesc, Table::New);

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
        // Get nr of rows in a tile.
        const int nrowTile = std::max(1u, bucketSize / (8*tileNcorr*tileNchan));
        TiledShapeStMan dataMan("TiledData",
                IPosition(3, tileNcorr, tileNchan, nrowTile));
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::DATA),
                dataMan);
        newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::FLAG),
                dataMan);
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

void copyPointing(const casa::MeasurementSet& source, casa::MeasurementSet& dest)
{
    const ROMSColumns srcMsc(source);
    const ROMSPointingColumns& sc = srcMsc.pointing();

    MSColumns destMsc(dest);
    MSPointingColumns& dc = destMsc.pointing();

    // Add new rows to the destination and copy the data
    dest.pointing().addRow(sc.nrow());

    // Improved way of copying the target & direction arrays
    // Need to copy the entire column (not just the first row), but
    // doing it with a single putColumn was found to be too slow (not
    // hanging, just taking a long time with the second call).
    ASKAPLOG_DEBUG_STR(logger, "Starting copy of direction & target columns");
    ASKAPCHECK(sc.direction().nrow()==sc.target().nrow(),
               "Different numbers of rows for POINTING table's DIRECTION & TARGET columns. Exiting.");
    for(unsigned int i=0;i<sc.direction().nrow();i++){
        dc.direction().put(i,sc.direction().get(i));
        dc.target().put(i,sc.target().get(i));
    }
    ASKAPLOG_DEBUG_STR(logger, "Finished copy of direction & target columns");
    
    dc.antennaId().putColumn(sc.antennaId());
    dc.interval().putColumn(sc.interval());
    dc.name().putColumn(sc.name());
    dc.numPoly().putColumn(sc.numPoly());
    dc.time().putColumn(sc.time());
    dc.timeOrigin().putColumn(sc.timeOrigin());
    dc.tracking().putColumn(sc.tracking());
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

void mergeMainTable(const std::vector< boost::shared_ptr<const ROMSColumns> >& srcMscs,
                    casa::MeasurementSet& dest)
{
    MSColumns dc(dest);

    // Add rows upfront
    const ROMSColumns& sc = *(srcMscs[0]);
    const casa::uInt nRows = sc.nrow();
    dest.addRow(nRows);

    // For each row
    for (uInt row = 0; row < nRows; ++row) {
        if (row % 10000 == 0) {
            ASKAPLOG_INFO_STR(logger,  "Merging row " << row << " of " << nRows);
        }

        // 1: Copy over the simple cells (i.e. those not needing merging)
        dc.scanNumber().put(row, sc.scanNumber()(row));
        dc.fieldId().put(row, sc.fieldId()(row));
        dc.dataDescId().put(row, 0);
        dc.time().put(row, sc.time()(row));
        dc.timeCentroid().put(row, sc.timeCentroid()(row));
        dc.arrayId().put(row, sc.arrayId()(row));
        dc.processorId().put(row, sc.processorId()(row));
        dc.exposure().put(row, sc.exposure()(row));
        dc.interval().put(row, sc.interval()(row));
        dc.observationId().put(row, sc.observationId()(row));
        dc.antenna1().put(row, sc.antenna1()(row));
        dc.antenna2().put(row, sc.antenna2()(row));
        dc.feed1().put(row, sc.feed1()(row));
        dc.feed2().put(row, sc.feed2()(row));
        dc.uvw().put(row, sc.uvw()(row));
        dc.flagRow().put(row, sc.flagRow()(row));
        dc.weight().put(row, sc.weight()(row));
        dc.sigma().put(row, sc.sigma()(row));

        // 2: Size the matrix for data and flag
        const uInt nPol = sc.data()(row).shape()(0);
        const uInt nChan = sc.data()(row).shape()(1);
        const uInt nChanTotal = nChan * srcMscs.size();

        casa::Matrix<casa::Complex> data(nPol, nChanTotal);
        casa::Matrix<casa::Bool> flag(nPol, nChanTotal);

        // 3: Copy the data from each input into the output matrix
        for (uInt i = 0; i < srcMscs.size(); ++i) {
            const casa::Matrix<casa::Complex> srcData = srcMscs[i]->data()(row);
            const casa::Matrix<casa::Bool> srcFlag = srcMscs[i]->flag()(row);
            for (uInt pol = 0; pol < nPol; ++pol) {
                for (uInt chan = 0; chan < nChan; ++chan) {
                    data(pol, chan + (nChan * i)) = srcData(pol, chan);
                    flag(pol, chan + (nChan * i)) = srcFlag(pol, chan);
                }
            }
        }

        // 4: Add those merged cells
        dc.data().put(row, data);
        dc.flag().put(row, flag);
    }
}

void merge(const std::vector<std::string>& inFiles, const std::string& outFile, casa::uInt tileNcorr = 4, casa::uInt tileNchan = 1)
{
    // Create the output measurement set
    ASKAPCHECK(!casa::File(outFile).exists(), "File or table "
            << outFile << " already exists!");
    boost::shared_ptr<casa::MeasurementSet> out(create(outFile,tileNcorr,tileNchan));

    // Open the input measurement sets
    std::vector< boost::shared_ptr<const casa::MeasurementSet> > in;
    std::vector< boost::shared_ptr<const ROMSColumns> > inColumns;
    std::vector<std::string>::const_iterator it;
    for (it = inFiles.begin(); it != inFiles.end(); ++it) {
        const boost::shared_ptr<const casa::MeasurementSet> p(new casa::MeasurementSet(*it));
        in.push_back(p);
        inColumns.push_back(boost::shared_ptr<const ROMSColumns>(new ROMSColumns(*p)));
    }

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

    // Merge main table
    ASKAPLOG_INFO_STR(logger,  "Merging main table");
    mergeMainTable(inColumns, *out);
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
        std::string outName;
        std::vector<std::string> inNamesVec;
        std::vector<std::string> inNames;
        namespace po = boost::program_options;
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message") 
        ("tileNcorr,x", po::value<int>(&tileNcorr)->default_value(4), "Number of correlations per tile") 
        ("tileNchan,c", po::value<int>(&tileNchan)->default_value(1), "Number of channels per tile") 
        ("input-file,i", po::value< vector<string> >(&inNames), "Input file(s) - you can also just list them after the other options")
        ("output-file,o",po::value<std::string>(&outName)->default_value("out.ms"),"Output filename");       
        
        po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

        if (vm.count("tileNcorr")) {
            cout << "Number of correlations per tile was set to " 
            << vm["tileNcorr"].as<int>() << ".\n";
        } 
        else 
        {
            cout << "tileNcorr was not setr, using default: " << tileNcorr << "\n";
        }

        if (vm.count("tileNchan")) {
            cout << "Number of channels per tile  was set to " 
            << vm["tileNchan"].as<int>() << ".\n";
        } 
        else 
        {
            cout << "tileNchan was not set, using default: " << tileNchan << "\n";
        }

        if (vm.count("input-file"))
        {
            cout << "Input files are: " 
            << vm["input-file"].as< vector<string> >() << "\n";
        }
 
        ASKAPLOG_INFO_STR(logger,
                "This program merges given measurement sets and writes the output into `"
                << outName << "`" << "with a tiling: " << tileNcorr << "," << tileNchan );

        // Turns inNames into vector<string>
        std::vector<std::string>::iterator it;
        for (it = inNames.begin(); it < inNames.end(); ++it) {
            // Check for file existence first - this array may be longer than required

            if (access((*it).c_str(),F_OK) != -1)
                inNamesVec.push_back(*it);
        }

        merge(inNamesVec, outName, tileNcorr, tileNchan); 

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
