/// @file Simulator.cc
///
/// @copyright (c) 2007 CSIRO
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

// Include own header file first
#include <simulation/Simulator.h>

// ASKAPsoft includes
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>

#include <ms/MeasurementSets/MSDerivedValues.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <casa/Logging/LogIO.h>
#include <casa/Containers/Record.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/StManAipsIO.h>
#include <tables/Tables/IncrementalStMan.h>
#include <tables/Tables/StandardStMan.h>
#include <tables/Tables/TiledShapeStMan.h>
#include <casa/BasicSL/Constants.h>
#include <casa/BasicMath/Random.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayIO.h>
#include <casa/Arrays/Slice.h>
#include <measures/Measures/Stokes.h>
#include <measures/Measures/MeasFrame.h>
#include <casa/Quanta/MVDirection.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Quanta/MVTime.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MeasData.h>
#include <measures/Measures.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Utilities/Assert.h>
#include <casa/Arrays/ArrayUtil.h>
#include <casa/Arrays/Slicer.h>
#include <casa/iostream.h>
#include <casa/fstream.h>
#include <casa/sstream.h>
#include <ms/MeasurementSets/MSTileLayout.h>
#include <scimath/Mathematics/RigidVector.h>
#include <scimath/Mathematics/SquareMatrix.h>

// temporary to get access to beam_offsets
#include <ms/MeasurementSets/MSIter.h>
//

using namespace casa;

namespace askap {
namespace synthesis {

/// a but ugly solution to use the feed table parser of MSIter
/// to extract antennaMounts and BeamOffsets.
/// @note (from MV) the same thing can probably be done via accessor feed table
/// handler (and it would work with time dependent tables correctly)
struct MSFeedParameterExtractor : protected MSIter {

    /// @brief constructor
    /// @param[in] ms measurement set
    MSFeedParameterExtractor(const casa::MeasurementSet &ms) {
        msc_p = new ROMSColumns(ms);
        msc_p->antenna().mount().getColumn(antennaMounts_p, True);
        checkFeed_p = True;
        setFeedInfo();
    }

    /// Return a string mount identifier for each antenna
    using MSIter::antennaMounts;

    /// Return a cube containing pairs of coordinate offset for each receptor
    /// of each feed (values are in radians, coordinate system is fixed with
    /// antenna and is the same as used to define the BEAM_OFFSET parameter
    /// in the feed table). The cube axes are receptor, antenna, feed.
    using MSIter::getBeamOffsets;

    /// True if all elements of the cube returned by getBeamOffsets are zero
    using MSIter::allBeamOffsetsZero;
};

void Simulator::defaults()
{
    fractionBlockageLimit_p = 1e-6;
    elevationLimit_p = Quantity(8.0, "deg");
    autoCorrelationWt_p = 1.0;
    telescope_p = "Unknown";
    qIntegrationTime_p = Quantity(10.0, "s");
    useHourAngle_p = True;
    Quantity today;
    MVTime::read(today, "today");
    mRefTime_p = MEpoch(today, MEpoch::UTC);
    itsNoiseRMS = 1.;
    itsRelAntennaWeight.assign(casa::Vector<double>());
}

Simulator::Simulator(const casa::String& MSName, int bucketSize,
                     int tileNcorr, int tileNchan) :
        ms_p(0), itsDishDiamForNoise(-1.), itsChanBandwidthForNoise(-100.), itsNoiseRMS(1.)
{
    try {
        defaults();

        // make MS with standard columns
        TableDesc msDesc(MS::requiredTableDesc());

        // Add the DATA column.
        MS::addColumnToDesc(msDesc, MS::DATA, 2);

        if (bucketSize < 2048) {
            bucketSize = 2048;
        }

        SetupNewTable newMS(MSName, msDesc, Table::New);

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
            if (tileNcorr <= 0) tileNcorr = 1;

            if (tileNchan <= 0) tileNchan = 1;

            int nrowTile = std::max(1, bucketSize / (8 * tileNcorr * tileNchan));
            TiledShapeStMan dataMan("TiledData",
                                    IPosition(3, tileNcorr,
                                              tileNchan, nrowTile));
            newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::DATA),
                             dataMan);
            newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::FLAG),
                             dataMan);
        }
        {
            int nrowTile = std::max(1, bucketSize / (4 * 8));
            TiledShapeStMan dataMan("TiledWeight",
                                    IPosition(2, 4, nrowTile));
            newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::SIGMA),
                             dataMan);
            newMS.bindColumn(MeasurementSet::columnName(MeasurementSet::WEIGHT),
                             dataMan);
        }

        // Now we can create the MeasurementSet and add the (empty) subtables
        ms_p = new MeasurementSet(newMS, 0);
        ms_p->createDefaultSubtables(Table::New);
        ms_p->flush();

        // Set the TableInfo
        {
            TableInfo& info(ms_p->tableInfo());
            info.setType(TableInfo::type(TableInfo::MEASUREMENTSET));
            info.setSubType(String("simulator"));
            info.readmeAddLine("This is a MeasurementSet Table holding simulated astronomical observations");
        }

        // We're done - wasn't that easy?

    } catch (std::exception& x) {
        ASKAPLOG_ERROR(logger, x.what());
        throw;
    }
}

Simulator::Simulator(casa::MeasurementSet& theMS) :
        ms_p(0), itsDishDiamForNoise(-1.), itsChanBandwidthForNoise(-100.), itsNoiseRMS(1.)
{
    defaults();

    ms_p = new MeasurementSet(theMS);

    ASKAPLOG_INFO_STR(logger, "Opening MeasurementSet " << ms_p->tableName() << " with "
                          << ms_p->nrow() << " rows");

    TableDesc td(ms_p->tableDesc());
    {
        MSColumns msc(*ms_p);
        MSSpWindowColumns& spwc = msc.spectralWindow();
        ASKAPLOG_INFO_STR(logger, "   last spectral window ID = " << spwc.nrow());
    }
}

Simulator::Simulator(const Simulator & mss)
{
    operator=(mss);
}

void Simulator::initAnt(const casa::String& telescope, const casa::Vector<double>& x,
                        const casa::Vector<double>& y, const casa::Vector<double>& z,
                        const casa::Vector<double>& dishDiameter, const casa::Vector<double>&,
                        const casa::Vector<casa::String>& mount, const casa::Vector<casa::String>& name,
                        const casa::String& coordsystem, const casa::MPosition& mRefLocation)
{
    telescope_p = telescope;

    Int nAnt = x.nelements();

    Vector<double> xx(x.nelements());
    Vector<double> yy(x.nelements());
    Vector<double> zz(x.nelements());

    if (coordsystem == "global") {
        xx = x;
        yy = y;
        zz = z;
        ASKAPLOG_INFO_STR(logger, "Using global coordinates for the antennas");
    } else if (coordsystem == "local") {

        MVAngle mvLong = mRefLocation.getAngle().getValue()(0);
        MVAngle mvLat = mRefLocation.getAngle().getValue()(1);

        ASKAPLOG_INFO_STR(logger, "Using local coordinates for the antennas: Reference position = "
                              << mvLong.string(MVAngle::ANGLE, 7) << " "
                              << mvLat.string(MVAngle::DIG2, 7));
        local2global(xx, yy, zz, mRefLocation, x, y, z);
    } else if (coordsystem == "longlat") {
        ASKAPLOG_INFO_STR(logger, "Using longitude-latitude coordinates for the antennas");
        longlat2global(xx, yy, zz, mRefLocation, x, y, z);
    } else {
        ASKAPLOG_INFO_STR(logger, "Unknown coordinate system type: " << coordsystem);
    }
    
    for (size_t i = 0; i<dishDiameter.size(); ++i) {
         if (i == 0) {
             itsDishDiamForNoise = dishDiameter[i];
         } else if (fabs(itsDishDiamForNoise - dishDiameter[i])>1e-6) {
             itsDishDiamForNoise = -1.;
             break;
         }
    }

    Vector<Int> antId(nAnt);
    Matrix<double> antXYZ(3, nAnt);

    for (Int i = 0; i < nAnt; i++) {
        antXYZ(0, i) = xx(i);
        antXYZ(1, i) = yy(i);
        antXYZ(2, i) = zz(i);
        antId(i) = i;
    }

    MSColumns msc(*ms_p);
    MSAntennaColumns& antc = msc.antenna();
    Int numOfAnt = antc.nrow();
    MSAntenna& ant = ms_p->antenna();

    ant.addRow(nAnt); // make nAnt rows
    Slicer antSlice(IPosition(1, numOfAnt),
                    IPosition(1, numOfAnt + nAnt - 1),
                    IPosition(1, 1), Slicer::endIsLast);
    antc.dishDiameter().putColumnRange(antSlice, dishDiameter);
    antc.mount().putColumnRange(antSlice, mount);
    antc.name().putColumnRange(antSlice, name);
    //  antc.offset().putColumnRange(antSlice,offset);
    antc.position().putColumnRange(antSlice, antXYZ);
    antc.station().fillColumn("");
    antc.flagRow().fillColumn(False);
    antc.type().fillColumn("GROUND-BASED");
    ASKAPLOG_INFO_STR(logger, "Added rows to ANTENNA table");
}

void Simulator::local2global(casa::Vector<double>& xGeo, casa::Vector<double>& yGeo,
                             casa::Vector<double>& zGeo, const casa::MPosition& mRefLocation,
                             const casa::Vector<double>& xLocal, const casa::Vector<double>& yLocal,
                             const casa::Vector<double>& zLocal)
{
    uInt nn = xLocal.nelements();
    xGeo.resize(nn);
    yGeo.resize(nn);
    zGeo.resize(nn);

    MPosition::Convert loc2(mRefLocation, MPosition::ITRF);
    MPosition locitrf(loc2());
    Vector<double> xyz = locitrf.get("m").getValue();

    Vector<double> ang = locitrf.getAngle("rad").getValue();
    double d1, d2;
    d1 = ang(0);
    d2 = ang(1);
    double cosLong = cos(d1);
    double sinLong = sin(d1);
    double cosLat = cos(d2);
    double sinLat = sin(d2);

    for (uInt i = 0; i < nn; i++) {

        double xG1 = -sinLat * yLocal(i) + cosLat * zLocal(i);
        double yG1 = xLocal(i);

        xGeo(i) = cosLong * xG1 - sinLong * yG1 + xyz(0);
        yGeo(i) = sinLong * xG1 + cosLong * yG1 + xyz(1);

        zGeo(i) = cosLat * yLocal(i) + sinLat * zLocal(i) + xyz(2);
    }

}

void Simulator::longlat2global(casa::Vector<double>& xReturned,
                               casa::Vector<double>& yReturned,
                               casa::Vector<double>& zReturned,
                               const casa::MPosition& mRefLocation,
                               const casa::Vector<double>& xIn,
                               const casa::Vector<double>& yIn,
                               const casa::Vector<double>& zIn)
{
    ASKAPLOG_INFO_STR(logger, "Simulator::longlat2global not yet implemented, passed parameters "<<xIn<<" "<<
                              yIn<<" "<<zIn<<" "<<mRefLocation<<", and will overwrite "<<xReturned<<" "<<yReturned<<" "<<zReturned);
}

void Simulator::initFields(const casa::String& sourceName,
                           const casa::MDirection& sourceDirection,
                           const casa::String& calCode)
{
    MSColumns msc(*ms_p);
    MSFieldColumns& fieldc = msc.field();
    Int baseFieldID = fieldc.nrow();

    ASKAPLOG_INFO_STR(logger, "Creating new field " << sourceName << ", ID " << baseFieldID
                      + 1);

    ms_p->field().addRow(1); //SINGLE DISH CASE
    fieldc.name().put(baseFieldID, sourceName);
    fieldc.code().put(baseFieldID, calCode);
    fieldc.time().put(baseFieldID, 0.0);
    fieldc.numPoly().put(baseFieldID, 0);
    fieldc.sourceId().put(baseFieldID, 0);
    Vector<MDirection> direction(1);
    direction(0) = sourceDirection;
    fieldc.delayDirMeasCol().put(baseFieldID, direction);
    fieldc.phaseDirMeasCol().put(baseFieldID, direction);
    fieldc.referenceDirMeasCol().put(baseFieldID, direction);

}

void Simulator::initSpWindows(const casa::String& spWindowName, const int& nChan,
                              const casa::Quantity& startFreq, const casa::Quantity& freqInc,
                              const casa::Quantity& freqRes, const casa::String& stokesString)
{
    ASKAPCHECK(fabs(freqInc.getValue("Hz")-freqRes.getValue("Hz"))<1, "freqInc="<<freqInc<<" and freqRes="<<
                    freqRes<<" look different");
    Vector<Int> stokesTypes(4);
    stokesTypes = Stokes::Undefined;
    String myStokesString = stokesString;
    Int nCorr = 0;

    for (Int j = 0; j < 4; j++) {
        while (myStokesString.at(0, 1) == " ") {
            myStokesString.del(0, 1);
        }

        if (myStokesString.length() == 0)
            break;

        stokesTypes(j) = Stokes::type(myStokesString.at(0, 2));
        myStokesString.del(0, 2);
        nCorr = j + 1;

        if (stokesTypes(j) == Stokes::Undefined) {
            ASKAPLOG_INFO_STR(logger, " Undefined polarization type in input");
        }
    }

    MSColumns msc(*ms_p);
    MSSpWindowColumns& spwc = msc.spectralWindow();
    MSDataDescColumns& ddc = msc.dataDescription();
    MSPolarizationColumns& polc = msc.polarization();
    Int baseSpWID = spwc.nrow();
    ASKAPLOG_DEBUG_STR(logger, "Creating new spectral window " << spWindowName << ", ID "
                          << baseSpWID + 1);
    // fill spectralWindow table
    ms_p->spectralWindow().addRow(1);
    ms_p->polarization().addRow(1);
    ms_p->dataDescription().addRow(1);
    spwc.numChan().put(baseSpWID, nChan);
    spwc.name().put(baseSpWID, spWindowName);
    spwc.netSideband().fillColumn(1);
    spwc.ifConvChain().fillColumn(0);
    spwc.freqGroup().fillColumn(0);
    spwc.freqGroupName().fillColumn("Group 1");
    spwc.flagRow().fillColumn(False);
    spwc.measFreqRef().fillColumn(MFrequency::TOPO);
    polc.flagRow().fillColumn(False);
    ddc.flagRow().fillColumn(False);
    polc.numCorr().put(baseSpWID, nCorr);
    Vector <double> freqs(nChan), bandwidth(nChan);
    bandwidth = freqInc.getValue("Hz");
    ddc.spectralWindowId().put(baseSpWID, baseSpWID);
    ddc.polarizationId().put(baseSpWID, baseSpWID);
    double vStartFreq(startFreq.getValue("Hz"));
    double vFreqInc(freqInc.getValue("Hz"));

    for (Int chan = 0; chan < nChan; chan++) {
        freqs(chan) = vStartFreq + chan * vFreqInc;
    }

    // translate stokesTypes into receptor products, catch invalid
    // fallibles.
    Matrix<Int> corrProduct(uInt(2), uInt(nCorr));
    Fallible<Int> fi;
    stokesTypes.resize(nCorr, True);

    for (Int j = 0; j < nCorr; j++) {
        fi = Stokes::receptor1(Stokes::type(stokesTypes(j)));
        corrProduct(0, j) = (fi.isValid() ? fi.value() : 0);
        fi = Stokes::receptor2(Stokes::type(stokesTypes(j)));
        corrProduct(1, j) = (fi.isValid() ? fi.value() : 0);
    }

    spwc.refFrequency().put(baseSpWID, vStartFreq);
    spwc.chanFreq().put(baseSpWID, freqs);
    spwc.chanWidth().put(baseSpWID, bandwidth);
    spwc.effectiveBW().put(baseSpWID, bandwidth);
    spwc.resolution().put(baseSpWID, bandwidth);
    spwc.totalBandwidth().put(baseSpWID, nChan*vFreqInc);
    polc.corrType().put(baseSpWID, stokesTypes);
    polc.corrProduct().put(baseSpWID, corrProduct);

    // store the bandwidth to be able to do noise estimate
    // we need only to take care of a single (e.g. the first defined)
    // spectral window. The consistency is checked during observations.
    // One limitation of this approach is that a particular spectral window
    // may not be used at all, but it will be checked for conformance
    if (itsChanBandwidthForNoise < -10) {
        // negative value less than -10 means that this field is undefined
        itsChanBandwidthForNoise = fabs(vFreqInc);
    }
}

// NOTE:  initAnt and initSpWindows must be called before this one!
void Simulator::initFeeds(const casa::String& mode, const casa::Vector<double>& x,
                          const casa::Vector<double>& y, const casa::Vector<casa::String>& pol)
{
    MSColumns msc(*ms_p);
    MSAntennaColumns& antc = msc.antenna();
    Int nAnt = antc.nrow();

    if (nAnt <= 0) {
        ASKAPLOG_INFO_STR(logger, "Simulator::initFeeds: must call initAnt() first");
    }

    Int nFeed = x.nelements();

    String feedPol0 = "R", feedPol1 = "L";
    Bool isList = False;

    if (nFeed > 0) {
        isList = True;

        if (x.nelements() != y.nelements()) {
            ASKAPLOG_FATAL_STR(logger, "Feed x and y must be the same length");
        }

        ASKAPCHECK(pol.nelements() == x.nelements(),
                   "Feed polarization list must be same length as the number of positions");
        ASKAPLOG_INFO_STR(logger, "Constructing FEED table from list");
    } else {
        nFeed = 1;

        // mode == "perfect R L" OR "perfect X Y"
        if (mode.contains("X", 0)) {
            feedPol0 = "X";
            feedPol1 = "Y";
        }
    }

    Int nRow = nFeed * nAnt;
    Vector<Int> feedAntId(nRow);
    Vector<Int> feedId(nRow);
    Vector<Int> feedSpWId(nRow);
    Vector<Int> feedBeamId(nRow);

    Vector<Int> feedNumRec(nRow);
    Cube<double> beamOffset(2, 2, nRow);

    Matrix<String> feedPol(2, nRow);
    Matrix<double> feedXYZ(3, nRow);
    Matrix<double> feedAngle(2, nRow);
    Cube<Complex> polResp(2, 2, nRow);

    Int iRow = 0;

    if (isList) {
        polResp = Complex(0.0, 0.0);

        for (Int i = 0; i < nAnt; i++) {
            for (Int j = 0; j < nFeed; j++) {
                feedAntId(iRow) = i;
                feedId(iRow) = j;
                feedSpWId(iRow) = -1;
                feedBeamId(iRow) = 0;
                feedNumRec(iRow) = 2;
                beamOffset(0, 0, iRow) = x(j);
                beamOffset(1, 0, iRow) = y(j);
                beamOffset(0, 1, iRow) = x(j);
                beamOffset(1, 1, iRow) = y(j);
                feedXYZ(0, iRow) = 0.0;
                feedXYZ(1, iRow) = 0.0;
                feedXYZ(2, iRow) = 0.0;
                feedAngle(0, iRow) = 0.0;
                feedAngle(1, iRow) = 0.0;

                if (pol(j).contains("X", 0)) {
                    feedPol(0, iRow) = "X";
                    feedPol(1, iRow) = "Y";
                } else {
                    feedPol(0, iRow) = "L";
                    feedPol(1, iRow) = "R";
                }

                polResp(0, 0, iRow) = polResp(1, 1, iRow) = Complex(1.0, 0.0);
                iRow++;
            }
        }
    } else {
        polResp = Complex(0.0, 0.0);

        for (Int i = 0; i < nAnt; i++) {
            feedAntId(iRow) = i;
            feedId(iRow) = 0;
            feedSpWId(iRow) = -1;
            feedBeamId(iRow) = 0;
            feedNumRec(iRow) = 2;
            beamOffset(0, 0, iRow) = 0.0;
            beamOffset(1, 0, iRow) = 0.0;
            beamOffset(0, 1, iRow) = 0.0;
            beamOffset(1, 1, iRow) = 0.0;
            feedXYZ(0, iRow) = 0.0;
            feedXYZ(1, iRow) = 0.0;
            feedXYZ(2, iRow) = 0.0;
            feedAngle(0, iRow) = 0.0;
            feedAngle(1, iRow) = 0.0;
            feedPol(0, iRow) = feedPol0;
            feedPol(1, iRow) = feedPol1;
            polResp(0, 0, iRow) = polResp(1, 1, iRow) = Complex(1.0, 0.0);
            iRow++;
        }
    }

    // fill Feed table - don't check to see if any of the positions match
    MSFeedColumns& feedc = msc.feed();
    Int numFeeds = feedc.nrow();
    Slicer feedSlice(IPosition(1, numFeeds), IPosition(1, nRow + numFeeds - 1),
                     IPosition(1, 1), Slicer::endIsLast);
    ms_p->feed().addRow(nRow);
    feedc.antennaId().putColumnRange(feedSlice, feedAntId);
    feedc.feedId().putColumnRange(feedSlice, feedId);
    feedc.spectralWindowId().putColumnRange(feedSlice, feedSpWId);
    feedc.beamId().putColumnRange(feedSlice, feedBeamId);
    feedc.numReceptors().putColumnRange(feedSlice, feedNumRec);
    feedc.position().putColumnRange(feedSlice, feedXYZ);
    const double forever = 1.e30;

    for (Int i = numFeeds; i < (nRow + numFeeds); i++) {
        feedc.beamOffset().put(i, beamOffset.xyPlane(i - numFeeds));
        feedc.polarizationType().put(i, feedPol.column(i - numFeeds));
        feedc.polResponse().put(i, polResp.xyPlane(i - numFeeds));
        feedc.receptorAngle().put(i, feedAngle.column(i - numFeeds));
        feedc.time().put(i, 0.0);
        feedc.interval().put(i, forever);
    }

    ASKAPLOG_INFO_STR(logger, "Added rows to FEED table");
}

Simulator::~Simulator()
{
}

Simulator & Simulator::operator=(const Simulator & other)
{
    if (this == &other) return *this;

    ASKAPTHROW(AskapError, "Copy constructor or assignment operator of Simulator are not supposed to be called");
    // copy state...
    return *this;
}

void Simulator::settimes(const casa::Quantity& qIntegrationTime,
                         const bool useHourAngle,
                         const casa::MEpoch& mRefTime)
{
    qIntegrationTime_p = qIntegrationTime;
    useHourAngle_p = useHourAngle;
    mRefTime_p = mRefTime;

    if (useHourAngle_p) {
        hourAngleDefined_p = False;
    }

    t_offset_p = 0.0;
}

void Simulator::observe(const casa::String& sourceName,
                        const casa::String& spWindowName,
                        const casa::Quantity& qStartTime,
                        const casa::Quantity& qStopTime)
{
    MSColumns msc(*ms_p);

    // Do we have antenna information?
    MSAntennaColumns& antc = msc.antenna();
    ASKAPCHECK(antc.nrow() > 0, "Antenna information not yet defined");

    Int nAnt = antc.nrow();
    Vector<double> antDiam;
    antc.dishDiameter().getColumn(antDiam);
    Matrix<double> antXYZ(3, nAnt);
    antc.position().getColumn(antXYZ);

    MSDerivedValues msd;
    msd.setAntennas(msc.antenna());

    // Do we have feed information?
    MSFeedColumns& feedc = msc.feed();
    ASKAPCHECK(feedc.nrow() > 0, "Feed information not yet defined");

    Int nFeed = feedc.nrow() / nAnt;

    // Spectral window
    MSSpWindowColumns& spwc = msc.spectralWindow();
    ASKAPCHECK(spwc.nrow() > 0, "Spectral window information not yet defined");

    Int baseSpWID = spwc.nrow();
    Int existingSpWID = -1;

    // Check for existing spectral window with correct name
    if (baseSpWID > 0) {
        Vector<String> spWindowNames;
        spwc.name().getColumn(spWindowNames);

        for (uInt i = 0; i < spWindowNames.nelements(); i++) {
            if (spWindowNames(i) == spWindowName) {
                existingSpWID = i;
                break;
            }
        }
    }

    ASKAPCHECK(existingSpWID > -1, "Spectral window named " + spWindowName + " not yet defined");

    MSPolarizationColumns& polc = msc.polarization();
    baseSpWID = existingSpWID;
    double startFreq, freqInc;
    Vector<double> resolution;
    spwc.refFrequency().get(baseSpWID, startFreq);
    spwc.resolution().get(baseSpWID, resolution);
    freqInc = resolution(0);
    Int nChan = resolution.nelements();
    Matrix<Int> corrProduct;
    polc.corrProduct().get(baseSpWID, corrProduct);
    Int nCorr = corrProduct.ncolumn();
    {
        ASKAPLOG_INFO_STR(logger, "Spectral window : " << spWindowName);
        ASKAPLOG_INFO_STR(logger, "   reference frequency : " << startFreq / 1.0e9 << "GHz");
        ASKAPLOG_INFO_STR(logger, "   number of channels : " << nChan);
        ASKAPLOG_INFO_STR(logger, "   total bandwidth : " << nChan*freqInc / 1.0e9 << "GHz");
        ASKAPLOG_INFO_STR(logger, "   number of correlations : " << nCorr);
    }
    
    if (itsChanBandwidthForNoise < -10) {
        // negative value less than -10 means that the first spectral window is processed
        itsChanBandwidthForNoise = fabs(freqInc);
    } else if (itsChanBandwidthForNoise > 0) {
        // check this spectral window for conformance
        if (fabs(itsChanBandwidthForNoise - fabs(freqInc))>0.1) {
            // this flag means that spectral windows which are observed had inhomogeneous resolutions
            itsChanBandwidthForNoise = -1;
        }
    }

    // Field
    MSFieldColumns& fieldc = msc.field();
    ASKAPCHECK(fieldc.nrow() > 0, "Field information not yet defined");

    Int baseFieldID = fieldc.nrow();
    Int existingFieldID = -1;

    // Check for existing field with correct name
    if (baseFieldID > 0) {
        Vector<String> fieldNames;
        fieldc.name().getColumn(fieldNames);

        for (uInt i = 0; i < fieldNames.nelements(); i++) {
            if (fieldNames(i) == sourceName) {
                existingFieldID = i;
                break;
            }
        }
    }

    ASKAPCHECK(existingFieldID > -1, "Field named " + sourceName + " not yet defined");

    baseFieldID = existingFieldID;
    Vector<MDirection> fcs(1);
    fieldc.phaseDirMeasCol().get(baseFieldID, fcs);
    msd.setFieldCenter(fcs(0));
    MDirection fieldCenter = fcs(0);
    {
        ASKAPLOG_INFO_STR(logger, "Observing source : " << sourceName
                              << "     direction : " << formatDirection(fieldCenter));
    }

    // A bit ugly solution to extract the information about beam offsets
    Cube<RigidVector<double, 2> > beam_offsets;
    Vector<String> antenna_mounts;
    { // to close MSIter, when the job is done
        MSFeedParameterExtractor msfpe_tmp(*ms_p);
        beam_offsets = msfpe_tmp.getBeamOffsets();
        antenna_mounts = msfpe_tmp.antennaMounts();
    }
    ASKAPCHECK(beam_offsets.nplane() == (uInt)nFeed && beam_offsets.ncolumn() == (uInt)nAnt,
               "Feed table format is incompatible with existing code of Simulator::observe");

    // Now we know where we are and where we are pointing, we can do the time calculations
    double Tstart, Tend, Tint;
    {
        Tint = qIntegrationTime_p.getValue("s");

        MEpoch::Ref tref(MEpoch::TAI);
        MEpoch::Convert tconvert(mRefTime_p, tref);
        MEpoch taiRefTime = tconvert();

        // until the qStartTime represents the starting Hour Angle
        if (useHourAngle_p && !hourAngleDefined_p) {
            msd.setEpoch(mRefTime_p);
            msd.setFieldCenter(fieldCenter);
            t_offset_p = - msd.hourAngle() * 3600.0 * 180.0 / C::pi / 15.0; // in seconds
            hourAngleDefined_p = True;
            ASKAPLOG_INFO_STR(logger, "Times specified are interpreted as hour angles for first source observed");
            ASKAPLOG_INFO_STR(logger, "     offset in time = " << t_offset_p / 3600.0 << " hours from "
                                  << formatTime(taiRefTime.get("s").getValue("s")));
        }

        Tstart = qStartTime.getValue("s") +
                 taiRefTime.get("s").getValue("s") + t_offset_p;
        Tend = qStopTime.getValue("s") +
               taiRefTime.get("s").getValue("s") + t_offset_p;
        ASKAPLOG_INFO_STR(logger, "Time range - start : " << formatTime(Tstart) << " stop  : " << formatTime(Tend));
    }

    // fill Observation Table for every call. Eventually we should fill
    // in the schedule information
    MSObservation& obs = ms_p->observation();
    MSObservationColumns& obsc = msc.observation();
    Int nobsrow = obsc.nrow();
    obs.addRow();
    obsc.telescopeName().put(nobsrow, telescope_p);
    Vector<double> timeRange(2);
    timeRange(0) = Tstart;
    timeRange(1) = Tend;
    obsc.timeRange().put(nobsrow, timeRange);
    obsc.observer().put(nobsrow, "ASKAP simulator");

    Int row = ms_p->nrow() - 1;
    Int maxObsId = -1;
    Int maxArrayId = 0;
    {
        Vector<Int> tmpids(row + 1);
        tmpids = msc.observationId().getColumn();

        if (tmpids.nelements() > 0) maxObsId = max(tmpids);

        tmpids = msc.arrayId().getColumn();

        if (tmpids.nelements() > 0) maxArrayId = max(tmpids);
    }

    double Time = Tstart;
    Bool firstTime = True;

    uInt nShadowed = 0;
    uInt nSubElevation = 0;

    // Start scan number from last one (if there was one)
    Int nMSRows = ms_p->nrow();

    // init counters past end
    Int scan = -1;

    if (nMSRows > 0) {
        msc.scanNumber().get(nMSRows - 1, scan);
    }

    // One call to observe corresponds to one scan
    scan++;

    // We can extend the ms just once
    Int nBaselines;

    if (autoCorrelationWt_p > 0.0) {
        nBaselines = nAnt * (nAnt + 1) / 2;
    } else {
        nBaselines = nAnt * (nAnt - 1) / 2;
    }

    Int nNewRows = nBaselines * nFeed;
    Int nIntegrations = max(1, Int(0.5 + (Tend - Tstart) / Tint));
    nNewRows *= nIntegrations;

    // We need to do addition in this order to get a new TSM file.

    // ... Next extend the table
    ASKAPLOG_INFO_STR(logger, "Adding " << nNewRows << " rows");
    ms_p->addRow(nNewRows);

    Matrix<Complex> data(nCorr, nChan);
    data.set(Complex(0.0));

    Matrix<Bool> flag(nCorr, nChan);
    flag = False;

    ASKAPLOG_INFO_STR(logger, "Calculating uvw coordinates for " << nIntegrations << " integrations");

    // Start of loop over time
    for (Int integration = 0; integration < nIntegrations; integration++) {
        MEpoch epUT1(Quantity(Time / C::day, "d"), MEpoch::UT1);
        MEpoch::Ref refGMST1(MEpoch::GMST1);
        MEpoch::Convert epGMST1(epUT1, refGMST1);
        double gmst = epGMST1().get("d").getValue("d");
        gmst = (gmst - Int(gmst)) * C::_2pi; // Into Radians

        MEpoch ep(Quantity((Time + Tint / 2), "s"));
        msd.setEpoch(ep);

        // current phase center for a beam without offset
        // For each individual beam pointing center always coincides
        // with the phase center

        // ???? May be we can use fcs defined earlier instead of fc ????
        MDirection fc = msc.field().phaseDirMeas(baseFieldID);
        msd.setFieldCenter(fc);
        msd.setAntenna(0); // assume for now that all par. angles are the same

        Vector<Bool> isShadowed(nAnt); isShadowed.set(False);
        Vector<Bool> isTooLow(nAnt); isTooLow.set(False);
        double fractionBlocked1 = 0.0, fractionBlocked2 = 0.0;
        Int startingRow = row;
        double diamMax2 = square(max(antDiam));

        // Start of loop over feed
        for (Int feed = 0; feed < nFeed; feed++) {
            // for now assume that all feeds have the same offsets w.r.t.
            // antenna frame for all antennas
            RigidVector<double, 2> beamOffset = beam_offsets(0, 0, feed);

            // fringe stopping center could be different for different feeds
            MDirection feed_phc = fc;

            // Do the first row outside the loop
            msc.scanNumber().put(row + 1, scan);
            msc.fieldId().put(row + 1, baseFieldID);
            msc.dataDescId().put(row + 1, baseSpWID);
            msc.time().put(row + 1, Time + Tint / 2);
            msc.arrayId().put(row + 1, maxArrayId);
            msc.processorId().put(row + 1, 0);
            msc.exposure().put(row + 1, Tint);
            msc.interval().put(row + 1, Tint);
            msc.observationId().put(row + 1, maxObsId + 1);
            msc.stateId().put(row + 1, -1);

            // assume also that all mounts are the same and posit. angle is the same
            if (antenna_mounts[0] == "ALT-AZ" || antenna_mounts[0] == "alt-az") {
                // parallactic angle rotation is necessary
                SquareMatrix<double, 2> xform(SquareMatrix<double, 2>::General);
                // SquareMatrix' default constructor is a bit strange, we probably
                // need to change it in the future


                const double pa = msd.parAngle();
                const double cpa = cos(pa);
                const double spa = sin(pa);
                xform(0, 0) = cpa;
                xform(1, 1) = cpa;
                xform(0, 1) = -spa;
                xform(1, 0) = spa;
                beamOffset *= xform;
            }

            // x direction is flipped to convert az-el type frame to ra-dec
            feed_phc.shift(-beamOffset(0), beamOffset(1), True);
            //            ASKAPLOG_DEBUG_STR(logger, "pointing/phase centre for beam="<<feed<<" is "<<
            //                               printDirection(feed_phc.getValue())<<" offsets: "<<beamOffset(0)/casa::C::pi*180<<" "<<
            //                               beamOffset(1)/casa::C::pi*180<<" mount="<<antenna_mounts[0]);

            double ra, dec; // current phase center
            ra = feed_phc.getAngle().getValue()(0);
            dec = feed_phc.getAngle().getValue()(1);

            // Transformation from antenna position difference (ant2-ant1) to uvw
            double H0 = gmst - ra, sH0 = sin(H0), cH0 = cos(H0), sd = sin(dec), cd = cos(dec);
            Matrix<double> trans(3, 3, 0);
            trans(0, 0) = -sH0; trans(0, 1) = -cH0;
            trans(1, 0) = sd * cH0; trans(1, 1) = -sd * sH0; trans(1, 2) = -cd;
            trans(2, 0) = -cd * cH0; trans(2, 1) = cd * sH0; trans(2, 2) = -sd;

            // Rotate antennas to correct frame
            Matrix<double> antUVW(3, nAnt);

            for (Int ant1 = 0; ant1 < nAnt; ant1++)
                antUVW.column(ant1) = product(trans, antXYZ.column(ant1));

            for (Int ant1 = 0; ant1 < nAnt; ant1++) {
                double x1 = antUVW(0, ant1), y1 = antUVW(1, ant1), z1 = antUVW(2, ant1);
                Int startAnt2 = ant1 + 1;

                if (autoCorrelationWt_p > 0.0) startAnt2 = ant1;

                for (Int ant2 = startAnt2; ant2 < nAnt; ant2++) {
                    row++;

                    msc.antenna1().put(row, ant1);
                    msc.antenna2().put(row, ant2);
                    msc.feed1().put(row, feed);
                    msc.feed2().put(row, feed);

                    double x2 = antUVW(0, ant2), y2 = antUVW(1, ant2), z2 = antUVW(2, ant2);
                    Vector<double> uvwvec(3);
                    uvwvec(0) = x2 - x1;
                    uvwvec(1) = y2 - y1;
                    uvwvec(2) = z2 - z1;
                    msc.uvw().put(row, uvwvec);

                    data.set(Complex(0., 0.));
                    msc.data().put(row, data);
                    msc.flag().put(row, flag);
                    msc.flagRow().put(row, False);

                    if (ant1 != ant2) {
                        blockage(fractionBlocked1, fractionBlocked2,
                                 uvwvec, antDiam(ant1), antDiam(ant2));

                        if (fractionBlocked1 > fractionBlockageLimit_p) {
                            isShadowed(ant1) = True;
                        }

                        if (fractionBlocked2 > fractionBlockageLimit_p) {
                            isShadowed(ant2) = True;
                        }
                    }
                    // case with variable Tsys/efficiency
                    double noiseRMS = itsNoiseRMS;
                    if (itsRelAntennaWeight.nelements() > 0) {
                        ASKAPCHECK(ant1 < casa::Int(itsRelAntennaWeight.nelements()), "encountered antenna index "<<ant1<<
                                   " which is beyond the array of Tsys and/or efficiencies (variable efficiency case)");
                        ASKAPCHECK(ant2 < casa::Int(itsRelAntennaWeight.nelements()), "encountered antenna index "<<ant2<<
                                   " which is beyond the array of Tsys and/or efficiencies (variable efficiency case)");
                        noiseRMS *= sqrt(itsRelAntennaWeight[ant1]*itsRelAntennaWeight[ant2]);
                    }

                    // Deal with differing diameter case
                    const Float sigma1 = diamMax2 / (antDiam(ant1) * antDiam(ant2)) * noiseRMS;
                    Float wt = 1 / square(sigma1);

                    if (ant1 == ant2) {
                        wt *= autoCorrelationWt_p;
                    }

                    Vector<Float> tmp(nCorr); tmp = wt;
                    msc.weight().put(row, tmp);
                    tmp = sigma1;
                    msc.sigma().put(row, tmp);
                }
            }

            // go back and flag weights based on shadowing
            // Future option: we could increase sigma based on
            // fraction shadowed.
            Matrix<Bool> trueFlag(nCorr, nChan);
            trueFlag = True;

            Int reRow = startingRow;

            for (Int ant1 = 0; ant1 < nAnt; ant1++) {
                Int startAnt2 = ant1 + 1;

                if (autoCorrelationWt_p > 0.0) startAnt2 = ant1;

                for (Int ant2 = startAnt2; ant2 < nAnt; ant2++) {
                    reRow++;

                    if (isShadowed(ant1) || isShadowed(ant2)) {
                        msc.flag().put(reRow, trueFlag);
                        msc.flagRow().put(reRow, True);
                        nShadowed++;
                    }
                }
            }

            // Find antennas pointing below the elevation limit
            Vector<double> azel(2);

            for (Int ant1 = 0; ant1 < nAnt; ant1++) {

                // We want to find elevation for each antenna separately (for VLBI)
                msd.setAntenna(ant1);
                azel = msd.azel().getAngle("rad").getValue("rad");

                if (azel(1) < elevationLimit_p.getValue("rad")) {
                    isTooLow(ant1) = True;
                }

                if (firstTime) {
                    firstTime = False;
                    double ha1 = msd.hourAngle() * 180.0 / C::pi / 15.0;
                    ASKAPLOG_INFO_STR(logger, "Starting conditions for antenna 1: ");
                    ASKAPLOG_INFO_STR(logger, "     time = " << formatTime(Time));
                    ASKAPLOG_INFO_STR(logger, "     scan = " << scan + 1);
                    ASKAPLOG_INFO_STR(logger, "     az   = " << azel(0) * 180.0 / C::pi << " deg");
                    ASKAPLOG_INFO_STR(logger, "     el   = " << azel(1) * 180.0 / C::pi << " deg");
                    ASKAPLOG_INFO_STR(logger, "     ha   = " << ha1 << " hours");
                }
            }

            // Now flag all antennas pointing below the elevation limit
            reRow = startingRow;

            for (Int ant1 = 0; ant1 < nAnt; ant1++) {
                Int startAnt2 = ant1 + 1;

                if (autoCorrelationWt_p > 0.0) startAnt2 = ant1;

                for (Int ant2 = startAnt2; ant2 < nAnt; ant2++) {
                    reRow++;

                    if (isTooLow(ant1) || isTooLow(ant2)) {
                        msc.flag().put(reRow, trueFlag);
                        msc.flagRow().put(reRow, True);
                        nSubElevation++;
                    }
                }
            }

            Int numpointrows = nAnt;
            MSPointingColumns& pointingc = msc.pointing();
            Int numPointing = pointingc.nrow();
            ms_p->pointing().addRow(numpointrows);
            numpointrows += numPointing;
            double Tint = qIntegrationTime_p.getValue("s");
            Vector<MDirection> direction(1);
            direction(0) = fieldCenter;

            for (Int m = numPointing; m < (numPointing + nAnt); m++) {
                pointingc.numPoly().put(m, 0);
                pointingc.interval().put(m, -1);
                pointingc.tracking().put(m, True);
                pointingc.time().put(m, Time);
                pointingc.timeOrigin().put(m, Tstart);
                pointingc.interval().put(m, Tint);
                pointingc.antennaId().put(m, m);
                pointingc.directionMeasCol().put(m, direction);
                pointingc.targetMeasCol().put(m, direction);
            }
        } // feeds

        Time += Tint;
    } // time ranges

    {
        msd.setAntenna(0);
        Vector<double> azel = msd.azel().getAngle("rad").getValue("rad");

        double ha1 = msd.hourAngle() * 180.0 / C::pi / 15.0;
        ASKAPLOG_INFO_STR(logger, "Stopping conditions for antenna 1: ");
        ASKAPLOG_INFO_STR(logger, "     time = " << formatTime(Time));
        ASKAPLOG_INFO_STR(logger, "     scan = " << scan + 1);
        ASKAPLOG_INFO_STR(logger, "     az   = " << azel(0) * 180.0 / C::pi << " deg");
        ASKAPLOG_INFO_STR(logger, "     el   = " << azel(1) * 180.0 / C::pi << " deg");
        ASKAPLOG_INFO_STR(logger, "     ha   = " << ha1 << " hours");
    }

    ASKAPLOG_INFO_STR(logger, (row + 1) << " visibilities simulated ");
    ASKAPLOG_INFO_STR(logger, nShadowed << " visibilities flagged due to shadowing ");
    ASKAPLOG_INFO_STR(logger, nSubElevation << " visibilities flagged due to elevation limit of " <<
                      elevationLimit_p.getValue("deg") << " degrees ");

}

// Calculates the fractional blockage of one antenna by another
// We will want to put this somewhere else eventually, but I don't yet know where!
// Till then.
// Stolen from Fred Schwab
void Simulator::blockage(double &fraction1, double &fraction2,
                         const casa::Vector<double>& uvw,
                         const double diam1,
                         const double diam2)
{
    double separation = sqrt(square(uvw(0)) + square(uvw(1)));
    double rmin = 0.5 * min(fabs(diam1), fabs(diam2));
    double rmax = 0.5 * max(fabs(diam1), fabs(diam2));

    if (separation >= (rmin + rmax)) {
        fraction1 = 0.0;
        fraction2 = 0.0;
    } else if ((separation + rmin) <= rmax) {
        fraction1 = min(1.0, square(fabs(diam2) / fabs(diam1)));
        fraction2 = min(1.0, square(fabs(diam1) / fabs(diam2)));
    } else {
        double c = separation / (0.5 * fabs(diam1));
        double s = fabs(diam2) / fabs(diam1);
        double sinb = sqrt(2.0 * (square(c * s) + square(c) + square(s)) - pow(c, 4.0) - pow(s, 4.0) - 1.0)
                      / (2.0 * c);
        double sina = sinb / s;
        //  Due to roundoff, sina or sinb might be ever so slightly larger than 1
        //  in the case of unequal radii, with the center of one antenna pattern
        //  inside the other:
        sinb = min(1.0, sinb);
        sina = min(1.0, sina);

        double b = asin(sinb);
        double a = asin(sina);
        double area = (square(s) * a + b) - (square(s) * sina * cos(a) + sinb * cos(b));
        fraction1 = area / C::pi;
        fraction2 = fraction1 / square(s);
    }

    // if antenna1 is in behind, w is > 0, 2 is NOT shadowed
    if (uvw(2) > 0.0) fraction2 = 0.0;

    // if antenna1 is in front, w is < 0, 1 is NOT shadowed
    if (uvw(2) < 0.0) fraction1 = 0.0;

    return;
}

String Simulator::formatDirection(const casa::MDirection& direction)
{
    MVAngle mvRa = direction.getAngle().getValue()(0);
    MVAngle mvDec = direction.getAngle().getValue()(1);
    ostringstream oss;
    oss.setf(ios::left, ios::adjustfield);
    oss.width(14);
    oss << mvRa(0.0).string(MVAngle::TIME, 8);
    oss.width(14);
    oss << mvDec.string(MVAngle::DIG2, 8);
    oss << "     " << MDirection::showType(direction.getRefPtr()->getType());
    return String(oss);
}

String Simulator::formatTime(const double time)
{
    MVTime mvtime(Quantity(time, "s"));
    return mvtime.string(MVTime::DMY, 7);
}

/// @brief return area times sqrt(bandwidth*int_time)
/// @details This quantity is used for automatic noise estimates. It is 
/// composed from itsChanBandwidthForNoise and itsDishDiamForNoise.
/// An exception is thrown if either array is inhomogeneous or 
/// multiple spectral resolutions are simulated. This method is supposed
/// to be called when the simulator is fully defined.
/// @return antenna area (m^2) multiplied by square root of the product of bandwidth(Hz) 
/// and integration time (s)
double Simulator::areaTimesSqrtBT() const
{
   ASKAPCHECK(itsDishDiamForNoise > 0., "Inhomogeneous antenna sizes have been detected (or sizes not defined at all). "
              "Automatic noise estimate is impossible. Please override with explicit value of 'rms' or 'variance'");
   ASKAPLOG_INFO_STR(logger, " Using antenna diameter of "<<itsDishDiamForNoise<<" m to estimate noise per visibility");
   ASKAPCHECK(itsChanBandwidthForNoise > 0., "Simulated observations contain multiple spectral resolutions "
              "(or they are not defined at all). "
              "Automatic noise estimate is impossible. Please override with explicit value of 'rms' or 'variance'");
   ASKAPLOG_INFO_STR(logger, " Using channel bandwidth of "<<itsChanBandwidthForNoise/1e3<<" kHz to estimate noise per visibility");

   const double inttime = qIntegrationTime_p.getValue("s");
   ASKAPCHECK(inttime > 0., "Integration time is supposed to be positive. You have "<<inttime<<" seconds");
   return casa::C::pi*casa::square(itsDishDiamForNoise)/4.*sqrt(itsChanBandwidthForNoise*inttime);
}


} // End namespace synthesis

} // End namespace askap
