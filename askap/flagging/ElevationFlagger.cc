/// @file ElevationFlagger.cc
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

// Include own header file first
#include "ElevationFlagger.h"

// Include package level header file
#include "askap/askap_synthesis.h"

// ASKAPsoft includes
#include "askap/askap/AskapLogging.h"
#include "askap/dataaccess/TableConstDataIterator.h"
#include "askap/dataaccess/IFlagDataAccessor.h"
#include "casacore/ms/MSOper/MSDerivedValues.h"

// Local package includes
#include "askap/flagging/FlaggingStats.h"

ASKAP_LOGGER(logger, ".ElevationFlagger");

using namespace askap;
using namespace casacore;
using namespace askap::synthesis;
using namespace askap::accessors;

vector<std::shared_ptr<IFlagger> > ElevationFlagger::build(
        const LOFAR::ParameterSet& parset)
{
    vector<std::shared_ptr<IFlagger> > flaggers;
    const string key = "elevation_flagger.enable";
    if (parset.isDefined(key) && parset.getBool(key)) {
        const LOFAR::ParameterSet subset = parset.makeSubset("elevation_flagger.");
        flaggers.push_back(std::shared_ptr<IFlagger>(new ElevationFlagger(subset)));
    }
    return flaggers;
}

ElevationFlagger:: ElevationFlagger(const LOFAR::ParameterSet& parset)
        : itsStats("ElevationFlagger"),
        itsHighLimit(parset.getFloat("high", 90.0), "deg"),
        itsLowLimit(parset.getFloat("low", 0.0), "deg"),
        itsTimeElevCalculated(0.0)
{
}

FlaggingStats ElevationFlagger::stats(void) const
{
    return itsStats;
}

casacore::Bool ElevationFlagger::processingRequired(const casacore::uInt pass) const
{
    return (pass==0);
}

void ElevationFlagger::updateElevations(const IDataSharedIter& di)
{
    // 1: Ensure the antenna elevation array is the correct size
    boost::shared_ptr<TableConstDataIterator> tdi =
        di.dynamicCast<TableConstDataIterator>();

    const casacore::uInt nAnt= tdi->subtableInfo().getAntenna().getNumberOfAntennas();
    if (itsAntennaElevations.size() != nAnt) {
        itsAntennaElevations.resize(nAnt);
    }

    // 2: Setup MSDerivedValues with antenna positions, field direction, and date/time
    MSDerivedValues msd;
    Vector<String> mount(nAnt);
    Vector<MPosition> pos(nAnt);
    for (uInt ant = 0; ant < nAnt; ant++) {
        mount(ant) = tdi->subtableInfo().getAntenna().getMount(ant);
        pos(ant) = tdi->subtableInfo().getAntenna().getPosition(ant);
    }
    //msd.setAntennaMount(mount);
    msd.setAntennaPositions(pos);
    msd.setEpoch(MEpoch(Quantity(di->time(),Unit("s"))));
    //ASKAPLOG_INFO_STR(logger,"pos="<<pos(0)<<" time="<<MEpoch(Quantity(di->time(),Unit("s"))));

    const casacore::Int fieldId = tdi->currentFieldID();
    msd.setFieldCenter(tdi->subtableInfo().getField().getReferenceDir(fieldId));
    // The following might be more accurate (feed pointing), but cflag version uses field dir.
    // msd.setFieldCenter(di->pointingDir1()(0));
    //ASKAPLOG_INFO_STR(logger,"fieldcenter="<<MDirection(di->pointingDir1()(0),MDirection::J2000).toString()<<" fieldid="<<fieldId<<" refdir="<<tdi->subtableInfo().getField().getReferenceDir(fieldId).toString());
    // 3: Calculate elevations for all antennas. Calculate each antenna
    // individually in case very long baselines exist.
    for (casacore::uInt i = 0 ; i < nAnt; ++i) {
        msd.setAntenna(i);
        const Vector<double> azel = msd.azel().getAngle("deg").getValue("deg");
        itsAntennaElevations(i) = Quantity(azel(1), "deg");
    }

    itsTimeElevCalculated = di->time();
}

void ElevationFlagger::processRows(const IDataSharedIter& di,
                         const casacore::Vector<bool>& rowFlag,
                         const casacore::uInt pass, const bool dryRun)
{
    updateElevations(di);
    casacore::uInt nRow = di->nRow();
    IFlagDataAccessor &fda=dynamic_cast<IFlagDataAccessor&>(*di);
    // 2: Do flagging
    for (casacore::uInt row = 0; row < nRow; row++) {
        if (!rowFlag(row)) {
            const int ant1 = di->antenna1()(row);
            const int ant2 = di->antenna2()(row);
            if (itsAntennaElevations(ant1) < itsLowLimit ||
                itsAntennaElevations(ant2) < itsLowLimit ||
                itsAntennaElevations(ant1) > itsHighLimit ||
                itsAntennaElevations(ant2) > itsHighLimit)
            {
                flagRow(fda.rwFlag(), row, dryRun);
            }
        }
    }
}

void ElevationFlagger::flagRow(casacore::Cube<casacore::Bool>& flag, const casacore::uInt row, const bool dryRun)
{

    //itsStats.visFlagged += flag.shape()(1) * flag.shape()(2);
    itsStats.visFlagged += flag.shape()(1) * flag.shape()(0);
    itsStats.rowsFlagged++;

    if (!dryRun) {
        //flag(casacore::Slice(row),casacore::Slice(),casacore::Slice()) = true;
        flag(casacore::Slice(),casacore::Slice(),casacore::Slice(row)) = true;
    }
}
