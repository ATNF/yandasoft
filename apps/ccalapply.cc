/// @file ccalapply.cc
///
/// @brief Calibration applicator
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

// Include package level header file
#include <askap/askap_synthesis.h>

// System include
#include <string>
#include <iostream>
#include <stdlib.h>

// ASKAPsoft includes
#include <askap/askap/Application.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/StatReporter.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/IConstDataSource.h>
#include <askap/dataaccess/IFlagDataAccessor.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/calibaccess/ICalSolutionConstSource.h>
#include <askap/calibaccess/CalibAccessFactory.h>
#include <askap/calibaccess/ServiceCalSolutionSourceStub.h>
#include <askap/calibaccess/ChanAdapterCalSolutionConstSource.h>

#ifdef USE_CAL_SERVICE
#include <askap/calserviceaccessor/ServiceCalSolutionSource.h>
#endif

#include <askap/dataaccess/OnDemandNoiseAndFlagDA.h>
#include <askap/askap/RangePartition.h>
#include <boost/shared_ptr.hpp>

// Local packages includes
#include <askap/measurementequation/ICalibrationApplicator.h>
#include <askap/measurementequation/CalibrationApplicatorME.h>
#include <askap/measurementequation/CalibrationIterator.h>
#include <askap/parallel/ParallelWriteIterator.h>
#include <askap/utils/TilingUtils.h>

// casacore includes
#include <casacore/casa/OS/Timer.h>


ASKAP_LOGGER(logger, ".ccalapply");

// Using
using namespace std;
using namespace askap;
using namespace askap::accessors;
using namespace askap::synthesis;

class CcalApplyApp : public askap::Application
{
    public:
        int run(int argc, char* argv[]) final
        {
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;
                LOFAR::ParameterSet subset(config().makeSubset("Ccalapply."));

                // distribute by channel with master writing if parallel except if
                // we have specified Tile based distribution instead
                itsDistribute = comms.isParallel() ? subset.getBool("distribute", !subset.isDefined("Tiles")) : false;

                if (itsDistribute) {
                    ASKAPLOG_INFO_STR(logger, "Data will be distributed between "<<(comms.nProcs() - 1)<<" workers, master will write data");
                }

                if (comms.isWorker() || !itsDistribute) {

                    // Get Measurement Set accessor
                    // could've passed optional parameters for the parallel iterator too, but it's fine tuning
                    IDataSharedIter it = itsDistribute ? IDataSharedIter(new ParallelWriteIterator(comms)) : getDataIterator(subset, comms);

                    // Setup calibration applicator
                    const boost::shared_ptr<ParallelWriteIterator> pwIt = it.dynamicCast<ParallelWriteIterator>();
                    boost::shared_ptr<ICalibrationApplicator> calME = buildCalApplicator(subset, pwIt ? pwIt->chanOffset() : 0u);

                    ASKAPDEBUGASSERT(it);
                    ASKAPDEBUGASSERT(calME);

                    casa::Timer timer;

                    // Apply calibration
                    uint64_t count = 1;
                    double calculationTime = 0.;
                    for (it.init(); it != it.end(); it.next()) {
                        if (count % 100 == 0) {
                            ASKAPLOG_DEBUG_STR(logger, "Progress - Chunk " << count<< " nRows = "<<it->nRow());
                        }
                        ++count;
                        timer.mark();

                        if (itsNoiseAndFlagDANeeded) {
                            // quick and dirty for now (mv: note, it will definitely ignore updates to noise and may ignore
                            // updates to flags as well, if the appropriate accessor doesn't support write operation!)
                            accessors::OnDemandNoiseAndFlagDA acc(*it);
                            acc.rwVisibility() = it->visibility();

                            calME->correct(acc);

                            it->rwVisibility() = acc.rwVisibility();

                            const boost::shared_ptr<IFlagDataAccessor> fda = boost::dynamic_pointer_cast<IFlagDataAccessor>(boost::shared_ptr<IDataAccessor>(it.operator->(), utility::NullDeleter()));
                            ASKAPCHECK(fda, "Data accessor is of type which does not support overwritting flag information");
                            fda->rwFlag() = acc.rwFlag();
                        } else {
                            calME->correct(*it);
                        }
                        calculationTime += timer.real();
                    }
                    ASKAPLOG_INFO_STR(logger, "Time spent in calculation and data movement (but excluding I/O): "<<calculationTime<<" seconds");
                } else if (itsDistribute) {
                    // server code. Note, noise is not propagated back (as for the serial case)
                    ParallelWriteIterator::masterIteration(comms, getDataIterator(subset, comms), ParallelWriteIterator::SYNCFLAG | ParallelWriteIterator::READ);
                }
                stats.logSummary();
            } catch (const askap::AskapError& x) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
                std::cerr << "Askap error in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            } catch (const std::exception& x) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            }

            return 0;
        }

    private:

        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
        /// @brief this flag indicates whether the underlying code needs to update flags or noise
        /// @details For now - quick and dirty fix to allow MRO tests to proceed. The ASKAP model
        /// is to apply calibration on the fly, so table accessor classes are unable to modify noise or flag
        /// information.
        bool itsNoiseAndFlagDANeeded;

        /// @brief this flag indicates that automatic distribution of data in channel space take space
        /// @details The logic in setting up adapters, etc depends on whether the data are already distributed or
        /// have to be split and distributed by this class. This flag controls the approach and is on by default.
        /// It's ignored in the serial mode.
        bool itsDistribute;

        static casa::MFrequency::Ref getFreqRefFrame(const LOFAR::ParameterSet& parset)
        {
            const string freqFrame = parset.getString("freqframe", "topo");
            if (freqFrame == "topo") {
                ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as topocentric");
                return casa::MFrequency::Ref(casa::MFrequency::TOPO);
            } else if (freqFrame == "lsrk") {
                ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as lsrk");
                return casa::MFrequency::Ref(casa::MFrequency::LSRK);
            } else if (freqFrame == "bary") {
                ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as barycentric");
                return casa::MFrequency::Ref(casa::MFrequency::BARY);
            } else {
                ASKAPTHROW(AskapError, "Unsupported frequency frame " << freqFrame);
            }
        }

        boost::shared_ptr<ICalibrationApplicator> buildCalApplicator(
                const LOFAR::ParameterSet& parset, const casa::uInt chanOffset = 0u)
        {
            // Create solution source
            ICalSolutionConstSource::ShPtr solutionSource =
                CalibAccessFactory::roCalSolutionSource(parset);
            ASKAPASSERT(solutionSource);

            // This is sloppy but I need to test whether this is likely to be a service
            // source as I need to reinstantiate the full implementation - as all we get from the factory
            // is a stub.

            const std::string calAccType = parset.getString("calibaccess","parset");

#ifdef USE_CAL_SERVICE
            if (calAccType == "service") {
              solutionSource.reset(new ServiceCalSolutionSource(parset));
              ASKAPLOG_INFO_STR(logger,"Yay I am a service source");

            }
            else {
              ASKAPLOG_INFO_STR(logger,"Boo I am not a service source!");

            }
#endif

            // wrap the solution source in an adapter, if necessary
            if (chanOffset != 0u) {
                ASKAPLOG_DEBUG_STR(logger, "Setup an adapter to adjust channels by "<<chanOffset);
                solutionSource.reset(new ChanAdapterCalSolutionConstSource(solutionSource, chanOffset));
            }

            // Create applicator
            boost::shared_ptr<ICalibrationApplicator> calME(new CalibrationApplicatorME(solutionSource));
            ASKAPASSERT(calME);
            const bool scaleNoise = parset.getBool("calibrate.scalenoise", false);
            const bool allowFlag = parset.getBool("calibrate.allowflag", false);
            itsNoiseAndFlagDANeeded = scaleNoise || allowFlag;
            calME->scaleNoise(scaleNoise);
            calME->allowFlag(allowFlag);
            calME->beamIndependent(parset.getBool("calibrate.ignorebeam", false));
            calME->channelIndependent(parset.getBool("calibrate.ignorechannel",false));
            calME->interpolateTime(parset.getBool("calibrate.interpolatetime",false));
            return calME;
        }

        IDataSharedIter getDataIterator(const LOFAR::ParameterSet& parset, const askap::askapparallel::AskapParallel &comms) const
        {
            // do substitutions
            LOFAR::ParameterSet subset(parset.makeSubset(""));
            if (subset.isDefined("Tiles")) {
                subset.replace("Tiles",comms.substitute(subset.getString("Tiles")));
            }
            const string ms = itsDistribute ? subset.getString("dataset") : comms.substitute(subset.getString("dataset"));

            TableDataSource ds(ms, TableDataSource::MEMORY_BUFFERS | TableDataSource::WRITE_DATA_ONLY);
            IDataSelectorPtr sel=ds.createSelector();
            sel << subset;
            IDataConverterPtr conv=ds.createConverter();
            conv->setFrequencyFrame(getFreqRefFrame(subset), "Hz");
            conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
            // ensure that time is counted in seconds since 0 MJD
            conv->setEpochFrame();

            // do automatic distribution over tiles if requested and
            //   if not distributing in other ways (by channel or ms)
            const bool distributeByTile = !itsDistribute && comms.isParallel() && subset.isDefined("Tiles") &&
                subset.getString("Tiles")=="auto" &&
                subset.getString("dataset")==comms.substitute(subset.getString("dataset"));
            if (distributeByTile) {
                utils::distributeByTile(sel, "DATA", comms.nProcs(), comms.rank());
            }

            if (subset.isDefined("maxchunkrows")) {
                const casa::uInt maxChunkSize = subset.getUint32("maxchunkrows");
                ASKAPCHECK(maxChunkSize > 0, "maxchunkrows parameter should be positive");
                ASKAPLOG_INFO_STR(logger, "Restricting the chunk size to at most "<<maxChunkSize<<" rows for each iteration");
                ds.configureMaxChunkSize(maxChunkSize);
            }

            return ds.createIterator(sel, conv);
        }
};

int main(int argc, char *argv[])
{
    CcalApplyApp app;
    return app.main(argc, argv);
}
