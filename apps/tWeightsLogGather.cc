/// @file tWeightsLogGather.cc
/// @details Run gather operation of the WeightsLog class in a pattern similar to that used in imager,
///          but in isolation. A number of parset parameters understood by imager are also understood by this tool
///          to setup communication classes the same way (whether it has any effect on it or not)
///
///    Example of the parset file:
///
///    # standard options which affect general setup of the distribution (as with imager)
///    # note, gathering is only done in the case of multiple writers and singleoutputfile=true, so
///    # the test becomes useless otherwise
///    nwriters = 2
///    nchanpercore = 1
///    singleoutputfile = true
///    nworkergroups = 1
///    # default delay before gather operation is called in seconds (set to 0 for no delay)
///    delay.default = 10
///    # delay for rank XX in seconds, overrides delay.default (set to 0 for no delay, omit the parameter to use the default one)
///    delay.rankXX = 0
///    # number of times the test is run
///    ncycles = 1
///    # if true, there is a barrier before the each cycle
///    cyclebarrier = false
///
/// @copyright (c) 2009,2016 CSIRO
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

// Include package level header file
#include <askap/askap_synthesis.h>

// System includes
#include <string>
#include <map>

// ASKAPsoft includes
#include <askap/askap/Application.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/StatReporter.h>
#include <Common/ParameterSet.h>
#include <askap/imageaccess/WeightsLog.h>

// Local Package includes
#include "askap/distributedimager/CubeComms.h"

using namespace askap;
using namespace askap::cp;

ASKAP_LOGGER(logger, ".tWeightsLogGather");

class WeightsLogGatherTestApp : public askap::Application
{
    public:
        void runActualTest(askap::cp::CubeComms &comms) {
             const int nwriters = config().getInt32("nwriters",1);
             ASKAPCHECK(nwriters>0,"Number of writers must be greater than 0");

             askap::accessors::WeightsLog weightslog;
             ASKAPLOG_INFO_STR(logger, "About to add weights list of size " << itsWeightsList.size() << " to the weights logger");
             weightslog.weightslist() = itsWeightsList;

             // MV: note, the code was essentially taken from imager to emulate the same distribution/same communication pattern
             // However, the gathering operation doesn't happen if nwriters == 1 or singleoutputfile is false in the parset rendering the whole
             // test useless.
             if (nwriters > 1 && config().getBool("singleoutputfile",false)) {
                 std::list<int> creators = comms.getCubeCreators();
                 ASKAPASSERT(creators.size() == 1);
                 const int creatorRank = creators.front();
                 ASKAPLOG_INFO_STR(logger, "Gathering all weights information, creator is rank " << creatorRank);
                 weightslog.gather(comms, creatorRank, false);
             }
             ASKAPLOG_INFO_STR(logger, "Finished the test, weights list has "<<weightslog.weightslist().size()<<" elements");
        }

        void prepareSomeWeightsList(askap::cp::CubeComms &comms) {

             // do the same math & logging as the imager worker
             // although we don't necessarily need this information, it may be handy
             // to be able to fill the list to test depending on this info
             const int nchanpercore = config().getInt32("nchanpercore", 1);

             // lets calculate a base
             const unsigned int nWorkers = comms.nProcs() - 1;
             const unsigned int nWorkersPerGroup = nWorkers / comms.nGroups();
             ASKAPCHECK(nWorkersPerGroup > 0, "Expect at least one worker per group, nWorkers="<<nWorkers<<" number of groups: "<<comms.nGroups());

             const unsigned int id = comms.rank();
             // e. g. rank 8, 3 per group should be pos. 1 (zero index)
             unsigned int posInGroup = (id % nWorkersPerGroup);

             if (posInGroup == 0) {
                 posInGroup = nWorkersPerGroup;
             }
             posInGroup = posInGroup - 1;

             const unsigned int baseChannel = posInGroup * nchanpercore;
             ASKAPLOG_INFO_STR(logger, "Distribution: Id " << id << " nWorkers " << nWorkers << " nGroups " << comms.nGroups()<<
                               " base channel " << baseChannel << " PosInGrp " << posInGroup);
             // now fill some weights
             for (unsigned int chan = 0; chan < nchanpercore; ++chan) {
                  itsWeightsList[chan+baseChannel] = 0.1 * id;
             }
             const unsigned int defaultDelay = config().getUint32("delay.default",10);
             const std::string delayParamStr = "delay.rank"+utility::toString(id);
             const unsigned int thisRankDelay = config().isDefined(delayParamStr) ? config().getUint32(delayParamStr) : defaultDelay;
             if (thisRankDelay > 0) {
                 ASKAPLOG_INFO_STR(logger, "Pausing for "<<thisRankDelay<<" seconds");
                 sleep(thisRankDelay);
             }
        }

        int run(int argc, char* argv[]) final
        {
            // Instantiate the comms class

            askap::cp::CubeComms comms(argc, const_cast<const char **>(argv));


            try {
                ASKAPCHECK(comms.isParallel(), "This test is only intended to be run as a parallel MPI job");
                // imager-specific configuration of the master/worker to allow groups of workers
                // MV: this is copied as is from the imager application, although I think it is a legacy of the old imager
                // and things are no longer rely on groups set up this way
                const int nWorkerGroups = config().getInt32("nworkergroups", 1);
                ASKAPCHECK(nWorkerGroups > 0, "nworkergroups is supposed to be greater than 0");
                if (nWorkerGroups > 1) {
                    ASKAPLOG_INFO_STR(logger, "Defining "<<nWorkerGroups<<" groups of workers");
                    comms.defineGroups(nWorkerGroups);
                } else {
                    ASKAPLOG_INFO_STR(logger, "All workers are treated as identical");
                }

                // this happens inside the ContinuumImager class
                comms.buildCommIndex();

                const unsigned int nwriters = config().getUint32("nwriters",1);
                const unsigned int nchanpercore = config().getUint32("nchanpercore", 1);
                // not sure where this is done in imager and whether this is necessary, but writers / cube creators are not
                // initialised and an exception is thrown if the following is not done
                comms.initWriters(nwriters,nchanpercore);
                comms.setSingleSink();

                // setup and run the pattern in workers
                const unsigned int ncycles = config().getUint32("ncycles",1);
                ASKAPCHECK(ncycles > 0, "Number of cycles should be at least one");
                const bool doBarrier = config().getBool("cyclebarrier",false);
                for (unsigned int cycle = 0; cycle < ncycles; ++cycle) {
                     ASKAPLOG_INFO_STR(logger, "Starting cycle "<<cycle + 1);
                     itsWeightsList.clear();
                     if (doBarrier) {
                         ASKAPLOG_INFO_STR(logger, "Barrier before proceeding to the test");
                         comms.barrier();
                     }
                     if (comms.isWorker()) {
                         prepareSomeWeightsList(comms);
                         runActualTest(comms);
                    }
                }

                ASKAPLOG_INFO_STR(logger, "Reached the barrier at the end (not sure if it is needed, but the imager has it)");
                comms.barrier();
                ASKAPLOG_INFO_STR(logger, "Passed the barrier at the end (not sure if it is needed, but the imager has it)");
            } catch (const askap::AskapError& e) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << e.what());
                std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
                comms.abort();
                return 1;
            } catch (const std::exception& e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << e.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what()
                    << std::endl;
                comms.abort();
                return 1;
            }
            return 0;
        }

    private:
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }

        /// @brief weights list to work with
        std::map<unsigned int, float> itsWeightsList;
};

int main(int argc, char *argv[])
{
    WeightsLogGatherTestApp app;
    return app.main(argc, argv);
}
