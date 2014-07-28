/// @file cimager.cc
///
/// @breif Synthesis imaging program
///
/// Performs synthesis imaging from a data source, using any of a number of
/// image solvers. Can run in serial or parallel (MPI) mode.
///
/// The data are accessed from the DataSource. This is and will probably remain
/// disk based. The images are kept purely in memory until the end.
///
/// Control parameters are passed in from a LOFAR ParameterSet file.
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>

// Package level header file
#include "askap_synthesis.h"

// System includes
#include <stdexcept>
#include <iostream>
#include <csignal>

// ASKAPsoft includes
#include "askap/Application.h"
#include "askap/AskapLogging.h"
#include "askap/AskapError.h"
#include "askap/SignalManagerSingleton.h"
#include "askap/SignalCounter.h"
#include "askap/StatReporter.h"
#include "boost/scoped_ptr.hpp"
#include <parallel/ImagerParallel.h>
#include <measurementequation/MEParsetInterface.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <fitting/Params.h>
#include <profile/AskapProfiler.h>


ASKAP_LOGGER(logger, ".cimager");

using namespace askap;
using namespace askap::synthesis;
using namespace askap::scimath;

class CimagerApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;
                LOFAR::ParameterSet subset(config().makeSubset("Cimager."));

                boost::scoped_ptr<askap::ProfileSingleton::Initialiser> profiler;
                if (parameterExists("profile")) {
                    std::string profileFileName("profile.cimager");
                    if (comms.isParallel()) {
                        profileFileName += ".rank"+utility::toString(comms.rank());
                    }
                    profiler.reset(new askap::ProfileSingleton::Initialiser(profileFileName));
                }

                // Put everything in scope to ensure that all destructors are called
                // before the final message
                {
                    // imager-specific configuration of the master/worker to allow groups of workers
                    const int nWorkerGroups = subset.getInt32("nworkergroups", 1);
                    ASKAPCHECK(nWorkerGroups > 0, "nworkergroups is supposed to be greater than 0");
                    if (nWorkerGroups > 1) {
                        ASKAPLOG_INFO_STR(logger, "Model parameters will be distributed between "<<nWorkerGroups<<
                                " groups of workers");
                        ASKAPCHECK(comms.isParallel(), "This option is only allowed in the parallel mode");
                        comms.defineGroups(nWorkerGroups);
                    } else {
                        ASKAPLOG_INFO_STR(logger, "All workers are treated as identical");
                    }   

                    // Perform %w substitutions for all keys.
                    // NOTE: This MUST happen after AskapParallel::defineGroups() is called
                    for (LOFAR::ParameterSet::iterator it = subset.begin();
                            it != subset.end(); ++it) {
                        it->second = LOFAR::ParameterValue(comms.substitute(it->second));
                    }

                    // new parset that includes any missing parameters that can be advised using metadata
                    // NOTE: This MUST happen after the %w substitutions. (Move both to act on makeSubset output above?)
                    LOFAR::ParameterSet fullset(ImagerParallel::autoSetParameters(comms, subset));

                    ImagerParallel imager(comms, fullset);
                    ASKAPLOG_INFO_STR(logger, "ASKAP synthesis imager " << ASKAP_PACKAGE_VERSION);

                    if (comms.isMaster()) {
                        LOFAR::ParameterSet tmpset = config();
                        tmpset.subtractSubset("Cimager."); // replace the original Cimager param set with the updated set.
                        tmpset.adoptCollection(fullset.makeSubset("","Cimager."));
                        ASKAPLOG_INFO_STR(logger, "Parset parameters:\n" << tmpset);
                    }

                    const double targetPeakResidual = SynthesisParamsHelper::convertQuantity(
                            fullset.getString("threshold.majorcycle", "-1Jy"), "Jy");
                    const bool writeAtMajorCycle = fullset.getBool("Images.writeAtMajorCycle", false);


                    const int nCycles = fullset.getInt32("ncycles", 0);
                    if (nCycles == 0) {
                        /// No cycling - just make a dirty image
                        imager.broadcastModel();
                        imager.receiveModel();
                        imager.calcNE();
                        imager.receiveNE();
                        //imager.solveNE();
                        //imager.zeroAllModelImages();
                    } else {
                        // Set up a new signal handler for SIGUSR1.
                        // This allows graceful exit from the major cycle loop
                        // upon receipt of a SIGUSR1.
                        SignalCounter sigcount;
                        SignalManagerSingleton::instance()->registerHandler(SIGUSR1, &sigcount);

                        /// Perform multiple major cycles
                        for (int cycle = 0; cycle < nCycles; ++cycle) {
                            imager.broadcastModel();
                            imager.receiveModel();
                            ASKAPLOG_INFO_STR(logger, "*** Starting major cycle " << cycle << " ***");
                            imager.calcNE();
                            imager.solveNE();

                            stats.logSummary();

                            if (comms.isMaster()) {
                                if (sigcount.getCount() > 0) {
                                    ASKAPLOG_INFO_STR(logger, "Signal SIGUSR1 receieved. Stopping.");
                                    break;
                                }

                                if (imager.params()->has("peak_residual")) {
                                    const double peak_residual = imager.params()->scalarValue("peak_residual");
                                    ASKAPLOG_INFO_STR(logger, "Reached peak residual of " << peak_residual);

                                    if (peak_residual < targetPeakResidual) {
                                        ASKAPLOG_INFO_STR(logger, "It is below the major cycle threshold of "
                                                << targetPeakResidual << " Jy. Stopping.");
                                        break;
                                    } else {
                                        if (targetPeakResidual < 0) {
                                            ASKAPLOG_INFO_STR(logger, "Major cycle flux threshold is not used.");
                                        } else {
                                            ASKAPLOG_INFO_STR(logger, "It is above the major cycle threshold of "
                                                    << targetPeakResidual << " Jy. Continuing.");
                                        }
                                    }
                                }
                            }

                            if (cycle + 1 >= nCycles) {
                                ASKAPLOG_INFO_STR(logger, "Reached " << nCycles << " cycle(s), the maximum number of major cycles. Stopping.");
                            }

                            if (writeAtMajorCycle) {
                                imager.writeModel(std::string(".majorcycle.") + utility::toString(cycle + 1));
                            }
                        }

                        imager.broadcastModel();
                        imager.receiveModel();
                        ASKAPLOG_INFO_STR(logger, "*** Finished major cycles ***");
                        imager.calcNE();
                        imager.receiveNE();
                    }
                    SignalManagerSingleton::instance()->removeHandler(SIGUSR1);

                    /// This is the final step - restore the image and write it out
                    imager.writeModel();
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
};

int main(int argc, char *argv[])
{
    CimagerApp app;
    app.addParameter("profile", "p", "Write profiling output files", false);
    return app.main(argc, argv);
}

