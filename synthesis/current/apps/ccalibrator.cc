/// @file ccalibrator.cc
///
/// @brief Perform calibration and write result in the parset file
/// @details This application performs calibration of a measurement set
/// and writes the solution to an external parset file
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

// Package level header file
#include <askap_synthesis.h>

// ASKAPsoft includes
#include <askap/Application.h>
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/StatReporter.h>
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <parallel/CalibratorParallel.h>

ASKAP_LOGGER(logger, ".ccalibrator");

using namespace std;
using namespace askap;
using namespace askap::synthesis;

class CcalibratorApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;
                LOFAR::ParameterSet subset(config().makeSubset("Ccalibrator."));

                // Perform %w substitutions for all keys
                for (LOFAR::ParameterSet::iterator it = subset.begin();
                        it != subset.end(); ++it) {
                    it->second = LOFAR::ParameterValue(comms.substitute(it->second));
                }

                CalibratorParallel calib(comms, subset);
                ASKAPLOG_INFO_STR(logger, "ASKAP synthesis calibrator " << ASKAP_PACKAGE_VERSION);

                if (comms.isMaster()) {
                    ASKAPLOG_INFO_STR(logger, "Parset file contents:\n" << config());
                }

                const int nCycles = subset.getInt32("ncycles", 1);
                ASKAPCHECK(nCycles >= 0, " Number of calibration iterations should be a non-negative number, you have " <<
                        nCycles);

                size_t solution = 0;
                for (bool continueFlag = true; continueFlag; ++solution) {
                    ASKAPLOG_INFO_STR(logger, "Calibration solution interval "<<solution + 1);
                    for (int cycle = 0; cycle < nCycles; ++cycle) {
                        ASKAPLOG_INFO_STR(logger, "*** Starting calibration iteration " << cycle + 1 << " ***");
                        calib.broadcastModel();
                        calib.receiveModel();
                        calib.calcNE();
                        calib.solveNE();
                        stats.logSummary();
                    }

                    ASKAPLOG_INFO_STR(logger,  "*** Finished calibration cycles ***");
                    calib.writeModel();

                    continueFlag = calib.getNextChunkFlag();
                    if (continueFlag) {
                        ASKAPLOG_INFO_STR(logger, "More data are available, continue to make solution for the next interval");
                        // initialise the model and measurement equation
                        calib.init(subset);
                    } else {
                        ASKAPLOG_INFO_STR(logger, "No more data are available, this was the last solution interval");
                    }

                    // Remove the next chunk flag since merge will not update the
                    // value if it already exists. This is important for the graph
                    // reduction of normal equations. This is really only needed for
                    // the master, but doesn't hurt at the worker.
                    calib.removeNextChunkFlag();
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
    CcalibratorApp app;
    return app.main(argc, argv);
}
