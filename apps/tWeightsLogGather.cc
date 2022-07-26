/// @file tWeightsLogGather.cc
/// @details Run gather operation of the WeightsLog class in a pattern similar to that used in imager,
///          but in isolation. A number of parset parameters understood by imager are also understood by this tool
///          to setup communication classes the same way (whether it has any effect on it or not)
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

// ASKAPsoft includes
#include <askap/askap/Application.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/StatReporter.h>
#include <Common/ParameterSet.h>

// Local Package includes
#include "askap/distributedimager/CubeComms.h"

using namespace askap;
using namespace askap::cp;

ASKAP_LOGGER(logger, ".tWeightsLogGather");

class WeightsLogGatherTestApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            // Instantiate the comms class

            askap::cp::CubeComms comms(argc, const_cast<const char **>(argv));


            try {
                ASKAPCHECK(comms.isParallel(), "This test is only intended to be run as a parallel MPI job");
                // imager-specific configuration of the master/worker to allow groups of workers
                const int nWorkerGroups = config().getInt32("nworkergroups", 1);
                ASKAPCHECK(nWorkerGroups > 0, "nworkergroups is supposed to be greater than 0");
                if (nWorkerGroups > 1) {
                    ASKAPLOG_INFO_STR(logger, "Defining "<<nWorkerGroups<<" groups of workers");
                    comms.defineGroups(nWorkerGroups);
                } else {
                    ASKAPLOG_INFO_STR(logger, "All workers are treated as identical");
                }

                // setup and run the pattern

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
        std::string getVersion() const override {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

int main(int argc, char *argv[])
{
    WeightsLogGatherTestApp app;
    return app.main(argc, argv);
}
