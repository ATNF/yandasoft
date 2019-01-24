/// @file imager.cc
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>
/// @author Stephen Ord <stephen.ord@csiro.au>

// Include package level header file
#include <askap_synthesis.h>

// System includes
#include <string>
#include <fstream>
#include <sstream>
// #include <mpi.h>

// Boost includes
#include <boost/scoped_ptr.hpp>

// ASKAPsoft includes
#include "askap/Application.h"
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/StatReporter.h>
#include <Common/ParameterSet.h>
#include <parallel/ImagerParallel.h>
#include <profile/AskapProfiler.h>

// Local Package includes
#include "distributedimager/ContinuumImager.h"
#include "distributedimager/CubeComms.h"

using namespace askap;
using namespace askap::cp;

ASKAP_LOGGER(logger, ".main");

class ImagerApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            // Instantiate the comms class

            askap::cp::CubeComms comms_p(argc, const_cast<const char **>(argv));


            try {

                StatReporter stats;

                // Create a subset

                LOFAR::ParameterSet subset(config().makeSubset("Cimager."));

                boost::scoped_ptr<askap::ProfileSingleton::Initialiser> profiler;
                if (parameterExists("profile")) {
                    std::string profileFileName("profile.imager");
                    if (subset.isDefined("Images.Names")){
                        profileFileName += "."+subset.getStringVector("Images.Names")[0];
                    }
                    if (comms_p.isParallel()) {
                        profileFileName += ".rank"+utility::toString(comms_p.rank());
                    }
                    profiler.reset(new askap::ProfileSingleton::Initialiser(profileFileName));
                }


                ASKAPCHECK(comms_p.isParallel(), "This imager can only be run as a parallel MPI job");
                // imager-specific configuration of the master/worker to allow groups of workers
                const int nWorkerGroups = subset.getInt32("nworkergroups", 1);
                ASKAPCHECK(nWorkerGroups > 0, "nworkergroups is supposed to be greater than 0");
                if (nWorkerGroups > 1) {
                    ASKAPLOG_INFO_STR(logger, "Model parameters will be distributed between "<<nWorkerGroups<<
                                      " groups of workers");
                    ASKAPCHECK(comms_p.isParallel(), "This option is only allowed in the parallel mode");
                    comms_p.defineGroups(nWorkerGroups);
                } else {
                    ASKAPLOG_INFO_STR(logger, "All workers are treated as identical");
                }
                // Instantiate the Distributed Imager
                ContinuumImager imager(subset, comms_p);

                // runit
                imager.run();
                comms_p.barrier(); // This is failing with craypat instrumentation on.
                stats.logSummary();

            } catch (const askap::AskapError& e) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << e.what());
                std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
                comms_p.abort();
                return 1;
            } catch (const std::exception& e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << e.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what()
                    << std::endl;
                comms_p.abort();
                return 1;
            }



            return 0;
        }
};

int main(int argc, char *argv[])
{
    ImagerApp app;
    app.addParameter("profile", "p", "Write profiling output files", false);
    return app.main(argc, argv);

}
