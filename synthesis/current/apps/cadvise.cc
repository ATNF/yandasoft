/// @file
///
/// Application to adivse parameters for synthesis imaging and calibration 
/// Control parameters are passed in from a LOFAR ParameterSet file. 
/// At least the measurement set (given by the dataset keyword) should be
/// defined. There is also an option to choose a user defined tangent point
/// as well as to enable snap-shot imaging by providing w-plane fitting threshold
/// (wtolerance).
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>

// Package level header file
#include <askap_synthesis.h>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/Application.h>
#include <askap/StatReporter.h>
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>

// Local package includes
#include <parallel/AdviseParallel.h>

ASKAP_LOGGER(logger, ".cadvise");

using namespace std;
using namespace askap;
using namespace askap::synthesis;

class CAdviseApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            StatReporter stats;

            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                
                // We cannot issue log messages until MPI is initialized!
                AdviseParallel cadvise(comms, config());

                ASKAPLOG_INFO_STR(logger, "ASKAP synthesis advise application " << ASKAP_PACKAGE_VERSION);

                if (comms.isMaster()) {
                    ASKAPLOG_INFO_STR(logger, "Parset file contents:\n" << config());
                }

                cadvise.estimate();
                cadvise.summary();
            } catch (const askap::AskapError& e) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << e.what());
                std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
                exit(1);
            } catch (const std::exception& e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << e.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what() << std::endl;
                exit(1);
            }
            stats.logSummary();
            return 0;
        }
};

// Main function
int main(int argc, char* argv[])
{
    CAdviseApp app;
    return app.main(argc, argv);
}
