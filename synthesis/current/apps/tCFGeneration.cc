/// @file
///
/// @breif performance tests for CF generation
///
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

// Package level header file
#include "askap_synthesis.h"

// ASKAPsoft includes
#include <fitting/Params.h>
#include <askap/StatReporter.h>
#include <casa/Logging/LogIO.h>
#include <askap/Log4cxxLogSink.h>
#include <CommandLineParser.h>
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <askap/AskapUtil.h>
#include <dataaccess/ParsetInterface.h>
#include <profile/AskapProfiler.h>

#include "askap/AskapLogging.h"
#include "askap/AskapError.h"

ASKAP_LOGGER(logger, ".tCFGeneration");

// this should be after logging
#include <gridding/TestCFGenPerformance.h>


using namespace askap;
using namespace askap::synthesis;
#ifdef _OPENMP
#include <omp.h>
#endif

// Main function
int main(int argc, const char** argv)
{
    // This class must have scope outside the main try/catch block
    askap::askapparallel::AskapParallel comms(argc, argv);

    try {
        // Ensure that CASA log messages are captured
        casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
        casa::LogSink::globalSink(globalSink);

        StatReporter stats;
        std::string profileFileName("profile.tCFGeneration");
        if (comms.isParallel()) {
            profileFileName += ".rank"+utility::toString(comms.rank());
        }
        ASKAP_INIT_PROFILING(profileFileName);

        // Put everything in scope to ensure that all destructors are called
        // before the final message
        {
                    cmdlineparser::Parser parser; // a command line parser
            // command line parameter
            cmdlineparser::FlaggedParameter<std::string> inputsPar("-inputs",
                    "tgridding.in");
            // this parameter is optional
            parser.add(inputsPar, cmdlineparser::Parser::return_default);

            parser.process(argc, argv);

            const std::string parsetFile = inputsPar;

            LOFAR::ParameterSet parset(parsetFile);
            LOFAR::ParameterSet subset(parset.isDefined("Cimager.gridder") ? parset.makeSubset("Cimager.") : parset);
            
            ASKAPLOG_INFO_STR(logger, "Setting up the gridder to test and the model");
            boost::shared_ptr<TestCFGenPerformance> tester = TestCFGenPerformance::createGridder(subset.makeSubset("gridder.AWProject."));
            ASKAPCHECK(tester, "Gridder is not defined");
            #ifdef _OPENMP
            const int nthreads = omp_get_max_threads();
            std::vector<boost::shared_ptr<TestCFGenPerformance> > testers(nthreads);
            ASKAPLOG_INFO_STR(logger, "Will attempt to run "<<nthreads<<" instances in parallel");
            for (int i=0; i<nthreads; ++i) {
                 testers[i] = tester->clone(); 
            }
            #endif
            { 
               scimath::Params model;
               boost::shared_ptr<scimath::Params> modelPtr(&model, utility::NullDeleter());
               SynthesisParamsHelper::setUpImages(modelPtr,subset.makeSubset("Images."));
               ASKAPLOG_INFO_STR(logger, "Model contains the following elements: "<<model);
               std::vector<std::string> names = model.freeNames();
               ASKAPCHECK(names.size() > 0, "At least one image should be defined!");
               ASKAPLOG_INFO_STR(logger, "Using "<<names[0]<<" to setup gridder");
               #ifdef _OPENMP
               for (int i=0; i<nthreads; ++i) {
                    testers[i]->initialiseGrid(model.axes(names[0]),model.value(names[0]).shape(),false);
               }
               #else
               tester->initialiseGrid(model.axes(names[0]),model.value(names[0]).shape(),false);            
               #endif
               // all additional memory used to hold the model is now released 
            }
            #ifdef _OPENMP
            #pragma omp parallel default(shared)
            { 
               #pragma omp for 
               for (int i=0; i<nthreads; ++i) {
                 testers[i]->run();
               }
            }
            #else
            tester->run();
            #endif
        }
        stats.logSummary();
        ///==============================================================================
    } catch (const cmdlineparser::XParser &ex) {
        ASKAPLOG_FATAL_STR(logger, "Command line parser error, wrong arguments " << argv[0]);
        std::cerr << "Usage: " << argv[0] << " [-inputs parsetFile]"
                      << std::endl;
    } catch (const askap::AskapError& x) {
        ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
        std::cerr << "Askap error in " << argv[0] << ": " << x.what()
                      << std::endl;
        exit(1);
    } catch (const std::exception& x) {
        ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
        std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what()
                      << std::endl;
        exit(1);
    }

    return 0;
}

