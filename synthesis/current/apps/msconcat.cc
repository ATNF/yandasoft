/// @file msconcat.cc
///
/// @brief
///
/// @copyright (c) 2012 CSIRO
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

// Package level header file
#include "askap_synthesis.h"

// System includes
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

// ASKAPsoft includes
#include "askap/AskapError.h"
#include "askap/AskapLogging.h"
#include "askap/StatReporter.h"
#include "askap/Log4cxxLogSink.h"
#include "boost/shared_ptr.hpp"
#include "CommandLineParser.h"
#include "casacore/casa/OS/File.h"
#include "casacore/casa/aips.h"
#include "casacore/casa/Quanta.h"
#include "casacore/casa/Arrays/Vector.h"
#include "casacore/casa/Arrays/MatrixMath.h"
#include "casacore/tables/Tables/TableDesc.h"
#include "casacore/tables/Tables/SetupNewTab.h"
#include "casacore/tables/Tables/IncrementalStMan.h"
#include "casacore/tables/Tables/StandardStMan.h"
#include "casacore/tables/Tables/TiledShapeStMan.h"

#include "casacore/ms/MSOper/MSConcat.h"
#include "casacore/ms/MeasurementSets/MeasurementSet.h"
#include "casacore/ms/MeasurementSets/MSColumns.h"

ASKAP_LOGGER(logger, ".msconcat");

using namespace askap;
using namespace casa;

void concat(const std::vector<std::string>& inFiles, const std::string& outFile)
{
    ASKAPCHECK(!casa::File(outFile).exists(),
            "File or table " << outFile << " already exists!");

    ASKAPCHECK(!inFiles.empty(), "No input measurement sets!");

    ASKAPLOG_INFO_STR(logger, "Concatenating " << inFiles.size() <<
            " measurement sets");

    // copy the first file to the output file
    const std::string &first = inFiles[0];

    ASKAPLOG_INFO_STR(logger, "Adding ms " << first);
    ASKAPCHECK(Table::isReadable(first), "ms "+first+" is not readable!");

    const string command = "cp -r "+first+" "+outFile;
    const int cp_status = system(command.c_str());
    ASKAPCHECK(cp_status==0, "Error copying "+first+" to "+outFile);
 
    // load the output ms and initialise MSConcat
    MeasurementSet ms_out(outFile, Table::Update);
    MSConcat mscat(ms_out);

    // check and concatenate the measurement sets
    std::vector<std::string>::const_iterator it;
    for (it = inFiles.begin(); it != inFiles.end(); ++it) {

        if (it == inFiles.begin()) {
            // already copied in
            continue;
        }

        ASKAPLOG_INFO_STR(logger, "Adding ms " << *it);
        ASKAPCHECK(Table::isReadable(*it), "ms "+*it+" is not readable!");

        // load the ms and concatenate
        MeasurementSet ms(*it, Table::Old);
        mscat.concatenate(ms);

    }

}

// Main function
int main(int argc, const char** argv)
{
    // Now we have to initialize the logger before we use it
    // If a log configuration exists in the current directory then
    // use it, otherwise try to use the programs default one
    std::ifstream config("askap.log_cfg", std::ifstream::in);
    if (config) {
        ASKAPLOG_INIT("askap.log_cfg");
    } else {
        std::ostringstream ss;
        ss << argv[0] << ".log_cfg";
        ASKAPLOG_INIT(ss.str().c_str());
    }

    // Ensure that CASA log messages are captured
    casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
    casa::LogSink::globalSink(globalSink);

    try {
        StatReporter stats;

        cmdlineparser::Parser parser; // a command line parser
        // command line parameter
        cmdlineparser::FlaggedParameter<std::string> outName("-o", "output.ms");
        // this parameter is required
        parser.add(outName, cmdlineparser::Parser::throw_exception);
        if (argc < 4) {
            throw cmdlineparser::XParser();
        }
        std::vector<cmdlineparser::GenericParameter<std::string> > inNames(argc-3);
        {
            std::vector<cmdlineparser::GenericParameter<std::string> >::iterator it;
            for (it = inNames.begin(); it < inNames.end(); ++it) {
                parser.add(*it);
            }
        }

        // Process command line options
        parser.process(argc, argv);
        if (!inNames.size()) {
            throw cmdlineparser::XParser();
        } 
        ASKAPLOG_INFO_STR(logger,
                "This program concatenates given measurement sets into "
                << outName.getValue());

        // Turns inNames into vector<string>
        std::vector<std::string> inNamesVec;
        std::vector<cmdlineparser::GenericParameter<std::string> >::iterator it;
        for (it = inNames.begin(); it < inNames.end(); ++it) {
            inNamesVec.push_back(it->getValue());
        }

        concat(inNamesVec, outName.getValue()); 

        stats.logSummary();
        ///==============================================================================
    } catch (const cmdlineparser::XParser &ex) {
        ASKAPLOG_FATAL_STR(logger, "Command line parser error, wrong arguments " << argv[0]);
        ASKAPLOG_FATAL_STR(logger, "Usage: " << argv[0] << " -o output.ms inMS1 ... inMSn");
        return 1;
    } catch (const askap::AskapError& x) {
        ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
        return 1;
    } catch (const std::exception& x) {
        ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
        return 1;
    }

    return 0;
}
