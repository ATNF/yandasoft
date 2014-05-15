///
/// @file : Test program for dynamic gridder load
///
/// @copyright (c) 2009 CSIRO
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
/// @author Ger van Diepen <diepen at astron dot nl>

#include <gridding/IVisGridder.h>
#include <gridding/VisGridderFactory.h>

#include <askapparallel/AskapParallel.h>
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/Log4cxxLogSink.h>
ASKAP_LOGGER(logger, ".testgridderload");

#include <stdexcept>
#include <iostream>


using namespace askap;
using namespace askap::synthesis;
using namespace LOFAR;

// Main function
int main (int argc, const char** argv)
{
  try {

    // Ensure that CASA log messages are captured
    casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
    casa::LogSink::globalSink (globalSink);
    // the following line is necessary just to initialise the logger
    askap::askapparallel::AskapParallel ap(argc, argv);
    ASKAPLOG_INFO_STR(logger, "Testing dynamic loading of the gridder");

    // Create a TestLoadGridder instance.
    ParameterSet parset;
    parset.add ("gridder", "TestLoadGridder");
    IVisGridder::ShPtr gridder = VisGridderFactory::make (parset);

  } catch (const std::exception& x) {
    std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what()
              << std::endl;
    ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": "
                       << x.what());
    return 1;
  }
  return 0;
}

