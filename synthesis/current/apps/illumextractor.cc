/// @file
/// 
/// @brief Utility to extract illumination pattern 
/// @details This program stores a given illumination pattern in an image
/// (any type supported by synthesis code). It allows to test the illumination 
/// pattern code standalone, but uses the same factory as the rest of the synthesis code.
/// 
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

// casa includes
#include <casa/Logging/LogIO.h>
#include <casa/OS/Timer.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>


// own includes
#include <askap/AskapLogging.h>
#include <askap/AskapUtil.h>
#include <askap/AskapError.h>
#include <Common/ParameterSet.h>
#include <askap/Log4cxxLogSink.h>
#include <gridding/UVPattern.h>
#include <utils/ImageUtils.h>
// just for logging
#include <askapparallel/AskapParallel.h>

// command line parser
#include <CommandLineParser.h>

// std includes
#include <stdexcept>

ASKAP_LOGGER(logger, ".illumextractor");

// this should be after logging
#include <gridding/AProjectGridderBase.h>

using namespace askap;
using namespace askap::synthesis;

// Main function
int main(int argc, const char** argv) { 

  // This class must have scope outside the main try/catch block
  // we need it to initialise the logging properly
  askapparallel::AskapParallel comms(argc, argv);
  try {
     casa::Timer timer;
     timer.mark();

     // Ensure that CASA log messages are captured
     casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
     casa::LogSink::globalSink(globalSink);
     {
        cmdlineparser::Parser parser; // a command line parser
        // command line parameter
        cmdlineparser::FlaggedParameter<std::string> inputsPar("-inputs",
                       "illumextractor.in");
        // this parameter is optional
        parser.add(inputsPar, cmdlineparser::Parser::return_default);

        parser.process(argc, argv);

        const std::string parsetFile = inputsPar;

        LOFAR::ParameterSet parset(parsetFile);
        if (parset.isDefined("Cimager.gridder")) {
            const std::string gridder = parset.getString("Cimager.gridder");
            ASKAPLOG_INFO_STR(logger,  "Using subset of the input parset file "<<parsetFile<<
                         " sliced at Cimager.gridder."<<gridder<<" to define illumination");
            SynthesisParamsHelper::setUpImageHandler(parset.makeSubset("Cimager."));             
            parset = parset.makeSubset("Cimager.gridder."+gridder+".");
        } else if (parset.isDefined("Csimulator.gridder")) {
            const std::string gridder = parset.getString("Csimulator.gridder");
            ASKAPLOG_INFO_STR(logger,  "Using subset of the input parset file "<<parsetFile<<
                         " sliced at Csimulator.gridder."<<gridder<<" to define illumination");
            SynthesisParamsHelper::setUpImageHandler(parset.makeSubset("Csimulator."));             
            parset = parset.makeSubset("Csimulator.gridder."+gridder+".");
        } else {
            SynthesisParamsHelper::setUpImageHandler(parset);
        }
            
        boost::shared_ptr<IBasicIllumination> illum = AProjectGridderBase::makeIllumination(parset);
        
        // hardcoded parameters at the moment
        UVPattern pattern(1024,1024, 10, 10, 4);
        
        ASKAPCHECK(illum, "No illumination pattern seems to be defined");
        illum->getPattern(1.4e9, pattern);
        
        casa::Array<casa::Float> buffer(pattern.pattern().shape());
        casa::convertArray<float,double>(buffer,casa::amplitude(pattern.pattern()));
        scimath::saveAsCasaImage("illum.img",buffer);
      }  
      ASKAPLOG_INFO_STR(logger,  "Total times - user:   " << timer.user()
               << " system: " << timer.system() << " real:   " << timer.real());
     
 ///==============================================================================
  } catch (const cmdlineparser::XParser &ex) {
        ASKAPLOG_FATAL_STR(logger, "Command line parser error, wrong arguments " << argv[0]);
        std::cerr << "Usage: " << argv[0] << " [-inputs parsetFile]" << std::endl;
        exit(1);
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
   
