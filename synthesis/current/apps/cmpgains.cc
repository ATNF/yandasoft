/// @file
/// Experiments to find a good way to compare two gain solutions
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

// casa includes
#include <casa/BasicSL/Complex.h>
#include <casa/Logging/LogIO.h>

// LOFAR
#include <Common/ParameterSet.h>

// own includes
#include <askap_synthesis.h>
#include <askap/AskapUtil.h>
#include <askap/AskapError.h>
#include <askap/AskapLogging.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <fitting/Params.h>
#include <askap/Log4cxxLogSink.h>
// just for logging
#include <askapparallel/AskapParallel.h>

// command line parser
#include <CommandLineParser.h>

// std includes
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

ASKAP_LOGGER(logger, ".cmpgains");

template<typename Cont>
bool exists(const Cont &cnt, const typename Cont::value_type &val) 
{
  const typename Cont::const_iterator ci = std::find(cnt.begin(),cnt.end(),val);
  return  ci != cnt.end();
}

int main(int argc, const char **argv)
{
   using namespace askap;
   using namespace askap::synthesis;

   // This class must have scope outside the main try/catch block
   // we need it to initialise the logging properly
   askapparallel::AskapParallel comms(argc, argv);

   try {
      // Ensure that CASA log messages are captured
      casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
      casa::LogSink::globalSink(globalSink);
        
      cmdlineparser::Parser parser; // a command line parser
      // command line parameters
      cmdlineparser::GenericParameter<std::string> gainsFileName1;
      cmdlineparser::GenericParameter<std::string> gainsFileName2;
      
      // required parameters
      parser.add(gainsFileName1);
      parser.add(gainsFileName2);
      
      parser.process(argc,argv);
      
      ASKAPLOG_INFO_STR(logger,"Loading gains from file "<<gainsFileName1.getValue());
      scimath::Params gains1;
      gains1 << LOFAR::ParameterSet(gainsFileName1.getValue());
      ASKAPLOG_INFO_STR(logger,"Loading gains from file "<<gainsFileName1.getValue());
      scimath::Params gains2;
      gains2 << LOFAR::ParameterSet(gainsFileName2.getValue());
      
      // build a union of all parameters
      const std::vector<std::string> names1 = gains1.names();
      const std::vector<std::string> names2 = gains2.names();
      std::set<std::string> allnames;
      std::set_union(names1.begin(),names1.end(),names2.begin(),names2.end(),std::inserter(allnames,allnames.begin()));
      std::ofstream os("out.dat");
      for (std::set<std::string>::const_iterator ci = allnames.begin(); ci!=allnames.end(); ++ci) {
           if (ci->find("g11") == std::string::npos) {
               continue;
           }
           if (!exists(names1,*ci) || !exists(names2,*ci)) {
               ASKAPLOG_INFO_STR(logger,"Gain parameter "<<*ci<<" is not present in both parameter sets");
               continue;
           }
           const casa::Complex g1 = gains1.complexValue(*ci);
           const casa::Complex g2 = gains2.complexValue(*ci);           
           //os<<*ci<<" "<<arg(g1)*180./casa::C::pi<<" "<<arg(g2)*180./casa::C::pi<<std::endl;
           os<<*ci<<" "<<abs(g1*conj(g2))<<" "<<arg(g1*conj(g2))<<std::endl;
      }
   }
   catch (const cmdlineparser::XParser &ex) {
      std::cerr<<"Usage: "<<argv[0]<<" gains1.par gains2.par"<<std::endl;
      std::cerr<<"gains1.par and gains2.par two parset files with gains"<<std::endl;
   }
   catch (const std::exception &ex) {
      std::cerr<<ex.what()<<std::endl;
   }       
}

