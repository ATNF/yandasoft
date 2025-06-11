/// @file
///
/// @breif performance and debugging tests for MultiScaleBasisFunction class
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
#include "askap/askap_synthesis.h"

// ASKAPsoft includes
#include <askap/askap/StatReporter.h>
#include <casacore/casa/Logging/LogIO.h>
#include <askap/askap/Log4cxxLogSink.h>

#include <askap/utils/CommandLineParser.h>

#include <askap/askap/Application.h>
#include <askap/askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <askap/askap/AskapUtil.h>
#include <askap/profile/AskapProfiler.h>
#include <askap/scimath/utils/SpheroidalFunction.h>

#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>

ASKAP_LOGGER(logger, ".tMultiScaleBasisFunction");

// this should be after logging
#include <askap/deconvolution/MultiScaleBasisFunction.h>

#include <fstream>


using namespace askap;
using namespace askap::synthesis;

class MultiScaleBasisFunctionTestApp : public askap::Application {
public:
   /// @brief do the actual test
   void doTest();

   /// @brief test spheroidal function calculator
   void doSphFuncTest();

   /// @brief run application
   /// @param[in] argc number of parameters
   /// @param[in] argv parameter vector
   /// @return exit code
   int run(int argc, char *argv[]) final;

protected:
   // test class to access protected method without making the whole app class a friend
   class SphFuncTestClass : public MultiScaleBasisFunction<float> {
   public:
       SphFuncTestClass() : MultiScaleBasisFunction<float>(casacore::Vector<float>()) {};

       using MultiScaleBasisFunction<float>::spheroidal;
       using MultiScaleBasisFunction<float>::spheroidalOld;
   };

private:
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

void MultiScaleBasisFunctionTestApp::doSphFuncTest() {
   scimath::SpheroidalFunction sphFunc(casacore::C::pi*3.,1);
   SphFuncTestClass tc;
   std::ofstream os("out.dat");
   const size_t nPoints = 100;
   for (size_t i = 0; i<nPoints; ++i) {
        const float nu = 1./nPoints * i;
        os<<i<<" "<<nu<<" "<<SphFuncTestClass::spheroidalOld(nu)<<" "<<sphFunc(nu)<<" "<<tc.spheroidal(nu)<<std::endl;
   }
}

void MultiScaleBasisFunctionTestApp::doTest() {
  const casa::uInt size = config().getUint32("size");
  const casacore::IPosition shape(2, size, size);
  MultiScaleBasisFunction<float> msbf(shape, casacore::Vector<float>(config().getFloatVector("scales")), config().getBool("orthogonal", false));
}

int MultiScaleBasisFunctionTestApp::run(int argc, char **argv) {
  // This class must have scope outside the main try/catch block
  askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

  try {
        StatReporter stats;
        std::string profileFileName("profile.tMultiScaleBasisFunction");
        if (comms.isParallel()) {
            profileFileName += ".rank"+utility::toString(comms.rank());
        }
        ASKAP_INIT_PROFILING(profileFileName);

        casa::Timer timer;
        timer.mark();
        // Put everything in scope to ensure that all destructors are called
        // before the final message
        {
           doTest();
           doSphFuncTest();
        }
        std::cerr<<"Job: "<<timer.real()<<std::endl;
        stats.logSummary();
        ///==============================================================================
  }
  catch(const AskapError &ce) {
     std::cerr<<"AskapError has been caught. "<<ce.what()<<std::endl;
     return -1;
  }
  catch(const std::exception &ex) {
     std::cerr<<"std::exception has been caught. "<<ex.what()<<std::endl;
     return -1;
  }
  return 0;
}

// Main function
int main(int argc, const char* argv[])
{
  MultiScaleBasisFunctionTestApp app;
  return app.main(argc,const_cast<char**>(argv));
}
