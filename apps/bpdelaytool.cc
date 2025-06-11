/// @file
/// @brief an utility to manipulate delays in bandpass solutions
/// @copyright (c) 2025 CSIRO
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

// package-level include
#include <askap/askap_synthesis.h>

// own includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/Application.h>
#include <askap/utils/BandpassDelayHelper.h>
#include <askap/calibaccess/CalibAccessFactory.h>
#include <askap/calibaccess/ICalSolutionConstSource.h>
#include <askap/calibaccess/ICalSolutionSource.h>

// casa includes
#include <casacore/casa/OS/Timer.h>

// LOFAR include
#include <Common/ParameterSet.h>

// std includes
#include <string>

// boost includes
#include <boost/shared_ptr.hpp>

ASKAP_LOGGER(logger, ".bpdelaytool");

using namespace askap;

class BPDelayToolApp : public askap::Application {
public:
   /// @brief run application
   /// @param[in] argc number of parameters
   /// @param[in] argv parameter vector
   /// @return exit code
   int run(int argc, char *argv[]) final;
private:
   /// @brief create bandpass delay helper class
   /// @details This method creates a new instance of the bandpass delay helper class ensuring
   /// the uniform treatment of parset parameters.
   /// @param[in] parset parameter set (with the application prefix, if any, removed)
   /// @param[in] prefix prefix added to keywords when they're looked for in the specified parset
   /// @return shared pointer to the new bandpass delay helper class
   /// @note The prefix is supplied through the parameters because we may setup a number of helper classes and
   /// accessors. But it is handy to default some parameters like the largest permissable number of antennas to
   /// that in the output section (as it needs to be specified there anyway). Therefore, the specified prefix is
   /// overriden by "output." to check if such keywords are present (but the application level prefix, if any, 
   /// is expected to be removed in all cases)
   static boost::shared_ptr<utils::BandpassDelayHelper> createBandpassDelayHelper(const LOFAR::ParameterSet &parset, const std::string &prefix);

   /// @brief initialise helper class with bandpass data, calculate delays
   /// @details This method initialises an instance of the bandpass delay helper class
   /// by setting up a calibration accessor and reading the bandpass from it. It does the delay calculation afterwards, so we can encapsulate
   /// getting the spectral resolution from the parset. This method is intended to encapsulate operations
   /// required for a given section in the parset (e.g. with "input, "add", "subtract" prefixes)
   /// @param[in] parset parameter set (with all possible prefixes removed)
   /// @param[in] bdh shared pointer to bandpass delay helper which does all the operations
   static void fill(const LOFAR::ParameterSet &parset, const boost::shared_ptr<utils::BandpassDelayHelper> &bdh);

   std::string getVersion() const final {
      const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
      return pkgVersion;
   }
};

/// @brief create bandpass delay helper class
/// @details This method creates a new instance of the bandpass delay helper class ensuring
/// the uniform treatment of parset parameters.
/// @param[in] parset parameter set (with the application prefix, if any, removed)
/// @param[in] prefix prefix added to keywords when they're looked for in the specified parset
/// @return shared pointer to the new bandpass delay helper class
/// @note The prefix is supplied through the parameters because we may setup a number of helper classes and
/// accessors. But it is handy to default some parameters like the largest permissable number of antennas to
/// that in the output section (as it needs to be specified there anyway). Therefore, the specified prefix is
/// overriden by "output." to check if such keywords are present (but the application level prefix, if any, 
/// is expected to be removed in all cases)
boost::shared_ptr<utils::BandpassDelayHelper> BPDelayToolApp::createBandpassDelayHelper(const LOFAR::ParameterSet &parset, const std::string &prefix)
{
   // take the number of antennas, etc from the output section, then own table section and finally just the prefix. If none present, take pre-defined default value.
   const casacore::uInt maxAnt = parset.getUint32(prefix+"maxant", parset.getUint32(prefix+"calibaccess.table.maxant", parset.getUint32("output.calibaccess.table.maxant", 36u)));
   const casacore::uInt maxBeam = parset.getUint32(prefix+"maxbeam", parset.getUint32(prefix+"calibaccess.table.maxbeam", parset.getUint32("output.calibaccess.table.maxbeam", 36u)));
   const casacore::uInt maxChan = parset.getUint32(prefix+"maxchan", parset.getUint32(prefix+"calibaccess.table.maxchan", parset.getUint32("output.calibaccess.table.maxchan", 15552u)));
   boost::shared_ptr<utils::BandpassDelayHelper> result(new utils::BandpassDelayHelper(maxAnt, maxBeam, maxChan));
   return result;
}

/// @brief initialise helper class with bandpass data, calculate delays
/// @details This method initialises an instance of the bandpass delay helper class
/// by setting up a calibration accessor and reading the bandpass from it. It does the delay calculation afterwards, so we can encapsulate
/// getting the spectral resolution from the parset. This method is intended to encapsulate operations
/// required for a given section in the parset (e.g. with "input, "add", "subtract" prefixes)
/// @param[in] parset parameter set (with all possible prefixes removed, including "input", "add", etc)
/// @param[in] bdh shared pointer to bandpass delay helper which does all the operations
void BPDelayToolApp::fill(const LOFAR::ParameterSet &parset, const boost::shared_ptr<utils::BandpassDelayHelper> &bdh)
{
   const boost::shared_ptr<accessors::ICalSolutionConstSource> solSrc = accessors::CalibAccessFactory::roCalSolutionSource(parset);
   ASKAPCHECK(solSrc, "Unable to get calibration solution source!");
   const long solId = solSrc->mostRecentSolution();
   boost::shared_ptr<accessors::ICalSolutionConstAccessor> acc = solSrc->roSolution(solId);
   const double resolution = asQuantity(parset.getString("resolution", "1MHz")).getValue("Hz");
   ASKAPASSERT(acc);
   ASKAPASSERT(bdh);
   bdh->loadBandpass(*acc, resolution);
   bdh->calcDelays();
}

int BPDelayToolApp::run(int, char **) {
   try {
      casa::Timer timer;

      boost::shared_ptr<utils::BandpassDelayHelper> input = createBandpassDelayHelper(config().makeSubset("bpdelaytool."),"input.");
      ASKAPASSERT(input);
      const std::string inputType = config().getString("bpdelaytool.input", "calibaccess");
      ASKAPCHECK(inputType == "calibaccess" || inputType == "ideal", "Invalid bpdelaytool.input keyword, only 'calibaccess' and 'ideal' are supported");
      if (inputType == "ideal") {
          ASKAPLOG_INFO_STR(logger, "Initialising input bandpass with zero phase and unity amplitude for all channels");
          // ideally, it would be nice to encapsulate reading the parameter with its defaults; it is also used in the fill method
          const double resolution = asQuantity(config().getString("bpdelaytool.input.resolution", "1MHz")).getValue("Hz");
          input->setIdealBandpass(resolution);
          // a bit of the waste doing the calculation, we could've set the delays to zero and raise the validity flag (may add such a method in the future)
          // the zeroDelays() method deliberately leaves the delay flags intact. This intidy way seems to be ok as the "ideal" option is really an add-on
          // just in case and is not expected to be used (much or at all)
          input->calcDelays();
      } else {
          ASKAPLOG_INFO_STR(logger, "Reading the input bandpass via the calibration solution accessor");
          fill(config().makeSubset("bpdelaytool.input."), input);
          ASKAPLOG_INFO_STR(logger, "Delays in the input bandpass table:");
          input->summary();
      }

      if (config().isDefined("bpdelaytool.subtract")) {
          ASKAPCHECK(config().getString("bpdelaytool.subtract") == "input", "'bpdelaytool.subtract' keyword should be either 'input' or left undefined");
          ASKAPCHECK(!config().isDefined("bpdelaytool.subtract.calibaccess"), "'bpdelaytool.subtract.calibaccess' is incompatible with 'bpdelaytool.subtract=input'");
          ASKAPLOG_INFO_STR(logger, "Delays found in the input bandpass table will be compensated");
          // this is a short cut case where we remove the delays found in the input bandpass table (to avoid reading it again)
          input->negateDelays();
      } else {
          // need to zero delays in the input helper class here preserving the validity flags!
          input->zeroDelays();
      }

      const std::vector<std::string> operations = {"add", "subtract"};
      for (std::string operation : operations) {
           if (config().isDefined("bpdelaytool."+operation+".calibaccess")) {
               ASKAPLOG_INFO_STR(logger, "Reading the bandpass via the calibration solution accessor to "+operation+" its delays to the input");
               boost::shared_ptr<utils::BandpassDelayHelper> bdh = createBandpassDelayHelper(config().makeSubset("bpdelaytool."),operation + ".");
               ASKAPASSERT(bdh);
               fill(config().makeSubset("bpdelaytool."+operation+"."), bdh);
               ASKAPLOG_INFO_STR(logger, "Delays in the '"+operation+" delay' bandpass table:");
               bdh->summary();
               if (operation == "subtract") {
                   bdh->negateDelays();
               }
               input->addDelays(*bdh);
           }
      }

      if (config().isDefined("bpdelaytool.output.calibaccess")) {
          ASKAPLOG_INFO_STR(logger, "Storing the bandpass into the output calibration solution accessor, applying the following delays:");
          input->summary();
          input->applyDelays();
          // writing part
          const boost::shared_ptr<accessors::ICalSolutionSource> solSrc = accessors::CalibAccessFactory::rwCalSolutionSource(config().makeSubset("bpdelaytool.output."));
          ASKAPCHECK(solSrc, "Unable to get output calibration solution source!");
          const long solId = solSrc->newSolutionID(0);
          boost::shared_ptr<accessors::ICalSolutionAccessor> acc = solSrc->rwSolution(solId);
          input->storeBandpass(*acc);
      }

      ASKAPLOG_DEBUG_STR(logger, "Job: "<<timer.real());
   }
   catch(const AskapError &ce) {
      ASKAPLOG_FATAL_STR(logger, "AskapError has been caught. "<<ce.what());
      return -1;
   }
   catch(const std::exception &ex) {
      ASKAPLOG_FATAL_STR(logger, "std::exception has been caught. "<<ex.what());
      return -1;
   }
   catch(...) {
      ASKAPLOG_FATAL_STR(logger, "An unexpected exception has been caught");
      return -1;
   }
   return 0;
}

int main(int argc, char *argv[]) {
  BPDelayToolApp app;
  return app.main(argc,argv);
}
