/// @file
/// @brief an utility to solve for antenna-based delays
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


#include <dataaccess/TableDataSource.h>
#include <askap/askap_synthesis.h>
#include <askap/AskapLogging.h>
#include <askap/AskapUtil.h>
ASKAP_LOGGER(logger, ".delaysolver");

#include <askap/AskapError.h>
#include <askap/Application.h>
#include <dataaccess/SharedIter.h>
#include <dataaccess/TableManager.h>
#include <dataaccess/IDataConverterImpl.h>
#include <utils/PolConverter.h>


#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/OS/Directory.h>

#include <askap/utils/DelaySolverImpl.h>

#include <Common/ParameterSet.h>

#include <iomanip>
#include <vector>
#include <string>

using namespace askap;
using namespace askap::accessors;

class DelaySolverApp : public askap::Application {
public:
   /// @brief process a single file
   /// @param[in] ds data source
   /// @param[in] currentDelays a vector with fixed delays (per antenna) used for observations
   void process(const IConstDataSource &ds, const std::vector<double>& currentDelays);
   
   /// @brief run application
   /// @param[in] argc number of parameters
   /// @param[in] argv parameter vector
   /// @return exit code
   virtual int run(int argc, char *argv[]);
protected:
   /// @brief helper method to fill antenna names
   /// @details This method extracts antenna names from the supplied parset of ingest pipeline
   /// (obtained via the SB number) and populates itsAntennaNames in the order of ingest indices.
   /// @param[in] parset parset used with ingest pipeline
   void fillAntennaNames(const LOFAR::ParameterSet &parset);

   /// @brief helper method to extract current delays from the ingest's parset
   /// @details
   /// @param[in] parset parset used with ingest pipeline
   /// @return vector with delays in ns
   std::vector<double> getCurrentDelays(const LOFAR::ParameterSet &parset);

private:
   /// @brief vector with antenna names
   /// @details Note, indices correspond to that internally used by ingest. If empty, the output
   /// parset (in the fcm format) will be in the form of ingest-specific vector (used in some test
   /// applications, although to be deprecated for normal operations). 
   /// @note If not empty, the size of the vector should match that of the delay vector
   std::vector<std::string> itsAntennaNames;
};

/// @brief helper method to fill antenna names
/// @details This method extracts antenna names from the supplied parset of ingest pipeline
/// (obtained via the SB number) and populates itsAntennaNames in the order of ingest indices.
void DelaySolverApp::fillAntennaNames(const LOFAR::ParameterSet &parset)
{
   if (!parset.isDefined("baselinemap.antennaidx") || !parset.isDefined("antennas")) {
       itsAntennaNames.resize(0);
       return;
   }

   // this is the name of antennas as specified in antenna.antXX.name keyword, we have to find
   // the name how it is called in the parset/fcm. In the ingest parset there is no cp.ingest prefix
   const std::vector<std::string> antennas = parset.getStringVector("baselinemap.antennaidx");

   // this is the list of all "parset/fcm" names of antennas, including those which may be unused at the moment
   // in the ingest parset there is no "common" prefix
   const std::vector<std::string> parsetNames = parset.getStringVector("antennas");

   std::map<std::string, std::string> ingestName2ParsetNameMap;
   for (size_t ant=0; ant<parsetNames.size(); ++ant) {
        const std::string curName = parset.getString("antenna." + parsetNames[ant] + ".name");
        ASKAPCHECK(ingestName2ParsetNameMap.find(curName) == ingestName2ParsetNameMap.end(), "Detected duplicated antenna name: "<<curName);
        ingestName2ParsetNameMap[curName] = parsetNames[ant];
   }

   
   itsAntennaNames.resize(antennas.size());
   for (size_t ant = 0; ant < itsAntennaNames.size(); ++ant) {
        const std::map<std::string, std::string>::const_iterator ci = ingestName2ParsetNameMap.find(antennas[ant]);
        ASKAPCHECK(ci != ingestName2ParsetNameMap.end(), "Antenna "<<antennas[ant]<<" is not defined, there is no antenna.XX.name keyword equal to "<<antennas[ant]);
        itsAntennaNames[ant] = ci->second;
   }
}

/// @brief helper method to extract current delays from the ingest's parset
/// @details
/// @param[in] parset parset used with ingest pipeline
/// @return vector with delays in ns
std::vector<double> DelaySolverApp::getCurrentDelays(const LOFAR::ParameterSet &parset)
{
   const std::string ingestSpecificFixedDelayKey = "tasks.FringeRotationTask.params.fixeddelays";
   if (parset.isDefined(ingestSpecificFixedDelayKey)) {
       ASKAPLOG_WARN_STR(logger, "Old-style fixed delay key ("<<ingestSpecificFixedDelayKey<<") is present in the ingest parset - ignoring");
       // uncomment the following line to use old-style fixed delay key instead of the new way to specify fixed delays
       // (may be handy for some commissioning experiments)
       //return parset.getDoubleVector("tasks.FringeRotationTask.params.fixeddelays");                                
   } 
   ASKAPCHECK(itsAntennaNames.size() > 0, "No antennas seem to be defined in the ingest parset, unable to load initial delays");
   // default delay - we use this value if antenna-specific keyword is not defined
   const std::string defaultDelay = parset.getString("antenna.ant.delay", "0s");
  
   std::vector<double> delays(itsAntennaNames.size(),0.);
   for (size_t ant=0; ant < delays.size(); ++ant) {
        const std::string delayKey = "antenna." + itsAntennaNames[ant] + ".delay";
        if (!parset.isDefined(delayKey)) {
            ASKAPLOG_WARN_STR(logger, "Antenna specific delay key ("<<delayKey<<") is not found in the ingest's parset, using the default: "<<defaultDelay);
        }
        const double delay = asQuantity(parset.getString(delayKey, defaultDelay)).getValue("ns");
        ASKAPLOG_INFO_STR(logger, "Initial delay for "<<itsAntennaNames[ant]<<" is "<<std::setprecision(9)<<delay<<" ns");
        delays[ant] = delay;
   }
   return delays;
}

void DelaySolverApp::process(const IConstDataSource &ds, const std::vector<double>& currentDelays) {
  IDataSelectorPtr sel=ds.createSelector();
  casa::uInt beam = config().getUint("beam",0);
  sel->chooseFeed(beam);
  sel->chooseCrossCorrelations();
  const int scan = config().getInt32("scan",-1);
  if (scan >= 0) {
      ASKAPLOG_INFO_STR(logger, "Process only scan "<<scan);
      sel->chooseUserDefinedIndex("SCAN_NUMBER",casa::uInt(scan));
  }
 
  IDataConverterPtr conv=ds.createConverter();  
  conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"Hz");
  conv->setEpochFrame(casa::MEpoch(casa::Quantity(55913.0,"d"),
                      casa::MEpoch::Ref(casa::MEpoch::UTC)),"s");
  conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));                    
 
  const double targetRes = config().getDouble("resolution",1e6);
  const std::string stokesStr = config().getString("stokes","XX");
  const casa::Vector<casa::Stokes::StokesTypes> stokesVector = scimath::PolConverter::fromString(stokesStr);
  ASKAPCHECK(stokesVector.nelements() == 1, "Exactly one stokes parameter should be defined, you have "<<stokesStr);
  const double ampCutoff = config().getDouble("cutoff",-1.);
  const casa::uInt refAnt = config().getUint("refant",1);
  const bool exclude13 = config().getBool("exclude13", false);
  utils::DelaySolverImpl solver(targetRes, stokesVector[0], ampCutoff, refAnt);
  if (exclude13) {
      solver.excludeBaselines(casa::Vector<std::pair<casa::uInt,casa::uInt> >(1,std::pair<casa::uInt,
             casa::uInt>(1,2)));
  }
  
  const bool estimateViaLags = config().getBool("uselags", false);
  
  if (estimateViaLags) {
      ASKAPLOG_INFO_STR(logger, "Initial delay will be estimated via lags using full resolution data");
      // the following means no averaging
      solver.setTargetResolution(1.);
      
      // initial pass over data
      for (IConstDataSharedIter it=ds.createConstIterator(sel,conv);it!=it.end();++it) {
           solver.process(*it);  
      }
      // solve via FFT
      const casa::Vector<double> delayApprox = solver.solve(true);
      // use FFT-based estimate as an approximation before averaging
      solver.setApproximateDelays(delayApprox);
      solver.init();
      solver.setTargetResolution(targetRes);      
  }
      
  for (IConstDataSharedIter it=ds.createConstIterator(sel,conv);it!=it.end();++it) {
       solver.process(*it);  
  }
  
  // corrections have the opposite sign from determined delays, hence the minus
  // the units in the fcm are in ns
  casa::Vector<double> delays = -solver.solve(false) * 1e9;
  ASKAPLOG_INFO_STR(logger, "Corrections (ns): "<<std::setprecision(9)<<delays);
  if (currentDelays.size() > 0) {
      ASKAPLOG_INFO_STR(logger, "Old delays (ns): "<< std::setprecision(9) << currentDelays);
      ASKAPCHECK(currentDelays.size() == delays.nelements(), "Number of antennas differ in fixeddelays (or antXX.delay) parameters and in the dataset");
      for (casa::uInt ant = 0; ant < delays.nelements(); ++ant) {
           delays[ant] += currentDelays[ant];
      }
      ASKAPLOG_INFO_STR(logger, "New delays (ns): "<< std::setprecision(9)<<delays);
      const std::string outParset = "corrected_fixeddelay.parset"; 
      {
          std::ofstream os(outParset.c_str());
          // write the file in the format directly understood by fcm put to simplify operations
          if ((itsAntennaNames.size() == 0) || config().getBool("oldfcmformat",false)) {
              ASKAPLOG_WARN_STR(logger, "Exporting delays in the old FCM format - use at your own risk");
              os << "cp.ingest.tasks.FringeRotationTask.params.fixeddelays = " << std::setprecision(9)<<delays << std::endl;
          } else {
              ASKAPCHECK(itsAntennaNames.size() == delays.size(), 
                    "Number of antennas defined in the ingest parset is different from the number of antennas delays are solved for");
              for (casa::uInt ant = 0; ant < delays.nelements(); ++ant) {
                   os << "common.antenna."<<itsAntennaNames[ant]<<".delay = " << std::setprecision(9)<<delays[ant] << "ns"<<std::endl;
                   ASKAPLOG_INFO_STR(logger, "      "<<std::setw(5)<<itsAntennaNames[ant]<< " (index "<<ant<<") ->  "<<std::setprecision(9)<<delays[ant]<<" ns");
              }
          }
      }
      ASKAPLOG_INFO_STR(logger, "The new delays are now stored in "<<outParset);       
  } else {
      ASKAPLOG_WARN_STR(logger, "No fixed delays specified in the parset -> no update");
  }
}

int DelaySolverApp::run(int, char **) {
  try {

     casa::Timer timer;
     std::string msName = parameterExists("ms") ? parameter("ms") : "";
     const std::string sbID = parameterExists("sb") ? parameter("sb") : "";

     // get current delays from the application's parset, this is only intended to be used if no scheduling block ID is given
     std::vector<double> currentDelays = config().getDoubleVector("cp.ingest.tasks.FringeRotationTask.params.fixeddelays", 
                                               std::vector<double>());
     if (config().isDefined("ms")) {
         ASKAPCHECK(msName == "", "Use either ms parset parameter or the command line argument, not both");
         msName = config().getString("ms");
     }
     
     if (sbID != "") {
         // scheduling block ID is specified, the file name will be taken from SB
         ASKAPCHECK(msName == "", "When the scheduling block ID is specified, the file name is taken from that SB. "
                    "Remove the -f command line parameter or ms keyword in the parset to continue.");
         casa::Path path2sb(config().getString("sbpath","./"));
         path2sb.append(sbID);
         const casa::Directory sbDir(path2sb);
         // do not follow symlinks, non-recursive
         const casa::Vector<casa::String> dirContent = sbDir.find(casa::Regex::fromPattern("*.ms"),casa::False, casa::False);
         ASKAPCHECK(dirContent.nelements() > 0, "Unable to find a measurement set file in "<<sbDir.path().absoluteName());
         int fileIndex = -1;
         if (dirContent.nelements() != 1) {
             ASKAPLOG_INFO_STR(logger, "Multiple MSs are found in "<<sbDir.path().absoluteName()<<" - assume one file per beam and index == beam");
             casa::uInt beam = config().getUint("beam",0);
             for (casa::uInt i = 0; i<dirContent.nelements(); ++i) {
                  const casa::String nameTemplate = "_"+utility::toString<casa::uInt>(beam)+".ms";
                  const size_t pos = dirContent[i].rfind(nameTemplate);
                  if ((pos != casa::String::npos) && (pos + nameTemplate.size() == dirContent[i].size())) {
                      ASKAPCHECK(fileIndex == -1, "Multiple measurement sets matching beam = "<<beam<<" are present in "<<sbDir.path().absoluteName());
                      ASKAPLOG_INFO_STR(logger, "Using "<<dirContent[i]);
                      fileIndex = static_cast<int>(i);
                  }
             }
             ASKAPCHECK(fileIndex >= 0, "Unable to find MS matching beam = "<<beam<<" selected in "<<sbDir.path().absoluteName());
         } else {
           fileIndex = 0;
         }
             
         ASKAPASSERT((fileIndex >= 0) && (fileIndex < static_cast<int>(dirContent.nelements())));
         casa::Path path2ms(path2sb);
         path2ms.append(dirContent[fileIndex]);
         msName = path2ms.absoluteName();
         // fixed delays will be taken from cpingest.in in the SB directory
         casa::Path path2cpingest(path2sb);
         path2cpingest.append("cpingest.in");
         ASKAPLOG_INFO_STR(logger, "Ingest parset: "<<path2cpingest.absoluteName());
         const LOFAR::ParameterSet ingestParset(path2cpingest.absoluteName());
         ASKAPCHECK(currentDelays.size() == 0, "When the scheduling block ID is specified, the current fixed delays are taken "
                    "from the ingest pipeline parset stored with that SB. Remove it from the application's parset to continue.");
         
         fillAntennaNames(ingestParset);
         // here we look at the actual ingest pipeline parset not the fcm, so there is no cp.ingest prefix
         currentDelays = getCurrentDelays(ingestParset);
     }
     timer.mark();
     ASKAPCHECK(msName != "", "Measurement set should be specified explicitly or the scheduling block should be given");
     ASKAPLOG_INFO_STR(logger, "Processing measurement set "<<msName);
     TableDataSource ds(msName,TableDataSource::MEMORY_BUFFERS);     
     std::cerr<<"Initialization: "<<timer.real()<<std::endl;
     timer.mark();
     process(ds, currentDelays);
     std::cerr<<"Job: "<<timer.real()<<std::endl;
  }
  catch(const AskapError &ce) {
     std::cerr<<"AskapError has been caught. "<<ce.what()<<std::endl;
     return -1;
  }
  catch(const std::exception &ex) {
     std::cerr<<"std::exception has been caught. "<<ex.what()<<std::endl;
     return -1;
  }
  catch(...) {
     std::cerr<<"An unexpected exception has been caught"<<std::endl;
     return -1;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  DelaySolverApp app;
  app.addParameter("ms","f", "Measurement set name (optional)","");
  app.addParameter("sb","s", "Scheduling block number (optional)","");  
  return app.main(argc,argv);
}

