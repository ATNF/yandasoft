/// @file ContinuumWorker.cc
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>

// Include own header file first
#include "ContinuumWorker.h"

// System includes
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <sys/stat.h>
#include <unistd.h>

#include "boost/shared_ptr.hpp"
#include "boost/filesystem.hpp"
// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/scimath/fitting/Equation.h>
#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/ImagingNormalEquations.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/gridding/IVisGridder.h>
#include <askap/gridding/VisGridderFactory.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/dataaccess/IConstDataSource.h>
#include <askap/dataaccess/TableConstDataSource.h>
#include <askap/dataaccess/IConstDataIterator.h>
#include <askap/dataaccess/IDataConverter.h>
#include <askap/dataaccess/IDataSelector.h>
#include <askap/dataaccess/IDataIterator.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/scimath/utils/PolConverter.h>
#include <Common/ParameterSet.h>
#include <Common/Exceptions.h>
#include <casacore/casa/OS/Timer.h>
#include <askap/parallel/ImagerParallel.h>
#include <askap/imageaccess/BeamLogger.h>
#include <askap/imagemath/linmos/LinmosAccumulator.h>
// CASA Includes

// Local includes
#include "askap/distributedimager/AdviseDI.h"
#include "askap/distributedimager/CalcCore.h"
#include "askap/distributedimager/MSSplitter.h"
#include "askap/messages/ContinuumWorkUnit.h"
#include "askap/messages/ContinuumWorkRequest.h"
#include "askap/distributedimager/CubeBuilder.h"
#include "askap/distributedimager/CubeComms.h"

using namespace std;
using namespace askap::cp;
using namespace askap;
using namespace askap::scimath;
using namespace askap::synthesis;
using namespace askap::accessors;

ASKAP_LOGGER(logger, ".ContinuumWorker");

ContinuumWorker::ContinuumWorker(LOFAR::ParameterSet& parset,
  CubeComms& comms)
  : itsParset(parset), itsComms(comms), itsBeamList()
{
    itsAdvisor = boost::shared_ptr<synthesis::AdviseDI> (new synthesis::AdviseDI(itsComms, itsParset));
    itsAdvisor->prepare();

    // lets properly size the storage
    const int nchanpercore = itsParset.getInt32("nchanpercore", 1);
    workUnits.resize(0);
    itsParsets.resize(0);

    // lets calculate a base
    unsigned int nWorkers = itsComms.nProcs() - 1;
    unsigned int nWorkersPerGroup = nWorkers / itsComms.nGroups();

    unsigned int id = itsComms.rank();
    // e. g. rank 8, 3 per group should be pos. 1 (zero index)
    unsigned int posInGroup = (id % nWorkersPerGroup);

    if (posInGroup == 0) {
      posInGroup = nWorkersPerGroup;
    }
    posInGroup = posInGroup - 1;

    this->baseChannel = posInGroup * nchanpercore;

    ASKAPLOG_INFO_STR(logger, "Distribution: Id " << id << " nWorkers " << nWorkers << " nGroups " << itsComms.nGroups());

    ASKAPLOG_INFO_STR(logger, "Distribution: Base channel " << this->baseChannel << " PosInGrp " << posInGroup);

    itsDoingPreconditioning = false;
    const vector<string> preconditioners = itsParset.getStringVector("preconditioner.Names", std::vector<std::string>());
    for (vector<string>::const_iterator pc = preconditioners.begin(); pc != preconditioners.end(); ++pc) {
      if ((*pc) == "Wiener" || (*pc) == "NormWiener" || (*pc) == "Robust" || (*pc) == "GaussianTaper") {
        itsDoingPreconditioning = true;
      }
    }

    itsGridderCanMosaick = false;
    std::string GridderStr = itsParset.getString("gridder",std::string());
    ASKAPLOG_INFO_STR(logger, "Gridder is " << GridderStr);

    if (GridderStr == "AWProject") {
      itsGridderCanMosaick = true;
      ASKAPLOG_INFO_STR(logger," Gridder <CAN> mosaick");
    }
    else {
      ASKAPLOG_INFO_STR(logger,"Gridder <CANNOT> mosaick");
    }


}

ContinuumWorker::~ContinuumWorker()
{



}

void ContinuumWorker::run(void)
{

  // Send the initial request for work
  ContinuumWorkRequest wrequest;

  ASKAPLOG_DEBUG_STR(logger, "Worker is sending request for work");

  wrequest.sendRequest(itsMaster, itsComms);


  while (1) {

    ContinuumWorkUnit wu;


    ASKAPLOG_DEBUG_STR(logger, "Worker is waiting for work allocation");
    wu.receiveUnitFrom(itsMaster, itsComms);
    if (wu.get_payloadType() == ContinuumWorkUnit::DONE) {
      ASKAPLOG_INFO_STR(logger, "Worker has received complete allocation");
      break;
    } else if (wu.get_payloadType() == ContinuumWorkUnit::NA) {
      ASKAPLOG_WARN_STR(logger, "Worker has received non applicable allocation");
      ASKAPLOG_WARN_STR(logger, "In new scheme we still process it ...");

    } else {

      ASKAPLOG_INFO_STR(logger, "Worker has received valid allocation");
    }
    const string ms = wu.get_dataset();
    ASKAPLOG_INFO_STR(logger, "Received Work Unit for dataset " << ms
      << ", local (topo) channel " << wu.get_localChannel()
      << ", global (topo) channel " << wu.get_globalChannel()
      << ", frequency " << wu.get_channelFrequency() / 1.e6 << " MHz"
      << ", width " << wu.get_channelWidth() / 1e3 << " kHz");
    try {
        ASKAPLOG_DEBUG_STR(logger, "Parset Reports (before): " << (itsParset.getStringVector("dataset", true)));
        processWorkUnit(wu);
        ASKAPLOG_DEBUG_STR(logger, "Parset Reports (after): " << (itsParset.getStringVector("dataset", true)));
    } catch (AskapError& e) {
        ASKAPLOG_WARN_STR(logger, "Failure processing workUnit");
        ASKAPLOG_WARN_STR(logger, "Exception detail: " << e.what());
    }


    wrequest.sendRequest(itsMaster, itsComms);

  } // while (1) // break when "DONE"
  ASKAPCHECK(workUnits.size() > 0, "No work at to do - something has broken in the setup");

  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " received data from master - waiting at barrier");
  itsComms.barrier(itsComms.theWorkers());
  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " passed barrier");

  const bool localSolver = itsParset.getBool("solverpercore", false);

  int nchanpercore = itsParset.getInt("nchanpercore", 1);

  int nWorkers = itsComms.nProcs() - 1;
  int nGroups = itsComms.nGroups();
  int nchanTotal = nWorkers * nchanpercore / nGroups;

  initialiseBeamLog(nchanTotal);

  if (localSolver) {
    ASKAPLOG_INFO_STR(logger, "In local solver mode - reprocessing allocations)");
    itsAdvisor->updateComms();
    int myMinClient = itsComms.rank();
    int myMaxClient = itsComms.rank();

    if (itsComms.isWriter()) {
      ASKAPLOG_DEBUG_STR(logger, "Getting client list for cube generation");
      std::list<int> myClients = itsComms.getClients();
      myClients.push_back(itsComms.rank());
      myClients.sort();
      myClients.unique();

      ASKAPLOG_DEBUG_STR(logger, "Client list " << myClients);
      if (myClients.size() > 0) {
        std::list<int>::iterator iter = std::min_element(myClients.begin(), myClients.end());

        myMinClient = *iter;
        iter = std::max_element(myClients.begin(), myClients.end());
        myMaxClient = *iter;

      }
      // these are in ranks
      // If a client is missing entirely from the list - the cube will be missing
      // channels - but they will be correctly labelled

      // e.g
      // bottom client rank is 4 - top client is 7
      // we have 4 chanpercore
      // 6*4 - 3*4
      // 3*4 = 12
      // (6 - 3 + 1) * 4
      if (!itsComms.isSingleSink()) {
        this->nchanCube = (myMaxClient - myMinClient + 1) * nchanpercore;
        this->baseCubeGlobalChannel = (myMinClient - 1) * nchanpercore;
        this->baseCubeFrequency = itsAdvisor->getBaseFrequencyAllocation((myMinClient - 1));
        ASKAPLOG_INFO_STR(logger, "MultiCube with multiple writers");
      } else {
        ASKAPLOG_INFO_STR(logger, "SingleCube with multiple writers");
        this->nchanCube = nchanTotal;
        this->baseCubeGlobalChannel = 0;
        this->baseCubeFrequency = itsAdvisor->getBaseFrequencyAllocation((0));

      }
      ASKAPLOG_INFO_STR(logger, "Number of channels in cube is: " << this->nchanCube);
      ASKAPLOG_INFO_STR(logger, "Base global channel of cube is " << this->baseCubeGlobalChannel);
    }
    this->baseFrequency = itsAdvisor->getBaseFrequencyAllocation(itsComms.rank() - 1);
  }
  ASKAPLOG_INFO_STR(logger, "Adding missing parameters");

  itsAdvisor->addMissingParameters();



  try {
    processChannels();
  } catch (AskapError& e) {
    ASKAPLOG_WARN_STR(logger, "Failure processing the channel allocation");
    ASKAPLOG_WARN_STR(logger, "Exception detail: " << e.what());
    throw;

  }

  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " finished");

  itsComms.barrier(itsComms.theWorkers());
  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " passed final barrier");
}
void ContinuumWorker::deleteWorkUnitFromCache(ContinuumWorkUnit& wu, LOFAR::ParameterSet& unitParset)
{

  const string ms = wu.get_dataset();

  struct stat buffer;

  if (stat(ms.c_str(), &buffer) == 0) {
    ASKAPLOG_WARN_STR(logger, "Split file " << ms << " exists - deleting");
    boost::filesystem::remove_all(ms.c_str());
  } else {
    ASKAPLOG_WARN_STR(logger, "Split file " << ms << " does not exist - nothing to do");
  }

}
void ContinuumWorker::clearWorkUnitCache()
{

  for (int cf = 0; cf < cached_files.size(); cf++) {
    const string ms = cached_files[cf];
    struct stat buffer;

    if (stat(ms.c_str(), &buffer) == 0) {
      ASKAPLOG_WARN_STR(logger, "Split file " << ms << " exists - deleting");
      boost::filesystem::remove_all(ms.c_str());
    } else {
      ASKAPLOG_WARN_STR(logger, "Split file " << ms << " does not exist - nothing to do");
    }

  }
}

void ContinuumWorker::cacheWorkUnit(ContinuumWorkUnit& wu, LOFAR::ParameterSet& unitParset)
{



  boost::filesystem::path mspath = boost::filesystem::path(wu.get_dataset());
  const string ms = mspath.filename().string();

  const string shm_root = unitParset.getString("tmpfs", "/dev/shm");

  std::ostringstream pstr;

  pstr << shm_root << "/" << ms << "_chan_" << wu.get_localChannel() + 1 << "_beam_" << wu.get_beam() << ".ms";

  const string outms = pstr.str();
  pstr << ".working";

  const string outms_flag = pstr.str();



  if (itsComms.inGroup(0)) {

    struct stat buffer;

    while (stat(outms_flag.c_str(), &buffer) == 0) {
      // flag file exists - someone is writing

      sleep(1);
    }
    if (stat(outms.c_str(), &buffer) == 0) {
      ASKAPLOG_WARN_STR(logger, "Split file already exists");
    } else if (stat(outms.c_str(), &buffer) != 0 && stat(outms_flag.c_str(), &buffer) != 0) {
      // file cannot be read

      // drop trigger
      ofstream trigger;
      trigger.open(outms_flag.c_str());
      trigger.close();
      MSSplitter mySplitter(unitParset);

      mySplitter.split(wu.get_dataset(), outms, wu.get_localChannel() + 1, wu.get_localChannel() + 1, 1, unitParset);
      unlink(outms_flag.c_str());
      this->cached_files.push_back(outms);

    }
    ///wait for all groups this rank to get here
    if (itsComms.nGroups() > 1) {
      ASKAPLOG_DEBUG_STR(logger, "Rank " << itsComms.rank() << " at barrier");
      itsComms.barrier(itsComms.interGroupCommIndex());
      ASKAPLOG_DEBUG_STR(logger, "Rank " << itsComms.rank() << " passed barrier");
    }
    wu.set_dataset(outms);

  }
  else {
    ASKAPLOG_WARN_STR(logger,"Cache MS requested but not done");
  }


}
void ContinuumWorker::processWorkUnit(ContinuumWorkUnit& wu)
{



  // This also needs to set the frequencies and directions for all the images
  ASKAPLOG_DEBUG_STR(logger, "In processWorkUnit");
  LOFAR::ParameterSet unitParset = itsParset;
  ASKAPLOG_DEBUG_STR(logger, "Parset Reports: (In process workunit)" << (itsParset.getStringVector("dataset", true)));

  char ChannelPar[64];

  sprintf(ChannelPar, "[1,%d]", wu.get_localChannel() + 1);
  bool perbeam = unitParset.getBool("perbeam", true);

  if (!perbeam) {
    string param = "beams";
    std::ostringstream bstr;

    bstr << "[" << wu.get_beam() << "]";

    unitParset.replace(param, bstr.str().c_str());
  }

  bool usetmpfs = unitParset.getBool("usetmpfs", false);
  bool localsolve = unitParset.getBool("solverpercore", false);

  if (usetmpfs && !localsolve) // only do this here if in continuum mode

  {
    cacheWorkUnit(wu, unitParset);
    sprintf(ChannelPar, "[1,1]");
  }

  unitParset.replace("Channels", ChannelPar);

  ASKAPLOG_DEBUG_STR(logger, "Getting advice on missing parameters");

  itsAdvisor->addMissingParameters(unitParset);

  ASKAPLOG_DEBUG_STR(logger, "Storing workUnit");
  workUnits.insert(workUnits.begin(),wu); //
  //workUnits.push_back(wu);
  ASKAPLOG_DEBUG_STR(logger, "Storing parset");
  itsParsets.insert(itsParsets.begin(),unitParset);
  // itsParsets.push_back(unitParset);
  ASKAPLOG_DEBUG_STR(logger, "Finished processWorkUnit");
  ASKAPLOG_DEBUG_STR(logger, "Parset Reports (leaving processWorkUnit): " << (itsParset.getStringVector("dataset", true)));

}


void ContinuumWorker::processSnapshot(LOFAR::ParameterSet& unitParset)
{
}
void ContinuumWorker::processChannels()
{
  ASKAPLOG_INFO_STR(logger, "Processing Channel Allocation");

  LOFAR::ParameterSet& unitParset = itsParset;

  if (workUnits.size() > 0) {
    unitParset = itsParsets[0];
  }

  const bool dumpgrids = unitParset.getBool("dumpgrids",false);

  if (dumpgrids) {
    ASKAPLOG_INFO_STR(logger,"Will output gridded visibilities");
  }

  const bool localSolver = unitParset.getBool("solverpercore", false);

  if (localSolver) {
    ASKAPLOG_INFO_STR(logger, "Processing multiple channels local solver mode");
  }
  else {
    ASKAPLOG_INFO_STR(logger, "Processing multiple channels in central solver mode");
  }



  const int nwriters = itsParset.getInt32("nwriters",1);
  ASKAPCHECK(nwriters>0,"Number of writers must be greater than 0");

  const bool updateDir = itsParset.getBool("updatedirection",false);

  ASKAPCHECK(!(updateDir && !localSolver), "Cannot <yet> Continuum image in on-the-fly mosaick mode - need to update the image parameter setup");

  // Define reference channel for giving restoring beam
  std::string reference = itsParset.getString("restore.beamReference", "mid");
  if (reference == "mid") {
    itsBeamReferenceChannel = nchanCube / 2;
  } else if (reference == "first") {
    itsBeamReferenceChannel = 0;
  } else if (reference == "last") {
    itsBeamReferenceChannel = nchanCube - 1;
  } else { // interpret reference as a 0-based channel nuumber
    unsigned int num = atoi(reference.c_str());
    if (num < nchanCube) {
      itsBeamReferenceChannel = num;
    } else {
      ASKAPLOG_WARN_STR(logger, "beamReference value (" << reference
      << ") not valid. Using middle value of " << nchanCube / 2);
      itsBeamReferenceChannel = nchanCube / 2;
    }
  }

  if (itsComms.isWriter()) {

    Quantity f0(this->baseCubeFrequency, "Hz");
    /// The width of a channel. THis does <NOT> take account of the variable width
    /// of Barycentric channels
    Quantity freqinc(workUnits[0].get_channelWidth(), "Hz");

    std::string root = "image";

    std::string img_name = root + std::string(".wr.") \
    + utility::toString(itsComms.rank());

    root = "psf";
    std::string psf_name = root + std::string(".wr.") \
    + utility::toString(itsComms.rank());

    root = "residual";

    std::string residual_name = root + std::string(".wr.") \
    + utility::toString(itsComms.rank());

    root = "weights";

    std::string weights_name = root + std::string(".wr.") \
    + utility::toString(itsComms.rank());

    root = "grid";

    std::string grid_name = root + std::string(".wr.") \
    + utility::toString(itsComms.rank());

    root = "pcf";
    std::string pcf_name = root + std::string(".wr.") \
    + utility::toString(itsComms.rank());
      
    if (itsComms.isSingleSink()) {
      // Need to reset the names to something eveyone knows
      img_name = "image";
      psf_name = "psf";
      residual_name = "residual";
      weights_name = "weights";
      grid_name = "grid";
      pcf_name = "pcf";
    }

    ASKAPLOG_DEBUG_STR(logger, "Configuring Spectral Cube");
    ASKAPLOG_DEBUG_STR(logger, "nchan: " << this->nchanCube << " base f0: " << f0.getValue("MHz")
    << " width: " << freqinc.getValue("MHz") << " (" << workUnits[0].get_channelWidth() << ")");


    if ( itsComms.isCubeCreator() ) {
      itsImageCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc,img_name));
      itsPSFCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, psf_name));
      itsResidualCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, residual_name));
      itsWeightsCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, weights_name));
      
      if ( dumpgrids ) {
        itsGriddedVis.reset(new CubeBuilder<casacore::Complex>(itsParset, this->nchanCube, f0, freqinc, grid_name));
        itsPCFCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, pcf_name));
      }
      

    }


    if (!itsComms.isCubeCreator()) {
      itsImageCube.reset(new CubeBuilder<casacore::Float>(itsParset, img_name));
      itsPSFCube.reset(new CubeBuilder<casacore::Float>(itsParset,  psf_name));
      itsResidualCube.reset(new CubeBuilder<casacore::Float>(itsParset,  residual_name));
      itsWeightsCube.reset(new CubeBuilder<casacore::Float>(itsParset, weights_name));

 
      if ( dumpgrids ) {
        itsGriddedVis.reset(new CubeBuilder<casacore::Complex>(itsParset, grid_name));
        itsPCFCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, pcf_name));

      }
      
    }

    if (itsParset.getBool("restore", false)) {
      root = "psf.image";
      std::string psf_image_name = root + std::string(".wr.") \
      + utility::toString(itsComms.rank());
      root = "image.restored";
      std::string restored_image_name = root + std::string(".wr.") \
      + utility::toString(itsComms.rank());
      if (itsComms.isSingleSink()) {
        // Need to reset the names to something eveyone knows
        psf_image_name = "psf.image";
        restored_image_name = "image.restored";

      }
      // Only create these if we are restoring, as that is when they get made
      if (itsComms.isCubeCreator()) {
        if (itsDoingPreconditioning) {
          itsPSFimageCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, psf_image_name));
        }
        itsRestoredCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, restored_image_name));
      }

      if (!itsComms.isCubeCreator()) {
        if (itsDoingPreconditioning) {
          itsPSFimageCube.reset(new CubeBuilder<casacore::Float>(itsParset,  psf_image_name));
        }
        itsRestoredCube.reset(new CubeBuilder<casacore::Float>(itsParset, restored_image_name));
      }
    }
  }

  ASKAPLOG_DEBUG_STR(logger, "You shall not pass. Waiting at a barrier for all ranks to have created the cubes ");
  itsComms.barrier(itsComms.theWorkers());
  ASKAPLOG_DEBUG_STR(logger, "Passed the barrier");

  if (workUnits.size() == 0) {
    ASKAPLOG_INFO_STR(logger,"No work todo");

    // write out the beam log
    ASKAPLOG_INFO_STR(logger, "About to log the full set of restoring beams");
    logBeamInfo();

    return;
  }

  /// What are the plans for the deconvolution?
  ASKAPLOG_DEBUG_STR(logger, "Ascertaining Cleaning Plan");
  const bool writeAtMajorCycle = itsParsets[0].getBool("Images.writeAtMajorCycle", false);
  const int nCycles = itsParset.getInt32("ncycles", 0);
  std::string majorcycle = itsParset.getString("threshold.majorcycle", "-1Jy");
  const double targetPeakResidual = SynthesisParamsHelper::convertQuantity(majorcycle, "Jy");

  const int uvwMachineCacheSize = itsParset.getInt32("nUVWMachines", 1);
  ASKAPCHECK(uvwMachineCacheSize > 0 ,
    "Cache size is supposed to be a positive number, you have "
    << uvwMachineCacheSize);

  const double uvwMachineCacheTolerance = SynthesisParamsHelper::convertQuantity(itsParset.getString("uvwMachineDirTolerance", "1e-6rad"), "rad");

  ASKAPLOG_DEBUG_STR(logger,
      "UVWMachine cache will store " << uvwMachineCacheSize << " machines");
      ASKAPLOG_DEBUG_STR(logger, "Tolerance on the directions is "
      << uvwMachineCacheTolerance / casacore::C::pi * 180. * 3600. << " arcsec");


  // the workUnits may include different epochs (for the same channel)
  // the order is strictly by channel - with multiple work units per channel.
  // so you can increment the workUnit until the frequency changes - then you know you
  // have all the workunits for that channel

  boost::shared_ptr<CalcCore> rootImagerPtr;
  bool gridder_initialized = false;

  for (int workUnitCount = 0; workUnitCount < workUnits.size();) {

    // NOTE:not all of these will have work
    // NOTE:this loop does not increment here.

    try {

      // spin for good workunit
      while (workUnitCount <= workUnits.size()) {
        if (workUnits[workUnitCount].get_payloadType() == ContinuumWorkUnit::DONE){
          workUnitCount++;
        }
        else if (workUnits[workUnitCount].get_payloadType() == ContinuumWorkUnit::NA) {
          if (itsComms.isWriter()) {
            // itsComms.removeChannelFromWriter(itsComms.rank());
            ASKAPLOG_WARN_STR(logger,"No longer removing whole channel from write as work allocation is bad. This may not work for multiple epochs");
          }
          workUnitCount++;
        }
        else {
          ASKAPLOG_INFO_STR(logger, "Good workUnit at number " << workUnitCount);
          break;
        }
      }
      if (workUnitCount >= workUnits.size()) {
        ASKAPLOG_INFO_STR(logger, "Out of work with workUnit " << workUnitCount);
        break;
      }

      ASKAPLOG_INFO_STR(logger, "Starting to process workunit " << workUnitCount+1 << " of " << workUnits.size());

      int initialChannelWorkUnit = workUnitCount;

      if (!updateDir) {

        // NOTE: this is because if we are mosaicking ON THE FLY. We do
        // not process the first workunit outside the imaging loop.
        // But for "normal" processing the first workunit is processed outside the loops
        // This adds alsorts of complications to the logic BTW.

        initialChannelWorkUnit = workUnitCount+1;
      }

      double frequency=workUnits[workUnitCount].get_channelFrequency();
      const string colName = itsParsets[workUnitCount].getString("datacolumn", "DATA");


      int localChannel;
      int globalChannel;

      bool usetmpfs = itsParsets[workUnitCount].getBool("usetmpfs", false);
      if (usetmpfs) {
        // probably in spectral line mode
        // copy the caching here ...
        cacheWorkUnit(workUnits[workUnitCount], itsParsets[workUnitCount]);

        localChannel = 0;

      } else {
        localChannel = workUnits[workUnitCount].get_localChannel();
      }

      const string ms = workUnits[workUnitCount].get_dataset();
      globalChannel = workUnits[workUnitCount].get_globalChannel();

      // MEMORY_BUFFERS mode opens the MS readonly
      TableDataSource ds(ms, TableDataSource::MEMORY_BUFFERS, colName);

      /// Need to set up the rootImager here
      if (updateDir == true) {
        itsAdvisor->updateDirectionFromWorkUnit(itsParsets[workUnitCount],workUnits[workUnitCount]);
      }
      if (updateDir || !gridder_initialized) {
          
        boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParsets[workUnitCount],itsComms,ds,localChannel));
        rootImagerPtr = tempIm;
        gridder_initialized = true;
      }
      else if (gridder_initialized){
        boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParsets[workUnitCount],itsComms,ds,rootImagerPtr->gridder(),localChannel));
        rootImagerPtr = tempIm;
      }
        
      CalcCore& rootImager = *rootImagerPtr; // just for the semantics
      //// CalcCore rootImager(itsParsets[workUnitCount], itsComms, ds, localChannel);
      /// set up the image for this channel
      /// this will actually build a full image for the first - it is not actually used tho.
      ///


      bool stopping = false;

      if (!localSolver) {
        // we need to wait for the first empty model.
        ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " at barrier");
        itsComms.barrier(itsComms.theWorkers());
        ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " passed barrier");


        ASKAPLOG_INFO_STR(logger, "Worker waiting to receive new model");
        rootImager.receiveModel();
        ASKAPLOG_INFO_STR(logger, "Worker received initial model for cycle 0");
      }
      else {
        // this assumes no subimage will be formed.
        setupImage(rootImager.params(), frequency,false);
      }

      ImagingNormalEquations &rootINERef =
      dynamic_cast<ImagingNormalEquations&>(*rootImager.getNE());


      try {

        rootImager.calcNE(); // why do this -
        // this essentially forces me to
        // image the full FOV for a single beam
        // but all I want is somthing to linmos into.
        // But I need this for the solver ....
        // I should find a away to get the NE initialised w/o regridding
        // which would be much better.
        // Why not just use a spheroidal for the PSF gridders / full FOV
        // FIXME
        if (updateDir == true) {
          rootINERef.weightType(FROM_WEIGHT_IMAGES);
          rootINERef.weightState(CORRECTED);
          rootImager.zero(); // then we delete all our work ....
        }

      }
      catch (const askap::AskapError& e) {
        ASKAPLOG_WARN_STR(logger,"Askap error in worker calcNE - rootImager failed");
        ASKAPLOG_WARN_STR(logger,"Incrementing workunit count as this one failed");
        workUnitCount++;

        throw;
      }

      /// need to put in the major and minor cycle loops
      /// If we are doing more than one major cycle I need to reset
      /// the workUnit count to permit a re-read of the input data.
      /// LOOP:

      /// For continuum we need to loop over epochs/beams and frequencies
      /// For "localSolver" or continuum we process each freuqency in turn.

      if (nCycles == 0) {
        stopping = true;
      }

      for (int majorCycleNumber = 0; majorCycleNumber <= nCycles; ++majorCycleNumber) {
        // NOTE: within this loop the workUnit is incremented.
        // so we need to check whether the frequency changes.
        // Perhaps something cleaner is needed.

        int tempWorkUnitCount = initialChannelWorkUnit;
        // clearer if it were called nextWorkUnit - but this is essentially the workunit we are starting this loop on.



        // now we are going to actually image this work unit
        // This loops over work units that are the same baseFrequency
        // but probably not the same epoch or beam ....


        while (tempWorkUnitCount < workUnits.size())   {

          /// need a working imager to allow a merge over epochs for this channel
          /// assuming subsequency workunits are the same channel but either different
          /// epochs or look directions.

          if (frequency != workUnits[tempWorkUnitCount].get_channelFrequency()) {
            if (localSolver) { // the frequencies should be the same.
              // THis is probably the normal spectral line or continuum cube mode.
              // each workunit is a different frequency
              ASKAPLOG_INFO_STR(logger,"Change of frequency for workunit");
              break;
            }
          }

          if (usetmpfs) {
            // probably in spectral line mode
            cacheWorkUnit(workUnits[tempWorkUnitCount], itsParsets[tempWorkUnitCount]);

            localChannel = 0;

          } else {
            localChannel = workUnits[tempWorkUnitCount].get_localChannel();
          }
          globalChannel = workUnits[tempWorkUnitCount].get_globalChannel();

          const string myMs = workUnits[tempWorkUnitCount].get_dataset();
          TableDataSource myDs(myMs, TableDataSource::DEFAULT, colName);
          myDs.configureUVWMachineCache(uvwMachineCacheSize, uvwMachineCacheTolerance);
          try {

            boost::shared_ptr<CalcCore> workingImagerPtr;

            if (updateDir) {
              itsAdvisor->updateDirectionFromWorkUnit(itsParsets[tempWorkUnitCount],workUnits[tempWorkUnitCount]);
              ///FIXME:
              // in updateDir mode I cannot cache the gridders as they have a tangent point.
              // So I have turned of caching for all modes. THis is performance hit on everyone for
              // a corner case .... FIX This!
              // FIXED: by just having 2 possible working imagers depending on the mode. ... easy really

              boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParsets[tempWorkUnitCount],itsComms,myDs,localChannel));
              workingImagerPtr = tempIm;
            }
            else {
              boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParsets[tempWorkUnitCount],itsComms,myDs,rootImager.gridder(),localChannel));
              workingImagerPtr = tempIm;
            }

            CalcCore& workingImager = *workingImagerPtr; // just for the semantics

            ///this loop does the calcNE and the merge of the residual images


            bool useSubSizedImages = false;

            if (updateDir) {

              useSubSizedImages = true;

              setupImage(workingImager.params(), frequency, useSubSizedImages);

              if (majorCycleNumber > 0) {
                copyModel(rootImager.params(),workingImager.params());
              }

            }
            else {
              workingImager.replaceModel(rootImager.params());
            }
            // grid and image
            try {
              workingImager.calcNE();
            }
            catch (const askap::AskapError& e) {
              ASKAPLOG_WARN_STR(logger,"Askap error in worker calcNE");
              // if this failed but the root did not one of two things may have happened
              // in continuum mode the gridding fails due to w projection errors - which
              // were not apparent in lower frequency observations - we have to just keep throwing
              // the exception up the tree in this case because we cannot recover.
              // in spectral line mode - this epoch/beam may have failed but other epochs succeeded.
              // what to do here. Do we continue with the accumulation or just fail ...
              throw;
            }

            // merge into root image if required.
            // this is required if there is more than one workunit per channel
            // either in time or by beam.

            ASKAPLOG_INFO_STR(logger,"About to merge into rootImager");
            ImagingNormalEquations &workingINERef =
            dynamic_cast<ImagingNormalEquations&>(*workingImager.getNE());
            if (updateDir) {
              workingINERef.weightType(FROM_WEIGHT_IMAGES);
              workingINERef.weightState(CORRECTED);
            }

            rootImager.getNE()->merge(*workingImager.getNE());
            ASKAPLOG_INFO_STR(logger,"Merged");


          }
          catch( const askap::AskapError& e) {
            ASKAPLOG_WARN_STR(logger, "Askap error in imaging - skipping accumulation: carrying on - this will result in a blank channel" << e.what());
            std::cerr << "Askap error in: " << e.what() << std::endl;
          }

          if (frequency == workUnits[tempWorkUnitCount].get_channelFrequency()) {
            tempWorkUnitCount++;
            // NOTE: here we increment the workunit count.
            // but the frequency in the same so this is just combining epochs or beams.
            // the accumulator does <not> have to be clean.
          }
          else {
            // the frequency has changed - which means for spectral line we break.
            // but for continuum we continue ...
            // this first condition has already been checked earlier in the loop.
            if (localSolver) {
              break;
            }
            else {
              // update the frequency
              frequency = workUnits[tempWorkUnitCount].get_channelFrequency();
              // we are now in the next channel
              // NOTE: we also need to increment the tempWorkUnitCount.
              tempWorkUnitCount++;

            }
          }

        }

        workUnitCount = tempWorkUnitCount; // this is to remember what finished on (important for localSolver).
        /// now if we are in spectral line mode we have a "full" set of NE we can SolveNE to update the model
        /// the solving is either done locally - or sent to a "master" for Solving
        /// IF dont locally then we solve - update the model and go again until we reach the majorcycle count.




        if (localSolver && (majorCycleNumber == nCycles)) { // done the last cycle
          stopping = true;
          break;
        }



        else if (!localSolver){ // probably continuum mode ....
          // If we are in continuum mode we have probaby ran through the whole allocation
          // lets send it to the master for processing.
          rootImager.sendNE();
          // now we have to wait for the model (solution) to come back.
          // we need to wait for the first empty model.
          ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " at barrier");
          itsComms.barrier(itsComms.theWorkers());
          ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " passed barrier");
          if (!stopping) { // if set then the master will not be sending a model
            ASKAPLOG_INFO_STR(logger, "Worker waiting to receive new model");
            rootImager.receiveModel();
            ASKAPLOG_INFO_STR(logger, "Worker received model for use in cycle " << majorCycleNumber+1);
          }
          else { // stopping == true.
            ASKAPLOG_INFO_STR(logger,"Worker stopping, the master will not be sending a new model");
            break;
          }

        }
        // check the model - have we reached a stopping threshold.

        if (rootImager.params()->has("peak_residual")) {
          const double peak_residual = rootImager.params()->scalarValue("peak_residual");
          ASKAPLOG_INFO_STR(logger, "Reached peak residual of " << peak_residual);
          if (peak_residual < targetPeakResidual) {
            ASKAPLOG_INFO_STR(logger, "It is below the major cycle threshold of "
            << targetPeakResidual << " Jy. Stopping.");
            stopping = true;
          }
          else {
            if (targetPeakResidual < 0) {
              ASKAPLOG_INFO_STR(logger, "Major cycle flux threshold is not used.");
            }
            else {
              ASKAPLOG_INFO_STR(logger, "It is above the major cycle threshold of "
              << targetPeakResidual << " Jy. Continuing.");
            }
          }

        }

        if (!localSolver && (majorCycleNumber == nCycles -1)) {
          stopping = true;
        }

        if (!stopping && localSolver) {
          try {
            rootImager.solveNE();
          } catch (const askap::AskapError& e) {
            ASKAPLOG_WARN_STR(logger, "Askap error in solver:" << e.what());

            throw;
          }
        }
        else if (stopping && localSolver) {
          break; // should be done if I am in local solver mode.
        }

        if (!stopping && updateDir){

          /// But we dont want to keep merging into the same NE
          /// so lets reset
          ASKAPLOG_INFO_STR(logger, "Continuuing - Reset normal equations");

          // this implies all workunits are processed independently including the first one - so I can completely
          // empty the NE

          // Actually I've found that I cannot completely empty the NE. As I need the full size PSF and this is stored in the NE
          // So this method pretty much only zeros the weights and the datavector(image)

          rootImager.zero();

          // the model is now updated but the NE are empty ... - lets go again
          // well they are not completely empty - the PSF is still there but the weights and image are zero
        }
        else if (!stopping && !updateDir) {
          // In this case the first workUnit is processed outside the workUnit loop.
          // So we need to calcNE again with the latest model before the major cycle starts.
          //
          // If we are using updateDir we reprocess all the workunits - so this is not needed.
          ASKAPLOG_INFO_STR(logger, "Continuuing - Reset normal equations");
          rootImager.getNE()->reset();

          // we have found that resetting the NE is causing some problems after r10290.

          try {
            rootImager.calcNE();
          }
          catch (const askap::AskapError& e) {
            ASKAPLOG_WARN_STR(logger, "Askap error in calcNE after majorcycle: " << e.what());

          }
        }
        else if (stopping && !localSolver) {
          ASKAPLOG_INFO_STR(logger, "Not local solver but last run - Reset normal equations");
          rootImager.getNE()->reset();

          if (!updateDir) {

            try {
              rootImager.calcNE();
            }
            catch (const askap::AskapError& e) {
              ASKAPLOG_WARN_STR(logger, "Askap error in calcNE after majorcycle: " << e.what());

            }

          }

        }



      }
      ASKAPLOG_INFO_STR(logger," Finished the major cycles");



      if (!localSolver) { // all my work is done - only continue if in local mode
        ASKAPLOG_INFO_STR(logger,"Finished imaging");
        ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " at barrier");
        itsComms.barrier(itsComms.theWorkers());
        ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " passed barrier");

        // write out the beam log
        ASKAPLOG_INFO_STR(logger, "About to log the full set of restoring beams");
        logBeamInfo();

        return;
      }

      if (localSolver) {
        rootImager.updateSolver();
      }

      // At this point we have finished our last major cycle. We have the "best" model from the
      // last minor cycle. Which should be in the archive - or full coordinate system
      // the residual image should be merged into the archive coordinated as well.

      ASKAPLOG_INFO_STR(logger,"Adding model.slice");
      ASKAPCHECK(rootImager.params()->has("image.slice"), "Params are missing image.slice parameter");
      rootImager.params()->add("model.slice", rootImager.params()->value("image.slice"));
      ASKAPCHECK(rootImager.params()->has("model.slice"), "Params are missing model.slice parameter");

      if (dumpgrids) {
        ASKAPLOG_INFO_STR(logger,"Adding grid.slice");
        casacore::Array<casacore::Complex> garr = rootImager.getGrid();
        casacore::Vector<casacore::Complex> garrVec(garr.reform(IPosition(1,garr.nelements())));
        rootImager.params()->addComplexVector("grid.slice",garrVec);
      } 

      rootImager.check();


      if (itsParsets[0].getBool("restore", false)) {
        ASKAPLOG_INFO_STR(logger, "Running restore");
        rootImager.restoreImage();
      }

      if (usetmpfs) {
        ASKAPLOG_INFO_STR(logger, "clearing cache");
        clearWorkUnitCache();
        ASKAPLOG_INFO_STR(logger, "done clearing cache");

      }

      ASKAPLOG_INFO_STR(logger, "writing channel into cube");

      if (itsComms.isWriter()) {

        ASKAPLOG_INFO_STR(logger, "I have (including my own) " << itsComms.getOutstanding() << " units to write");
        ASKAPLOG_INFO_STR(logger, "I have " << itsComms.getClients().size() << " clients with work");
        int cubeChannel = workUnits[workUnitCount - 1].get_globalChannel() - this->baseCubeGlobalChannel;
        ASKAPLOG_INFO_STR(logger, "Attempting to write channel " << cubeChannel << " of " << this->nchanCube);
        ASKAPCHECK((cubeChannel >= 0 || cubeChannel < this->nchanCube), "cubeChannel outside range of cube slice");
        handleImageParams(rootImager.params(), cubeChannel);
        ASKAPLOG_INFO_STR(logger, "Written channel " << cubeChannel);

        itsComms.removeChannelFromWriter(itsComms.rank());

        itsComms.removeChannelFromWorker(itsComms.rank());

        /// write everyone elses

        /// one per client ... I dont care what order they come in at

        int targetOutstanding = itsComms.getOutstanding() - itsComms.getClients().size();
        if (targetOutstanding < 0) {
          targetOutstanding = 0;
        }
        ASKAPLOG_INFO_STR(logger, "this iteration target is " << targetOutstanding);
        ASKAPLOG_INFO_STR(logger, "iteration count is " << itsComms.getOutstanding());

        while (itsComms.getOutstanding() > targetOutstanding) {
          if (itsComms.getOutstanding() <= (workUnits.size() - workUnitCount)) {
            ASKAPLOG_INFO_STR(logger, "local remaining count is " << (workUnits.size() - workUnitCount)) ;

            break;
          }


          ContinuumWorkRequest result;

          int id;
          /// this is a blocking receive
          ASKAPLOG_INFO_STR(logger, "Waiting for a write request");
          result.receiveRequest(id, itsComms);
          ASKAPLOG_INFO_STR(logger, "Received a request to write from rank " << id);
          int cubeChannel = result.get_globalChannel() - this->baseCubeGlobalChannel;

          try {
            ASKAPLOG_INFO_STR(logger, "Attempting to write channel " << cubeChannel << " of " << this->nchanCube);
            ASKAPCHECK((cubeChannel >= 0 || cubeChannel < this->nchanCube), "cubeChannel outside range of cube slice");

            handleImageParams(result.get_params(), cubeChannel);

            ASKAPLOG_INFO_STR(logger, "Written the slice from rank" << id);
          }

          catch (const askap::AskapError& e) {
            ASKAPLOG_WARN_STR(logger, "Failed to write a channel to the cube: " << e.what());
          }

          itsComms.removeChannelFromWriter(itsComms.rank());
          ASKAPLOG_INFO_STR(logger, "this iteration target is " << targetOutstanding);
          ASKAPLOG_INFO_STR(logger, "iteration count is " << itsComms.getOutstanding());
        }

      } else {

        ContinuumWorkRequest result;
        result.set_params(rootImager.params());
        result.set_globalChannel(workUnits[workUnitCount - 1].get_globalChannel());
        /// send the work to the writer with a blocking send
        result.sendRequest(workUnits[workUnitCount - 1].get_writer(), itsComms);
        itsComms.removeChannelFromWorker(itsComms.rank());

      }

      /// outside the clean-loop write out the slice
    }

    catch (const askap::AskapError& e) {

      if (!localSolver) {
        ASKAPLOG_WARN_STR(logger, "Askap error processing a channel in continuum mode");
        throw;
      }

      ASKAPLOG_WARN_STR(logger, "Askap error in channel processing skipping: " << e.what());
      std::cerr << "Askap error in: " << e.what() << std::endl;

      // Need to either send an empty map - or
      if (itsComms.isWriter()) {
        ASKAPLOG_INFO_STR(logger, "Marking bad channel as processed in count for writer\n");
        itsComms.removeChannelFromWriter(itsComms.rank());
      } else {
        int goodUnitCount = workUnitCount - 1; // last good one - needed for the correct freq label and writer
        ASKAPLOG_INFO_STR(logger, "Failed on count " << goodUnitCount);
        ASKAPLOG_INFO_STR(logger, "Sending blankparams to writer " << workUnits[goodUnitCount].get_writer());
        askap::scimath::Params::ShPtr blankParams;

        blankParams.reset(new Params);
        ASKAPCHECK(blankParams, "blank parameters (images) not initialised");

        setupImage(blankParams, workUnits[goodUnitCount].get_channelFrequency());


        ContinuumWorkRequest result;
        result.set_params(blankParams);
        result.set_globalChannel(workUnits[goodUnitCount].get_globalChannel());
        /// send the work to the writer with a blocking send
        result.sendRequest(workUnits[goodUnitCount].get_writer(), itsComms);
        ASKAPLOG_INFO_STR(logger, "Sent\n");
      }
      // No need to increment workunit. Although this assumes that we are here becuase we failed the solveNE not the calcNE


    }

    catch (const std::exception& e) {
      ASKAPLOG_WARN_STR(logger, "Unexpected exception in: " << e.what());
      std::cerr << "Unexpected exception in: " << e.what();
      // I need to repeat the bookkeeping here as errors other than AskapErrors are thrown by solveNE
      if (!localSolver) {
        /// this is MFS/continuum mode
        /// throw this further up - this avoids a failure in continuum mode generating bogus - or furphy-like
        /// error messages
        throw e;
      }

      // Need to either send an empty map - or
      if (itsComms.isWriter()) {
        ASKAPLOG_INFO_STR(logger, "Marking bad channel as processed in count for writer\n");
        itsComms.removeChannelFromWriter(itsComms.rank());
      } else {
        int goodUnitCount = workUnitCount - 1; // last good one - needed for the correct freq label and writer
        ASKAPLOG_INFO_STR(logger, "Failed on count " << goodUnitCount);
        ASKAPLOG_INFO_STR(logger, "Sending blankparams to writer " << workUnits[goodUnitCount].get_writer());
        askap::scimath::Params::ShPtr blankParams;

        blankParams.reset(new Params);
        ASKAPCHECK(blankParams, "blank parameters (images) not initialised");

        setupImage(blankParams, workUnits[goodUnitCount].get_channelFrequency());


        ContinuumWorkRequest result;
        result.set_params(blankParams);
        result.set_globalChannel(workUnits[goodUnitCount].get_globalChannel());
        /// send the work to the writer with a blocking send
        result.sendRequest(workUnits[goodUnitCount].get_writer(), itsComms);
        ASKAPLOG_INFO_STR(logger, "Sent\n");
      }
      // No need to increment workunit. Although this assumes that we are here becuase we failed the solveNE not the calcNE

    }

  } // next workunit if required.

  // cleanup
  if (itsComms.isWriter()) {

    while (itsComms.getOutstanding() > 0) {
      ASKAPLOG_INFO_STR(logger, "I have " << itsComms.getOutstanding() << "outstanding work units");
      ContinuumWorkRequest result;
      int id;
      result.receiveRequest(id, itsComms);
      ASKAPLOG_INFO_STR(logger, "Received a request to write from rank " << id);
      int cubeChannel = result.get_globalChannel() - this->baseCubeGlobalChannel;
      try {


        ASKAPLOG_INFO_STR(logger, "Attempting to write channel " << cubeChannel << " of " << this->nchanCube);
        ASKAPCHECK((cubeChannel >= 0 || cubeChannel < this->nchanCube), "cubeChannel outside range of cube slice");
        handleImageParams(result.get_params(), cubeChannel);
        ASKAPLOG_INFO_STR(logger, "Written the slice from rank" << id);

      } catch (const askap::AskapError& e) {
        ASKAPLOG_WARN_STR(logger, "Failed to write a channel to the cube: " << e.what());
      }

      itsComms.removeChannelFromWriter(itsComms.rank());
    }
  }

  // write out the beam log
  ASKAPLOG_INFO_STR(logger, "About to log the full set of restoring beams");
  logBeamInfo();

}

void ContinuumWorker::copyModel(askap::scimath::Params::ShPtr SourceParams, askap::scimath::Params::ShPtr SinkParams)
{
  askap::scimath::Params& src = *SourceParams;
  askap::scimath::Params& dest = *SinkParams;
  // ASKAPLOG_WARN_STR(logger, "Names are " << src.names());
  // before the restore the image is the model ....
  SynthesisParamsHelper::copyImageParameter(src, dest,"image.slice");

}
void ContinuumWorker::handleImageParams(askap::scimath::Params::ShPtr params, unsigned int chan)
{


  // Pre-conditions
  if (!params->has("model.slice")) {
    ASKAPLOG_WARN_STR(logger, "Params are missing model parameter");

  } else // Write image
  {
    ASKAPLOG_INFO_STR(logger, "Writing model for (local) channel " << chan);
    const casacore::Array<double> imagePixels(params->value("model.slice"));
    casacore::Array<float> floatImagePixels(imagePixels.shape());
    casacore::convertArray<float, double>(floatImagePixels, imagePixels);
    itsImageCube->writeSlice(floatImagePixels, chan);
  }



  if (!params->has("psf.slice")) {
    ASKAPLOG_WARN_STR(logger,  "Params are missing psf parameter");
  }
  else {
    // Write PSF

    ASKAPLOG_INFO_STR(logger, "Writing PSF");
    const casacore::Array<double> imagePixels(params->value("psf.slice"));
    casacore::Array<float> floatImagePixels(imagePixels.shape());
    casacore::convertArray<float, double>(floatImagePixels, imagePixels);
    itsPSFCube->writeSlice(floatImagePixels, chan);

  }
  if (!params->has("residual.slice")) {
    ASKAPLOG_WARN_STR(logger,  "Params are missing residual parameter");

  } else
  {
    // Write residual
    ASKAPLOG_INFO_STR(logger, "Writing Residual");
    const casacore::Array<double> imagePixels(params->value("residual.slice"));
    casacore::Array<float> floatImagePixels(imagePixels.shape());
    casacore::convertArray<float, double>(floatImagePixels, imagePixels);
    itsResidualCube->writeSlice(floatImagePixels, chan);
  }

  if (!params->has("weights.slice")) {
    ASKAPLOG_WARN_STR(logger, "Params are missing weights parameter");

  }
  else
  {
    ASKAPLOG_INFO_STR(logger, "Writing Weights");
    const casacore::Array<double> imagePixels(params->value("weights.slice"));
    casacore::Array<float> floatImagePixels(imagePixels.shape());
    casacore::convertArray<float, double>(floatImagePixels, imagePixels);
    itsWeightsCube->writeSlice(floatImagePixels, chan);
  }
  // Write the grids

  if (params->has("grid.slice")) {
    ASKAPLOG_INFO_STR(logger, "Writing Grid");
    const casacore::Vector<casacore::Complex> gr(params->complexVectorValue("grid.slice"));
    casacore::Array<casacore::Complex> grid(gr.reform(params->value("psf.slice").shape()));
    itsGriddedVis->writeSlice(grid,chan);
  }

  if (itsParset.getBool("restore", false)) {
    ASKAPCHECK(params->has("image.slice"), "Params are missing image parameter");
    if (itsDoingPreconditioning) {
      ASKAPCHECK(params->has("psf.image.slice"), "Params are missing psf.image parameter");
    }
  }

  if (itsParset.getBool("restore", false)) {
    // Record the restoring beam
    const askap::scimath::Axes &axes = params->axes("image.slice");
    recordBeam(axes, chan);
  }

  if (itsParset.getBool("restore", false)) {

    if (itsDoingPreconditioning) {
      // Write preconditioned PSF image
      {
        ASKAPLOG_INFO_STR(logger, "Writing preconditioned PSF");
        const casacore::Array<double> imagePixels(params->value("psf.image.slice"));
        casacore::Array<float> floatImagePixels(imagePixels.shape());
        casacore::convertArray<float, double>(floatImagePixels, imagePixels);
        itsPSFimageCube->writeSlice(floatImagePixels, chan);
      }
    }

    // Write Restored image

    ASKAPLOG_INFO_STR(logger, "Writing Restored Image");
    const casacore::Array<double> imagePixels(params->value("image.slice"));
    casacore::Array<float> floatImagePixels(imagePixels.shape());
    casacore::convertArray<float, double>(floatImagePixels, imagePixels);
    itsRestoredCube->writeSlice(floatImagePixels, chan);

  }

}

void ContinuumWorker::initialiseBeamLog(const unsigned int numChannels)
{

    casa::Vector<casa::Quantum<double> > beamVec(3);
    beamVec[0] = casa::Quantum<double>(0., "rad");
    beamVec[1] = casa::Quantum<double>(0., "rad");
    beamVec[2] = casa::Quantum<double>(0., "deg");

    for(unsigned int i=0;i<numChannels;i++) {
        itsBeamList[i] = beamVec;
    }

}

void ContinuumWorker::recordBeam(const askap::scimath::Axes &axes, const unsigned int cubeChannel)
{

  if (axes.has("MAJMIN")) {
    // this is a restored image with beam parameters set
    ASKAPCHECK(axes.has("PA"), "PA axis should always accompany MAJMIN");
    ASKAPLOG_DEBUG_STR(logger, "Found beam for image.slice, channel " <<
                       cubeChannel << ", with shape " <<
                       axes.start("MAJMIN") * 180. / M_PI * 3600. << "x" <<
                       axes.end("MAJMIN") * 180. / M_PI * 3600. << ", " <<
                       axes.start("PA") * 180. / M_PI);

    casacore::Vector<casacore::Quantum<double> > beamVec(3, 0.);
    beamVec[0] = casacore::Quantum<double>(axes.start("MAJMIN"), "rad");
    beamVec[1] = casacore::Quantum<double>(axes.end("MAJMIN"), "rad");
    beamVec[2] = casacore::Quantum<double>(axes.start("PA"), "rad");

    itsBeamList[cubeChannel] = beamVec;

  }

}


void ContinuumWorker::storeBeam(const unsigned int cubeChannel)
{
  if (cubeChannel == itsBeamReferenceChannel) {
    itsRestoredCube->addBeam(itsBeamList[cubeChannel]);
  }
}

void ContinuumWorker::logBeamInfo()
{

  if (itsParset.getBool("restore", false)) {
    askap::accessors::BeamLogger beamlog;
    ASKAPLOG_INFO_STR(logger, "Channel-dependent restoring beams will be written to log file " << beamlog.filename());
    ASKAPLOG_DEBUG_STR(logger, "About to add beam list of size " << itsBeamList.size() << " to the beam logger");
    beamlog.beamlist() = itsBeamList;

    if (itsParset.getInt32("nwriters",1)>1 && itsParset.getBool("singleoutputfile",false)) {
      std::list<int> creators = itsComms.getCubeCreators();
      ASKAPASSERT(creators.size() == 1);
      int creatorRank = creators.front();
      ASKAPLOG_DEBUG_STR(logger, "Gathering all beam information beam creator is rank " << creatorRank);
      beamlog.gather(itsComms, creatorRank,false);
    }
    if (itsComms.isCubeCreator()) {

        ASKAPLOG_DEBUG_STR(logger, "Writing list of individual channel beams to beam log "
                           << beamlog.filename());
        beamlog.setFilename("beamlog." + itsRestoredCube->filename() + ".txt");
        beamlog.write();

        ASKAPLOG_DEBUG_STR(logger, "Writing restoring beam to header of restored cube");
        casa::Vector<casa::Quantum<double> > refbeam = beamlog.beam(itsBeamReferenceChannel);
        itsRestoredCube->addBeam(refbeam);

    }
  }

}


void ContinuumWorker::setupImage(const askap::scimath::Params::ShPtr& params,double channelFrequency, bool shapeOverride)
{
  try {
    ASKAPLOG_DEBUG_STR(logger, "Setting up image");
    const LOFAR::ParameterSet parset = itsParset.makeSubset("Images.");

    const int nfacets = parset.getInt32("nfacets", 1);
    const string name("image.slice");
    const vector<string> direction = parset.getStringVector("direction");
    const vector<string> cellsize = parset.getStringVector("cellsize");
    vector<int> shape = parset.getInt32Vector("shape");
    //const vector<double> freq = parset.getDoubleVector("frequency");
    const int nchan = 1;

    if (shapeOverride == true) {
      string param = "subshape";
      if (parset.isDefined(param)) {
        ASKAPLOG_INFO_STR(logger,"Over-riding image shape from parset");
        shape = parset.getInt32Vector("subshape");
        ASKAPLOG_INFO_STR(logger,"Image shape now " << shape);
      }
      else {
        ASKAPLOG_WARN_STR(logger,"Shape over-ride requested but no subshape parameter in parset");
      }
    }


    if (!parset.isDefined("polarisation")) {
      ASKAPLOG_DEBUG_STR(logger, "Polarisation frame is not defined, "
      << "only stokes I will be generated");
    }
    const vector<string> stokesVec = parset.getStringVector("polarisation",
    vector<string>(1, "I"));

    // there could be many ways to define stokes, e.g. ["XX YY"] or ["XX","YY"] or "XX,YY"
    // to allow some flexibility we have to concatenate all elements first and then
    // allow the parser from PolConverter to take care of extracting the products.
    string stokesStr;
    for (size_t i = 0; i < stokesVec.size(); ++i) {
      stokesStr += stokesVec[i];
    }
    const casacore::Vector<casacore::Stokes::StokesTypes>
    stokes = scimath::PolConverter::fromString(stokesStr);

    const bool ewProj = parset.getBool("ewprojection", false);
    if (ewProj) {
      ASKAPLOG_DEBUG_STR(logger, "Image will have SCP/NCP projection");
    } else {
      ASKAPLOG_DEBUG_STR(logger, "Image will have plain SIN projection");
    }

    ASKAPCHECK(nfacets > 0, "Number of facets is supposed to be a positive number, you gave " << nfacets);
    ASKAPCHECK(shape.size() >= 2, "Image is supposed to be at least two dimensional. " << "check shape parameter, you gave " << shape);

    if (nfacets == 1) {
      SynthesisParamsHelper::add(*params, name, direction, cellsize, shape, ewProj,
        channelFrequency, channelFrequency, nchan, stokes);
        // SynthesisParamsHelper::add(*params, name, direction, cellsize, shape, ewProj,
        //                            freq[0], freq[1], nchan, stokes);
    } else {
        // this is a multi-facet case
        const int facetstep = parset.getInt32("facetstep", casacore::min(shape[0], shape[1]));
        ASKAPCHECK(facetstep > 0,"facetstep parameter is supposed to be positive, you have " << facetstep);
        ASKAPLOG_DEBUG_STR(logger, "Facet centers will be " << facetstep << " pixels apart, each facet size will be " << shape[0] << " x " << shape[1]);
        // SynthesisParamsHelper::add(*params, name, direction, cellsize, shape, ewProj,
        //                            freq[0], freq[1], nchan, stokes, nfacets, facetstep);
        SynthesisParamsHelper::add(*params, name, direction, cellsize, shape, ewProj,channelFrequency, channelFrequency, nchan, stokes, nfacets, facetstep);
    }


  } catch (const LOFAR::APSException &ex) {
    throw AskapError(ex.what());
  }
}
