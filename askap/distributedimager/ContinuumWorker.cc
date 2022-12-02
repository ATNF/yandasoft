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
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/scimath/fitting/Equation.h>
#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/ImagingNormalEquations.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fft/FFTWrapper.h>
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
#include <askap/imageaccess/WeightsLog.h>
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
using utility::toString;

ASKAP_LOGGER(logger, ".ContinuumWorker");

ContinuumWorker::ContinuumWorker(LOFAR::ParameterSet& parset,
  CubeComms& comms, StatReporter& stats)
  : itsParset(parset), itsComms(comms), itsStats(stats), itsBeamList(),itsWeightsList()
{


    itsAdvisor = boost::shared_ptr<synthesis::AdviseDI> (new synthesis::AdviseDI(itsComms, itsParset));
    itsAdvisor->prepare();

    // lets properly size the storage
    const int nchanpercore = itsParset.getInt32("nchanpercore", 1);
    workUnits.resize(0);

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

    if (GridderStr == "AWProject" || GridderStr == "AProjectWStack") {
      itsGridderCanMosaick = true;
      ASKAPLOG_INFO_STR(logger," Gridder <CAN> mosaick");
    }
    else {
      ASKAPLOG_INFO_STR(logger,"Gridder <CANNOT> mosaick");
    }

    itsRestore = itsParset.getBool("restore", false); // do restore and write restored image
    itsWriteResidual = itsParset.getBool("residuals",false); // write residual image
    itsWriteResidual = itsParset.getBool("write.residualimage",itsWriteResidual); // alternative param name
    itsWritePsfRaw = itsParset.getBool("write.psfrawimage", false); // write unnormalised, natural wt psf
    itsWritePsfImage = itsParset.getBool("write.psfimage", true); // write normalised, preconditioned psf
    itsWriteWtLog = itsParset.getBool("write.weightslog", false); // write weights log file
    itsWriteWtImage = itsParset.getBool("write.weightsimage", false); // write weights image
    itsWriteModelImage = itsParset.getBool("write.modelimage", !itsRestore); // clean model
    itsWriteGrids = itsParset.getBool("dumpgrids", false); // write (dump) the gridded data, psf and pcf
    itsWriteGrids = itsParset.getBool("write.grids",itsWriteGrids); // new name
    itsGridType = itsParset.getString("imagetype","casa");
    itsGridCoordUV = itsParset.getBool("write.grids.uvcoord", itsGridType=="casa"); // label grid with UV coordinates
    itsGridFFT = itsParset.getBool("write.grids.fft",false); // write fft of grid (i.e. dirty image, psf)
    const int nwriters = itsParset.getInt32("nwriters",1);
    ASKAPCHECK(nwriters>0,"Number of writers must be greater than 0");
    if (itsGridType == "casa" && itsParset.getBool("singleoutputfile",false) && nwriters > 1){
      ASKAPLOG_WARN_STR(logger,"Reducing number of writers to 1 because we are writing a single casa image cube");
      itsNumWriters = 1;
    } else {
      itsNumWriters = nwriters;
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

      ASKAPLOG_DEBUG_STR(logger, "Worker has received valid allocation");
    }
    const string ms = wu.get_dataset();
    ASKAPLOG_DEBUG_STR(logger, "Received Work Unit for dataset " << ms
      << ", local (topo) channel " << wu.get_localChannel()
      << ", global (topo) channel " << wu.get_globalChannel()
      << ", frequency " << wu.get_channelFrequency() / 1.e6 << " MHz"
      << ", width " << wu.get_channelWidth() / 1e3 << " kHz");
    try {
        ASKAPLOG_DEBUG_STR(logger, "Parset Reports (before): " << (itsParset.getStringVector("dataset", true)));
        preProcessWorkUnit(wu);
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
        ASKAPLOG_INFO_STR(logger, "MultiCube with multiple writers");
        this->nchanCube = (myMaxClient - myMinClient + 1) * nchanpercore;
        this->baseCubeGlobalChannel = (myMinClient - 1) * nchanpercore;
        this->baseCubeFrequency = itsAdvisor->getBaseFrequencyAllocation((myMinClient - 1));
      } else {
        ASKAPLOG_INFO_STR(logger, "SingleCube with multiple writers");
        this->nchanCube = nchanTotal;
        this->baseCubeGlobalChannel = 0;
        this->baseCubeFrequency = itsAdvisor->getBaseFrequencyAllocation((0));
      }
      initialiseBeamLog(this->nchanCube);
      initialiseWeightsLog(this->nchanCube);

      ASKAPLOG_INFO_STR(logger, "Number of channels in cube is: " << this->nchanCube);
      ASKAPLOG_INFO_STR(logger, "Base global channel of cube is " << this->baseCubeGlobalChannel);
    }
    this->baseFrequency = itsAdvisor->getBaseFrequencyAllocation(itsComms.rank() - 1);
  }
  else {
      bool combineChannels = itsParset.getBool("combinechannels", false);
      if (combineChannels) {
          ASKAPLOG_INFO_STR(logger, "Not in localsolver (spectral line) mode - and combine channels is set so compressing channel allocations)");
          compressWorkUnits();
      }
      initialiseBeamLog(nchanTotal);
      initialiseWeightsLog(nchanTotal);

  }
  ASKAPLOG_INFO_STR(logger, "Adding all missing parameters");

  itsAdvisor->addMissingParameters(true);

  try {
    processChannels();
  } catch (AskapError& e) {
    ASKAPLOG_WARN_STR(logger, "Failure processing the channel allocation");
    ASKAPLOG_WARN_STR(logger, "Exception detail: " << e.what());
    throw;

  }

  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " finished");

  itsComms.barrier(itsComms.theWorkers());
  const bool singleoutputfile = itsParset.getBool("singleoutputfile", false);
  const bool calcstats = itsParset.getBool("calcstats", false);
  if ( singleoutputfile && calcstats ) {
    writeCubeStatistics();
    itsComms.barrier(itsComms.theWorkers());
  }
  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " passed final barrier");
}

// unused
void ContinuumWorker::deleteWorkUnitFromCache(ContinuumWorkUnit& wu)
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

void ContinuumWorker::cacheWorkUnit(ContinuumWorkUnit& wu)
{



  boost::filesystem::path mspath = boost::filesystem::path(wu.get_dataset());
  const string ms = mspath.filename().string();

  const string shm_root = itsParset.getString("tmpfs", "/dev/shm");

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
      MSSplitter mySplitter(itsParset);

      mySplitter.split(wu.get_dataset(), outms, wu.get_localChannel() + 1, wu.get_localChannel() + 1, 1, itsParset);
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
void ContinuumWorker::compressWorkUnits() {

    // This takes the list of workunits and reprocesses them so that all the contiguous
    // channels are compressed into single workUnits for multiple channels
    // this is not applicable for the spectral line experiment but can markedly reduce
    // the number of FFT required for the continuum processesing mode

    // In preProcessWorkUnit we made a list of all the channels in the allocation
    // but the workunit may contain different measurement sets so I suppose it is
    // globalChannel that is more important for the sake of allocation ... but
    // the selector only works on one measurement set.

    // So the upshot is this simple scheme cannot combine channels from different
    // measurement sets into the same grid as we are using the MS accessor as the vehicle to
    // provide the integration.

    // So we need to loop through our workunit list and make a new list that just contains a
    // single workunit for each contiguous group of channels.

    // First lets loop through our workunits

    vector<ContinuumWorkUnit> compressedList; // probably easier to generate a new list

    ContinuumWorkUnit startUnit = workUnits[0];

    unsigned int contiguousCount = 1;
    if (workUnits.size() == 1) {
        ASKAPLOG_WARN_STR(logger,"Asked to compress channels but workunit count 1");
    }
    for ( int count = 1; count < workUnits.size(); count++) {

        ContinuumWorkUnit nextUnit = workUnits[count];

        std::string startDataset = startUnit.get_dataset();
        int startChannel = startUnit.get_localChannel();
        std::string nextDataset = nextUnit.get_dataset();
        int nextChannel = nextUnit.get_localChannel();

        if ( startDataset.compare(nextDataset) == 0 ) { // same dataset
            if (nextChannel == (startChannel + contiguousCount)) { // next channel is contiguous to previous
                contiguousCount++;
                ASKAPLOG_DEBUG_STR(logger, "contiguous channel detected: count " << contiguousCount);
                startUnit.set_nchan(contiguousCount); // update the nchan count for this workunit
                // Now need to update the parset details
                string ChannelParam = "["+toString(contiguousCount)+","+toString(startUnit.get_localChannel())+"]";
                ASKAPLOG_DEBUG_STR(logger, "compressWorkUnit: ChannelParam = "<<ChannelParam);
                itsParset.replace("Channels",ChannelParam);
            }
            else { // no longer contiguous channels reset the count
                contiguousCount = 0;
            }
        }
        else { // different dataset reset the count
            ASKAPLOG_DEBUG_STR(logger, "Datasets differ resetting count");
            contiguousCount = 0;
        }
        if (count == (workUnits.size()-1) || contiguousCount == 0) { // last unit
            ASKAPLOG_DEBUG_STR(logger, "Adding unit to compressed list");
            compressedList.insert(compressedList.end(),startUnit);
            startUnit = nextUnit;
        }

    }
    if (compressedList.size() > 0) {
        ASKAPLOG_INFO_STR(logger, "Replacing workUnit list of size " << workUnits.size() << " with compressed list of size " << compressedList.size());
        ASKAPLOG_INFO_STR(logger,"A corresponding change has been made to the parset");
        workUnits = compressedList;
    }
    else {
        ASKAPLOG_WARN_STR(logger,"No compression performed");
    }
    ASKAPCHECK(compressedList.size() < 2, "The number of compressed workunits is greater than one. Channel parameters may be incorrect - see AXA-1004 and associated technical debt tickets");
}
void ContinuumWorker::preProcessWorkUnit(ContinuumWorkUnit& wu)
{
  // This also needs to set the frequencies and directions for all the images
  ASKAPLOG_DEBUG_STR(logger, "In preProcessWorkUnit");
  ASKAPLOG_DEBUG_STR(logger, "Parset Reports: (In preProcess workunit)" << (itsParset.getStringVector("dataset", true)));

  const bool localsolve = itsParset.getBool("solverpercore", false);

  // We're processing spectral data one channel at a time, but needs the stats for all, try setting this here
  // For continuum this is done in compressWorkUnits (if combinechannels is set, which it should be)
  // Channel numbers are zero based
  const int n = (localsolve ? itsParset.getInt("nchanpercore", 1) : 1);
  string ChannelPar = "["+toString(n)+","+toString(wu.get_localChannel())+"]";
  int last_beam = -1;

  // AXA-1004 this will not be unique as the parsets are passed by reference
  // if we are expecting multiple beams in the work units then these will be clobbered
  // unless the accessor gets the beam information some other way

  const bool perbeam = itsParset.getBool("perbeam", true);
  if (!perbeam) {
    string param = "beams";
    string bstr = "[" + toString(wu.get_beam()) + "]";
    if (last_beam != -1) {
      ASKAPCHECK(last_beam == wu.get_beam(), "beam index changed in perbeam processing parset - clearly in the expectation that this will do something but the parset is stored by reference AXA-1004");
    }
    itsParset.replace(param, bstr);
  }

  const bool usetmpfs = itsParset.getBool("usetmpfs", false);

  if (usetmpfs && !localsolve) {
    // only do this here if in continuum mode
    cacheWorkUnit(wu);
    ChannelPar="[1,0]";
  }

  // only add channel selection for valid workunits and topo frame
  // other frames have shifted channel allocations which can't be handled this way
  // if combinechannels==false the Channel parameter is only used for advise
  // Ord AXA-1004 removing this as it is not used by subsequent code and can break advise
  // in some corner cases.
  // ASKAPLOG_DEBUG_STR(logger, "In preProcessWorkUnit - replacing Channels parameter "<<
  // itsParset.getString("Channels","none")<<" with "<<ChannelPar<<" if topo="<<
  // unitParset.getString("freqframe","topo")<< " and "<< (wu.get_dataset()!=""));
  // if (wu.get_dataset()!="" && unitParset.getString("freqframe","topo")=="topo") {
  //     unitParset.replace("Channels", ChannelPar);
  //}

  ASKAPLOG_DEBUG_STR(logger, "Getting advice on missing parameters");

  itsAdvisor->addMissingParameters();

  ASKAPLOG_DEBUG_STR(logger, "Storing workUnit");
  workUnits.insert(workUnits.begin(),wu); //
  //workUnits.push_back(wu);
  ASKAPLOG_DEBUG_STR(logger, "Finished preProcessWorkUnit");
  ASKAPLOG_DEBUG_STR(logger, "Parset Reports (leaving preProcessWorkUnit): " << (itsParset.getStringVector("dataset", true)));
}

void ContinuumWorker::processSnapshot()
{
}
void ContinuumWorker::processChannels()
{
  ASKAPLOG_INFO_STR(logger, "Processing Channel Allocation");

  if (itsWriteGrids) {
    ASKAPLOG_INFO_STR(logger,"Will output gridded visibilities");
  }

  const bool localSolver = itsParset.getBool("solverpercore", false);

  if (localSolver) {
    ASKAPLOG_INFO_STR(logger, "Processing multiple channels local solver mode");
  }
  else {
    ASKAPLOG_INFO_STR(logger, "Processing multiple channels in central solver mode");
  }

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

    // This code is only used in the spectral line/local solver case -
    //   continuum images are written from ImagerParallel::writeModel in ContinuumMaster

    Quantity f0(this->baseCubeFrequency, "Hz");
    /// The width of a channel. THis does <NOT> take account of the variable width
    /// of Barycentric channels
    Quantity freqinc(workUnits[0].get_channelWidth(), "Hz");

    // add rank based postfix if we're writing to multiple cubes
    std::string postfix = (itsComms.isSingleSink() ? "" : std::string(".wr.") + utility::toString(itsComms.rank()));

    std::string img_name = "image" + postfix;
    std::string psf_name = "psf" + postfix;
    std::string residual_name = "residual" + postfix;
    // may need this name for the weightslog
    std::string weights_name = "weights" + postfix;
    std::string visgrid_name = "visgrid" + postfix;
    std::string pcfgrid_name = "pcfgrid" + postfix;
    std::string psfgrid_name = "psfgrid" + postfix;
    std::string psf_image_name = "psf.image" + postfix;
    std::string restored_image_name = "image.restored" + postfix;

    ASKAPLOG_DEBUG_STR(logger, "Configuring Spectral Cube");
    ASKAPLOG_DEBUG_STR(logger, "nchan: " << this->nchanCube << " base f0: " << f0.getValue("MHz")
    << " width: " << freqinc.getValue("MHz") << " (" << workUnits[0].get_channelWidth() << ")");

    if (itsWriteWtLog) {
        itsWeightsName = CubeBuilder<casacore::Float>::makeImageName(itsParset,weights_name);
    }

    LOFAR::ParameterSet gridParset = itsParset.makeSubset("");
    gridParset.remove("Images.extraoversampling");

    if ( itsComms.isCubeCreator() ) {

      // Get keywords to write to the image header
      if (!itsParset.isDefined("header.DATE-OBS")) {
        // We want the start of observations stored in the image keywords
        // The velocity calculations use the first MS for this, so we'll do that too
        casacore::MVEpoch dateObs = itsAdvisor->getEpoch(0);
        String date, timesys;
        casacore::FITSDateUtil::toFITS(date, timesys, casacore::MVTime(dateObs));
        // replace adds if non-existant
        itsParset.replace("header.DATE-OBS","["+date+",Start of observation]");
        itsParset.replace("header.TIMESYS","["+timesys+",Time System]");
      }

      if (itsWriteModelImage) {
        itsImageCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, img_name));
      }
      if (itsWritePsfRaw) {
        itsPSFCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, psf_name));
      }
      if (itsWriteResidual) {
        itsResidualCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, residual_name));
        itsResidualStatsAndMask.reset(new askap::utils::StatsAndMask(itsComms,itsResidualCube->filename(),itsResidualCube->imageHandler()));
        ASKAPLOG_INFO_STR(logger,"Created StatsAndMask object for residual cube");
      }
      if (itsWriteWtImage) {
        itsWeightsCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, weights_name));
      }
      if (itsWriteGrids) {
        if (itsGridFFT) {
          itsVisGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, visgrid_name));
          itsPCFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, pcfgrid_name));
          itsPSFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, psfgrid_name));
        } else {
          if (itsGridType == "casa") {
              itsVisGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, this->nchanCube, f0, freqinc, visgrid_name, true));
              itsPCFGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, this->nchanCube, f0, freqinc, pcfgrid_name, true));
              itsPSFGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, this->nchanCube, f0, freqinc, psfgrid_name, true));
          } else {
              itsVisGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, visgrid_name+".real", itsGridCoordUV));
              itsPCFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, pcfgrid_name+".real", itsGridCoordUV));
              itsPSFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, psfgrid_name+".real", itsGridCoordUV));
              itsVisGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, visgrid_name+".imag", itsGridCoordUV));
              itsPCFGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, pcfgrid_name+".imag", itsGridCoordUV));
              itsPSFGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, this->nchanCube, f0, freqinc, psfgrid_name+".imag", itsGridCoordUV));
          }
        }
      }
      if (itsRestore) {
        // Only create these if we are restoring, as that is when they get made
          if (itsDoingPreconditioning) {
            if (itsWritePsfImage) {
              itsPSFimageCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, psf_image_name));
            }
          }
          itsRestoredCube.reset(new CubeBuilder<casacore::Float>(itsParset, this->nchanCube, f0, freqinc, restored_image_name));
          // we are only interested to collect statistics for the restored image cube
          itsRestoredStatsAndMask.reset(new askap::utils::StatsAndMask(itsComms,itsRestoredCube->filename(),itsRestoredCube->imageHandler()));
          ASKAPLOG_INFO_STR(logger,"Created StatsAndMask object");
      }

    } else {

      if (itsWriteModelImage) {
        itsImageCube.reset(new CubeBuilder<casacore::Float>(itsParset, img_name));
      }
      if (itsWritePsfRaw) {
        itsPSFCube.reset(new CubeBuilder<casacore::Float>(itsParset, psf_name));
      }
      if (itsWriteResidual) {
        itsResidualCube.reset(new CubeBuilder<casacore::Float>(itsParset,  residual_name));
        itsResidualStatsAndMask.reset(new askap::utils::StatsAndMask(itsComms,itsResidualCube->filename(),itsResidualCube->imageHandler()));
      }
      if (itsWriteWtImage) {
        itsWeightsCube.reset(new CubeBuilder<casacore::Float>(itsParset,  weights_name));
      }

      if (itsWriteGrids) {
        if (itsGridFFT) {
          itsVisGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, visgrid_name));
          itsPCFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, pcfgrid_name));
          itsPSFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, psfgrid_name));
        } else {
          if (itsGridType == "casa") {
              itsVisGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, visgrid_name));
              itsPCFGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, pcfgrid_name));
              itsPSFGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, psfgrid_name));
          } else {
              itsVisGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, visgrid_name+".real"));
              itsPCFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, pcfgrid_name+".real"));
              itsPSFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, psfgrid_name+".real"));
              itsVisGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, visgrid_name+".imag"));
              itsPCFGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, pcfgrid_name+".imag"));
              itsPSFGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, psfgrid_name+".imag"));
          }
        }
      }
      if (itsRestore) {
        // Only create these if we are restoring, as that is when they get made
          if (itsDoingPreconditioning) {
            if (itsWritePsfImage) {
              itsPSFimageCube.reset(new CubeBuilder<casacore::Float>(itsParset, psf_image_name));
            }
          }
          itsRestoredCube.reset(new CubeBuilder<casacore::Float>(itsParset, restored_image_name));
          // we are only interested to collect statistics for the restored image cube
          itsRestoredStatsAndMask.reset(new askap::utils::StatsAndMask(itsComms,itsRestoredCube->filename(),itsRestoredCube->imageHandler()));
          ASKAPLOG_INFO_STR(logger,"Created StatsAndMask object");
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
    logWeightsInfo();

    return;
  }

  /// What are the plans for the deconvolution?
  ASKAPLOG_DEBUG_STR(logger, "Ascertaining Cleaning Plan");
  const bool writeAtMajorCycle = itsParset.getBool("Images.writeAtMajorCycle", false);
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
      itsStats.logSummary();
      ASKAPLOG_INFO_STR(logger, "Starting to process workunit " << workUnitCount+1 << " of " << workUnits.size());

      int initialChannelWorkUnit = workUnitCount;

      if (!updateDir) {

        // NOTE: this is because if we are mosaicking ON THE FLY. We do
        // not process the first workunit outside the imaging loop.
        // But for "normal" processing the first workunit is processed outside the loops
        // This adds all sorts of complications to the logic BTW.

        initialChannelWorkUnit = workUnitCount+1;
      }

      double frequency=workUnits[workUnitCount].get_channelFrequency();
      const string colName = itsParset.getString("datacolumn", "DATA");


      int localChannel;
      int globalChannel;

      bool usetmpfs = itsParset.getBool("usetmpfs", false);
      bool clearcache = itsParset.getBool("clearcache", false);

      if (usetmpfs) {
        // probably in spectral line mode
        // copy the caching here ...
        cacheWorkUnit(workUnits[workUnitCount]);

        localChannel = 0;

      } else {
        localChannel = workUnits[workUnitCount].get_localChannel();
        if (clearcache) {
            cached_files.push_back(workUnits[workUnitCount].get_dataset());
        }
      }

      double globalFrequency = workUnits[workUnitCount].get_channelFrequency();

      const string ms = workUnits[workUnitCount].get_dataset();
      globalChannel = workUnits[workUnitCount].get_globalChannel();

      // MEMORY_BUFFERS mode opens the MS readonly
      TableDataSource ds(ms, TableDataSource::MEMORY_BUFFERS, colName);

      /// Need to set up the rootImager here
      if (updateDir) {
        itsAdvisor->updateDirectionFromWorkUnit(workUnits[workUnitCount]);
      }
      if (updateDir) {
            // change gridder for initial calcNE in updateDir mode
            LOFAR::ParameterSet tmpParset = itsParset.makeSubset("");
            tmpParset.replace("gridder","SphFunc");
            boost::shared_ptr<CalcCore> tempIm(new CalcCore(tmpParset,itsComms,ds,localChannel,globalFrequency));
            rootImagerPtr = tempIm;
      } else if (!gridder_initialized) {
            boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParset,itsComms,ds,localChannel,globalFrequency));
            rootImagerPtr = tempIm;
            gridder_initialized = true;
      } else {
        boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParset,itsComms,ds,rootImagerPtr->gridder(),localChannel,globalFrequency));
        rootImagerPtr = tempIm;
      }

      CalcCore& rootImager = *rootImagerPtr; // just for the semantics
      //// CalcCore rootImager(itsParsets[workUnitCount], itsComms, ds, localChannel);
      /// set up the image for this channel
      /// this will actually build a full image for the first - it is not actually used tho.
      ///
      ASKAPLOG_INFO_STR(logger, "Initialised imager & gridder");

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
        setupImage(rootImager.params(), frequency, false);
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
          rootINERef.weightState(WEIGHTED);
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
            cacheWorkUnit(workUnits[tempWorkUnitCount]);

            localChannel = 0;

          } else {
            localChannel = workUnits[tempWorkUnitCount].get_localChannel();
            if (clearcache) {
                cached_files.push_back(workUnits[tempWorkUnitCount].get_dataset());
            }
          }
          globalChannel = workUnits[tempWorkUnitCount].get_globalChannel();
          double globalFrequency = workUnits[workUnitCount].get_channelFrequency();

          const string myMs = workUnits[tempWorkUnitCount].get_dataset();
          TableDataSource myDs(myMs, TableDataSource::MEMORY_BUFFERS, colName);
          myDs.configureUVWMachineCache(uvwMachineCacheSize, uvwMachineCacheTolerance);
          try {

            boost::shared_ptr<CalcCore> workingImagerPtr;

            if (updateDir) {
              itsAdvisor->updateDirectionFromWorkUnit(workUnits[tempWorkUnitCount]);
              // in updateDir mode I cannot cache the gridders as they have a tangent point.
              // FIXED: by just having 2 possible working imagers depending on the mode. ... easy really

              boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParset,itsComms,myDs,localChannel,globalFrequency));
              workingImagerPtr = tempIm;
            }
            else {
              boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParset,itsComms,myDs,rootImager.gridder(),localChannel,globalFrequency));
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
            itsStats.logSummary();

            // merge into root image if required.
            // this is required if there is more than one workunit per channel
            // either in time or by beam.

            ASKAPLOG_INFO_STR(logger,"About to merge into rootImager");
            ImagingNormalEquations &workingINERef =
            dynamic_cast<ImagingNormalEquations&>(*workingImager.getNE());
            if (updateDir) {
              workingINERef.weightType(FROM_WEIGHT_IMAGES);
              workingINERef.weightState(WEIGHTED);
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
          ASKAPLOG_INFO_STR(logger, "Reached peak residual of " << abs(peak_residual));
          if (peak_residual < targetPeakResidual) {
            if (peak_residual < 0) {
              ASKAPLOG_WARN_STR(logger, "Clean diverging, did not reach the major cycle threshold of "
                              << targetPeakResidual << " Jy. Stopping.");
            } else {
              ASKAPLOG_INFO_STR(logger, "It is below the major cycle threshold of "
              << targetPeakResidual << " Jy. Stopping.");
            }
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
            itsStats.logSummary();

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
        itsStats.logSummary();


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
        logWeightsInfo();

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
      // before archiving "image.slice" as the model, check if a high-resolution "fullres.slice" has been set up
      if (rootImager.params()->has("fullres.slice")) {
        rootImager.params()->add("model.slice", rootImager.params()->valueT("fullres.slice"));
      }
      else {
        rootImager.params()->add("model.slice", rootImager.params()->valueT("image.slice"));
      }
      ASKAPCHECK(rootImager.params()->has("model.slice"), "Params are missing model.slice parameter");

      if (itsWriteGrids) {
        ASKAPLOG_INFO_STR(logger,"Adding grid.slice");
        casacore::Array<casacore::Complex> garr = rootImager.getGrid();
        casacore::Vector<casacore::Complex> garrVec(garr.reform(IPosition(1,garr.nelements())));
        rootImager.params()->addComplexVector("grid.slice",garrVec);
        ASKAPLOG_INFO_STR(logger,"Adding pcf.slice");
        casacore::Array<casacore::Complex> pcfarr = rootImager.getPCFGrid();
        casacore::Vector<casacore::Complex> pcfVec(pcfarr.reform(IPosition(1,pcfarr.nelements())));
        rootImager.params()->addComplexVector("pcf.slice",pcfVec);
        ASKAPLOG_INFO_STR(logger,"Adding psfgrid.slice");
        casacore::Array<casacore::Complex> psfarr = rootImager.getPSFGrid();
        casacore::Vector<casacore::Complex> psfVec(psfarr.reform(IPosition(1,psfarr.nelements())));
        rootImager.params()->addComplexVector("psfgrid.slice",psfVec);
      }

      rootImager.check();


      if (itsRestore) {
        ASKAPLOG_INFO_STR(logger, "Running restore");
        rootImager.restoreImage();
      }

      if (usetmpfs) {
        ASKAPLOG_INFO_STR(logger, "clearing cache");
        clearWorkUnitCache();
        ASKAPLOG_INFO_STR(logger, "done clearing cache");

      } else if (clearcache) {
        // clear the hypercube caches (with 1 channel tiles we won't use it again)
        static int count = 0;
        for (string fileName : cached_files) {
            ROTiledStManAccessor tsm(Table(fileName),colName,True);
            ASKAPLOG_INFO_STR(logger, "Clearing Table cache for " <<colName<< " column");
            tsm.clearCaches();
            // Not sure we should clear the FLAG cache everytime, flags are normally stored in tile with 8 channels
            if (count == 16) {
                ROTiledStManAccessor tsm2(Table(fileName),"FLAG",True);
                ASKAPLOG_INFO_STR(logger, "Clearing Table cache for FLAG column");
                tsm2.clearCaches();
            }
        }
        cached_files.clear();
        if (++count > 16) {
            count = 0;
        }
      }

      itsStats.logSummary();

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

        blankParams.reset(new Params(true));
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

        blankParams.reset(new Params(true));
        ASKAPCHECK(blankParams, "blank parameters (images) not initialised");

        setupImage(blankParams, workUnits[goodUnitCount].get_channelFrequency());


        ContinuumWorkRequest result;
        result.set_params(blankParams);
        result.set_globalChannel(workUnits[goodUnitCount].get_globalChannel());
        /// send the work to the writer with a blocking send
        result.sendRequest(workUnits[goodUnitCount].get_writer(), itsComms);
        ASKAPLOG_INFO_STR(logger, "Sent\n");
      }
      // No need to increment workunit. Although this assumes that we are here because we failed the solveNE not the calcNE

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
  logWeightsInfo();

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

  // Write image
  if (itsImageCube) {
      if (!params->has("model.slice")) {
        ASKAPLOG_WARN_STR(logger, "Params are missing model parameter");
      }
      else {
      ASKAPLOG_INFO_STR(logger, "Writing model for (local) channel " << chan);
      if (params->has("fullres.slice")) {
        // model has been set at a higher resolution
        // If the restored image has full resolution, so will the model. Avoid further oversampling
        ASKAPDEBUGASSERT(params->shape("model.slice").isEqual(params->shape("fullres.slice")));
        itsImageCube->writeRigidSlice(params->valueF("model.slice"), chan);
      }
      else {
        itsImageCube->writeFlexibleSlice(params->valueF("model.slice"), chan);
      }
    }
  }

  // Write PSF
  if (itsPSFCube) {
      if (!params->has("psf.slice")) {
        ASKAPLOG_WARN_STR(logger,  "Params are missing psf parameter");
      }
      else {
          if (!itsWriteGrids) {
            itsPSFCube->writeFlexibleSlice(params->valueF("psf.slice"), chan);
          }
      }
  }

  // Write residual
  if (itsResidualCube) {
      if (!params->has("residual.slice")) {
          ASKAPLOG_WARN_STR(logger,  "Params are missing residual parameter");
      }
      else {
          ASKAPLOG_INFO_STR(logger, "Writing Residual");
          itsResidualCube->writeFlexibleSlice(params->valueF("residual.slice"), chan);
          itsResidualStatsAndMask->calculate(itsResidualCube->filename(),chan,params->valueF("residual.slice"));
      }
  }

  // Write weights
  if (!params->has("weights.slice")) {
      ASKAPLOG_WARN_STR(logger, "Params are missing weights parameter");
  } else {
      if (itsWeightsCube) {
          ASKAPLOG_INFO_STR(logger, "Writing Weights");
          itsWeightsCube->writeFlexibleSlice(params->valueF("weights.slice"), chan);
      } else {
          Array<float> wts = params->valueF("weights.slice");
          float wt = wts.data()[0];
          if (allEQ(wts,wt)) {
            recordWeight(wt, chan);
            ASKAPLOG_INFO_STR(logger,"Writing Weights " << (itsWriteWtLog ? "log" : "extension"));
          } else {
            ASKAPLOG_WARN_STR(logger,"Weights are not identical across image, disabling weights "<< (itsWriteWtLog ? "log" : "extension"));
            recordWeight(-1.0, chan);
          }
      }
  }

  // Write the grids
  if (params->has("grid.slice") && (itsVisGridCube||itsVisGridCubeReal)) {
    const casacore::Vector<casacore::Complex> gr(params->complexVectorValue("grid.slice"));
    casacore::Array<casacore::Complex> grid(gr.reform(params->shape("psf.slice")));
    if (itsGridFFT) {
      ASKAPLOG_INFO_STR(logger, "FFTing Vis Grid and writing it as a real image");
      askap::scimath::fft2d(grid,false);
      itsVisGridCubeReal->writeRigidSlice(casacore::real(grid),chan);
    } else {
      if (itsGridType == "casa") {
        ASKAPLOG_INFO_STR(logger, "Writing Vis Grid");
        itsVisGridCube->writeRigidSlice(grid,chan);
      } else {
        ASKAPLOG_INFO_STR(logger, "Writing Vis Grid as real & imag FITS images");
        itsVisGridCubeReal->writeRigidSlice(casacore::real(grid),chan);
        itsVisGridCubeImag->writeRigidSlice(casacore::imag(grid),chan);
      }
    }
  }
  if (params->has("pcf.slice") && (itsPCFGridCube||itsPCFGridCubeReal)) {
    const casacore::Vector<casacore::Complex> gr(params->complexVectorValue("pcf.slice"));
    casacore::Array<casacore::Complex> grid(gr.reform(params->shape("psf.slice")));
    if (itsGridFFT) {
      ASKAPLOG_INFO_STR(logger, "FFTing PCF Grid and writing it as a real image");
      askap::scimath::fft2d(grid,false);
      itsPCFGridCubeReal->writeRigidSlice(casacore::real(grid),chan);

    } else {
      if (itsGridType == "casa") {
        ASKAPLOG_INFO_STR(logger, "Writing PCF Grid");
        itsPCFGridCube->writeRigidSlice(grid,chan);
      } else {
        ASKAPLOG_INFO_STR(logger, "Writing PCF Grid as real & imag FITS images");
        itsPCFGridCubeReal->writeRigidSlice(casacore::real(grid),chan);
        itsPCFGridCubeImag->writeRigidSlice(casacore::imag(grid),chan);
      }
    }
  }
  if (params->has("psfgrid.slice") && (itsPSFGridCube||itsPSFGridCubeReal)) {
    const casacore::Vector<casacore::Complex> gr(params->complexVectorValue("psfgrid.slice"));
    casacore::Array<casacore::Complex> grid(gr.reform(params->shape("psf.slice")));
    if (itsGridFFT) {
      ASKAPLOG_INFO_STR(logger, "FFTing PSF Grid and writing it as a real image");
      askap::scimath::fft2d(grid,false);
      itsPSFGridCubeReal->writeRigidSlice(casacore::real(grid),chan);
    } else {
      if (itsGridType == "casa") {
        ASKAPLOG_INFO_STR(logger, "Writing PSF Grid");
        itsPSFGridCube->writeRigidSlice(grid,chan);
      } else {
        ASKAPLOG_INFO_STR(logger, "Writing PSF Grid as real & imag FITS images");
        itsPSFGridCubeReal->writeRigidSlice(casacore::real(grid),chan);
        itsPSFGridCubeImag->writeRigidSlice(casacore::imag(grid),chan);
      }
    }
  }
  if (params->has("psf.raw.slice") && itsPSFCube) {
    if (itsWriteGrids) {
      ASKAPLOG_INFO_STR(logger, "Writing un-normalised PSF");
      itsPSFCube->writeFlexibleSlice(params->valueF("psf.raw.slice"), chan);
    }
  }

  // Restored images
  if (itsRestore) {
    if (itsDoingPreconditioning) {
      // Write preconditioned PSF image
      if (itsPSFimageCube) {
          ASKAPCHECK(params->has("psf.image.slice"), "Params are missing psf.image parameter");
          ASKAPLOG_INFO_STR(logger, "Writing preconditioned PSF");
          itsPSFimageCube->writeFlexibleSlice(params->valueF("psf.image.slice"), chan);
      }
    }

    // Record the restoring beam
    const askap::scimath::Axes &axes = params->axes("image.slice");
    recordBeam(axes, chan);

    // Write Restored image
    if (itsRestoredCube) {
        ASKAPLOG_INFO_STR(logger, "Writing Restored Image");
        if (params->has("fullres.slice")) {
          // Restored image has been generated at full resolution, so avoid further oversampling
          ASKAPLOG_INFO_STR(logger, "Writing fullres.slice");
          itsRestoredCube->writeRigidSlice(params->valueF("fullres.slice"), chan);
          calculateImageStats(itsRestoredStatsAndMask,itsRestoredCube,chan,params->valueF("fullres.slice"));
        }
        else {
          ASKAPCHECK(params->has("image.slice"), "Params are missing image parameter");
          ASKAPLOG_INFO_STR(logger, "Writing image.slice");
          itsRestoredCube->writeFlexibleSlice(params->valueF("image.slice"), chan);
          itsRestoredStatsAndMask->calculate(itsRestoredCube->filename(),chan,params->valueF("image.slice"));
        }
    }

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

void ContinuumWorker::initialiseWeightsLog(const unsigned int numChannels)
{
  for(unsigned int i=0;i<numChannels;i++) {
      itsWeightsList[i] = 0.0;
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

void ContinuumWorker::recordWeight(float wt, const unsigned int cubeChannel)
{
  itsWeightsList[cubeChannel] = wt;
}

// void ContinuumWorker::storeBeam(const unsigned int cubeChannel)
// {
//   if (cubeChannel == itsBeamReferenceChannel) {
//     itsRestoredCube->addBeam(itsBeamList[cubeChannel]);
//   }
// }

void ContinuumWorker::logBeamInfo()
{
    bool beamlogAsFile = itsParset.getBool("write.beamlog",true);
    askap::accessors::BeamLogger beamlog;
    if (beamlogAsFile) {
      ASKAPLOG_INFO_STR(logger, "Channel-dependent restoring beams will be written to log file " << beamlog.filename());
    } else if (itsRestoredCube) {
      ASKAPLOG_INFO_STR(logger, "Channel-dependent restoring beams will be written to image " << itsRestoredCube->filename());
    }
    ASKAPLOG_DEBUG_STR(logger, "About to add beam list of size " << itsBeamList.size() << " to the beam logger");
    beamlog.beamlist() = itsBeamList;

    if (itsNumWriters > 1 && itsParset.getBool("singleoutputfile",false)) {
      std::list<int> creators = itsComms.getCubeCreators();
      ASKAPASSERT(creators.size() == 1);
      int creatorRank = creators.front();
      ASKAPLOG_DEBUG_STR(logger, "Gathering all beam information, beam creator is rank " << creatorRank);
      beamlog.gather(itsComms, creatorRank,false);
    }
    if (itsComms.isCubeCreator()) {
        if (itsRestoredCube) {
            if (beamlogAsFile) {
              ASKAPLOG_DEBUG_STR(logger, "Writing list of individual channel beams to beam log");
              beamlog.setFilename("beamlog." + itsRestoredCube->filename() + ".txt");
              beamlog.write();
            } else {
              ASKAPLOG_DEBUG_STR(logger, "Writing list of individual channel beams to image file");
              itsRestoredCube->addBeamList(beamlog.beamlist());
            }

            if (beamlogAsFile || itsParset.getString("imagetype") == "fits") {
              // can't write ref beam to casa image if per channel beams are stored
              ASKAPLOG_DEBUG_STR(logger, "Writing reference restoring beam to header of restored cube");
              casa::Vector<casa::Quantum<double> > refbeam = beamlog.beam(itsBeamReferenceChannel);
              itsRestoredCube->addBeam(refbeam);
            }
        }
    }

}

void ContinuumWorker::logWeightsInfo()
{

  if (!itsWriteWtImage) {
    const string wtLogExt = (itsWriteWtLog ? "log file" : "extension");
    askap::accessors::WeightsLog weightslog;
    ASKAPLOG_INFO_STR(logger, "Channel-dependent weights will be written to "<<wtLogExt);
    ASKAPLOG_DEBUG_STR(logger, "About to add weights list of size " << itsWeightsList.size() << " to the weights logger");
    weightslog.weightslist() = itsWeightsList;

    if (itsNumWriters > 1 && itsParset.getBool("singleoutputfile",false)) {
      std::list<int> creators = itsComms.getCubeCreators();
      ASKAPASSERT(creators.size() == 1);
      int creatorRank = creators.front();
      ASKAPLOG_DEBUG_STR(logger, "Gathering all weights information, creator is rank " << creatorRank);
      weightslog.gather(itsComms, creatorRank,false);
    }
    if (itsComms.isCubeCreator()) {

        // First check weightslog is valid
        for(const auto & wt : itsWeightsList) {
            if (wt.second < 0) {
                ASKAPLOG_WARN_STR(logger, "Weights log invalid - not writing out the channel weights "<<wtLogExt);
                return;
            }
        }
        if (itsWriteWtLog) {
          weightslog.setFilename(itsWeightsName + ".txt");
          ASKAPLOG_INFO_STR(logger, "Writing list of individual channel weights to weights log "
              << weightslog.filename());
          weightslog.write();
        } else {
          ASKAPLOG_INFO_STR(logger, "Writing list of individual channel weights to image extension");
          casacore::Record wtInfo = weightslog.toRecord();
          if (itsRestoredCube) {
            itsRestoredCube->setInfo(wtInfo);
          }
          if (itsResidualCube) {
            itsResidualCube->setInfo(wtInfo);
          }
          // TODO: which other images should have the weights?
        }
    }
  }

}


void ContinuumWorker::setupImage(const askap::scimath::Params::ShPtr& params,
                                 double channelFrequency, bool shapeOverride)
{
  try {
    const LOFAR::ParameterSet imParset = itsParset.makeSubset("Images.");
    ASKAPLOG_DEBUG_STR(logger, "Setting up image");

    const int nfacets = imParset.getInt32("nfacets", 1);
    const string name("image.slice");
    vector<string> direction = imParset.getStringVector("direction");

    const vector<string> cellsize = imParset.getStringVector("cellsize");
    vector<int> shape = imParset.getInt32Vector("shape");
    const int nchan = 1;

    if (shapeOverride == true) {
      string param = "subshape";
      if (imParset.isDefined(param)) {
        ASKAPLOG_INFO_STR(logger,"Over-riding image shape from parset");
        shape = imParset.getInt32Vector("subshape");
        ASKAPLOG_INFO_STR(logger,"Image shape now " << shape);
      }
      else {
        ASKAPLOG_WARN_STR(logger,"Shape over-ride requested but no subshape parameter in parset");
      }
    } else if (itsParset.getBool("updatedirection",false)) {
          // override with image specific direction if present - for mosaic case - combined image direction
          vector<string> names = imParset.getStringVector("Names",{},false);
          if (names.size()>0) {
              if (imParset.isDefined(names[0]+".direction")) {
                  ASKAPLOG_INFO_STR(logger,"Using image direction from parset instead of tangent point from advise");
                  direction = imParset.getStringVector(names[0]+".direction");
              }
          }
    }


    if (!imParset.isDefined("polarisation")) {
      ASKAPLOG_DEBUG_STR(logger, "Polarisation frame is not defined, "
      << "only stokes I will be generated");
    }
    const vector<string> stokesVec = imParset.getStringVector("polarisation",
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

    const bool ewProj = imParset.getBool("ewprojection", false);
    if (ewProj) {
      ASKAPLOG_DEBUG_STR(logger, "Image will have SCP/NCP projection");
    } else {
      ASKAPLOG_DEBUG_STR(logger, "Image will have plain SIN projection");
    }

    ASKAPCHECK(nfacets > 0, "Number of facets is supposed to be a positive number, you gave " << nfacets);
    ASKAPCHECK(shape.size() >= 2, "Image is supposed to be at least two dimensional. " << "check shape parameter, you gave " << shape);

    ASKAPLOG_DEBUG_STR(logger,"setupImage : direction = "<<direction<< " shape = "<< shape);

    if (nfacets == 1) {
      SynthesisParamsHelper::add(*params, name, direction, cellsize, shape, ewProj,
        channelFrequency, channelFrequency, nchan, stokes);
        // SynthesisParamsHelper::add(*params, name, direction, cellsize, shape, ewProj,
        //                            freq[0], freq[1], nchan, stokes);
    } else {
        // this is a multi-facet case
        const int facetstep = imParset.getInt32("facetstep", casacore::min(shape[0], shape[1]));
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

void ContinuumWorker::calculateImageStats(boost::shared_ptr<askap::utils::StatsAndMask> statsAndMask,
                                          boost::shared_ptr<CubeBuilder<casacore::Float> > imgCube,
                                          int channel, const casacore::Array<float>& arr)
{
    boost::optional<float> oversamplingFactor = imgCube->oversamplingFactor();
    if ( oversamplingFactor ) {
        // Image param is stored at a lower resolution, so increase to desired resolution before writing
        casacore::Array<float> fullresarr(scimath::PaddingUtils::paddedShape(arr.shape(),*oversamplingFactor));
        scimath::PaddingUtils::fftPad(arr,fullresarr);
        statsAndMask->calculate(imgCube->filename(),channel,fullresarr);
    } else {
        statsAndMask->calculate(imgCube->filename(),channel,arr);
    }
}

void ContinuumWorker::writeCubeStatistics()
{
  const std::string imageType = itsParset.getString("imagetype","casa");
  const std::string statsFile = itsParset.getString("outputStats","");

  // get one of the ranks that create the cube
  std::list<int> creators = itsComms.getCubeCreators();
  creators.sort();

  if ( creators.size() > 0 ) {
    // make the first rank in the cube creator to be the receiver of the cube statistics
    // and the rest of the workers send stats to it
    // NOTE: the compiler wont let me do this :
    // unsigned int statsCollectorRank; //  = reinterpret_cast<unsigned int> (creators.front());
    int temp = creators.front();
    unsigned int statsCollectorRank = static_cast<unsigned int> (temp);
    ASKAPLOG_INFO_STR(logger, "statsCollectorRank = " << statsCollectorRank);
    if ( itsRestore ) {
      if ( itsComms.rank() == statsCollectorRank ) {
        if ( itsRestoredStatsAndMask ) {
          std::string fullFilename = itsRestoredCube->filename();
          if ( itsNumWriters > 1 ) {
            // if there are more than one writers, then the statistics is distributed over more than one
            // ranks so we designate one of the ranks (statsCollectorRank) to be the collector of the
            // stats and the other ranks send their stats to it
            ASKAPLOG_INFO_STR(logger, "Collecting cube statistics with itsNumWriters = " << itsNumWriters);
            std::set<unsigned int> excludedRanks {0,statsCollectorRank};
            itsRestoredStatsAndMask->receiveStats(excludedRanks);
          }
          ASKAPLOG_INFO_STR(logger,"Writing statistic to restored image: " << fullFilename);
          if ( statsFile != "" ) {
            const std::string restoredStatsFile = std::string("Restored_") + statsFile;
            itsRestoredStatsAndMask->writeStatsToFile(restoredStatsFile);
          }
          itsRestoredStatsAndMask->writeStatsToImageTable(fullFilename);
        } else {
            ASKAPLOG_INFO_STR(logger, "itsRestoredStatsAndMask of statsCollectorRank: " << statsCollectorRank << " is null");
        }
      } else {
        // not all the ranks that send the stats to the receiver rank have stats to send (i.e
        // their StatsAndMask is null) so we create a dummy stats for them and force them to send
        // null (0) stats to the receiver.
        if ( itsRestoredStatsAndMask ) {
          ASKAPLOG_INFO_STR(logger, "Rank: " << itsComms.rank() << " sends stats to rank " << statsCollectorRank);
          itsRestoredStatsAndMask->sendStats(statsCollectorRank);
        } else {
          askap::utils::StatsAndMask dummy {itsComms};
          ASKAPLOG_INFO_STR(logger, "Rank: " << itsComms.rank() << " sends dummy stats to rank " << statsCollectorRank);
          dummy.sendStats(statsCollectorRank);
        }
      }
    }
    ASKAPLOG_INFO_STR(logger,"Waiting for all ranks to finish");
    itsComms.barrier(itsComms.theWorkers());
    // write the stats to the residual cube
    if ( itsWriteResidual ) {
      if ( itsComms.rank() == statsCollectorRank ) {
        if ( itsResidualStatsAndMask ) {
          std::string fullFilename = itsResidualCube->filename();
          if ( itsNumWriters > 1 ) {
            // if there are more than one writers, then the statistics is distributed over more than one
            // ranks so we designate one of the ranks (statsCollectorRank) to be the collector of the
            // stats and the other ranks send their stats to it
            ASKAPLOG_INFO_STR(logger, "Collecting cube statistics with itsNumWriters = " << itsNumWriters);
            std::set<unsigned int> excludedRanks {0,statsCollectorRank};
            itsResidualStatsAndMask->receiveStats(excludedRanks);
          }
          if ( statsFile != "" ) {
            const std::string residualStatsFile = std::string("Residual_") + statsFile;
            itsResidualStatsAndMask->writeStatsToFile(residualStatsFile);
          }
          ASKAPLOG_INFO_STR(logger,"Writing statistic to residual image: " << fullFilename);
          itsResidualStatsAndMask->writeStatsToImageTable(fullFilename);
        }
      } else {
        // not all the ranks that send the stats to the receiver rank have stats to send (i.e
        // their StatsAndMask is null) so we create a dummy stats for them and force them to send
        // null (0) stats to the receiver.
        if ( itsResidualStatsAndMask ) {
          ASKAPLOG_INFO_STR(logger, "Rank: " << itsComms.rank() << " sends stats to rank " << statsCollectorRank);
          itsResidualStatsAndMask->sendStats(statsCollectorRank);
        } else {
          askap::utils::StatsAndMask dummy {itsComms};
          ASKAPLOG_INFO_STR(logger, "Rank: " << itsComms.rank() << " sends dummy stats to rank " << statsCollectorRank);
          dummy.sendStats(statsCollectorRank);
        }
      }
    }
  } else {
    ASKAPLOG_INFO_STR(logger,"creators.size() < 0");
  }
}
