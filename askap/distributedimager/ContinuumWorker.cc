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
/// Refactoring by Max Voronkov

// Include own header file first
#include "ContinuumWorker.h"

// System includes
#include <string>
#include <stdexcept>
#include <vector>

// boost includes
#include "boost/shared_ptr.hpp"

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/profile/AskapProfiler.h>
#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/ImageParamsHelper.h>
#include <askap/scimath/utils/PolConverter.h>
#include <Common/Exceptions.h>
#include <casacore/casa/OS/Timer.h>
#include <askap/parallel/ImagerParallel.h>
#include <askap/imageaccess/BeamLogger.h>
#include <askap/imageaccess/WeightsLog.h>
#include <askap/gridding/UVWeightParamsHelper.h>

// Local includes
#include "askap/distributedimager/AdviseDI.h"
#include "askap/messages/ContinuumWorkRequest.h"

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
  : itsParset(parset), itsComms(comms), itsStats(stats),
    itsDoingPreconditioning(doingPreconditioning(parset)),
    // setup whether we solve locally (spectral line mode) or on the master (continuum mode)
    itsLocalSolver(parset.getBool("solverpercore", false)),
    // read an initial model cube (standard naming), spectral line case
    itsReadStartingModelCube(parset.getBool("Images.reuse",false)),
    // setup whether we restore and write restored image
    itsRestore(parset.getBool("restore", false)),
    // setup whether to write the residual image + support of an alternative parameter name
    itsWriteResidual(parset.getBool("write.residualimage", parset.getBool("residuals",false))),
    // write unnormalised, natural wt psf
    itsWritePsfRaw(parset.getBool("write.psfrawimage", false)),
    // write normalised, preconditioned psf
    itsWritePsfImage(parset.getBool("write.psfimage", true)),
    // write weights image
    itsWriteWtImage(parset.getBool("write.weightsimage", false)),
    // write weights log file
    itsWriteWtLog(parset.getBool("write.weightslog", false)),
    // clean model (itsRestore is already initialised by this point and can be used for default value)
    // if we're reading in the model cube, we also write it out again
    itsWriteModelImage(parset.getBool("write.modelimage", !itsRestore)||itsReadStartingModelCube),
    // write (dump) the gridded data, psf and pcf + support of an alternative parameter name
    itsWriteGrids(parset.getBool("write.grids", parset.getBool("dumpgrids", false))),
    // grid image type
    itsGridType(parset.getString("imagetype","casa")),
    // label grid with UV coordinates
    itsGridCoordUV(parset.getBool("write.grids.uvcoord", itsGridType=="casa")),
    // write fft of grid (i.e. dirty image, psf)
    itsGridFFT(parset.getBool("write.grids.fft",false)),
    // number of cube writers (itsGridType and itsParset are already defined and could be used)
    itsNumWriters(configureNumberOfWriters()),
    // joint gridding of multiple beams/directions
    itsUpdateDir(parset.getBool("updatedirection",false)),
    // flag to mask mosaic output with weights image & tolerance
    itsMaskOutput(parset.getBool("maskmosaic",true) && itsUpdateDir),
    // masking threshold for mosaic output
    itsMaskLevel(parset.getFloat("solver.Clean.tolerance",0.1)),
    // use MFS starting model (for spectral mode)
    itsMFSStartingModel(parset.getBool("mfsstartingmodel",false)),
    // flag that we do traditional weighting (note, this is somewhat ugly to setup calculators only to check whether 
    // the shared pointer is not empty. But this is cheap. We can clear this up later)
    // uv-weight calculator object; at the moment, it is an empty pointer if no traditional weighting is done
    itsUVWeightCalculator(ImagerParallel::createUVWeightCalculator(parset))
{
    ASKAPTRACE("ContinuumWorker::constructor");

    ASKAPCHECK(!(itsUpdateDir && !itsLocalSolver), "Cannot <yet> Continuum image in on-the-fly mosaick mode - need to update the image parameter setup");

    itsAdvisor = boost::shared_ptr<synthesis::AdviseDI> (new synthesis::AdviseDI(itsComms, itsParset));
    itsAdvisor->prepare();

    // lets properly size the storage
    const int nchanpercore = itsParset.getInt32("nchanpercore", 1);

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

    itsBaseChannel = posInGroup * nchanpercore;

    ASKAPLOG_INFO_STR(logger, "Distribution: Id " << id << " nWorkers " << nWorkers << " nGroups " << itsComms.nGroups());

    ASKAPLOG_INFO_STR(logger, "Distribution: Base channel " << itsBaseChannel << " PosInGrp " << posInGroup);

    ASKAPCHECK(!itsParset.getBool("usetmpfs", false), "usetmpfs option is no longer supported by the code");

    const bool dopplerTracking = itsParset.getBool("dopplertracking",false);
    if (dopplerTracking) {
        std::vector<string> direction = itsParset.getStringVector("dopplertracking.direction",{},false);
        if (direction.size() == 3) {
            casacore::MDirection itsVelRefDir(asMDirection(direction));
            ASKAPLOG_INFO_STR(logger,"Velocity reference direction =  "<<printDirection(itsVelRefDir.getValue()));
        } else {
            ASKAPLOG_INFO_STR(logger,"Using default velocity reference direction ( = pointing/field centre)");
        }
    }

    if (isSampleDensityGridNeeded() && itsParset.getBool("sphfuncforpsf", false) && !itsUpdateDir) {
        ASKAPLOG_WARN_STR(logger, "Using together sphfuncforpsf = true and the traditional weighting mode may lead to unexpected results. Consider using updatedirection=true. See AXA-2792 for details.");
    }

    // check if we can load in the model cube
    if (itsLocalSolver && itsReadStartingModelCube) {
      const int nWriters = itsParset.getInt32("nwriters",1);
      const string imageType = itsParset.getString("imagetype","casa");
      ASKAPCHECK(imageType=="fits" && nWriters == nWorkers,"The Images.reuse option in spectral mode requires imagetype=fits and nwriters==nworkers");
    }
}

/// @brief figure out if preconditioning is to be done
/// @details This method encapsulates checks of the parset indicating that preconditioning is going to be done.
/// It is necessary to configure writing of additional data products, although preconditioning itself
/// is enabled and done by the appropriate solver class.
/// @param[in] parset parset to use (this method is expected to be used in the constructor, so it is handy not
/// to rely on the itsParset data field).
/// @return true if preconditioning is to be done, false otherwise
bool ContinuumWorker::doingPreconditioning(const LOFAR::ParameterSet &parset)
{
   const vector<string> preconditioners = parset.getStringVector("preconditioner.Names", std::vector<std::string>());
   for (vector<string>::const_iterator pc = preconditioners.begin(); pc != preconditioners.end(); ++pc) {
        if ((*pc) == "Wiener" || (*pc) == "NormWiener" || (*pc) == "Robust" || (*pc) == "GaussianTaper") {
            return true;
        }
   }
   return false;
}

/// @brief helper method to obtain the number of writers for the cube
/// @details This method obtains the number of writers from the parset and adjusts it if necessary.
/// It is intended to be used in the constructor to fill itsNumWriters data field and requires
/// itsGridType and itsParset to be valid.
/// @return number of writer ranks for grid export
int ContinuumWorker::configureNumberOfWriters()
{
   const int nwriters = itsParset.getInt32("nwriters",1);
   ASKAPCHECK(nwriters>0,"Number of writers must be greater than 0");
   if (itsGridType == "casa" && itsParset.getBool("singleoutputfile",false) && nwriters > 1){
       ASKAPLOG_WARN_STR(logger,"Reducing number of writers to 1 because we are writing a single casa image cube");
       return 1;
   }
   return nwriters;
}

void ContinuumWorker::run(void)
{
  ASKAPTRACE("ContinuumWorker::run");

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
      // MV: in the original code before refactoring NA allocations would be processed and skipped later
      // in principle, this allows us to do other actions (e.g. related to writing). However, I am not sure it
      // ever worked correctly before the refactoring (would've written a blank image into wrong channel). After the
      // refactoring it triggered some asserts and, therefore, the decision was made to skip such work units here.
      // If such work units are found necessary later, the logic needs to be altered in processChannels to deal with them.
      // And probably the logic of writing needs to be cleared up at some point (there is still some technical debt which
      // is causing issues here)
      ASKAPLOG_WARN_STR(logger, "Worker has received non applicable allocation - ignoring");
    } else {

      ASKAPLOG_DEBUG_STR(logger, "Worker has received valid allocation");

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
    }

    wrequest.sendRequest(itsMaster, itsComms);

  } // while (1) // break when "DONE"

  // MV: see AXA-3083, we still find some corner cases. For now disable this check but it can cause issues somewhere else
  // (probably need refactoring / rewrite of the logic of the part generating work units and doing the write - it was untouched
  // during major ContinuumWorker refactoring)
  //ASKAPCHECK(itsWorkUnits.size() > 0, "No work at to do - something has broken in the setup");

  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " received data from master - waiting at barrier");
  itsComms.barrier(itsComms.theWorkers());
  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " passed barrier");

  configureChannelAllocation();

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
  const bool calcstats = itsParset.getBool("calcstats", true);
  if ( singleoutputfile && calcstats ) {
    writeCubeStatistics();
    itsComms.barrier(itsComms.theWorkers());
  }
  ASKAPLOG_INFO_STR(logger, "Rank " << itsComms.rank() << " passed final barrier");
}

/// @brief configure allocation in channels
/// @details This method sets up channel allocation for cube writing (in local solver mode)
/// or combines channels in the global solver mode if configured in the parset.
void ContinuumWorker::configureChannelAllocation()
{
   const int nchanpercore = itsParset.getInt("nchanpercore", 1);

   const int nWorkers = itsComms.nProcs() - 1;
   const int nGroups = itsComms.nGroups();
   const int nchanTotal = nWorkers * nchanpercore / nGroups;
   ASKAPDEBUGASSERT(itsAdvisor);


   if (itsLocalSolver) {
       ASKAPLOG_INFO_STR(logger, "In local solver mode - reprocessing allocations)");
       itsAdvisor->updateComms();
       itsComms.buildWriterIndex(itsComms.theWorkers());
       int myMinClient = itsComms.rank();
       int myMaxClient = itsComms.rank();

       if (itsComms.isWriter()) {
           ASKAPLOG_DEBUG_STR(logger, "Getting client list for cube generation");
           // MV - probably should've used std::set here (no need to sort + unique by default)
           std::list<int> myClients = itsComms.getClients();
           myClients.push_back(itsComms.rank());
           myClients.sort();
           myClients.unique();

           ASKAPLOG_DEBUG_STR(logger, "Client list " << myClients);
           if (myClients.size() > 0) {
               typedef std::list<int>::const_iterator IterType;
               const std::pair<IterType, IterType> extrema = std::minmax_element(myClients.begin(), myClients.end());

               myMinClient = *(extrema.first);
               myMaxClient = *(extrema.second);
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
           if (!itsComms.isSingleSink() && (itsGridType != "adios")) {
               ASKAPLOG_INFO_STR(logger, "MultiCube with multiple writers");
               itsNChanCube = (myMaxClient - myMinClient + 1) * nchanpercore;
               itsBaseCubeGlobalChannel = (myMinClient - 1) * nchanpercore;
               itsBaseCubeFrequency = itsAdvisor->getBaseFrequencyAllocation((myMinClient - 1));
           } else {
               ASKAPLOG_INFO_STR(logger, "SingleCube with multiple writers");
               itsNChanCube = nchanTotal;
               itsBaseCubeGlobalChannel = 0;
               itsBaseCubeFrequency = itsAdvisor->getBaseFrequencyAllocation((0));
           }
           initialiseBeamLog(itsNChanCube);
           initialiseWeightsLog(itsNChanCube);

           ASKAPLOG_INFO_STR(logger, "Number of channels in cube is: " << itsNChanCube);
           ASKAPLOG_INFO_STR(logger, "Base global channel of cube is " << itsBaseCubeGlobalChannel);
       }
       itsBaseFrequency = itsAdvisor->getBaseFrequencyAllocation(itsComms.rank() - 1);
   } else {
       const bool combineChannels = itsParset.getBool("combinechannels", false);
       if (combineChannels) {
           ASKAPLOG_INFO_STR(logger, "Not in localsolver (spectral line) mode - and combine channels is set so compressing channel allocations");
           compressWorkUnits();
       }
       initialiseBeamLog(nchanTotal);
       initialiseWeightsLog(nchanTotal);
       itsComms.buildWriterIndex(itsComms.theWorkers());
   }
}

void ContinuumWorker::compressWorkUnits() {

   itsWorkUnits.mergeAdjacentChannels();
   ASKAPCHECK(itsWorkUnits.size() < 2, "The number of compressed workunits is greater than one. Channel parameters may be incorrect - see AXA-1004 and associated technical debt tickets");
   if (itsWorkUnits.size() > 0) {
       // Now need to update the parset details
       const cp::ContinuumWorkUnit& wu = *itsWorkUnits.begin();
       const std::string ChannelParam = "["+toString(wu.get_nchan())+","+
                    toString(wu.get_localChannel())+"]";
       ASKAPLOG_DEBUG_STR(logger, "compressWorkUnit: ChannelParam = "<<ChannelParam);
       itsParset.replace("Channels",ChannelParam);
   }
}

void ContinuumWorker::preProcessWorkUnit(ContinuumWorkUnit& wu)
{
  // This also needs to set the frequencies and directions for all the images
  ASKAPLOG_DEBUG_STR(logger, "In preProcessWorkUnit");
  ASKAPLOG_DEBUG_STR(logger, "Parset Reports: (In preProcess workunit)" << (itsParset.getStringVector("dataset", true)));

  ASKAPLOG_DEBUG_STR(logger, "Getting advice on missing parameters");

  itsAdvisor->addMissingParameters();

  ASKAPLOG_DEBUG_STR(logger, "Storing workUnit");
  itsWorkUnits.add(wu); //
  ASKAPLOG_DEBUG_STR(logger, "Finished preProcessWorkUnit");
  ASKAPLOG_DEBUG_STR(logger, "Parset Reports (leaving preProcessWorkUnit): " << (itsParset.getStringVector("dataset", true)));
}

/// @brief configure reference channel used for the restoring beam
/// @details This method populates itsBeamReferenceChannel based on the parset and
/// the number of channels allocated to the cube handled by this rank
void ContinuumWorker::configureReferenceChannel()
{
   // Define reference channel for giving restoring beam
   std::string reference = itsParset.getString("restore.beamReference", "mid");
   if (reference == "mid") {
       itsBeamReferenceChannel = itsNChanCube / 2;
   } else if (reference == "first") {
        itsBeamReferenceChannel = 0;
   } else if (reference == "last") {
         itsBeamReferenceChannel = itsNChanCube - 1;
   } else { // interpret reference as a 0-based channel nuumber
         const unsigned int num = utility::fromString<unsigned int>(reference);
         if (num < itsNChanCube) {
             itsBeamReferenceChannel = num;
         } else {
             itsBeamReferenceChannel = itsNChanCube / 2;
             ASKAPLOG_WARN_STR(logger, "beamReference value (" << reference
                   << ") not valid. Using middle value of " << itsBeamReferenceChannel);
         }
   }
}

/// @brief initialise cube writing if in local solver mode
/// @details This method encapsulates the code which handles cube writing in the local solver mode
/// (i.e. when it is done from the worker). Safe to call in continuum mode too (as itsComms.isWriter
/// would return false in this case)
void ContinuumWorker::initialiseCubeWritingIfNecessary()
{
   // MV: I factored out this code with little modification (largely constness and the like), mainly to avoid
   // having a giant processChannels method. I feel that more restructuring can be done to this code
   if (itsComms.isWriter()) {

       // This code is only used in the spectral line/local solver case -
       //   continuum images are written from ImagerParallel::writeModel in ContinuumMaster
       ASKAPDEBUGASSERT(itsLocalSolver);

       const Quantity f0(itsBaseCubeFrequency, "Hz");
       /// The width of a channel. THis does <NOT> take account of the variable width
       /// of Barycentric channels
       const Quantity freqinc(itsWorkUnits[0].get_channelWidth(), "Hz");

       // add rank based postfix if we're writing to multiple cubes
       const std::string postfix = ((itsComms.isSingleSink() || (itsParset.getString("imagetype") == "adios")) ? "" : std::string(".wr.") + utility::toString(itsComms.rank()));

       const std::string img_name = "image" + postfix;
       const std::string psf_name = "psf" + postfix;
       const std::string residual_name = "residual" + postfix;
       // may need this name for the weightslog
       const std::string weights_name = "weights" + postfix;
       const std::string visgrid_name = "visgrid" + postfix;
       const std::string pcfgrid_name = "pcfgrid" + postfix;
       const std::string psfgrid_name = "psfgrid" + postfix;
       const std::string psf_image_name = "psf.image" + postfix;
       const std::string restored_image_name = "image.restored" + postfix;

       ASKAPLOG_DEBUG_STR(logger, "Configuring Spectral Cube");
       ASKAPLOG_DEBUG_STR(logger, "nchan: " << itsNChanCube << " base f0: " << f0.getValue("MHz")
            << " width: " << freqinc.getValue("MHz") << " (" << itsWorkUnits[0].get_channelWidth() << ")");

       if (itsWriteWtLog) {
           itsWeightsName = CubeBuilder<casacore::Float>::makeImageName(itsParset,weights_name);
       }

       LOFAR::ParameterSet gridParset = itsParset.makeSubset("");
       gridParset.remove("Images.extraoversampling");
       size_t comm_index = -1;

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

            if (itsReadStartingModelCube) {
                itsImageCube.reset(new CubeBuilder<casacore::Float>(itsParset, img_name));             
            } else if (itsWriteModelImage) {
                if ((itsGridType == "adios") && (itsParset.getString("imageaccess", "individual") == "collective")) {
                    if (comm_index == -1) {
                        comm_index = itsComms.theWriters();
                    }
                    itsImageCube.reset(new CubeBuilder<casacore::Float>(itsComms, comm_index, itsParset, itsNChanCube, f0, freqinc, img_name));
                } else {
                    itsImageCube.reset(new CubeBuilder<casacore::Float>(itsParset, itsNChanCube, f0, freqinc, img_name));
                }
            }
            if (itsWritePsfRaw) {
                itsPSFCube.reset(new CubeBuilder<casacore::Float>(itsParset, itsNChanCube, f0, freqinc, psf_name));
            }
            if (itsWriteResidual) {
                itsResidualCube.reset(new CubeBuilder<casacore::Float>(itsParset, itsNChanCube, f0, freqinc, residual_name));
                itsResidualStatsAndMask.reset(new askap::utils::StatsAndMask(itsComms,itsResidualCube->filename(),itsResidualCube->imageHandler()));
                ASKAPLOG_INFO_STR(logger,"Created StatsAndMask object for residual cube");
            }
            if (itsWriteWtImage) {
                itsWeightsCube.reset(new CubeBuilder<casacore::Float>(itsParset, itsNChanCube, f0, freqinc, weights_name));
            }
            if (itsWriteGrids) {
                if (itsGridFFT) {
                    itsVisGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, visgrid_name));
                    itsPCFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, pcfgrid_name));
                    itsPSFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, psfgrid_name));
                } else {
                    if (itsGridType == "casa") {
                        itsVisGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, itsNChanCube, f0, freqinc, visgrid_name, true));
                        itsPCFGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, itsNChanCube, f0, freqinc, pcfgrid_name, true));
                        itsPSFGridCube.reset(new CubeBuilder<casacore::Complex>(gridParset, itsNChanCube, f0, freqinc, psfgrid_name, true));
                    } else if ((itsGridType == "adios") && (itsParset.getString("imageaccess", "individual") == "collective")) {
                        if (comm_index == -1) {
                            comm_index = itsComms.theWriters();
                        }
                        itsVisGridCubeReal.reset(new CubeBuilder<casacore::Float>(itsComms, comm_index, gridParset, itsNChanCube, f0, freqinc, visgrid_name+".real", itsGridCoordUV));
                        itsPCFGridCubeReal.reset(new CubeBuilder<casacore::Float>(itsComms, comm_index, gridParset, itsNChanCube, f0, freqinc, pcfgrid_name+".real", itsGridCoordUV));
                        itsPSFGridCubeReal.reset(new CubeBuilder<casacore::Float>(itsComms, comm_index, gridParset, itsNChanCube, f0, freqinc, psfgrid_name+".real", itsGridCoordUV));
                        itsVisGridCubeImag.reset(new CubeBuilder<casacore::Float>(itsComms, comm_index, gridParset, itsNChanCube, f0, freqinc, visgrid_name+".imag", itsGridCoordUV));
                        itsPCFGridCubeImag.reset(new CubeBuilder<casacore::Float>(itsComms, comm_index, gridParset, itsNChanCube, f0, freqinc, pcfgrid_name+".imag", itsGridCoordUV));
                        itsPSFGridCubeImag.reset(new CubeBuilder<casacore::Float>(itsComms, comm_index, gridParset, itsNChanCube, f0, freqinc, psfgrid_name+".imag", itsGridCoordUV));
                    } else {
                        itsVisGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, visgrid_name+".real", itsGridCoordUV));
                        itsPCFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, pcfgrid_name+".real", itsGridCoordUV));
                        itsPSFGridCubeReal.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, psfgrid_name+".real", itsGridCoordUV));
                        itsVisGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, visgrid_name+".imag", itsGridCoordUV));
                        itsPCFGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, pcfgrid_name+".imag", itsGridCoordUV));
                        itsPSFGridCubeImag.reset(new CubeBuilder<casacore::Float>(gridParset, itsNChanCube, f0, freqinc, psfgrid_name+".imag", itsGridCoordUV));
                    }
                }
            }
            if (itsRestore) {
                // Only create these if we are restoring, as that is when they get made
                if (itsDoingPreconditioning) {
                    if (itsWritePsfImage) {
                        itsPSFimageCube.reset(new CubeBuilder<casacore::Float>(itsParset, itsNChanCube, f0, freqinc, psf_image_name));
                    }
                }
                itsRestoredCube.reset(new CubeBuilder<casacore::Float>(itsParset, itsNChanCube, f0, freqinc, restored_image_name));
                // we are only interested to collect statistics for the restored image cube
                itsRestoredStatsAndMask.reset(new askap::utils::StatsAndMask(itsComms,itsRestoredCube->filename(),itsRestoredCube->imageHandler()));
                ASKAPLOG_INFO_STR(logger,"Created StatsAndMask object for restored cube");
            }

       } else {
            // this is a cube writer rather than creator
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
            }
       }
   }

   // MV - moved code pretty much as it was with only minor polishing during refactoring, but perhaps this barrier could be moved inside
   // the if-statement above as presumably it is not needed in the continuum case
   ASKAPLOG_DEBUG_STR(logger, "You shall not pass. Waiting at a barrier for all ranks to have created the cubes ");
   itsComms.barrier(itsComms.theWorkers());
   ASKAPLOG_DEBUG_STR(logger, "Passed the barrier");
}

/// @brief helper method to create and configure work and (optionally) root imagers
/// @details This method encapsulates the part of single work unit processing where the work and root imagers are created. 
/// Using two imager objects is a bit of the technical debt - ideally, one has to merge normal equations or models directly.
/// But this is deeply in the design of this application and left as is for now. Normally, all gridding of data is taken place
/// in the 'work imager' and the results are merged into 'root imager' when ready. If the root imager is not defined, the work imager
/// becomes one for the subsequent data merge. However, the logic of creating these imagers depends on the mode (e.g. itsUpdateDir, itsLocalSolver).
/// This method encapsulates all the logic, so it can be repeated easily for both normal gridding and sample grid calculation for traditional weighting.
/// @note in the case of central solver (and no joint deconvolution), this method expects to receive model from the master rank if root imager is not
/// defined (i.e. when the work imager created inside this method would become a new root imager for future accumulation)
/// @param[in] wu work unit to work with (parameters like frequency, channel and the dataset may be used). Note, access to data currently happens
/// in the joint imaging mode where the full image is created through dummy iteration hack.
/// @param[inout] rootImagerPtr shared pointer to CalcCore object to use as the root imager (if empty, it is created in the joint imaging mode)
/// @return shared pointer to CalcCore object to be used as a work imager for the given work unit
boost::shared_ptr<CalcCore> ContinuumWorker::createImagers(const cp::ContinuumWorkUnit &wu, boost::shared_ptr<CalcCore> &rootImagerPtr) const
{
   ASKAPASSERT(itsDSM);
   const int localChannel = wu.get_localChannel();
   const double globalFrequency = wu.get_channelFrequency();
   const uInt globalChannel = wu.get_globalChannel();
   TableDataSource& ds = itsDSM->dataSource(wu.get_dataset());

   if (itsUpdateDir) {
       // note, this can update the parset which is then used to construct CalcCore objects
       itsAdvisor->updateDirectionFromWorkUnit(wu);
   }

   // to accumulate data from the current work unit
   boost::shared_ptr<CalcCore> workingImagerPtr;
   if (itsUpdateDir || !rootImagerPtr) {
       // for itsUpdateDir mode gridders cannot be cached as they have a tangent point
       // so use the appropriate CalcCore constructor (and same case if the cache is empty
       // solver is initialised if root imager has not yet been defined and it is not joint deconvolution mode
       // (i.e. if this working imager will eventually become the root imager for subsequent data accumulation)
       workingImagerPtr.reset(new CalcCore(itsParset,itsComms,ds,localChannel,globalFrequency, !rootImagerPtr && !itsUpdateDir));
   } else {
       // in this case we can reuse gridder from the rootImager (created previously)
       // also, inhibit solver initialisation because we only run solver in the root imager
       workingImagerPtr.reset(new CalcCore(itsParset,itsComms,ds,rootImagerPtr->gridder(),localChannel,globalFrequency,false));
   }

   CalcCore& workingImager = *workingImagerPtr; // just for the semantics
   if (itsUpdateDir) {
       // MV: we could've had more intelligent checking whether we need a new subpatch here
       // (i.e. there could be multiple epochs for the same beam). Leave it as it was before
       // the refactoring.
       const bool useSubSizedImages = true;
       setupImage(workingImager.params(), globalFrequency, useSubSizedImages);
       // Note, the original code prior to refactoring had "if (majorCycleNumber > 0)" condition for model copy
       // the following may need to be adjusted if I (MV) didn't understand the logic correctly
       if (rootImagerPtr && rootImagerPtr->params()) {
           copyModel(rootImagerPtr->params(),workingImager.params());
       } else {
           // MV: use the same approach with a dummy iteration as in the code prior to refactoring
           // (this is to setup a full size NE to linmos into. This condition indicates that we have
           // the first work unit processed. On subsequent major cycles rootImager will contain the new model and
           // for subsequent work units copyModel shouldn't do any harm
           // I added the following assert condition as I don't expect rootImager defined without a model. Another
           // wrapping if-statement will be necessary if this assumption doesn't hold
           ASKAPASSERT(!rootImagerPtr);
           // change gridder for initial calcNE in itsUpdateDir mode
           LOFAR::ParameterSet tmpParset = itsParset.makeSubset("");
           tmpParset.replace("gridder","SphFunc");
           // the following will initialise the solver (i.e. default parameter)
           rootImagerPtr.reset(new CalcCore(tmpParset,itsComms,ds,localChannel,globalFrequency));
           // setup full size image
           setupImage(rootImagerPtr->params(), globalFrequency, false);

           if (isSampleDensityGridNeeded()) {
               // note, this effectively emulates the same approach as used in the traditional weighting hack, i.e.
               // it would generate separate weight for the dummy pass done with the full size grid using the root imager
               // (this is only done in the joint deconvolution mode). This approach is fundamentally flawed (see AXA-2792)
               // because it uses PSF built with the first work unit only. But if the data are only split by beams it may be
               // a viable short-cut. So the intention is to have this approach available as an option (the only choice until
               // other approaches are implemented), at least in the short term. Besides, more accurate ways to construct PSF
               // which take into account all work units (and, therefore, need a way to select beams and fields) would also
               // require some form of casting or reprojection between grids of different shape. Therefore, it would unlikely
               // produce the same PSF as one could get using the traditional weighting hack (but PSF is approximate anyway in
               // the joint deconvolution case).

               ASKAPLOG_DEBUG_STR(logger, "Making separate UV weights for the PSF obtained using the full size grid in the joint imaging mode");
               CalcCore& rootImager = *rootImagerPtr; // just for the semantics
               rootImager.setUVWeightCalculator(itsUVWeightCalculator);
               rootImager.setupUVWeightBuilder();
               try {
                  rootImager.accumulateUVWeights();
               }
               catch (const askap::AskapError& e) {
                  ASKAPLOG_WARN_STR(logger,"Askap error in uv weight accumulation - dummy run for rootImager in updatedirection mode, ignoring the current work unit");
                  rootImagerPtr.reset();
                  throw;
               }
               // this will compute weights and add them to the model held by root imager
               rootImager.computeUVWeights();
               rootImager.recreateNormalEquations();
           }
           try {
              rootImagerPtr->calcNE(); // dummy pass (but unlike the code prior to refactoring it is not used for anything else
              rootImagerPtr->configureNormalEquationsForMosaicing();
              rootImagerPtr->zero(); // then we delete all our work ....
              // load a starting model if required
              if (loadStartingModel(rootImagerPtr->params(), globalChannel, globalFrequency)) {
                  copyModel(rootImagerPtr->params(),workingImager.params());
              };
           }
           catch (const askap::AskapError& e) {
                  ASKAPLOG_WARN_STR(logger,"Askap error in worker calcNE - dummy run for rootImager in updatedirection mode");
                  ASKAPLOG_WARN_STR(logger,"Ignoring the current work unit");
                  // reset rootImagerPtr to ensure we attempt full size NE generation with the next work unit
                  rootImagerPtr.reset();
                  throw;
           }

           // this block of code matches similar if-statement above (see notes up there) -  we need to remove all UV weights now (PSF has been generated)
           // to avoid conflict with the weights generated for individual pointings/beams (as those grids can have different size)
           if (isSampleDensityGridNeeded()) {
               UVWeightParamsHelper hlp(rootImagerPtr->params());
               // we could've written a general method (e.g. in ImagerParallel) to remove all uv-weight related parameters, but the following should be sufficient for now
               ASKAPDEBUGASSERT(hlp.exists("slice"));
               hlp.remove("slice");
           }
       }
   } else {
       if (rootImagerPtr) {
           workingImager.replaceModelByReference(rootImagerPtr->params());
       } else {
           if (itsLocalSolver) {
              // load starting model or setup empty model image
              loadStartingModel(workingImager.params(), globalChannel, globalFrequency);
           } else {
               // need to receve the model from master
               // we may need an option to force this behaviour, although alternatively if rootImagerPtr is defined, we
               // can receive the model outside of this method
               ASKAPLOG_DEBUG_STR(logger, "Worker waiting to receive new model");
               workingImager.receiveModel();
               ASKAPLOG_DEBUG_STR(logger, "Worker received initial model for either uv-weights calculation or cycle 0");
           }
       }
   }
   return workingImagerPtr;
}

bool ContinuumWorker::loadStartingModel(const askap::scimath::Params::ShPtr& params, uInt channel, double frequency) const {
    if (itsReadStartingModelCube) {
        // load plane from model cube as starting model
        loadImage(params, channel);
    } else if (itsMFSStartingModel) {
        // Calculate the spectral plane model and set that as the first model
        // subsequent cleaning would be cumulative in this model and it would be restored - problem if
        // models are not on the same grid
        loadImageFromMFSModel(params, frequency, channel);
    } else {
          if (!itsUpdateDir) {
              // set up a zero valued model image
              setupImage(params, frequency, false);
          }
        return false;
    }
    return true;
}

/// @brief helper method to accumulate uv-weights for a single work unit
/// @details It is alalogous to processOneWorkUnit, but is used in the section responsible for traditional weighting.
/// The resulting sample density grid (wrapped into an NE-like object) is either added to the root imager passed as a parameter or
/// a new CalcCore object is created and returned if the passed shared pointer is null.
/// @param[inout] rootImagerPtr shared pointer to CalcCore object to update or create (if empty shared pointer is passed)
/// @param[in] wu work unit to process
void ContinuumWorker::accumulateUVWeightsForOneWorkUnit(boost::shared_ptr<CalcCore> &rootImagerPtr, const cp::ContinuumWorkUnit &wu) const
{
   try {
        // model is received from master or created from scratch inside createImagers taking operating mode into account
        const boost::shared_ptr<CalcCore> workingImagerPtr = createImagers(wu, rootImagerPtr);
        ASKAPDEBUGASSERT(workingImagerPtr);
        if (rootImagerPtr && rootImagerPtr->notStashedNormalEquations()) {
            // this block of code is only executed in the joint deconvolution mode (itsUpdateDir==true) on the first pass, when
            // rootImagerPtr is created. We then stash the normal equations to be able to restore it later before the imaging can begin.
            // This allows us to avoid extra operation to create the full image again.

            rootImagerPtr->stashNormalEquations();

            // just a cross-check, because we don't exect this code to be used in other modes
            ASKAPDEBUGASSERT(itsUpdateDir);

            // the following is required because in the joint deconvolution mode the root imager is not a copy of one of the previous work imagers.

            rootImagerPtr->setUVWeightCalculator(itsUVWeightCalculator);
            // this is necessary to get the right type of NE as we merge data into rootImager later on
            rootImagerPtr->setupUVWeightBuilder();
        }
        CalcCore& workingImager = *workingImagerPtr; // just for the semantics
        // just to get a consistent picture w.r.t. uv-weight calculation mode everywhere
        workingImager.setUVWeightCalculator(itsUVWeightCalculator);
        workingImager.setupUVWeightBuilder();
        workingImager.accumulateUVWeights();
        // if rootImager is empty (i.e. this is the first work unit contributing to it) we just copy
        // the shared pointer
        if (rootImagerPtr) {
            rootImagerPtr->mergeNormalEquations(workingImager);
            ASKAPLOG_DEBUG_STR(logger,"Merged sample density grid for uv-weights");
        } else {
            rootImagerPtr = workingImagerPtr;
        }

        ASKAPDEBUGASSERT(rootImagerPtr);
   }
   catch(const askap::AskapError& e) {
         if (itsLocalSolver) {
             ASKAPLOG_ERROR_STR(logger, "Askap error in uv-weight accumulation - skipping accumulation and carrying on - this may adversely affect PSF (dataset " <<
                                wu.get_dataset()<<" global channel "<<wu.get_globalChannel()<<"(: "<< e.what());
         } else {
             // MV: it is not clear to me whether we should ignore this error, but keep the same behaviour as it was prior to the refactoring
             // for normal imaging 
             ASKAPLOG_ERROR_STR(logger, "Askap error in uv-weight accumulation - skipping accumulation of some or all data in "<<
                               wu.get_dataset()<<" and carrying on: " << e.what());
         }
   }
}

/// @brief helper method to process a single work unit
/// @details This method encapsulates a part of the old processChannel calculating and merging NE for a single work
/// unit. The resulting NE is either added to the root imager passed as the parameter or a new CalcCore object is
/// created and returned if the passed shared pointer is empty.
/// @param[inout] rootImagerPtr shared pointer to CalcCore object to update or create (if empty shared pointer is passed)
/// @param[in] wu work unit to process
/// @param[in] lastcycle if this parameter is true and itsWriteGrids is true as well, the grids are extracted into root imager
/// for writing later on. We only do this in the last major cycle, hence the name. In the central solver case this option has
/// no effect.
void ContinuumWorker::processOneWorkUnit(boost::shared_ptr<CalcCore> &rootImagerPtr, const cp::ContinuumWorkUnit &wu, bool lastcycle) const
{
   try {
        // to accumulate data from the current work unit
        const boost::shared_ptr<CalcCore> workingImagerPtr = createImagers(wu, rootImagerPtr);
        ASKAPDEBUGASSERT(workingImagerPtr);
        CalcCore& workingImager = *workingImagerPtr; // just for the semantics

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
        // MV: I'm not sure we want to log the summary again (it is done just before the processng of this work unit,
        // but this matches the old behaviour prior to refactoring as far as I understand it
        itsStats.logSummary();

        // merge into root image if required.
        // this is required if there is more than one workunit per channel
        // either in time or by beam.
        ASKAPLOG_DEBUG_STR(logger,"About to merge into rootImager");
        if (itsUpdateDir) {
            workingImager.configureNormalEquationsForMosaicing();
        }

        // if rootImager is empty (i.e. this is the first work unit contributing to it) we just copy
        // the shared pointer
        if (rootImagerPtr) {
            rootImagerPtr->mergeNormalEquations(workingImager);
        } else {
            rootImagerPtr = workingImagerPtr;
        }
        ASKAPLOG_DEBUG_STR(logger,"Merged");

        ASKAPDEBUGASSERT(rootImagerPtr);
        if (lastcycle == false) {
            lastcycle = checkStoppingThresholds(rootImagerPtr->params());
        }
        if (itsWriteGrids && lastcycle && itsLocalSolver) {
            ASKAPLOG_INFO_STR(logger, "Extracting grids and summing them in the root imager");
            // the following would work regarless whether root imager and working imager are the same object or not
            workingImager.addGridsToModel(rootImagerPtr->params());
        }
        // MV: it would be nice to think about proper reuse of measurement equation (and associated gridders), currently they are recreated
        // for every work unit which may not be ideal. For now we can dispose the measurement equation in the root imager as it is not
        // needed there. If it doesn't exist there the following call is very cheap.
        rootImagerPtr->resetMeasurementEquation();
   }
   catch(const askap::AskapError& e) {
         if (itsLocalSolver) {
             ASKAPLOG_ERROR_STR(logger, "Askap error in imaging - skipping accumulation and carrying on - this will result in a blank channel (dataset " <<
                                wu.get_dataset()<<" global channel "<<wu.get_globalChannel()<<"(: "<< e.what());
         } else {
             // MV: it is not clear to me whether we should ignore this error, but keep the same behaviour as it was prior to the refactoring
             // (just change the message to better reflect the situation)
             ASKAPLOG_ERROR_STR(logger, "Askap error in imaging - skipping accumulation of some or all data in "<<
                               wu.get_dataset()<<" and carrying on: " << e.what());
         }
   }
}

// this is a rewrite of processChannel method
void ContinuumWorker::processChannels()
{
   ASKAPTRACE("ContinuumWorker::processChannels");

   ASKAPLOG_INFO_STR(logger, "Processing Channel Allocation");

   if (itsWriteGrids) {
       ASKAPLOG_INFO_STR(logger,"Will output gridded visibilities");
   }

   if (itsLocalSolver) {
       ASKAPLOG_INFO_STR(logger, "Processing channels in local solver mode");
   } else {
       ASKAPLOG_INFO_STR(logger, "Processing channels in central solver mode");
   }

   // MV: we can remove this - it is now checked in the constructor (but the code is not too far off from supporting both modes, may look into this after immediate needs are satisfied)
   ASKAPCHECK(!(itsUpdateDir && !itsLocalSolver), "Cannot <yet> Continuum image in on-the-fly mosaick mode - need to update the image parameter setup");

   configureReferenceChannel();

   initialiseCubeWritingIfNecessary();

   if (itsWorkUnits.size() == 0) {
      ASKAPLOG_INFO_STR(logger,"No work todo");

      // write out the beam log
      ASKAPLOG_DEBUG_STR(logger, "About to log the full set of restoring beams");

      logBeamInfo();
      logWeightsInfo();

      return;
   }

   /// What are the plans for the deconvolution?
   ASKAPLOG_DEBUG_STR(logger, "Ascertaining Cleaning Plan");
   const bool writeAtMajorCycle = itsParset.getBool("Images.writeAtMajorCycle", false);
   if (writeAtMajorCycle && itsLocalSolver) {
       ASKAPLOG_WARN_STR(logger, "Images.writeAtMajorCycle is not supported in the local solver case - ignoring");
   }
   const int nCycles = itsParset.getInt32("ncycles", 0);

   const int uvwMachineCacheSize = itsParset.getInt32("nUVWMachines", 1);
   ASKAPCHECK(uvwMachineCacheSize > 0 ,
       "Cache size is supposed to be a positive number, you have "
       << uvwMachineCacheSize);

   const double uvwMachineCacheTolerance = SynthesisParamsHelper::convertQuantity(itsParset.getString("uvwMachineDirTolerance", "1e-6rad"), "rad");

   ASKAPLOG_DEBUG_STR(logger,
       "UVWMachine cache will store " << uvwMachineCacheSize << " machines");
   ASKAPLOG_DEBUG_STR(logger, "Tolerance on the directions is "
       << uvwMachineCacheTolerance / casacore::C::pi * 180. * 3600. << " arcsec");

   const string colName = itsParset.getString("datacolumn", "DATA");
   const bool clearcache = itsParset.getBool("clearcache", false);

   // setup data source manager
   itsDSM.reset(new DataSourceManager(colName, clearcache, static_cast<size_t>(uvwMachineCacheSize), uvwMachineCacheTolerance));

   // itsWorkUnits may include different epochs (for the same channel)
   // the order is strictly by channel - with multiple work units per channel.
   // we use appropriately configured iterators to iterate over all work units with the same
   // frequency if necessary.

   // for continuum process all allocated frequency channels together, so just one block (note - we need to grab the right iterator later on)
   const size_t numberOfFrequencyBlocks = itsLocalSolver ? itsWorkUnits.numberOfFrequencyBlocks() : 1;

   for (size_t frequencyBlock = 0; frequencyBlock < numberOfFrequencyBlocks; ++frequencyBlock) {
        ASKAPLOG_DEBUG_STR(logger, "Processing frequency block "<<frequencyBlock + 1<<" out of "<<numberOfFrequencyBlocks);

        // this will be used to accumulate normal equations across multiple work units
        // define it inside the frequency block loop because each frequency block is considered to be independent
        // (for the continuum case there is only one frequency block comprising all data)
        boost::shared_ptr<CalcCore> rootImagerPtr;

        // begin and end iterators for the part of work units we need (all of them in continuum, given frequency
        // for the spectral line mode)
        const WorkUnitContainer::const_iterator wuBeginIt = itsLocalSolver ? itsWorkUnits.begin(frequencyBlock) : itsWorkUnits.begin();
        const WorkUnitContainer::const_iterator wuEndIt = itsLocalSolver ? itsWorkUnits.end(frequencyBlock) : itsWorkUnits.end();

        if (wuBeginIt == wuEndIt) {
            ASKAPLOG_WARN_STR(logger, "No work assigned for this rank for frequency block "<<frequencyBlock + 1);
            continue;
        }

        // iterator to the last processed work unit (to use in writing, outside the work unit loop)
        // initialise it with end iterator to double check that it got updated within the loop and no wrong data
        // will get through by chance
        WorkUnitContainer::const_iterator wuLastProcessedIt = wuEndIt;

        // the counter is for reporting and to size the writing job in some circumstances
        // define it in the frequency loop as opposed to inner loops for that latter cause
        size_t workUnitCounter = 0u;
        try {
             if (isSampleDensityGridNeeded()) {
                 ASKAPLOG_DEBUG_STR(logger, "Worker rank "<<itsComms.rank()<<" is about to compute weight grid for its portion of the data");
                 workUnitCounter = 0u;
                 for (WorkUnitContainer::const_iterator wuIt = wuBeginIt; wuIt != wuEndIt; ++wuIt) {
                      // counter just for reporting and to spread writing jobs more evenly (although this can probably
                      // be done differently - use the same approach as in the code prior to refactoring for now)
                      ++workUnitCounter;
                      wuLastProcessedIt = wuIt;
                      ASKAPLOG_DEBUG_STR(logger, "Processing work unit "<<workUnitCounter<<" to obtain sample density grid");
                      itsStats.logSummary();
                      accumulateUVWeightsForOneWorkUnit(rootImagerPtr, *wuIt);
                 }
                 ASKAPCHECK(rootImagerPtr, "It is expected that rootImager should be defined by this point in the process of the sample density grid accumulation");
                 if (itsLocalSolver) {
                     // this will compute weights and add them to the model
                     rootImagerPtr->computeUVWeights();
                     ASKAPLOG_DEBUG_STR(logger, "uv-weight has been added to the model");
                 } else {
                     // the following call sends the weight grid back to the master for merging and processing,
                     // the result will be sent back along with the model
                     rootImagerPtr->sendNE();

                     // need to receive the model explicitly (as it is only done inside createImagers call (executed from
                     // accumulateUVWeightsFromOneWorkUnit and processOneWorkUnit) if rootImagerPtr is empty
                     ASKAPLOG_DEBUG_STR(logger, "Worker waiting to receive new model");
                     rootImagerPtr->receiveModel();
                     ASKAPLOG_DEBUG_STR(logger, "Worker received initial model for cycle 0 with uv-weights");
                 }
                 // if needed, set the mosaicing flag to the newly created imaging normal equations to avoid type mismatch
                 if (itsUpdateDir) {
                     // in the joint deconvolution mode restore normal equations (stashed inside accumulateUVWeightsForOneWorkUnit
                     // when root imager has been created) from the buffer because they have some coordinate system information for the full image
                     rootImagerPtr->popNormalEquations();
                 } else {
                     // in all other modes simply recreate normal equations, this will revert it back to the type suitable for imaging
                     // (although, in principle, we could have stash-pop approach for all cases)
                     rootImagerPtr->recreateNormalEquations();
                 }
             }
             for (int majorCycleNumber = 0; majorCycleNumber <= nCycles; ++majorCycleNumber) {
                  if (majorCycleNumber != nCycles) {
                      ASKAPLOG_INFO_STR(logger, "Starting major cycle "<<majorCycleNumber + 1<<" out of "<<nCycles<<
                                            " (frequency block "<<frequencyBlock + 1<<" out of "<<numberOfFrequencyBlocks<<")");
                  } else {
                      ASKAPLOG_INFO_STR(logger, "Concluding major cycles (last pass over data, frequency block "<<frequencyBlock + 1<<" out of "<<numberOfFrequencyBlocks<<")");
                  }
                  workUnitCounter = 0u;
                  for (WorkUnitContainer::const_iterator wuIt = wuBeginIt; wuIt != wuEndIt; ++wuIt) {
                       // MV: the original code prior to refactoring had the codition skipping DONE and NA workunits here
                       // Currently, we don't add them, hence the condition has been removed

                       // counter just for reporting and to spread writing jobs more evenly (although this can probably
                       // be done differently - use the same approach as in the code prior to refactoring for now)
                       ++workUnitCounter;
                       wuLastProcessedIt = wuIt;
                       itsStats.logSummary();

                       ASKAPLOG_DEBUG_STR(logger, "Processing work unit "<<workUnitCounter);

                       processOneWorkUnit(rootImagerPtr, *wuIt, majorCycleNumber == nCycles);

                  } // for loop over workunits of the given frequency block (or all of them in continuum mode)

                  // minor cycle would fail if rootImagerPtr is void (exception handler later on would write a blank image if necessary)
                  if (runMinorCycleSolver(rootImagerPtr, majorCycleNumber < nCycles)) {
                      // we're here if the early termination condition has been triggered (i.e. on thresholds)
                      if (itsLocalSolver || majorCycleNumber == nCycles) {
                          break;
                      }
                      // we need to do one more round to send the final residuals
                      majorCycleNumber = nCycles - 1;
                  }
             } // for loop over major cycles

             ASKAPLOG_INFO_STR(logger," Finished the major cycles");
             // no work is left in the continuum mode except the final logging which is done outside the frequency loop
             if (itsLocalSolver) {
                 ASKAPASSERT(rootImagerPtr);
                 CalcCore& rootImager = *rootImagerPtr; // just for semantics;
                 rootImager.updateSolver();

                 // At this point we have finished our last major cycle. We have the "best" model from the
                 // last minor cycle. Which should be in the archive - or full coordinate system
                 // the residual image should be merged into the archive coordinated as well.
                 addImageAsModel(rootImager.params());

                 rootImager.check();


                 if (itsRestore) {
                     ASKAPLOG_DEBUG_STR(logger, "Running restore");
                     rootImager.restoreImage();
                 }

                 // force cache clearing here to match the code behaviour prior to refactoring.
                 // It will be no operation if clearcache is false
                 if (itsDSM) {
                     itsDSM->reset();
                 }

                 itsStats.logSummary();

                 ASKAPLOG_DEBUG_STR(logger, "writing channel into cube");
                 ASKAPCHECK(wuLastProcessedIt != wuEndIt, "No work units seem to be processed - logic error or bad configuration?");

                 if (itsComms.isWriter()) {

                     // write own portion first
                     performOwnWriteJob(wuLastProcessedIt->get_globalChannel(), rootImager.params());

                     /// write everyone elses

                     /// one per client ... I dont care what order they come in at

                     performOutstandingWriteJobs(itsComms.getOutstanding() > itsComms.getClients().size() ?
                                  itsComms.getOutstanding() - itsComms.getClients().size() : 0,
                             itsWorkUnits.size() - workUnitCounter);

                 } else {
                     // this rank is not the writer, send the result elsewhere
                     ContinuumWorkRequest result;
                     result.set_params(rootImager.params());
                     result.set_globalChannel(wuLastProcessedIt->get_globalChannel());
                     /// send the work to the writer with a blocking send
                     result.sendRequest(wuLastProcessedIt->get_writer(), itsComms);
                     itsComms.removeChannelFromWorker(itsComms.rank());
                 }
             } // if in the spectral line mode
        } // try-block to be able to fail processing of the single spectral channel gracefully
        catch (const std::exception& e) {
               if (!itsLocalSolver) {
                   /// this is MFS/continuum mode
                   /// throw this further up - this avoids a failure in continuum mode generating bogus - or furphy-like
                   /// error messages
                   ASKAPLOG_WARN_STR(logger, "Error processing a channel in continuum mode");
                   throw;
               }
               ASKAPLOG_WARN_STR(logger, "Error in channel processing, skipping: " << e.what());
               std::cerr << "Skipping channel due to error and continuing: " << e.what() << std::endl;

               // Need to either send an empty map - or remove the channel from the list to write
               if (itsComms.isWriter()) {
                   ASKAPLOG_DEBUG_STR(logger, "Marking bad channel as processed in count for writer");
                   itsComms.removeChannelFromWriter(itsComms.rank());
               } else {
                   // MV: I am not sure we have to report the counter here as the failure can happen during solve
                   ASKAPLOG_DEBUG_STR(logger, "Failed on count " << workUnitCounter);
                   ASKAPASSERT(wuLastProcessedIt != wuEndIt);
                   sendBlankImageToWriter(*wuLastProcessedIt);
               }
        } // catch block bypassing channel write in spectral line mode or passing the exception in continuum
   } // for loop over frequency blocks (just one pass for the continuum case)

   ASKAPLOG_DEBUG_STR(logger,"Finished imaging");

   if (itsLocalSolver) {
       // cleanup
       performOutstandingWriteJobs();
   }
   /*
   // Allowing ranks with no workunits (See AXA-3083) caused problems with this barrier, so removing it for now as
   // I don't think it is necessary anyway - see my notes below from the time of refactoring
   //
   // MV: I don't think the barrier is necesary here. If it is needed to ensure all channel write operations are
   // performed before going further, then we have to move it into performOutstandingWriteJobs.
   // Anyway, leave it as it was prior to refactoring for now
   ASKAPLOG_DEBUG_STR(logger, "Rank " << itsComms.rank() << " at processChannels barrier");
   itsComms.barrier(itsComms.theWorkers());
   ASKAPLOG_DEBUG_STR(logger, "Rank " << itsComms.rank() << " passed processChannels barrier");
   */

   // write out the beam log
   ASKAPLOG_DEBUG_STR(logger, "About to log the full set of restoring beams");
   logBeamInfo();
   logWeightsInfo();
}

/// @brief helper method to perform minor cycle activities
/// @details This method encapsulates running the solver at the conclusion of each major cycle and
/// associated data transfers, if necessary. It can be viewed at the place where the minor cycle is
/// performed in the case of the local solver (or the interface part, if the master perofms it as
/// it happens in the continuum mode).
/// @param[in] rootImagerPtr shared pointer to the CalcCore object containing the result of the current
///                          major cycle for the given worker.
/// @param[in] haveMoreMajorCycles flag that more major cycles are to be done (subject to thresholds). In the
///                          central solver mode we expect to receive the new model in this flag is true or
///                          perform new solution ourselves in the case of the local solver.
/// @return true if major cycles have to be terminated due to thresholds being reached
/// @note empty rootImagerPtr causes an exception. This may happen if no data are processed.
bool ContinuumWorker::runMinorCycleSolver(const boost::shared_ptr<CalcCore> &rootImagerPtr, bool haveMoreMajorCycles) const
{
   ASKAPCHECK(rootImagerPtr, "Root imager is empty - no data processed?");
   if (!itsLocalSolver) {
       // we run through the whole allocation
       // lets send it to the master for processing.
       rootImagerPtr->sendNE();

       // MV: for now leave the original barrier in place although it is not required
       ASKAPLOG_DEBUG_STR(logger, "Rank " << itsComms.rank() << " at runMinorCycleSolver barrier");
       itsComms.barrier(itsComms.theWorkers());
       ASKAPLOG_DEBUG_STR(logger, "Rank " << itsComms.rank() << " passed runMinorCycleSolver barrier");

       // now we have to wait for the model (solution) to come back.
       // MV: except on the very last major cycle
       if (haveMoreMajorCycles) {
           ASKAPLOG_DEBUG_STR(logger, "Worker waiting to receive new model");
           rootImagerPtr->receiveModel();
           ASKAPLOG_DEBUG_STR(logger, "Worker received model for use in the next major cycle");
       }
   }
   const bool forcedStopping = checkStoppingThresholds(rootImagerPtr->params());
   // MV: it would be nice to check if continuum and spectral line mode do the same thing in terms of the number of major cycles 
   // (in both cases of stopping on thresholds and on reaching the limit of major cycles)
   const bool lastCycle = forcedStopping || !haveMoreMajorCycles;
   if (itsLocalSolver && !lastCycle) {
       try {
          rootImagerPtr->solveNE();
          itsStats.logSummary();
       }
       catch (const askap::AskapError& e) {
              ASKAPLOG_WARN_STR(logger, "Askap error in solver:" << e.what());
              throw;
       }
   }
   // MV: although fine for now, the following if-statement conceptually is not a part of the minor cycle. It may be better to have it
   // in a separate method making this clean up action more explicit.
   if (haveMoreMajorCycles && !(forcedStopping && itsLocalSolver)) {
       ASKAPLOG_DEBUG_STR(logger, "Preparing for one more iteration - Reset normal equations");
       if (itsUpdateDir) {
           // MV: the original code has the following comment:
           // Actually I've found that I cannot completely empty the NE. As I need the full size PSF and this is stored in the NE
           // So this method pretty much only zeros the weights and the datavector(image)
           rootImagerPtr->zero();
       } else {
           rootImagerPtr->reset();
       }
   }
   itsStats.logSummary();
   return forcedStopping;
}


// this is the old version of processChannels method prior to refactoring. Kept in place for now just in case we need to flip between the
// old and new versions quickly for debugging. It is no longer used.
void ContinuumWorker::processChannelsOld()
{
  ASKAPTRACE("ContinuumWorker::processChannels");

  ASKAPLOG_INFO_STR(logger, "Processing Channel Allocation");

  if (itsWriteGrids) {
    ASKAPLOG_INFO_STR(logger,"Will output gridded visibilities");
  }

  if (itsLocalSolver) {
    ASKAPLOG_INFO_STR(logger, "Processing multiple channels in local solver mode");
  }
  else {
    ASKAPLOG_INFO_STR(logger, "Processing multiple channels in central solver mode");
  }

  ASKAPCHECK(!(itsUpdateDir && !itsLocalSolver), "Cannot <yet> Continuum image in on-the-fly mosaick mode - need to update the image parameter setup");

  configureReferenceChannel();

  initialiseCubeWritingIfNecessary();

  if (itsWorkUnits.size() == 0) {
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

  const int uvwMachineCacheSize = itsParset.getInt32("nUVWMachines", 1);
  ASKAPCHECK(uvwMachineCacheSize > 0 ,
    "Cache size is supposed to be a positive number, you have "
    << uvwMachineCacheSize);

  const double uvwMachineCacheTolerance = SynthesisParamsHelper::convertQuantity(itsParset.getString("uvwMachineDirTolerance", "1e-6rad"), "rad");

  ASKAPLOG_DEBUG_STR(logger,
      "UVWMachine cache will store " << uvwMachineCacheSize << " machines");
  ASKAPLOG_DEBUG_STR(logger, "Tolerance on the directions is "
      << uvwMachineCacheTolerance / casacore::C::pi * 180. * 3600. << " arcsec");

  const string colName = itsParset.getString("datacolumn", "DATA");
  const bool clearcache = itsParset.getBool("clearcache", false);

  itsDSM.reset(new DataSourceManager(colName, clearcache, static_cast<size_t>(uvwMachineCacheSize), uvwMachineCacheTolerance));

  // the itsWorkUnits may include different epochs (for the same channel)
  // the order is strictly by channel - with multiple work units per channel.
  // so you can increment the workUnit until the frequency changes - then you know you
  // have all the workunits for that channel

  boost::shared_ptr<CalcCore> rootImagerPtr;
  bool gridder_initialized = false;

  for (int workUnitCount = 0; workUnitCount < itsWorkUnits.size();) {

    // NOTE:not all of these will have work
    // NOTE:this loop does not increment here.

    try {

      // spin for good workunit
      while (workUnitCount <= itsWorkUnits.size()) {
        if (itsWorkUnits[workUnitCount].get_payloadType() == ContinuumWorkUnit::DONE){
          workUnitCount++;
        }
        else if (itsWorkUnits[workUnitCount].get_payloadType() == ContinuumWorkUnit::NA) {
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
      if (workUnitCount >= itsWorkUnits.size()) {
        ASKAPLOG_INFO_STR(logger, "Out of work with workUnit " << workUnitCount);
        break;
      }
      itsStats.logSummary();
      ASKAPLOG_INFO_STR(logger, "Starting to process workunit " << workUnitCount+1 << " of " << itsWorkUnits.size());

      int initialChannelWorkUnit = workUnitCount;

      if (!itsUpdateDir) {

        // NOTE: this is because if we are mosaicking ON THE FLY. We do
        // not process the first workunit outside the imaging loop.
        // But for "normal" processing the first workunit is processed outside the loops
        // This adds all sorts of complications to the logic BTW.

        initialChannelWorkUnit = workUnitCount+1;
      }

      const cp::ContinuumWorkUnit& currentWorkUnit = itsWorkUnits[workUnitCount];

      double frequency=currentWorkUnit.get_channelFrequency();

      int localChannel = currentWorkUnit.get_localChannel();

      double globalFrequency = currentWorkUnit.get_channelFrequency();
      int globalChannel = currentWorkUnit.get_globalChannel();

      TableDataSource& ds = itsDSM->dataSource(currentWorkUnit.get_dataset());

      /// Need to set up the rootImager here
      if (itsUpdateDir) {
            itsAdvisor->updateDirectionFromWorkUnit(currentWorkUnit);
            // change gridder for initial calcNE in itsUpdateDir mode
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
      /// set up the image for this channel
      /// this will actually build a full image for the first - it is not actually used tho.
      ///
      ASKAPLOG_INFO_STR(logger, "Initialised imager & gridder");
      bool stopping = false;

      if (!itsUpdateDir) {
          // this method just sets up weight calculator if traditional weighting is done or a null shared pointer if not
          // for itsUpdateDir option we have to do weighting in the working imager, root imager just handles the linmos
          // (although this is probably a bit of the technical debt)
          rootImager.createUVWeightCalculator();
      }
      if (!itsLocalSolver) {
        // for central solver weight grid computation happens here, if it is required
        if (rootImager.isSampleDensityGridNeeded()) {
            // the code below is expected to be called in normal continuum case, incompatible with updatedir
            ASKAPASSERT(!itsUpdateDir);
            // MV: a bit of the technical debt here, we don't need the whole model for weights, but we need coordinate systems, shapes and names distributed the right way
            ASKAPLOG_INFO_STR(logger, "Worker waiting to receive new model (just for uv-weight calculation)");
            rootImager.receiveModel();
            ASKAPLOG_DEBUG_STR(logger, "Worker rank "<<itsComms.rank()<<" is about to compute weight grid for its portion of the data");
            rootImager.setupUVWeightBuilder();
            rootImager.accumulateUVWeights();
            // the following call sends the weight grid back to the master for merging and processing,
            // the result will be sent back along with the model
            rootImager.sendNE();
            // revert normal equations back to the type suitable for imaging
            rootImager.recreateNormalEquations();
        }

        //
        // MV: technically, this barrier should be redundant as we'd wait in receiveModel anyway
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

        // for local solver build weights locally too without interrank communication
        // Note, the check for !itsUpdateDir is technically redundant here as traditional weighting will only
        // be setup for rootImager if itsUpdateDir is false. But add it here for clarify.
        if (rootImager.isSampleDensityGridNeeded() && !itsUpdateDir) {
            ASKAPLOG_DEBUG_STR(logger, "Worker rank "<<itsComms.rank()<<" is about to compute weight grid for its portion of the data");
            rootImager.setupUVWeightBuilder();
            rootImager.accumulateUVWeights();
            // this will compute weights and add them to the model
            rootImager.computeUVWeights();
            ASKAPLOG_DEBUG_STR(logger, "uv-weight has been added to the model");
            // revert normal equations back to the type suitable for imaging
            rootImager.recreateNormalEquations();
        }
      }


      try {

        rootImager.calcNE(); // why do this -
        // this essentially forces me to
        // image the full FOV for a single beam
        // but all I want is something to linmos into.
        // But I need this for the solver ....
        // I should find a away to get the NE initialised w/o regridding
        // which would be much better.
        // Why not just use a spheroidal for the PSF gridders (use sphfuncforpsf)/ full FOV (done)
        // FIXME
        if (itsUpdateDir) {
            rootImager.configureNormalEquationsForMosaicing();
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
        // This loops over work units that are the same itsBaseFrequency
        // but probably not the same epoch or beam ....

        while (tempWorkUnitCount < itsWorkUnits.size())   {

          /// need a working imager to allow a merge over epochs for this channel
          /// assuming subsequent workunits are the same channel but either different
          /// epochs or look directions.

          const cp::ContinuumWorkUnit& tempWorkUnit = itsWorkUnits[tempWorkUnitCount];

          if (frequency != tempWorkUnit.get_channelFrequency()) {
            if (itsLocalSolver) { // the frequencies should be the same.
              // THis is probably the normal spectral line or continuum cube mode.
              // each workunit is a different frequency
              ASKAPLOG_INFO_STR(logger,"Change of frequency for workunit");
              break;
            }
          }

          localChannel = tempWorkUnit.get_localChannel();

          globalFrequency = tempWorkUnit.get_channelFrequency();
          TableDataSource& myDs = itsDSM->dataSource(tempWorkUnit.get_dataset());
          try {

            boost::shared_ptr<CalcCore> workingImagerPtr;

            if (itsUpdateDir) {
              itsAdvisor->updateDirectionFromWorkUnit(tempWorkUnit);
              // in itsUpdateDir mode I cannot cache the gridders as they have a tangent point.
              // FIXED: by just having 2 possible working imagers depending on the mode. ... easy really

              boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParset,itsComms,myDs,localChannel,globalFrequency));
              workingImagerPtr = tempIm;

              // this method just sets up weight calculator if traditional weighting is done or a null shared pointer if not
              // (it is used as a flag indicating whether to do traditional weighting)
              // MV: for now set this up only for itsUpdateDir=true and follow the old logic for all other cases,
              // however it may be worth while to be able to regenerate weight in workingImager instead of reusing what has
              // been done in rootImager in the case of itsUpdateDir=false. There is a bit of untidy design / technical debt here.
              workingImagerPtr->createUVWeightCalculator();
            }
            else {

              boost::shared_ptr<CalcCore> tempIm(new CalcCore(itsParset,itsComms,myDs,rootImager.gridder(),localChannel,globalFrequency));
              workingImagerPtr = tempIm;
            }

            CalcCore& workingImager = *workingImagerPtr; // just for the semantics

            ///this loop does the calcNE and the merge of the residual images

            if (itsUpdateDir) {

              const bool useSubSizedImages = true;
              setupImage(workingImager.params(), frequency, useSubSizedImages);

              // if traditional weighting is enabled compute the weight. The code below assumes local
              // computation without interrank communication
              if (workingImager.isSampleDensityGridNeeded()) {
                  ASKAPLOG_DEBUG_STR(logger, "Worker rank "<<itsComms.rank()<<" is about to compute weight grid for its portion of the data");
                  ASKAPASSERT(itsLocalSolver);
                  workingImager.setupUVWeightBuilder();
                  workingImager.accumulateUVWeights();
                  // this will compute weights and add them to the model
                  workingImager.computeUVWeights();
                  ASKAPLOG_DEBUG_STR(logger, "uv-weight has been added to the model");
                  // revert normal equations back to the type suitable for imaging
                  workingImager.recreateNormalEquations();
              }
              ASKAPLOG_DEBUG_STR(logger, "model after uv-weight generation workingImager = "<<*workingImager.params()<<" rootImager = "<<*rootImager.params());
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

            ASKAPLOG_DEBUG_STR(logger,"About to merge into rootImager");
            if (itsUpdateDir) {
              workingImager.configureNormalEquationsForMosaicing();
            }

            rootImager.mergeNormalEquations(workingImager);
            ASKAPLOG_DEBUG_STR(logger,"Merged");
          }
          catch( const askap::AskapError& e) {
            ASKAPLOG_WARN_STR(logger, "Askap error in imaging - skipping accumulation: carrying on - this will result in a blank channel" << e.what());
            std::cerr << "Askap error in: " << e.what() << std::endl;
          }

          if (frequency == tempWorkUnit.get_channelFrequency()) {
            tempWorkUnitCount++;
            // NOTE: here we increment the workunit count.
            // but the frequency is the same so this is just combining epochs or beams.
            // the accumulator does <not> have to be clean.
          }
          else {
            // the frequency has changed - which means for spectral line we break.
            // but for continuum we continue ...
            // this first condition has already been checked earlier in the loop.
            if (itsLocalSolver) {
              break;
            }
            else {
              // update the frequency
              frequency = tempWorkUnit.get_channelFrequency();
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




        if (itsLocalSolver && (majorCycleNumber == nCycles)) { // done the last cycle
          stopping = true;
          break;
        }



        else if (!itsLocalSolver){ // probably continuum mode ....
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
        stopping |= checkStoppingThresholds(rootImager.params());

        if (!itsLocalSolver && (majorCycleNumber == nCycles -1)) {
          stopping = true;
        }

        if (!stopping && itsLocalSolver) {
          try {
            rootImager.solveNE();
            itsStats.logSummary();
          } catch (const askap::AskapError& e) {
            ASKAPLOG_WARN_STR(logger, "Askap error in solver:" << e.what());

            throw;
          }
        }
        else if (stopping && itsLocalSolver) {
          break; // should be done if I am in local solver mode.
        }

        if (!stopping && itsUpdateDir){

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
        else if (!stopping && !itsUpdateDir) {
          // In this case the first workUnit is processed outside the workUnit loop.
          // So we need to calcNE again with the latest model before the major cycle starts.
          //
          // If we are using itsUpdateDir we reprocess all the workunits - so this is not needed.
          ASKAPLOG_INFO_STR(logger, "Continuuing - Reset normal equations");
          rootImager.reset();

          // we have found that resetting the NE is causing some problems after r10290.

          try {
            rootImager.calcNE();
          }
          catch (const askap::AskapError& e) {
            ASKAPLOG_WARN_STR(logger, "Askap error in calcNE after majorcycle: " << e.what());
          }
        }
        else if (stopping && !itsLocalSolver) {
          ASKAPLOG_INFO_STR(logger, "Not local solver but last run - Reset normal equations");
          rootImager.reset();

          if (!itsUpdateDir) {

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



      if (!itsLocalSolver) { // all my work is done - only continue if in local mode
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

      rootImager.updateSolver();

      // At this point we have finished our last major cycle. We have the "best" model from the
      // last minor cycle. Which should be in the archive - or full coordinate system
      // the residual image should be merged into the archive coordinated as well.
      addImageAsModel(rootImager.params());

      if (itsWriteGrids) {
          rootImager.addGridsToModel(rootImager.params());
      }

      rootImager.check();


      if (itsRestore) {
        ASKAPLOG_INFO_STR(logger, "Running restore");
        rootImager.restoreImage();
      }

      // force cache clearing here (although it would be done automatically at the end of the method) to match the
      // code behaviour prior to refactoring. It will be no operation if clearcache is false
      itsDSM->reset();

      itsStats.logSummary();

      ASKAPLOG_INFO_STR(logger, "writing channel into cube");

      if (itsComms.isWriter()) {

        // write own portion first
        performOwnWriteJob(itsWorkUnits[workUnitCount - 1].get_globalChannel(), rootImager.params());

        /// write everyone elses

        /// one per client ... I dont care what order they come in at

        performOutstandingWriteJobs(itsComms.getOutstanding() > itsComms.getClients().size() ?
                                    itsComms.getOutstanding() - itsComms.getClients().size() : 0,
                                    itsWorkUnits.size() - workUnitCount);

      } else {

        ContinuumWorkRequest result;
        result.set_params(rootImager.params());
        result.set_globalChannel(itsWorkUnits[workUnitCount - 1].get_globalChannel());
        /// send the work to the writer with a blocking send
        result.sendRequest(itsWorkUnits[workUnitCount - 1].get_writer(), itsComms);
        itsComms.removeChannelFromWorker(itsComms.rank());

      }

      /// outside the clean-loop write out the slice
    }

    catch (const std::exception& e) {

      if (!itsLocalSolver) {
        /// this is MFS/continuum mode
        /// throw this further up - this avoids a failure in continuum mode generating bogus - or furphy-like
        /// error messages
        ASKAPLOG_WARN_STR(logger, "Error processing a channel in continuum mode");
        throw;
      }

      ASKAPLOG_WARN_STR(logger, "Error in channel processing, skipping: " << e.what());
      std::cerr << "Skipping channel due to error and continuing: " << e.what() << std::endl;

      // Need to either send an empty map - or
      if (itsComms.isWriter()) {
        ASKAPLOG_INFO_STR(logger, "Marking bad channel as processed in count for writer\n");
        itsComms.removeChannelFromWriter(itsComms.rank());
      } else {
        const int goodUnitCount = workUnitCount - 1; // last good one - needed for the correct freq label and writer
        ASKAPLOG_INFO_STR(logger, "Failed on count " << goodUnitCount);
        sendBlankImageToWriter(itsWorkUnits[goodUnitCount]);
      }
      // No need to increment workunit. Although this assumes that we are here because we failed the solveNE not the calcNE

    }

  } // next workunit if required.

  // cleanup
  performOutstandingWriteJobs();

  // write out the beam log
  ASKAPLOG_INFO_STR(logger, "About to log the full set of restoring beams");
  itsComms.barrier(itsComms.theWorkers());
  logBeamInfo();
  logWeightsInfo();

}

/// @brief send blank image to writer
/// @details This method is expected to be used when calculation of a spectral plane is failed for some reason,
/// but some other rank is responsible for writing it. Essentially it sends a blank image with parameters
/// (like frequency and channel) filled from the work unit.
/// @param[in] wu work unit to take the information from
/// @note (MV:) This doesn't seem like a good design, but the behaviour is left the same as it was prior to
/// the refactoring.
void ContinuumWorker::sendBlankImageToWriter(const cp::ContinuumWorkUnit &wu) const
{
   ASKAPLOG_DEBUG_STR(logger, "Sending blankparams to writer " << wu.get_writer());
   askap::scimath::Params::ShPtr blankParams;

   blankParams.reset(new Params(true));
   ASKAPCHECK(blankParams, "blank parameters (images) not initialised");
   setupImage(blankParams, wu.get_channelFrequency());

   ContinuumWorkRequest result;
   result.set_params(blankParams);
   result.set_globalChannel(wu.get_globalChannel());
   /// send the work to the writer with a blocking send
   result.sendRequest(wu.get_writer(), itsComms);
   ASKAPLOG_DEBUG_STR(logger, "Sent");
}

/// @brief perform write job allocated to this rank
/// @details unlike performOutstandingWriteJobs or performSingleWriteJob this method deals with the write
/// job handled entirely by this rank (i.e. its own write job) and, hence, provided explicitly rather than
/// received from another rank.
/// @param[in] globalChannel global channel (i.e. channel in the whole cube) to write
/// @param[in] params shared pointer to the model with required info (should not be empty)
void ContinuumWorker::performOwnWriteJob(unsigned int globalChannel, const boost::shared_ptr<scimath::Params> &params)
{
   ASKAPLOG_INFO_STR(logger, "I have (including my own) " << itsComms.getOutstanding() << " units to write");
   ASKAPLOG_INFO_STR(logger, "I have " << itsComms.getClients().size() << " clients with work");
   ASKAPDEBUGASSERT(params);
   const int cubeChannel = globalChannel - itsBaseCubeGlobalChannel;
   ASKAPLOG_DEBUG_STR(logger, "Attempting to write channel " << cubeChannel << " of " << itsNChanCube);
   ASKAPCHECK((cubeChannel >= 0 || cubeChannel < itsNChanCube), "cubeChannel outside range of cube slice");
   handleImageParams(params, cubeChannel);
   ASKAPLOG_INFO_STR(logger, "Written channel " << cubeChannel);

   // MV: I'm not sure why we remove the channel from the list twice, but it has been copied from the original code prior to refactoring
   itsComms.removeChannelFromWriter(itsComms.rank());

   itsComms.removeChannelFromWorker(itsComms.rank());
}

/// @brief perform one write job for a remote client
/// @details This method is expected to be used for cube writing ranks only. It receives a single
/// write job and performs it.
void ContinuumWorker::performSingleWriteJob()
{
   ASKAPDEBUGASSERT(itsComms.isWriter());
   ContinuumWorkRequest result;
   int id;
   result.receiveRequest(id, itsComms);
   ASKAPLOG_INFO_STR(logger, "Received a request to write from rank " << id);
   const int cubeChannel = result.get_globalChannel() - itsBaseCubeGlobalChannel;
   try {
        ASKAPLOG_INFO_STR(logger, "Attempting to write channel " << cubeChannel << " of " << itsNChanCube);
        ASKAPCHECK((cubeChannel >= 0 || cubeChannel < itsNChanCube), "cubeChannel outside range of cube slice");
        handleImageParams(result.get_params(), cubeChannel);
        ASKAPLOG_INFO_STR(logger, "Written the slice from rank" << id);

   } catch (const askap::AskapError& e) {
        ASKAPLOG_WARN_STR(logger, "Failed to write a channel to the cube: " << e.what());
   }

   itsComms.removeChannelFromWriter(itsComms.rank());
}

/// work units and performs write operation assigned to this rank.
/// @param[in] targetOutstanding desired number of outstanding write jobs at the end of execution
///                              (to spread writing across the iteration), default is all jobs
/// @param[in] minOutstanding    minimal number of outstanding jobs to remain (default - none)
/// @note (MV) I didn't fully understand the logic behind targetOutstanding and minOutstanding (one should be
/// sufficient), the same behaviour as we had prior to refactoring has been implemented.
void ContinuumWorker::performOutstandingWriteJobs(int targetOutstanding, int minOutstanding)
{
   if (itsComms.isWriter()) {
       ASKAPLOG_DEBUG_STR(logger, "this iteration target is " << targetOutstanding);
       ASKAPLOG_DEBUG_STR(logger, "iteration count is " << itsComms.getOutstanding());

       while (itsComms.getOutstanding() > targetOutstanding) {
              if (itsComms.getOutstanding() <= minOutstanding) {
                  ASKAPLOG_DEBUG_STR(logger, "local remaining count is " << minOutstanding);
                  break;
              }
              ASKAPLOG_DEBUG_STR(logger, "I have " << itsComms.getOutstanding() << "outstanding work units");
              performSingleWriteJob();
       }
   }
}

/// @brief check stopping thresholds in the model
/// @details This method is used at the end of minor cycle deconvolution to check whether to continue iterations.
/// @param[in] model shared pointer to the scimath::Params object with the model
/// @return true if stopping is required
/// @note We do similar checks in both the master and in workers. So this method can be moved somewhere else to be shared.
bool ContinuumWorker::checkStoppingThresholds(const boost::shared_ptr<scimath::Params> &model) const
{
   const std::string majorcycle = itsParset.getString("threshold.majorcycle", "-1Jy");
   const double targetPeakResidual = SynthesisParamsHelper::convertQuantity(majorcycle, "Jy");
   if (model && model->has("peak_residual")) {
       const double peak_residual = model->scalarValue("peak_residual");
       ASKAPLOG_INFO_STR(logger, "Reached peak residual of " << abs(peak_residual));
       if (peak_residual < targetPeakResidual) {
           if (peak_residual < 0) {
               ASKAPLOG_WARN_STR(logger, "Clean diverging, did not reach the major cycle threshold of "
                                         << targetPeakResidual << " Jy. Stopping.");
           } else {
              ASKAPLOG_INFO_STR(logger, "It is below the major cycle threshold of "
                                        << targetPeakResidual << " Jy. Stopping.");
           }
           return true;
       } else {
           if (targetPeakResidual < 0) {
               ASKAPLOG_INFO_STR(logger, "Major cycle flux threshold is not used.");
           } else {
               if (model->has("noise_threshold_reached") &&
                   model->scalarValue("noise_threshold_reached")>0) {
                   ASKAPLOG_INFO_STR(logger, "It is below the noise threshold. Stopping.");
                   return true;
               } else {
                 ASKAPLOG_INFO_STR(logger, "It is above the major cycle threshold of "
                                           << targetPeakResidual << " Jy. Continuing.");
               }
          }
       }
   }
   return false;
}

void ContinuumWorker::copyModel(askap::scimath::Params::ShPtr SourceParams, askap::scimath::Params::ShPtr SinkParams) const
{
  ASKAPCHECK(SourceParams && SinkParams, "Either input or output models are not defined");
  askap::scimath::Params& src = *SourceParams;
  askap::scimath::Params& dest = *SinkParams;
  // ASKAPLOG_WARN_STR(logger, "Names are " << src.names());
  // before the restore the image is the model ....
  SynthesisParamsHelper::copyImageParameter(src, dest,"image.slice");

  // uv-weight related parameters are stored as part of the model. If present, the corresponding parameter name 
  // as accepted by UVWeightParamsHelper would be without the leading "image". 
  UVWeightParamsHelper hlp(src);
  hlp.copyTo(dest, "slice");
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
          const casacore::Array<float> arr = itsResidualCube->writeFlexibleSlice(params->valueF("residual.slice"), chan);
          ASKAPLOG_INFO_STR(logger, "Calculating residual stats");
          itsResidualStatsAndMask->calculate(chan,arr);
      }
  }

  // Write weights
  // Keep weights around in case we need them for masking restored image
  casacore::Array<float> wts;
  if (!params->has("weights.slice")) {
      ASKAPLOG_WARN_STR(logger, "Params are missing weights parameter");
  } else {
      if (itsWeightsCube) {
          ASKAPLOG_INFO_STR(logger, "Writing Weights");
          wts = itsWeightsCube->writeFlexibleSlice(params->valueF("weights.slice"), chan);
      } else {
          wts = params->valueF("weights.slice");
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

  // Write the grids - we write all or none, so only need to set gridShape once
  IPosition gridShape;
  // Limit number of fft threads to 8 (more is slower for our fft sizes)
  scimath::FFT2DWrapper<casacore::Complex> fft2d(true,8);

  if (params->has("grid.slice") && (itsVisGridCube||itsVisGridCubeReal)) {
    const casacore::Vector<casacore::Complex> gr(params->complexVectorValue("grid.slice"));
    ASKAPLOG_INFO_STR(logger,"grid.slice shape= "<<params->shape("grid.slice"));
    // turns out we don't always have psf.slice, so need some alternative way to get shape
    if (params->has("psf.slice")) {
        gridShape = params->shape("psf.slice").getFirst(2);
    } else if (itsVisGridCubeReal) {
        gridShape = itsVisGridCubeReal->imageHandler()->shape(itsVisGridCubeReal->filename()).getFirst(2);
    } else {
        gridShape = itsVisGridCube->imageHandler()->shape(itsVisGridCube->filename()).getFirst(2);
    }
    casacore::Matrix<casacore::Complex> grid(gr.reform(gridShape));
    if (itsGridFFT) {
      ASKAPLOG_INFO_STR(logger, "FFTing Vis Grid and writing it as a real image");
      fft2d(grid,false);
      grid *= static_cast<casacore::Float>(grid.nelements());
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
    casacore::Matrix<casacore::Complex> grid(gr.reform(gridShape));
    if (itsGridFFT) {
      ASKAPLOG_INFO_STR(logger, "FFTing PCF Grid and writing it as a real image");
      fft2d(grid,false);
      grid *= static_cast<casacore::Float>(grid.nelements());
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
    casacore::Matrix<casacore::Complex> grid(gr.reform(gridShape));
    if (itsGridFFT) {
      ASKAPLOG_INFO_STR(logger, "FFTing PSF Grid and writing it as a real image");
      fft2d(grid,false);
      grid *= static_cast<casacore::Float>(grid.nelements());
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
      ASKAPCHECK(params->has("fullres.slice")||params->has("image.slice"), "Params are missing image parameter");
      // use fullres if available, otherwise oversample if needed
      casacore::Array<float> arr = (params->has("fullres.slice") ? params->valueF("fullres.slice") :
        itsRestoredCube->createFlexibleSlice(params->valueF("image.slice")));
      if (itsMaskOutput) {
        // mask mosaic output before writing it out
        maskOutput(arr, wts);
      }
      itsRestoredCube->writeRigidSlice(arr, chan);
      itsRestoredStatsAndMask->calculate(chan,arr);
    }
  }
}

void ContinuumWorker::maskOutput(casacore::Array<float>& arr, const casacore::Array<float>& wts)
{
    const double maxWt(casacore::max(wts));

    ASKAPCHECK(maxWt>0., "Maximum wts element is supposed to be positive, check that at least some data were gridded, maxWt="
        <<maxWt);
    const double cutoff=itsMaskLevel*maxWt;

    ASKAPLOG_DEBUG_STR(logger, "Maximum weight " <<maxWt<<
        ", cutoff weight is "<<itsMaskLevel*100<<"\% of the largest diagonal element");
    ASKAPDEBUGASSERT(arr.contiguousStorage() && wts.contiguousStorage());
    float *pArr = arr.data();
    const float * pWts = wts.data();
    for (size_t elem=0; elem<wts.nelements(); ++elem) {
        if (pWts[elem]<=cutoff) {
            casacore::setNaN(pArr[elem]);
        }
    }
}


/// @brief add current image as a model
/// @details This method adds fullres (if present) or ordinary image as model.slice in the given params object.
/// It is expected that model.slice will be absent and either fullres or Nyquist resolution image should be
/// present.
/// @param[in] params shared pointer to the params object to work with (should be non-empty)
/// @note I (MV) think there could be untidy design here - we probably make an extra copy which could be avoided
void ContinuumWorker::addImageAsModel(const boost::shared_ptr<scimath::Params> &params)
{
   ASKAPLOG_INFO_STR(logger,"Adding model.slice");
   ASKAPDEBUGASSERT(params);

   ASKAPCHECK(params->has("image.slice"), "Params are missing image.slice parameter");
   // before archiving "image.slice" as the model, check if a high-resolution "fullres.slice" has been set up
   if (params->has("fullres.slice")) {
       params->add("model.slice", params->valueT("fullres.slice"));
   } else {
       params->add("model.slice", params->valueT("image.slice"));
   }
   ASKAPCHECK(params->has("model.slice"), "Params are missing model.slice parameter");
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

void ContinuumWorker::logBeamInfo() const
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

void ContinuumWorker::logWeightsInfo() const
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
        }
    }
  }
}

void ContinuumWorker::loadImage(const askap::scimath::Params::ShPtr& params, int channel) const
{
  casacore::Array<float> imagePixels = itsImageCube->readRigidSlice(channel);
  const std::string imageName = itsImageCube->filename();
  const casacore::CoordinateSystem imageCoords = itsImageCube->imageHandler()->coordSys(imageName);
  const string name("image.slice");
  boost::optional<float> extraOversampleFactor = itsImageCube->oversamplingFactor();
  SynthesisParamsHelper::loadImageParameter(*params, name, imagePixels, imageCoords,
    extraOversampleFactor, channel);
}

void ContinuumWorker::loadImageFromMFSModel(const askap::scimath::Params::ShPtr& params, double freq, int channel) const
{
  const std::vector<std::string> sources = itsParset.getStringVector("sources.names");
  ASKAPCHECK(sources.size()==1,"Only a single source model image is allowed when preloading an MFS model");
  const std::string modelPar = std::string("sources.")+sources[0]+".model";
  ASKAPCHECK(itsParset.isDefined(modelPar),"No MFS model image specified for "<< modelPar);
  const std::vector<std::string> modelNames = itsParset.getStringVector(modelPar);
  const int nTaylorTerms = itsParset.getInt32(std::string("sources.")+sources[0]+".nterms",1);
  ASKAPCHECK(nTaylorTerms>0, "Number of Taylor terms is supposed to be a positive number, you gave "<<
                      nTaylorTerms);
  ASKAPCHECK(modelNames.size()==1 || modelNames.size()==nTaylorTerms,
    "Model images should be specified as single base name or all taylor terms");
  // convert taylor model to channel model
  casacore::Array<float> imagePixels;
  casacore::DirectionCoordinate inputCoords;
  ImageParamsHelper iph("image."+sources[0]);

  for (int order = nTaylorTerms-1; order >=0; --order) {
      if (nTaylorTerms > 1) {
          // this is an MFS case, setup Taylor terms
          iph.makeTaylorTerm(order);
          ASKAPLOG_DEBUG_STR(logger,"Processing Taylor term "<<order);
      }
      std::string model;
      if (modelNames.size()==1) {
        // only base name is given, need to add taylor suffix
        model = modelNames[0] + iph.suffix();
      } else {
        model = modelNames[order];
      }
      const std::string name = iph.paramName();

      ASKAPLOG_DEBUG_STR(logger, "Adding image " << model << " as model for "<< sources[0]
                          << ", parameter name: "<< name);
      SynthesisParamsHelper::loadImageParameter(*params, name, model);
      const Axes axes = params->axes(name);
      ASKAPDEBUGASSERT(axes.has("FREQUENCY") && axes.start("FREQUENCY")>0);
      const double f0 = axes.start("FREQUENCY");
      const float w = (freq - f0)/f0;
      if (imagePixels.size()==0) {
        imagePixels = params->valueF(name);
        inputCoords = axes.directionAxis();
      } else {
        imagePixels *= w;
        imagePixels += params->valueF(name);
      }
      // clean up
      params->remove(name);
  }
  const std::string imageName = itsImageCube->filename();
  const casacore::CoordinateSystem imageCoords = itsImageCube->imageHandler()->coordSys(imageName);
  const string name("image.slice");
  const boost::optional<float> extraOversampleFactor = itsImageCube->oversamplingFactor();
  IPosition inShape = imagePixels.shape().getFirst(2); 
  const IPosition outShape = itsImageCube->imageHandler()->shape(imageName).getFirst(2);
  ASKAPCHECK(inputCoords.hasSquarePixels()&&imageCoords.directionCoordinate().hasSquarePixels(),"Can't deal with non square pixels yet");
  const double inputInc = abs(inputCoords.increment()(0));
  const double outputInc = abs(imageCoords.directionCoordinate().increment()(0));
  ASKAPCHECK(inputInc>0 && outputInc>0 && inputInc/outputInc<1.00001,"starting model should have the same or smaller cellsize: input="<<inputInc<<", output="<<outputInc);
  ASKAPCHECK(inShape(0)==inShape(1) && outShape(0)==outShape(1),"Can only deal with square images for now");
  // adjust the cellsize of the input to match the output
  SynthesisParamsHelper::adjustCellsize(imagePixels,inputInc, outputInc, inShape(0), outShape(0));
  // get input size again, as it may have changed
  inShape = imagePixels.shape().getFirst(2); 
  // Option to subset the input to the size of the output
  ASKAPCHECK(inShape(0) >= outShape(0), "Model MFS image should be the same size or larger than output image");
  Array<float> pixels;
  if (inShape(0) == outShape(0)) {
      pixels.reference(imagePixels);
  } else {
      ASKAPLOG_DEBUG_STR(logger,"Shape mismatch of MFS model image and output cube, using central part");
      IPosition blc(4,inShape(0)/2 - outShape(0)/2, inShape(1)/2 - outShape(1)/2, 0, 0);
      IPosition trc(4,inShape(0)/2 + outShape(0)/2 - 1, inShape(1)/2 + outShape(1)/2 -1, 0, 0);
      pixels = imagePixels(blc, trc);
  }

  SynthesisParamsHelper::loadImageParameter(*params, name, pixels, imageCoords,
    extraOversampleFactor, channel);
}



void ContinuumWorker::setupImage(const askap::scimath::Params::ShPtr& params,
                                 double channelFrequency, bool shapeOverride) const
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
