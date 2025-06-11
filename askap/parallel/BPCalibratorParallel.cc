/// @file
///
/// BPCalibratorParallel: part of the specialised tool to do optimised bandpass & leakage calibration with
/// limited functionality. Unlike CalibratorParallel, this class
///
///      * solves for bandpass & leakage only
///      * works only with preaveraging calibration approach
///      * does not support multiple chunks in time (i.e. only one solution is made for the whole dataset)
///      * does not support data distribution except per beam
///      * does not support a distributed model (e.h. with individual workers dealing with individual Taylor terms)
///      * does not require exact match between number of workers and number of channel chunks, data are dealt with
///        serially by each worker with multiple iterations over data, if required.
///      * solves normal equations at the worker level in the parallel case
///
/// This specialised tool matches closely BETA needs and will be used for BETA initially (at least until we converge
/// on the best approach to do bandpass calibration). The lifetime of this tool is uncertain at present. In many
/// instances the code is quick and dirty, just to suit our immediate needs.
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
///
/// last change: Wasim Raja <Wasim.Raja@csiro.au>
///      -> XX and YY-visibility phases are referenced independently to XX and YY
///         visibilities of the reference antenna.

// Include own header file first
#include <askap/parallel/BPCalibratorParallel.h>

// logging stuff
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

// own includes
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askapparallel/BlobIBufMW.h>
#include <askap/askapparallel/BlobOBufMW.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobOStream.h>
#include <askap/profile/AskapProfiler.h>

#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/dataaccess/DataIteratorStub.h>

#include <askap/scimath/fitting/LinearSolver.h>
#include <askap/scimath/fitting/GenericNormalEquations.h>
#include <askap/scimath/fitting/Params.h>

#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/MEParsetInterface.h>
#include <askap/measurementequation/CalibrationME.h>
#include <askap/measurementequation/PreAvgCalMEBase.h>
#include <askap/measurementequation/ComponentEquation.h>
#include <askap/measurementequation/NoXPolGain.h>
#include <askap/measurementequation/KeepXPolOnly.h>
#include <askap/measurementequation/LeakageTerm.h>
#include <askap/measurementequation/Product.h>
#include <askap/measurementequation/ImagingEquationAdapter.h>
#include <askap/gridding/VisGridderFactory.h>
#include <askap/askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <askap/calibaccess/CalParamNameHelper.h>
#include <askap/calibaccess/CalibAccessFactory.h>

#include <askap/parallel/CalibratorParallel.h>

#ifdef USE_CAL_SERVICE
#include <askap/calibaccess/ServiceCalSolutionSourceStub.h>
#include <calserviceaccessor/ServiceCalSolutionSource.h>
#endif

// casa includes
#include <casacore/casa/aips.h>
#include <casacore/casa/OS/Timer.h>

#include <algorithm>


namespace askap {

namespace synthesis {

using askap::operator<<;

/// @brief Constructor from ParameterSet
/// @details The parset is used to construct the internal state. We could
/// also support construction from a python dictionary (for example).
/// @param[in] comms communication object
/// @param[in] parset ParameterSet for inputs
BPCalibratorParallel::BPCalibratorParallel(askap::askapparallel::AskapParallel& comms,
          const LOFAR::ParameterSet& parset) : MEParallelApp(comms,emptyDatasetKeyword(parset),false),
      itsPerfectModel(new scimath::Params()), itsRefAntenna(-1), itsSolutionID(-1),
      itsSolveLeakage(false), itsSolveBandpass(false), itsStoreLeakage(false), itsStoreBandpass(false),
      itsPerformGainThresholding(parset.getBool("threshold.gain.enable", false)),
      itsExpectedGainAmplitude(parset.getDouble("threshold.gain.expected", 1.)),
      itsPerformLeakageThresholding(parset.getBool("threshold.leakage.enable", false)),
      itsUseXPolOnly(parset.getBool("xpol_only", false))
{
  ASKAPLOG_INFO_STR(logger, "Bandpass or Leakage will be solved for using a specialised pipeline");
  const std::string what2solve = parset.getString("solve","bandpass");

  const std::string what2store = parset.getString("store", what2solve);
  if (what2store.find("leakages") != std::string::npos) {
      itsStoreLeakage = true;
  }
  if (what2store.find("bandpass") != std::string::npos) {
      itsStoreBandpass = true;
  }

  if (what2solve.find("leakages") != std::string::npos) {
      if (itsStoreLeakage) {
          ASKAPLOG_INFO_STR(logger, "Leakages will be solved for (solve='"<<what2solve<<"')");
      } else {
          ASKAPLOG_INFO_STR(logger, "Leakages will be solved for (solve='"<<what2solve<<"'), but not stored (store='"<<what2store<<"')");
      }
      itsSolveLeakage = true;
  }

  if (what2solve.find("bandpass") != std::string::npos) {
      if (itsStoreBandpass) {
          ASKAPLOG_INFO_STR(logger, "Bandpass will be solved for (solve='"<<what2solve<<"')");
      } else {
          ASKAPLOG_INFO_STR(logger, "Bandpass will be solved for (solve='"<<what2solve<<"'), but not stored (store='"<<what2store<<"')");
      }
      itsSolveBandpass = true;
  }
  ASKAPCHECK(itsSolveLeakage || itsSolveBandpass,"Need to specify solve=leakages or solve=bandpass");
  if (itsUseXPolOnly) {
      ASKAPCHECK(itsSolveLeakage && !itsSolveBandpass, "xpol_only = true option is allowed for leakage-only solutions");
      ASKAPLOG_INFO_STR(logger, "Parallel-hand products will be ignored when making equations");
  }

  if ((itsStoreLeakage != itsSolveLeakage) && (itsStoreBandpass != itsSolveBandpass)) {
      ASKAPLOG_INFO_STR(logger, "Essentially a dry run will occur (store='"<<what2store<<"', solve='"<<what2solve<<"')");
  }

  if (itsPerformGainThresholding) {
      ASKAPCHECK(itsSolveBandpass, "Unable to perform gain thresholding on the solution as the solver is not configured to solve for bandpass!");
      const std::string tolerancePar("threshold.gain.tolerance");
      ASKAPCHECK(parset.isDefined(tolerancePar), "Tolerance for bandpass gains is not defined - unable to perform thresholding");
      itsGainAmplitudeTolerance = parset.getDouble(tolerancePar);
      ASKAPLOG_INFO_STR(logger, "For each beam and channel, the whole solution will be flagged invalid for antennas with gain amplitude deviating by more than "<<
           itsGainAmplitudeTolerance<<" from the expected value of "<<itsExpectedGainAmplitude);
  }
  if (itsPerformLeakageThresholding) {
      ASKAPCHECK(itsSolveLeakage, "Unable to perform leakage thresholding on the solution as the solver is not configured to solve for leakages!");
      const std::string tolerancePar("threshold.leakage.tolerance");
      ASKAPCHECK(parset.isDefined(tolerancePar), "Tolerance for bandpass leakages is not defined - unable to perform thresholding");
      itsLeakageTolerance = parset.getDouble(tolerancePar);
      ASKAPLOG_INFO_STR(logger, "For each beam and channel, the whole solution will be flagged invalid for antennas with leakages deviating by more than "<<
           itsLeakageTolerance<<" from zero");
  }
  const std::string stOverridePar("solution_time");
  if (parset.isDefined(stOverridePar)) {
      itsSolutionTimeOverride = parset.getDouble(stOverridePar);
      ASKAPLOG_INFO_STR(logger, "Solution time will be set to "<<*itsSolutionTimeOverride<<" MJD discarding data time stamp (user override)");
  }

  if (itsComms.isMaster()) {
      // setup solution source (or sink to be exact, because we're writing the solution here)
      itsSolutionSource = accessors::CalibAccessFactory::rwCalSolutionSource(parset);
      ASKAPASSERT(itsSolutionSource);

      if (itsComms.isParallel()) {
          ASKAPLOG_INFO_STR(logger, "The work will be distributed between "<<itsComms.nProcs() - 1<<" workers");
      } else {
          ASKAPLOG_INFO_STR(logger, "The work will be done in serial by the current process");
      }

      // This is sloppy but I need to test whether this is likely to be a service
      // source as I need to reinstantiate the full implementation - as all we get from the factory
      // is a stub.

      const std::string calAccType = parset.getString("calibaccess","parset");
#ifdef USE_CAL_SERVICE
      if (calAccType == "service") {
        itsSolutionSource.reset(new ServiceCalSolutionSource(parset));
        ASKAPLOG_INFO_STR(logger,"Obtaining calibration information from service source");
        ASKAPLOG_INFO_STR(logger,"SolutionID determined by ServiceSource");
        // we are only solving for bandpass
        boost::shared_ptr<ServiceCalSolutionSource> src = boost::dynamic_pointer_cast<ServiceCalSolutionSource>(itsSolutionSource);
        src->solveBandpass();

      }
#endif

  }
  if (itsComms.isWorker()) {
      // set datasets (we cannot rely on the code in base classes because we don't distribute by node here
      setMeasurementSets(parset.getStringVector("dataset"));

      /// Create solver in workers
      itsSolver.reset(new scimath::LinearSolver(1e3));
      ASKAPCHECK(itsSolver, "Solver not defined correctly");
      ASKAPCHECK(!parset.isDefined("refgain"), "usage of refgain is deprecated, define reference antenna instead");
      itsRefAntenna = parset.getInt32("refantenna",-1);
      if (itsRefAntenna >= 0) {
          ASKAPLOG_INFO_STR(logger, "Phases will be rotated, so antenna "<<itsRefAntenna<<" has zero phase for all channels and beams in the first polarisation");
          ASKAPCHECK(itsRefAntenna < static_cast<int>(nAnt()), "Requested reference antenna doesn't exist, nAnt="<<nAnt());
      } else {
          ASKAPLOG_INFO_STR(logger, "No phase rotation will be done between iterations");
      }

      // load sky model, populate itsPerfectModel
      readModels();
      if (itsComms.isParallel()) {
          // setup work units in the parallel case, make beams the first (fastest to change) parameter to achieve
          // greater benefits if multiple measurement sets are present (more likely to be scheduled for different ranks)
          ASKAPLOG_INFO_STR(logger, "Work for "<<nBeam()<<" beams and "<<nChan()<<" channels will be split between "<<
                   (itsComms.nProcs() - 1)<<" ranks, this one handles chunk "<<(itsComms.rank() - 1));
          itsWorkUnitIterator.init(casacore::IPosition(2, nBeam(), nChan()), itsComms.nProcs() - 1, itsComms.rank() - 1);
      }

      ASKAPCHECK((measurementSets().size() == 1) || (measurementSets().size() == nBeam()),
          "Number of measurement sets given in the parset ("<<measurementSets().size()<<
          ") should be either 1 or equal the number of beams solved for ("<<nBeam()<<")");

      // optional beam index mapping
      const std::vector<uint> beamIndices = parset.getUintVector("beamindices", std::vector<uint>());
      if (beamIndices.size() != 0) {
          ASKAPCHECK(beamIndices.size() == nBeam(),
                "The number of explicit beam indices is different from the number of beams to solve for ("<<nBeam()<<")");
          ASKAPLOG_INFO_STR(logger, "Will solve for beams: "<<beamIndices);
          for (size_t index = 0; index < beamIndices.size(); ++index) {
               itsBeamIndexConverter.add(static_cast<int>(index), static_cast<int>(beamIndices[index]));
          }
      }
  }
  if (!itsComms.isParallel()) {
      // setup work units in the serial case - all work to be done here
      ASKAPLOG_INFO_STR(logger, "All work for "<<nBeam()<<" beams and "<<nChan()<<" channels will be handled by this rank");
      itsWorkUnitIterator.init(casacore::IPosition(2, nBeam(), nChan()));
  }

}

/// @brief helper method to remove the dataset name from a parset
/// @details We deal with multiple measurement sets in a dit different
/// way from the other synthesis applications (they are not per worker
/// here). This method allows to remove the string with measurement sets
/// in the parset passed to base classes and replace it by empty string
/// @param[in] parset input parset
/// @return a copy without the dataset keyword
LOFAR::ParameterSet BPCalibratorParallel::emptyDatasetKeyword(const LOFAR::ParameterSet &parset)
{
  LOFAR::ParameterSet result(parset.makeSubset(""));
  result.replace("dataset",std::string());
  return result;
}

/// @brief method which does the main job
/// @details it iterates over all channels/beams and writes the result.
/// In the parallel mode each worker iterates over their own portion of work and
/// then sends the result to master for writing.
void BPCalibratorParallel::run()
{
  if (itsComms.isWorker()) {
      ASKAPDEBUGASSERT(itsModel);
      const int nCycles = parset().getInt32("ncycles", 1);
      ASKAPCHECK(nCycles >= 0, " Number of calibration iterations should be a non-negative number, you have " <<
                       nCycles);
      for (itsWorkUnitIterator.origin(); itsWorkUnitIterator.hasMore(); itsWorkUnitIterator.next()) {
           // this will force creation of the new measurement equation for this beam/channel pair
           itsEquation.reset();
           itsModel->reset();

           const std::pair<casacore::uInt, casacore::uInt> indices = currentBeamAndChannel();
           if (itsSolveBandpass) {
               ASKAPLOG_INFO_STR(logger, "Initialise bandpass (unknowns) for "<<nAnt()<<" antennas for beam="<<indices.first<<
                                 " and channel="<<indices.second);
               for (casacore::uInt ant = 0; ant<nAnt(); ++ant) {
                    itsModel->add(accessors::CalParamNameHelper::paramName(ant, indices.first, casacore::Stokes::XX), casacore::Complex(1.,0.));
                    itsModel->add(accessors::CalParamNameHelper::paramName(ant, indices.first, casacore::Stokes::YY), casacore::Complex(1.,0.));
               }
           }

           if (itsSolveLeakage) {
               ASKAPLOG_INFO_STR(logger, "Initialise leakages (unknowns) for "<<nAnt()<<" antennas for beam="<<indices.first<<
                " and channel="<<indices.second);
               for (casa::uInt ant = 0; ant<nAnt(); ++ant) {
                         itsModel->add(accessors::CalParamNameHelper::paramName(ant, indices.first, casa::Stokes::XY),casa::Complex(0.,0.));
                         itsModel->add(accessors::CalParamNameHelper::paramName(ant, indices.first, casa::Stokes::YX),casa::Complex(0.,0.));
               }
           }

           // setup reference gain, if needed
           if (itsRefAntenna >= 0) {
               //itsRefGain = accessors::CalParamNameHelper::paramName(itsRefAntenna, indices.first, casacore::Stokes::XX);
	       //wasim was here
               itsRefGainXX = accessors::CalParamNameHelper::paramName(itsRefAntenna, indices.first, casacore::Stokes::XX);
               itsRefGainYY = accessors::CalParamNameHelper::paramName(itsRefAntenna, indices.first, casacore::Stokes::YY);
           } else {
               itsRefGainXX = "";
               itsRefGainYY = "";
           }

	   ASKAPLOG_INFO_STR(logger, "*** Starting iterations for beam="<< indices.first<<" and channel="<<indices.second<<" ***");
           for (int cycle = 0; (cycle < nCycles) && validSolution(); ++cycle) {
                ASKAPLOG_DEBUG_STR(logger, "*** Starting calibration iteration " << cycle + 1 << " for beam="<<
                              indices.first<<" and channel="<<indices.second<<" ***");
                // iterator is used to access the current work unit inside calcNE
                calcNE();
                solveNE();
           }
           if (itsComms.isParallel()) {
               // send the model to the master, add beam and channel tags first
               itsModel->add("beam",static_cast<double>(indices.first));
               itsModel->add("channel",static_cast<double>(indices.second));
               itsModel->fix("beam");
               itsModel->fix("channel");
               // unlike for ccalibtator we don't send normal equations to master. Therefore, we need to
               // extract solution time on the worker explicitly
               itsModel->add("solution_time", solutionTime());
               itsModel->fix("solution_time");
               sendModelToMaster();
           } else {
               // serial operation, just write the result
               if (validSolution()) {
                   writeModel();
               }
           }
      }
  }
  if (itsComms.isMaster() && itsComms.isParallel()) {
      const casacore::uInt numberOfWorkUnits = nBeam() * nChan();
      for (casacore::uInt chunk = 0; chunk < numberOfWorkUnits; ++chunk) {
           // asynchronously receive result from workers
           receiveModelFromWorker();
           if (validSolution()) {
               writeModel();
           }
      }
  }
  // Destroy the accessor, which should call syncCache and write the table out.
  ASKAPLOG_DEBUG_STR(logger, "Syncing the cached bandpass table to disk");
  itsSolAcc.reset();
}

/// @brief verify that the current solution is valid
/// @details We use a special keywork 'invalid' in the model to
/// signal that a particular solution failed. for whatever reason.
/// This flag is checked to avoid writing the solution (which would
/// automatically set validity flag
/// @return true, if the current solution is valid
bool BPCalibratorParallel::validSolution() const
{
   if (itsModel) {
       return !itsModel->has("invalid");
   }
   return false;
}

/// @brief extract current beam/channel pair from the iterator
/// @details This method encapsulates interpretation of the output of itsWorkUnitIterator.cursor() for workers and
/// in the serial mode. However, it extracts the current beam and channel info out of the model for the master
/// in the parallel case. This is done because calibration data are sent to the master asynchronously and there is no
/// way of knowing what iteration in the worker they correspond to without looking at the data.
/// @return pair of beam (first) and channel (second) indices
std::pair<casacore::uInt, casacore::uInt> BPCalibratorParallel::currentBeamAndChannel() const
{
  if (itsComms.isMaster() && itsComms.isParallel()) {
      ASKAPDEBUGASSERT(itsModel);
      ASKAPDEBUGASSERT(itsModel->has("beam") && itsModel->has("channel"));
      const double beam = itsModel->scalarValue("beam");
      const double channel = itsModel->scalarValue("channel");
      ASKAPDEBUGASSERT((beam >= 0.) && (channel >= 0.));
      const std::pair<casacore::uInt,casacore::uInt> result(static_cast<casacore::uInt>(beam), static_cast<casacore::uInt>(channel));
      ASKAPDEBUGASSERT(result.first < nBeam());
      ASKAPDEBUGASSERT(result.second < nChan());
      return result;
  } else {
      const casacore::IPosition cursor = itsWorkUnitIterator.cursor();
      ASKAPDEBUGASSERT(cursor.nelements() == 2);
      ASKAPDEBUGASSERT((cursor[0] >= 0) && (cursor[1] >= 0));
      ASKAPDEBUGASSERT(static_cast<casacore::uInt>(cursor[0]) < nBeam());
      const std::pair<casacore::uInt,casacore::uInt> result(static_cast<casacore::uInt>(itsBeamIndexConverter(cursor[0])), static_cast<casacore::uInt>(cursor[1]));
      ASKAPDEBUGASSERT(result.first < nBeam());
      ASKAPDEBUGASSERT(result.second < nChan());
      return result;
  }
}

#define BPCALIBRATOR_PARALLEL_BLOB_STREAM_VERSION 1

/// @brief send current model to the master
/// @details This method is supposed to be called from workers in the parallel mode and
/// sends the current results to the master rank
void BPCalibratorParallel::sendModelToMaster() const
{
   ASKAPDEBUGTRACE("BPCalibratorParallel::sendModelToMaster");
   ASKAPLOG_DEBUG_STR(logger, "Sending results to the master");
   itsComms.notifyMaster();
   ASKAPDEBUGASSERT(itsModel);

   askapparallel::BlobOBufMW bobmw(itsComms, 0);
   LOFAR::BlobOStream out(bobmw);
   out.putStart("calmodel", BPCALIBRATOR_PARALLEL_BLOB_STREAM_VERSION);
   out << *itsModel;
   out.putEnd();
   bobmw.flush();
}

/// @brief asynchronously receive model from one of the workers
/// @details This method is supposed to be used in the master rank in the parallel mode. It
/// waits until the result becomes available from any of the workers and then stores it
/// in itsModel.
void BPCalibratorParallel::receiveModelFromWorker()
{
   ASKAPDEBUGTRACE("BPCalibratorParallel::receiveModelFromWorker");
   itsModel.reset(new scimath::Params);

   // wait for the notification
   const int source = itsComms.waitForNotification().first;
   ASKAPLOG_DEBUG_STR(logger, "Receiving results from rank "<<source);

   askapparallel::BlobIBufMW bibmw(itsComms, source);
   LOFAR::BlobIStream in(bibmw);
   const int version = in.getStart("calmodel");
   ASKAPASSERT(version == BPCALIBRATOR_PARALLEL_BLOB_STREAM_VERSION);
   in >> *itsModel;
   in.getEnd();
}



/// @brief Calculate the normal equations (runs in workers)
/// @details Model, either image-based or component-based, is used in conjunction with
/// CalibrationME to calculate the generic normal equations.
void BPCalibratorParallel::calcNE()
{
  ASKAPDEBUGASSERT(itsComms.isWorker());

  // create a new instance of the normal equations class
  boost::shared_ptr<scimath::GenericNormalEquations> gne(new scimath::GenericNormalEquations);
  itsNe = gne;

  ASKAPDEBUGASSERT(itsNe);

  // obtain details on the current iteration, i.e. beam and channel
  ASKAPDEBUGASSERT(itsWorkUnitIterator.hasMore());

  // first is beam, second is channel
  const std::pair<casacore::uInt, casacore::uInt> indices = currentBeamAndChannel();

  ASKAPDEBUGASSERT((measurementSets().size() == 1) || (indices.first < measurementSets().size()));

  const std::string ms = (measurementSets().size() == 1 ? measurementSets()[0] : measurementSets()[indices.first]);

  // actual computation
  calcOne(ms, indices.second, indices.first);
}

/// @brief helper method to invalidate current solution
void BPCalibratorParallel::invalidateSolution() {
   ASKAPDEBUGASSERT(itsModel);
   itsModel->add("invalid",1.);
   itsModel->fix("invalid");
}


/// @brief Solve the normal equations (runs in workers)
/// @details Parameters of the calibration problem are solved for here
void BPCalibratorParallel::solveNE()
{
  if (itsComms.isWorker()) {
      ASKAPLOG_DEBUG_STR(logger, "Solving normal equations");
      ASKAPDEBUGASSERT(itsNe);
      const std::vector<std::string> unknowns = itsNe->unknowns();
      if (unknowns.size() == 0) {
          ASKAPLOG_WARN_STR(logger, "Normal equations are empty - no valid data found, flagging the solution as bad");
          invalidateSolution();
          return;
      }
      casacore::Timer timer;
      timer.mark();
      scimath::Quality q;
      ASKAPDEBUGASSERT(itsSolver);
      ASKAPDEBUGASSERT(itsModel);
      // if additional selectors are used, the shape of the MS may be such that some antennas/beams are not present at all
      // this class uses automatic resizing of buffers and therefore may attempt solving for parameter which is not in the
      // normal equations. The code below fixes such parameters.
      const std::vector<std::string> freeNames = itsModel->freeNames();
      for (std::vector<std::string>::const_iterator ci = freeNames.begin(); ci != freeNames.end(); ++ci) {
           if (std::find(unknowns.begin(), unknowns.end(), *ci) == unknowns.end()) {
               ASKAPLOG_DEBUG_STR(logger, "Parameter "<<*ci<<" is missing in the normal equations - no data");
               itsModel->fix(*ci);
           }
      }
      // now all missing parameters should be fixed

      itsSolver->init();
      itsSolver->addNormalEquations(*itsNe);

      const std::string solverType = parset().getString("solver", "SVD");
      itsSolver->setAlgorithm(solverType);

      if (solverType == "LSQR") {
          std::map<std::string, std::string> params = CalibratorParallel::getLSQRSolverParameters(parset());
          itsSolver->setParameters(params);
      }

      itsSolver->solveNormalEquations(*itsModel,q);
      ASKAPLOG_DEBUG_STR(logger, "Solved normal equations in "<< timer.real() << " seconds ");
      ASKAPLOG_DEBUG_STR(logger, "Solution quality: "<<q);

      const unsigned int minRank = parset().getUint32("minrank",15u);
      if (q.rank() < minRank) {
          ASKAPLOG_WARN_STR(logger, "Solution failed - minimum rank is "<<minRank<<", normal matrix has rank = "<<q.rank());
          invalidateSolution();
          return;
      }

      //wasim was here
      /*
      if (itsRefGain != "") {
          ASKAPLOG_DEBUG_STR(logger, "Rotating phases to have that of "<<itsRefGain<<" equal to 0");
          rotatePhases();
      }
      */
      if (itsRefGainXX != "") {
	      if (itsRefGainXX == itsRefGainYY){
		       ASKAPLOG_DEBUG_STR(logger, "Rotating both XX and YY phases to have that of "<<
                  itsRefGainXX<<" equal to 0");
	      } else{
		       ASKAPLOG_DEBUG_STR(logger, "Rotating XX phases to have that of "<<
                  itsRefGainXX<<" equal to 0 and YY phases to have that of "<<
                  itsRefGainYY<<" equal to 0");
	      }
          rotatePhases();
      }
  }
}

/// @brief Write the results (runs in master)
/// @details The solution (calibration parameters) is reported via solution accessor
void BPCalibratorParallel::writeModel(const std::string &)
{
  ASKAPDEBUGASSERT(itsComms.isMaster());

  const std::pair<casacore::uInt, casacore::uInt> indices = currentBeamAndChannel();

  ASKAPLOG_DEBUG_STR(logger, "Writing results of the calibration for beam="<<indices.first<<" channel="<<indices.second);

  ASKAPCHECK(itsSolutionSource, "Solution source has to be defined by this stage");

  ASKAPDEBUGASSERT(itsModel);

  size_t nDiscardedIntentionally = 0;
  size_t nDiscarded = 0;

  std::vector<std::string> parlist = itsModel->freeNames();
  for (std::vector<std::string>::const_iterator it = parlist.begin(); it != parlist.end(); ++it) {
       const casacore::Complex val = itsModel->complexValue(*it);
       const std::pair<accessors::JonesIndex, casacore::Stokes::StokesTypes> paramType =
             accessors::CalParamNameHelper::parseParam(*it);
       // beam is also coded in the parameters, although we don't need it because the data are partitioned
       // just cross-check it
       ASKAPDEBUGASSERT(static_cast<casacore::uInt>(paramType.first.beam()) == indices.first);
       const bool isBandpass = paramType.second == casacore::Stokes::XX || paramType.second == casacore::Stokes::YY;
       const bool isLeakage = paramType.second == casacore::Stokes::XY || paramType.second == casacore::Stokes::YX;
       bool toBeStored = (isBandpass && itsStoreBandpass) || (isLeakage && itsStoreLeakage);

       // it is not clear whether it is better to build the list up front which would result in another loop over
       // all parameters + parsing or, as we do here, build it on-the-fly which results in repeats for all products
       // but the number of parameters is small anyway as we partition by beam and channel
       if (toBeStored) {
           if (itsPerformGainThresholding) {
               const std::string xGainPar = accessors::CalParamNameHelper::paramName(paramType.first.antenna(), paramType.first.beam(), casacore::Stokes::XX);
               const std::string yGainPar = accessors::CalParamNameHelper::paramName(paramType.first.antenna(), paramType.first.beam(), casacore::Stokes::YY);
               if ((casacore::abs(casacore::abs(itsModel->complexValue(xGainPar)) - itsExpectedGainAmplitude) >= itsGainAmplitudeTolerance) ||
                   (casacore::abs(casacore::abs(itsModel->complexValue(yGainPar)) - itsExpectedGainAmplitude) >= itsGainAmplitudeTolerance)) {
                    toBeStored = false;
               }
           }
       } else {
         ++nDiscardedIntentionally;
       }
       if (toBeStored) {
           if (itsPerformLeakageThresholding) {
               const std::string leakagePar1 = accessors::CalParamNameHelper::paramName(paramType.first.antenna(), paramType.first.beam(), casacore::Stokes::XY);
               const std::string leakagePar2 = accessors::CalParamNameHelper::paramName(paramType.first.antenna(), paramType.first.beam(), casacore::Stokes::YX);
               if ((casacore::abs(itsModel->complexValue(leakagePar1)) >= itsLeakageTolerance) ||
                   (casacore::abs(itsModel->complexValue(leakagePar2)) >= itsLeakageTolerance)) {
                    toBeStored = false;
               }
           }
       }
       if (toBeStored) {
           if (!itsSolAcc) {
               // this is the first attempt to write calibration information - set up the accessor
               // solution accessor shared pointer acts as a flag that we need to set everything up
               ASKAPLOG_DEBUG_STR(logger, "About to set the solution accessor");
               ASKAPDEBUGASSERT(itsComms.isMaster());
               // obtain solution ID only once, the results can come in random order and the
               // accessor is responsible for aggregating all of them together. This is done based on this ID.
               itsSolutionID = itsSolutionSource->newSolutionID(solutionTime());

               itsSolAcc = itsSolutionSource->rwSolution(itsSolutionID);
               ASKAPLOG_DEBUG_STR(logger, "Have set up the solution accessor with id = "<<itsSolutionID);
               ASKAPDEBUGASSERT(itsSolAcc);
           }
           itsSolAcc->setBandpassElement(paramType.first, paramType.second, indices.second, val);
       } else {
         ++nDiscarded;
       }
  }
  ASKAPLOG_DEBUG_STR(logger, "Total number of calibration parameters: "<<parlist.size()<<" intentionally discarded: "<<nDiscardedIntentionally<<
       " discarded due to thresholding: "<<nDiscarded - nDiscardedIntentionally);
}

/// @brief create measurement equation
/// @details This method initialises itsEquation with shared pointer to a proper type.
/// It uses internal flags to create a correct type (i.e. polarisation calibration or
/// just antenna-based gains). Parameters are passed directly to the constructor of
/// CalibrationME template.
/// @param[in] dsi data shared iterator
/// @param[in] perfectME uncorrupted measurement equation
void BPCalibratorParallel::createCalibrationME(const accessors::IConstDataSharedIter &dsi,
                const boost::shared_ptr<IMeasurementEquation const> &perfectME)
{
   ASKAPDEBUGASSERT(itsModel);
   ASKAPDEBUGASSERT(perfectME);

   // it is handy to have a shared pointer to the base type because it is
   // not templated
   boost::shared_ptr<PreAvgCalMEBase> preAvgME;
   // solve as normal gains (rather than bandpass) because only one channel is supposed to be selected
   // this also opens a possibility to use several (e.g. 54 = coarse resolution) channels to get one gain
   // solution which is then replicated to all channels involved.
   if (itsSolveBandpass && !itsSolveLeakage) {
          preAvgME.reset(new CalibrationME<NoXPolGain, PreAvgCalMEBase>());
   } else if (itsSolveLeakage && !itsSolveBandpass) {
          if (itsUseXPolOnly) {
              preAvgME.reset(new CalibrationME<Product<KeepXPolOnly,LeakageTerm>, PreAvgCalMEBase>());
          } else {
              preAvgME.reset(new CalibrationME<LeakageTerm, PreAvgCalMEBase>());
          }
   } else if (itsSolveLeakage && itsSolveBandpass) {
          preAvgME.reset(new CalibrationME<Product<NoXPolGain,LeakageTerm>, PreAvgCalMEBase>());
   }
   ASKAPDEBUGASSERT(preAvgME);

   ASKAPDEBUGASSERT(dsi.hasMore());
   preAvgME->accumulate(dsi,perfectME);
   itsEquation = preAvgME;
   // after a call to accumulate the buffer will be setup appropriately, so we can query the stokes vector
   const casa::Vector<casa::Stokes::StokesTypes> stokes = preAvgME->stokes();

   // go through parameters and fix them if there is no data
   const std::vector<std::string> params(itsModel->freeNames());
   for (std::vector<std::string>::const_iterator ci = params.begin(); ci != params.end(); ++ci) {
        const std::pair<accessors::JonesIndex, casa::Stokes::StokesTypes> parsed = accessors::CalParamNameHelper::parseParam(*ci);
        casa::uInt pol = 0;
        for (; pol < stokes.nelements(); ++pol) {
             if (stokes[pol] == parsed.second) {
                 break;
             }
        }
        if (pol < stokes.nelements()) {
            if (preAvgME->hasDataAccumulated(parsed.first.antenna(), parsed.first.beam(), pol)) {
                continue;
            }
        }
        // no data for the given parameter - fix it
        itsModel->fix(*ci);
   }
   //

   // this is just because we bypass setting the model for the first major cycle
   // in the case without pre-averaging
   itsEquation->setParameters(*itsModel);
}

/// @brief helper method to rotate all phases
/// @details This method rotates the phases of all gains in itsModel
/// to have the phase of itsRefGain exactly 0. This operation does
/// not seem to be necessary for SVD solvers, however it simplifies
/// "human eye" analysis of the results (otherwise the phase degeneracy
/// would make the solution different from the simulated gains).
/// @note The method throws exception if itsRefGain is not among
/// the parameters of itsModel
void BPCalibratorParallel::rotatePhases()
{
  // the intention is to rotate phases in worker (for this class)
  ASKAPDEBUGASSERT(itsComms.isWorker());
  ASKAPDEBUGASSERT(itsModel);
  //wasim was here
 /*
  ASKAPCHECK(itsModel->has(itsRefGain), "phase rotation to `"<<itsRefGain<<
             "` is impossible because this parameter is not present in the model");
  casacore::Complex  refPhaseTerm = casacore::polar(1.f,-arg(itsModel->complexValue(itsRefGain)));
 */
  ASKAPCHECK(itsModel->has(itsRefGainXX), "phase rotation to `"<<itsRefGainXX<<
             "` is impossible because this parameter is not present in the model");
  ASKAPCHECK(itsModel->has(itsRefGainYY), "phase rotation to `"<<itsRefGainYY<<
             "` is impossible because this parameter is not present in the model");
  casacore::Complex  refPhaseTermXX = casacore::polar(1.f,-arg(itsModel->complexValue(itsRefGainXX)));
  casacore::Complex  refPhaseTermYY = casacore::polar(1.f,-arg(itsModel->complexValue(itsRefGainYY)));
  std::vector<std::string> names(itsModel->freeNames());
  for (std::vector<std::string>::const_iterator it=names.begin(); it!=names.end();++it)  {
       const std::string parname = *it;
       //wasim was here
       /*if (parname.find("gain") != std::string::npos) {
           itsModel->update(parname, itsModel->complexValue(parname) * refPhaseTerm);
       } */
       if (parname.find("gain.g11") != std::string::npos) {
           itsModel->update(parname, itsModel->complexValue(parname) * refPhaseTermXX);
       }
       else if (parname.find("gain.g22") != std::string::npos) {
           itsModel->update(parname, itsModel->complexValue(parname) * refPhaseTermYY);
       }
  }
}

/// @brief helper method to extract solution time from NE or model
/// @details To be able to time tag the calibration solutions we add
/// start and stop times extracted from the dataset as metadata to normal
/// equations. This method extracts start time from NE. Because NEs are not
/// shipped to the master for bp-calibrator, if called on master the method looks
/// for a different fixed keyword in the model (which should be created on the worker
/// before the model is sent to master).
/// @return solution time (seconds since 0 MJD)
/// @note if no start/stop time metadata are present in the normal equations or model
/// this method returns 0. In addition, if itsSolutionTimeOverride field is defined,
/// it will be returned instead.
double BPCalibratorParallel::solutionTime() const
{
  if (itsSolutionTimeOverride) {
      // override field is in days for user convenience
      return (*itsSolutionTimeOverride) * 86400;
  }
  // for bp-calibrator normal equations with actual solution time are available only on workers (in the parallel case)
  // for master we rely on explicit copy shipped together with the solution
  if (itsComms.isParallel() && itsComms.isMaster()) {
      if (itsModel) {
          if (itsModel->has("solution_time")) {
              return itsModel->scalarValue("solution_time");
          }
      }
      return 0;
  }

  // use the earliest time corresponding to the data used to make this calibration solution
  // to tag the solution. A request for any latest time than this would automatically
  // extract this solution as most recent.
  ASKAPASSERT(itsNe);

  boost::shared_ptr<scimath::GenericNormalEquations> gne = boost::dynamic_pointer_cast<scimath::GenericNormalEquations>(itsNe);
  if (gne) {
      const scimath::Params& metadata = gne->metadata();
      if (metadata.has("min_time")) {
          return metadata.scalarValue("min_time");
      }
  }
  return 0.;
}

/// Calculate normal equations for one data set, channel and beam
/// @param[in] ms Name of data set
/// @param[in] chan channel to work with
/// @param[in] beam beam to work with
void BPCalibratorParallel::calcOne(const std::string& ms, const casacore::uInt chan, const casacore::uInt beam)
{
  casacore::Timer timer;
  timer.mark();
  ASKAPLOG_DEBUG_STR(logger, "Calculating normal equations for " << ms <<" channel "<<chan<<" beam "<<beam);
  // First time around we need to generate the equation
  if (!itsEquation) {
      ASKAPLOG_DEBUG_STR(logger, "Creating measurement equation" );
      accessors::TableDataSource ds(ms, accessors::TableDataSource::MEMORY_BUFFERS, dataColumn());
      ds.configureUVWMachineCache(uvwMachineCacheSize(),uvwMachineCacheTolerance());
      accessors::IDataSelectorPtr sel=ds.createSelector();
      sel << parset();
      sel->chooseChannels(1,chan);
      sel->chooseFeed(beam);
      accessors::IDataConverterPtr conv=ds.createConverter();
      conv->setFrequencyFrame(getFreqRefFrame(), "Hz");
      conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));
      // ensure that time is counted in seconds since 0 MJD
      conv->setEpochFrame();
      accessors::IConstDataSharedIter it=ds.createConstIterator(sel, conv);
      ASKAPCHECK(it.hasMore(), "No data seem to be available for channel "<<chan<<" and beam "<<beam);

      ASKAPCHECK(itsModel, "Initial assumption of parameters is not defined");

      if (!itsPerfectME) {
          ASKAPLOG_DEBUG_STR(logger, "Constructing measurement equation corresponding to the uncorrupted model");
          ASKAPCHECK(itsPerfectModel, "Uncorrupted model not defined");
          if (SynthesisParamsHelper::hasImage(itsPerfectModel)) {
              ASKAPCHECK(!SynthesisParamsHelper::hasComponent(itsPerfectModel),
                         "Image + component case has not yet been implemented");
              // have to create an image-specific equation
              boost::shared_ptr<ImagingEquationAdapter> ieAdapter(new ImagingEquationAdapter);
              ASKAPCHECK(gridder(), "Gridder not defined");
              ieAdapter->assign<ImageFFTEquation>(*itsPerfectModel, gridder());
              itsPerfectME = ieAdapter;
          } else {
              // model is a number of components, don't need an adapter here

              // it doesn't matter which iterator is passed below. It is not used
              accessors::IDataSharedIter stubIter(new accessors::DataIteratorStub(1));

              boost::shared_ptr<ComponentEquation>
                  compEq(new ComponentEquation(*itsPerfectModel,stubIter));
              itsPerfectME = compEq;
          }
      }
      // now we could've used class data members directly instead of passing them to createCalibrationME
      createCalibrationME(it,itsPerfectME);
      ASKAPCHECK(itsEquation, "Equation is not defined");
  } else {
      ASKAPLOG_DEBUG_STR(logger, "Reusing measurement equation" );
      // we need to update the model held by measurement equation
      // because it has been cloned at construction
      ASKAPCHECK(itsEquation, "Equation is not defined");
      ASKAPCHECK(itsModel, "Model is not defined");
      itsEquation->setParameters(*itsModel);
  }
  ASKAPCHECK(itsNe, "NormalEquations are not defined");
  itsEquation->calcEquations(*itsNe);
  ASKAPLOG_DEBUG_STR(logger, "Calculated normal equations for "<< ms << " channel "<<chan<<" beam " <<beam<<" in "<< timer.real()
                     << " seconds ");
}


} // namespace synthesis

} // namespace askap
