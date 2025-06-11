/// @file ContinuumMaster.cc
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
#include "ContinuumMaster.h"

// System includes
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>
#include <askap/askapparallel/AskapParallel.h>

#include <Common/ParameterSet.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/Axes.h>
/*
#include <askap/dataaccess/IConstDataSource.h>
#include <askap/dataaccess/TableConstDataSource.h>
#include <askap/dataaccess/IConstDataIterator.h>
#include <askap/dataaccess/IDataConverter.h>
#include <askap/dataaccess/IDataSelector.h>
#include <askap/dataaccess/IDataIterator.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/TableInfoAccessor.h>
*/
#include <casacore/casa/Quanta.h>
#include <askap/imageaccess/BeamLogger.h>
#include <askap/parallel/ImagerParallel.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>

// Local includes
#include "askap/distributedimager/AdviseDI.h"
#include "askap/distributedimager/CalcCore.h"
#include "askap/distributedimager/CubeComms.h"
#include "askap/messages/ContinuumWorkUnit.h"
#include "askap/messages/ContinuumWorkRequest.h"


//casacore includes
#include "casacore/ms/MeasurementSets/MeasurementSet.h"
#include "casacore/ms/MeasurementSets/MSColumns.h"

using namespace std;
using namespace askap::cp;
using namespace askap;

ASKAP_LOGGER(logger, ".ContinuumMaster");

ContinuumMaster::ContinuumMaster(LOFAR::ParameterSet& parset,
                                       CubeComms& comms, StatReporter& stats)
    : itsParset(parset), itsComms(comms), itsStats(stats), itsBeamList()
{
}

ContinuumMaster::~ContinuumMaster()
{
}

void ContinuumMaster::run(void)
{
    ASKAPTRACE("ContinuumMaster::run");
    // print out the parset
    ASKAPLOG_INFO_STR(logger,"Parset: \n"<<itsParset);
    // Read from the configuration the list of datasets to process
    const vector<string> ms = getDatasets();
    if (ms.size() == 0) {
        ASKAPTHROW(std::runtime_error, "No datasets specified in the parameter set file");
    }
    // Need to break these measurement sets into groups
    // there are three posibilties:
    // 1 - the different measurement sets have the same epoch - but different
    //      frequencies
    // 2 - they have different epochs but the same TOPO centric frequencies

    vector<int> theBeams = getBeams();

    const double targetPeakResidual = synthesis::SynthesisParamsHelper::convertQuantity(
                itsParset.getString("threshold.majorcycle", "-1Jy"), "Jy");


    const bool writeAtMajorCycle = itsParset.getBool("Images.writeAtMajorCycle", false);
    const int nCycles = itsParset.getInt32("ncycles", 0);
    const bool localSolver = itsParset.getBool("solverpercore",false);
    synthesis::AdviseDI diadvise(itsComms,itsParset);

    try {

        diadvise.prepare();
        diadvise.addMissingParameters(false);

        ASKAPLOG_DEBUG_STR(logger,"*****");
        ASKAPLOG_DEBUG_STR(logger,"Parset" << itsParset);
        ASKAPLOG_DEBUG_STR(logger,"*****");

        const int totalChannels = diadvise.getTopoFrequencies().size();

        ASKAPLOG_INFO_STR(logger,"AdviseDI reports " << totalChannels << " channels to process");
        ASKAPLOG_INFO_STR(logger,"AdviseDI reports " << diadvise.getWorkUnitCount() << " work units to allocate");
    }

    catch (AskapError& e) {
        ASKAPLOG_WARN_STR(logger, "Failure adding extra params");
        ASKAPLOG_WARN_STR(logger, "Exception detail: " << e.what());
    }
    catch (...) {
        ASKAPLOG_WARN_STR(logger, "Unknown exeption thrown in diadvise");
    }
    size_t beam = theBeams[0];
    // Iterate over all measurement sets
    // Lets sort out the output frames ...
    // iterate over the measurement sets and lets look at the
    // channels
    int id; // incoming rank ID
    int remainingWorkers = itsComms.nProcs() - 1;
    while(diadvise.getWorkUnitCount()>0 || remainingWorkers>0) {

        ContinuumWorkRequest wrequest;
        ASKAPLOG_DEBUG_STR(logger,"Waiting for a request " << diadvise.getWorkUnitCount() \
        << " units remaining");
        wrequest.receiveRequest(id, itsComms);
        ASKAPLOG_DEBUG_STR(logger,"Received a request from " << id);
        /// Now we can just pop a work allocation off the stack for this rank
        ContinuumWorkUnit wu = diadvise.getAllocation(id-1);
        ASKAPLOG_DEBUG_STR(logger,"Sending Allocation to  " << id);
        wu.sendUnit(id,itsComms);
        ASKAPLOG_DEBUG_STR(logger,"Sent Allocation to " << id);
        if (wu.get_payloadType() == ContinuumWorkUnit::DONE) {
            ASKAPLOG_INFO_STR(logger,"Sent DONE to " << id);
            remainingWorkers--;
        }

    }
    // all the work units allocated - lets send the DONEs
    // now finish the advice for remaining parameters
    diadvise.addMissingParameters(true);

    itsStats.logSummary();


    if (localSolver) {
        ASKAPLOG_INFO_STR(logger, "Master no longer required");
        return;
    }
    // this parset need to know direction and frequency for the final maps/models
    // But I dont want to run Cadvise as it is too specific to the old imaging requirements

    synthesis::ImagerParallel imager(itsComms, itsParset);
    // do a separate loop to build weights (in workers) if we are doing traditional weighting and build new weight grid
    imager.createUVWeightCalculator();
    if (imager.isSampleDensityGridNeeded()) {
        // MV: a bit of technical debt / waste of resources here. Technically we don't need the initial model for weight calculation
        // what we need is just the appropriate image names, their coordinate systems and shapes. Perhaps, in the future another communication
        // pattern can be added to send just the required information and don't bother with broadcasting the pixel array
        ASKAPLOG_DEBUG_STR(logger, "Master is about to broadcast initial model for uv-weight calculation");
        imager.broadcastModel();
        ASKAPLOG_DEBUG_STR(logger, "Master will merge weight grids from workers and compute final weights");
        imager.setupUVWeightBuilder();
        // cannot use receiveNE below, because it automatically assigns NE to the solver and this causes issues later
        // (wrong type of NE will be stuck inside the solver - caused by some technical debt in the design, we should've
        //  separated receiving NE and assigning it to the solver)
        imager.reduceNE(imager.getNE());
        // this will compute weights and add them to the model (which is distributed back to workers later on)
        imager.computeUVWeights();
        ASKAPLOG_DEBUG_STR(logger, "uv-weight has been added to the model");
        // reset normal equations back to the state suitable for imaging
        // (this is actually redundant in the case of the master, because calcNE later on would do this, but the code is
        // cleaner this way when we do it explicitly)
        imager.recreateNormalEquations();
        ASKAPLOG_DEBUG_STR(logger, "Master is about to broadcast the first model with uv-weights");
    } else {
        ASKAPLOG_DEBUG_STR(logger, "Master is about to broadcast the first <empty> model");
    }

    if (nCycles == 0) { // no solve if ncycles is 0
        ASKAPLOG_DEBUG_STR(logger, "Master beginning single - empty model");
        imager.broadcastModel(); // initially empty model

        imager.calcNE(); // Needed here because it resets the itsNE
        imager.receiveNE();
        imager.writeModel();
        itsStats.logSummary();

    }
    else {
        for (int cycle = 0; cycle < nCycles; ++cycle) {
            ASKAPLOG_DEBUG_STR(logger, "Master beginning major cycle ** " << cycle+1);

            if (cycle==0) {
                imager.broadcastModel(); // initially empty model (with or without uv-weights)
            }
            /// Minor Cycle

            imager.calcNE(); // Needed here because it resets the itsNE as Master
                            // Nothing else is done
            imager.solveNE(); /// Implicit receiveNE in here


            if (imager.params()->has("peak_residual")) {
                const double peak_residual = imager.params()->scalarValue("peak_residual");
                ASKAPLOG_INFO_STR(logger, "Major Cycle " << cycle+1 << " Reached peak residual of " << abs(peak_residual) << " after solve");

                if (peak_residual < targetPeakResidual) {

                    if (peak_residual < 0) {
                      ASKAPLOG_WARN_STR(logger, "Clean diverging, did not reach the major cycle threshold of "
                                      << targetPeakResidual << " Jy. Stopping.");
                    } else {
                      ASKAPLOG_INFO_STR(logger, "It is below the major cycle threshold of "
                                      << targetPeakResidual << " Jy. Stopping.");

                    }
                    ASKAPLOG_INFO_STR(logger, "Broadcasting final model");
                    imager.broadcastModel();
                    ASKAPLOG_INFO_STR(logger, "Broadcasting final model - done");
                    break;

                    // we have reached a peak residual after the

                } else {
                    if (targetPeakResidual < 0) {
                        ASKAPLOG_INFO_STR(logger, "Major cycle flux threshold is not used.");
                    } else {
                        if (imager.params()->has("noise_threshold_reached") &&
                            imager.params()->scalarValue("noise_threshold_reached")>0) {
                            ASKAPLOG_INFO_STR(logger, "It is below the noise threshold. Stopping.");
                            ASKAPLOG_INFO_STR(logger, "Broadcasting final model");
                            imager.broadcastModel();
                            ASKAPLOG_INFO_STR(logger, "Broadcasting final model - done");
                            break;
                        } else {
                        ASKAPLOG_INFO_STR(logger, "It is above the major cycle threshold of "
                                          << targetPeakResidual << " Jy. Continuing.");
                        }
                    }
                }
            }
            ASKAPLOG_INFO_STR(logger, "Broadcasting latest model");
            imager.broadcastModel();
            ASKAPLOG_INFO_STR(logger, "Broadcasting latest model - done");

            if (writeAtMajorCycle && (cycle != nCycles-1) ) {
                ASKAPLOG_INFO_STR(logger, "Writing out model");
                imager.writeModel(std::string(".beam") + utility::toString(beam) + \
                std::string(".majorcycle.") + utility::toString(cycle));
            }

            else {
                ASKAPLOG_DEBUG_STR(logger, "Not writing out model");
            }
            itsStats.logSummary();

        }
        ASKAPLOG_INFO_STR(logger, "Cycles complete - Receiving residuals for latest model");
        imager.calcNE(); // Needed here because it resets the itsNE as Master
                        // Nothing else is done
        imager.receiveNE(); // updates the residuals from workers
        ASKAPLOG_INFO_STR(logger, "Writing out model");
        imager.writeModel();
        itsStats.logSummary();

    }



}

// Utility function to get dataset names from parset.
std::vector<std::string> ContinuumMaster::getDatasets()
{
    if (itsParset.isDefined("dataset") && itsParset.isDefined("dataset0")) {
        ASKAPTHROW(std::runtime_error,
                   "Both dataset and dataset0 are specified in the parset");
    }

    // First look for "dataset" and if that does not exist try "dataset0"
    vector<string> ms;
    if (itsParset.isDefined("dataset")) {
        ms = itsParset.getStringVector("dataset", true);
    } else {
        string key = "dataset0";   // First key to look for
        long idx = 0;
        while (itsParset.isDefined(key)) {
            const string value = itsParset.getString(key);
            ms.push_back(value);

            ostringstream ss;
            ss << "dataset" << idx + 1;
            key = ss.str();
            ++idx;
        }
    }

    return ms;
}
// Utility function to get beam names from parset.
std::vector<int> ContinuumMaster::getBeams()
{
    std::vector<int> bs;

    if (itsParset.isDefined("beams")) {
        bs = itsParset.getInt32Vector("beams",bs);

    }
    else {
        bs.push_back(0);
    }
    return bs;
}
