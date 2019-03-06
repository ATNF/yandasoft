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
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askapparallel/AskapParallel.h>

#include <Common/ParameterSet.h>
#include <fitting/Params.h>
#include <fitting/Axes.h>
#include <dataaccess/IConstDataSource.h>
#include <dataaccess/TableConstDataSource.h>
#include <dataaccess/IConstDataIterator.h>
#include <dataaccess/IDataConverter.h>
#include <dataaccess/IDataSelector.h>
#include <dataaccess/IDataIterator.h>
#include <dataaccess/SharedIter.h>
#include <dataaccess/TableInfoAccessor.h>
#include <casacore/casa/Quanta.h>
#include <imageaccess/BeamLogger.h>
#include <parallel/ImagerParallel.h>
#include <measurementequation/SynthesisParamsHelper.h>

// Local includes
#include "distributedimager/AdviseDI.h"
#include "distributedimager/CalcCore.h"
#include "distributedimager/CubeComms.h"
#include "messages/ContinuumWorkUnit.h"
#include "messages/ContinuumWorkRequest.h"


//casacore includes
#include "casacore/ms/MeasurementSets/MeasurementSet.h"
#include "casacore/ms/MeasurementSets/MSColumns.h"

using namespace std;
using namespace askap::cp;
using namespace askap;

ASKAP_LOGGER(logger, ".ContinuumMaster");

ContinuumMaster::ContinuumMaster(LOFAR::ParameterSet& parset,
                                       CubeComms& comms)
    : itsParset(parset), itsComms(comms), itsBeamList()
{
}

ContinuumMaster::~ContinuumMaster()
{
}

void ContinuumMaster::run(void)
{
    // Read from the configruation the list of datasets to process
    const vector<string> ms = getDatasets(itsParset);
    if (ms.size() == 0) {
        ASKAPTHROW(std::runtime_error, "No datasets specified in the parameter set file");
    }
    // Need to break these measurement sets into groups
    // there are three posibilties:
    // 1 - the different measurement sets have the same epoch - but different
    //      frequencies
    // 2 - they have different epochs but the same TOPO centric frequencies

    vector<int> theBeams = getBeams();
    int totalChannels = 0;

    const double targetPeakResidual = synthesis::SynthesisParamsHelper::convertQuantity(
                itsParset.getString("threshold.majorcycle", "-1Jy"), "Jy");







    LOFAR::ParameterSet unitParset = itsParset;



    const bool writeAtMajorCycle = unitParset.getBool("Images.writeAtMajorCycle", false);
    const int nCycles = unitParset.getInt32("ncycles", 0);
    const bool localSolver = unitParset.getBool("solverpercore",false);
    synthesis::AdviseDI diadvise(itsComms,unitParset);

    try {

        diadvise.prepare();
        diadvise.addMissingParameters();

        ASKAPLOG_DEBUG_STR(logger,"*****");
        ASKAPLOG_DEBUG_STR(logger,"Parset" << diadvise.getParset());
        ASKAPLOG_DEBUG_STR(logger,"*****");

        totalChannels = diadvise.getBaryFrequencies().size();

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

    while(diadvise.getWorkUnitCount()) {

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
    }


    if (localSolver) {
        ASKAPLOG_INFO_STR(logger, "Master no longer required");
        return;
    }
    ASKAPLOG_DEBUG_STR(logger, "Master is about to broadcast first <empty> model");

    // this parset need to know direction and frequency for the final maps/models
    // But I dont want to run Cadvise as it is too specific to the old imaging requirements


    if (nCycles == 0) { // no solve if ncycles is 0
        synthesis::ImagerParallel imager(itsComms, diadvise.getParset());
        ASKAPLOG_DEBUG_STR(logger, "Master beginning single - empty model");
        imager.broadcastModel(); // initially empty model

        imager.calcNE(); // Needed here becuase it resets the itsNE
        imager.receiveNE();
        imager.writeModel();

    }
    else {
        synthesis::ImagerParallel imager(itsComms, diadvise.getParset());
        for (int cycle = 0; cycle < nCycles; ++cycle) {
            ASKAPLOG_DEBUG_STR(logger, "Master beginning major cycle ** " << cycle);

            if (cycle==0) {
                imager.broadcastModel(); // initially empty model
            }
            /// Minor Cycle

            imager.calcNE(); // Needed here becuase it resets the itsNE as Master
                            // Nothing else is done
            imager.solveNE(); /// Implicit receiveNE in here


            if (imager.params()->has("peak_residual")) {
                const double peak_residual = imager.params()->scalarValue("peak_residual");
                ASKAPLOG_INFO_STR(logger, "Major Cycle " << cycle << " Reached peak residual of " << peak_residual << " after solve");

                if (peak_residual < targetPeakResidual) {

                    ASKAPLOG_INFO_STR(logger, "It is below the major cycle threshold of "
                                      << targetPeakResidual << " Jy. Stopping.");

                    ASKAPLOG_INFO_STR(logger, "Broadcasting final model");
                    imager.broadcastModel();
                    ASKAPLOG_INFO_STR(logger, "Broadcasting final model - done");
                    break;

                    // we have reached a peak residual after the

                } else {
                    if (targetPeakResidual < 0) {
                        ASKAPLOG_INFO_STR(logger, "Major cycle flux threshold is not used.");
                    } else {
                        ASKAPLOG_INFO_STR(logger, "It is above the major cycle threshold of "
                                          << targetPeakResidual << " Jy. Continuing.");
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

        }
        ASKAPLOG_INFO_STR(logger, "Cycles complete - Receiving residuals for latest model");
        imager.calcNE(); // Needed here becuase it resets the itsNE as Master
                        // Nothing else is done
        imager.receiveNE(); // updates the residuals from workers
        ASKAPLOG_INFO_STR(logger, "Writing out model");
        imager.writeModel();

    }



}

// Utility function to get dataset names from parset.
std::vector<std::string> ContinuumMaster::getDatasets(const LOFAR::ParameterSet& parset)
{
    if (parset.isDefined("dataset") && parset.isDefined("dataset0")) {
        ASKAPTHROW(std::runtime_error,
                   "Both dataset and dataset0 are specified in the parset");
    }

    // First look for "dataset" and if that does not exist try "dataset0"
    vector<string> ms;
    if (parset.isDefined("dataset")) {
        ms = itsParset.getStringVector("dataset", true);
    } else {
        string key = "dataset0";   // First key to look for
        long idx = 0;
        while (parset.isDefined(key)) {
            const string value = parset.getString(key);
            ms.push_back(value);

            ostringstream ss;
            ss << "dataset" << idx + 1;
            key = ss.str();
            ++idx;
        }
    }

    return ms;
}
// Utility function to get dataset names from parset.
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
