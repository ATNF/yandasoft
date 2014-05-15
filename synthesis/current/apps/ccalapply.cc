/// @file ccalapply.cc
///
/// @brief Calibration applicator
///
/// @copyright (c) 2013 CSIRO
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

// Include package level header file
#include "askap_synthesis.h"

// System include
#include <string>
#include <iostream>
#include <stdlib.h>

// ASKAPsoft includes
#include "askap/Application.h"
#include "askap/AskapError.h"
#include "askap/AskapLogging.h"
#include "askap/StatReporter.h"
#include "dataaccess/TableDataSource.h"
#include "dataaccess/IConstDataSource.h"
#include "dataaccess/SharedIter.h"
#include "dataaccess/ParsetInterface.h"
#include "calibaccess/ICalSolutionConstSource.h"
#include "calibaccess/CalibAccessFactory.h"
#include "dataaccess/OnDemandNoiseAndFlagDA.h"
#include "boost/shared_ptr.hpp"

// Local packages includes
#include "measurementequation/ICalibrationApplicator.h"
#include "measurementequation/CalibrationApplicatorME.h"
#include "measurementequation/CalibrationIterator.h"

ASKAP_LOGGER(logger, ".ccalapply");

// Using
using namespace std;
using namespace askap;
using namespace askap::accessors;
using namespace askap::synthesis;

class CcalApplyApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            try {
                StatReporter stats;
                LOFAR::ParameterSet subset(config().makeSubset("Ccalapply."));

                // Setup calibration applicator
                boost::shared_ptr<ICalibrationApplicator> calME = buildCalApplicator(subset);

                // Get Measurement Set accessor
                IDataSharedIter it = getDataIterator(subset);
                
                ASKAPDEBUGASSERT(it);
                ASKAPDEBUGASSERT(calME);

                // Apply calibration
                uint64_t count = 1;
                for (it.init(); it != it.end(); it.next()) {
                    if (count % 100 == 0) {
                        ASKAPLOG_DEBUG_STR(logger, "Progress - Chunk " << count);
                    }
                    ++count;
                    if (itsNoiseAndFlagDANeeded) {
                        // quick and dirty for now
                        accessors::OnDemandNoiseAndFlagDA acc(*it);
                        acc.rwVisibility() = it->visibility();

                        calME->correct(acc);

                        it->rwVisibility() = acc.rwVisibility().copy();
                    } else {
                        calME->correct(*it);
                    }
                }

                stats.logSummary();
            } catch (const askap::AskapError& x) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
                std::cerr << "Askap error in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            } catch (const std::exception& x) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            }

            return 0;
        }

    private:
        
        /// @brief this flag indicates whether the underlying code needs to update flags or noise
        /// @details For now - quick and dirty fix to allow MRO tests to proceed. The ASKAP model
        /// is to apply calibration on the fly, so table accessor classes are unable to modify noise or flag
        /// information. 
        bool itsNoiseAndFlagDANeeded;
    
        static casa::MFrequency::Ref getFreqRefFrame(const LOFAR::ParameterSet& parset)
        {
            const string freqFrame = parset.getString("freqframe", "topo");
            if (freqFrame == "topo") {
                ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as topocentric");
                return casa::MFrequency::Ref(casa::MFrequency::TOPO);
            } else if (freqFrame == "lsrk") {
                ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as lsrk");
                return casa::MFrequency::Ref(casa::MFrequency::LSRK);
            } else if (freqFrame == "bary") {
                ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as barycentric");
                return casa::MFrequency::Ref(casa::MFrequency::BARY);
            } else {
                ASKAPTHROW(AskapError, "Unsupported frequency frame " << freqFrame);
            }
        }

        boost::shared_ptr<ICalibrationApplicator> buildCalApplicator(
                const LOFAR::ParameterSet& parset)
        {
            // Create solution source
            ICalSolutionConstSource::ShPtr solutionSource =
                CalibAccessFactory::roCalSolutionSource(parset);
            ASKAPASSERT(solutionSource);

            // Create applicator
            boost::shared_ptr<ICalibrationApplicator> calME(new CalibrationApplicatorME(solutionSource));
            ASKAPASSERT(calME);
            const bool scaleNoise = parset.getBool("calibrate.scalenoise", false);
            const bool allowFlag = parset.getBool("calibrate.allowflag", false);
            itsNoiseAndFlagDANeeded = scaleNoise || allowFlag;  
            calME->scaleNoise(scaleNoise);
            calME->allowFlag(allowFlag);
            calME->beamIndependent(parset.getBool("calibrate.ignorebeam", false));
            return calME;
        }

        static IDataSharedIter getDataIterator(const LOFAR::ParameterSet& parset)
        {
            const string ms = parset.getString("dataset");
            TableDataSource ds(ms);

            IDataSelectorPtr sel=ds.createSelector();
            sel << parset;
            IDataConverterPtr conv=ds.createConverter();
            conv->setFrequencyFrame(getFreqRefFrame(parset), "Hz");
            conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
            // ensure that time is counted in seconds since 0 MJD
            conv->setEpochFrame();

            return ds.createIterator(sel, conv);
        }
};

int main(int argc, char *argv[])
{
    CcalApplyApp app;
    return app.main(argc, argv);
}
