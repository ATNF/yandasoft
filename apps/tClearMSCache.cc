/// @file tClearMSCache.cc
///
/// @breif A tool to exercise features of the DataSourceManager class in parallel / at scale
/// @details
/// Although this class largely exists to encapsulate the lifecycle of DataSource class (i.e. 
/// access to the measurement sets), one non-trivial bit of functionality it provides is related to
/// clearing caches of data and flag column storage managers via calling the appropriate low-level methods. 
/// It is not clear to me (MV) if this is needed at all (there could be a reason for this earlier, but now
/// the DataSource objects are not left unnecessarily, so the cleanup should happen automatically). Anyway, 
/// this test can be used to exercise this functionality (clearing caches is enabled by default) along with
/// other benchmarking. To ensure the test is realistic, the code performs reading of the measurement sets
/// and calculates the sum of all unflagged data. The data selection via parset is supported to allow 
/// simulation of parallel access to the same or different measurement sets.
///
/// Control parameters are passed in from a LOFAR ParameterSet file in a standard way for all yandasoft apps.
///
/// @copyright (c) 2024 CSIRO
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
#include <askap/askap_synthesis.h>

// System includes
#include <stdexcept>
#include <iostream>

// boost includes
#include <boost/shared_ptr.hpp>

// ASKAPsoft includes
#include <askap/askap/Application.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/StatReporter.h>
#include <askap/distributedimager/DataSourceManager.h>
#include <askap/askapparallel/AskapParallel.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/dataaccess/SharedIter.h>


//#include <askap/measurementequation/MEParsetInterface.h>


ASKAP_LOGGER(logger, ".tClearMSCache");

using namespace askap;
using namespace askap::synthesis;

class DataSourceManagerTestApp : public askap::Application
{
    public:

        void doTest(const std::string& dataset, const LOFAR::ParameterSet& parset, StatReporter &stats)
        {
            ASKAPDEBUGASSERT(itsDSM);
            ASKAPLOG_INFO_STR(logger, "About to perform the test with "<<dataset);
            stats.logSummary();
            accessors::TableDataSource& ds = itsDSM->dataSource(dataset);
            accessors::IDataSelectorPtr sel = ds.createSelector();
            sel << parset;
            accessors::IDataConverterPtr conv = ds.createConverter();
            conv->setFrequencyFrame(casacore::MFrequency::Ref(casa::MFrequency::TOPO),"Hz");
            conv->setEpochFrame(casacore::MEpoch(casa::Quantity(53635.5,"d"),
                                casacore::MEpoch::Ref(casacore::MEpoch::UTC)),"s");
            conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));                    
            casacore::Complex sumBuffer(0.f);
            size_t counter = 0u;
            for (accessors::IConstDataSharedIter it=ds.createConstIterator(sel,conv);it!=it.end();++it) {  
                 const accessors::IConstDataAccessor& acc = *it;
                 const casacore::Cube<casacore::Complex>& vis = acc.visibility();
                 const casacore::Cube<casacore::Bool>& flag = acc.flag();
                 for (casacore::uInt row = 0; row < acc.nRow(); ++row) {
                      for (casacore::uInt chan = 0; chan < acc.nChannel(); ++chan) {
                           for (casacore::uInt pol = 0; pol < acc.nPol(); ++pol) {
                                if (!flag(pol, chan, row)) {
                                    sumBuffer += vis(pol, chan, row);
                                    ++counter;
                                }
                           }
                      }
                 }
            }
            ASKAPLOG_INFO_STR(logger, "Finished iteration, sum = "<<sumBuffer<<" number of unflagged visibilities: "<<counter);
            if (counter > 0u) {
                ASKAPLOG_INFO_STR(logger, "      mean visibility: "<<sumBuffer / float(counter));
            }
            stats.logSummary();
            // explicit disposal of the old data source, if configured
            if (parset.getBool("cleardsfirst", false)) {
                ASKAPLOG_INFO_STR(logger, "Disposing of the old data source object");
                itsDSM->forceNewDataSourceNextTime();
            }
            // explcit reset of the manager would cause the cleanup of caches (if enabled) and disposal of the old data source
            // otherwise it will happen at the very end when this class is destroyed
            if (parset.getBool("dsmreset", true)) {
                ASKAPLOG_INFO_STR(logger, "Resetting data source manager");
                itsDSM->reset();
            }
            stats.logSummary();
        }

        int run(int argc, char* argv[]) final
        {
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;
                // in principle, we can remove Cimager prefix, but having it simplifies copying
                // from real imager parsets
                //LOFAR::ParameterSet subset(config().makeSubset("Cimager."));
                LOFAR::ParameterSet subset(config());

                // imager-specific configuration of the master/worker to allow groups of workers
                // in principle, we don't need to set aside a master rank in the parallel mode but it is handy to
                // reuse existing imager infrastructure and do substitutions in a consistent way
                const int nWorkerGroups = subset.getInt32("nworkergroups", 1);
                ASKAPCHECK(nWorkerGroups > 0, "nworkergroups is supposed to be greater than 0");
                if (nWorkerGroups > 1) {
                    ASKAPLOG_INFO_STR(logger, "There are "<<nWorkerGroups<<
                            " groups of workers, (each) measurement set will be read multiple times");
                    ASKAPCHECK(comms.isParallel(), "This option is only allowed in the parallel mode");
                    comms.defineGroups(nWorkerGroups);
                } else {
                    ASKAPLOG_INFO_STR(logger, "All workers are treated as identical");
                }

                if (comms.isWorker()) {

                    // Perform %w substitutions for all keys.
                    // NOTE: This MUST happen after AskapParallel::defineGroups() is called
                    for (LOFAR::ParameterSet::iterator it = subset.begin();
                            it != subset.end(); ++it) {
                         it->second = LOFAR::ParameterValue(comms.substitute(it->second));
                    }
                    // get effective worker number via substitution to avoid reimplementing that logic,
                    // although string to integer conversion looks ugly
                    const int workerSeqNumber = comms.isParallel() ? utility::fromString<int>(comms.substitute("\%w")) : 0;
                    ASKAPCHECK(workerSeqNumber >= 0, "For some reason, substitution of \%w returned a negative number");
     
                    const std::vector<std::string> datasets = subset.getStringVector("dataset");
                    ASKAPCHECK(datasets.size() > 0, "At least one dataset should be given, not an empty vector");

                    // setup data source manager, get parameters from parset which are required (column name, caching, etc)
                    const int uvwMachineCacheSize = subset.getInt32("nUVWMachines", 1);
                    ASKAPCHECK(uvwMachineCacheSize > 0 ,
                               "Cache size is supposed to be a positive number, you have "
                               << uvwMachineCacheSize);

                    const double uvwMachineCacheTolerance = SynthesisParamsHelper::convertQuantity(subset.getString("uvwMachineDirTolerance", "1e-6rad"), "rad");

                    ASKAPLOG_DEBUG_STR(logger,
                             "UVWMachine cache will store " << uvwMachineCacheSize << " machines");
                    ASKAPLOG_DEBUG_STR(logger, "Tolerance on the directions is "
                             << uvwMachineCacheTolerance / casacore::C::pi * 180. * 3600. << " arcsec");

                    const string colName = subset.getString("datacolumn", "DATA");
                    const bool clearcache = subset.getBool("clearcache", true);

                    itsDSM.reset(new DataSourceManager(colName, clearcache, static_cast<size_t>(uvwMachineCacheSize), uvwMachineCacheTolerance));
                    //

                    const bool distributeDatasets = subset.getBool("distributedatasets", true);
                    if (distributeDatasets && (datasets.size() > 1)) {
                        ASKAPCHECK(logger, "Datasets will be distributed between worker ranks");
                        const std::string dataset = datasets[workerSeqNumber % datasets.size()];
                        doTest(dataset, subset, stats);
                    } else {
                        ASKAPCHECK(logger, "Each worker will read all datasets in a sequence");
                        for (size_t index = 0; index < datasets.size(); ++index) {
                             doTest(datasets[index], subset, stats);
                        }
                    }

                    ASKAPLOG_INFO_STR(logger, "ASKAP tClearMSCache test utility " << ASKAP_PACKAGE_VERSION);
                } else {
                    ASKAPLOG_INFO_STR(logger, "Nothing to do for master in the parallel mode - finishing");
                }
                itsDSM.reset();
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
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
        
        /// @brief shared pointer to the manager object 
        boost::shared_ptr<DataSourceManager> itsDSM;
};

int main(int argc, char *argv[])
{
    DataSourceManagerTestApp app;
    return app.main(argc, argv);
}
