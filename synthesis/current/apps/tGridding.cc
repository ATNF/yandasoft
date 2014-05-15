/// @file
///
/// @breif tests for gridding run in parallel
///
/// Performs synthesis imaging from a data source, using any of a number of
/// image solvers. Can run in serial or parallel (MPI) mode.
///
/// The data are accessed from the DataSource. This is and will probably remain
/// disk based. The images are kept purely in memory until the end.
///
/// Control parameters are passed in from a LOFAR ParameterSet file.
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

// Package level header file
#include "askap_synthesis.h"

// ASKAPsoft includes
#include "askap/AskapLogging.h"
#include "askap/AskapError.h"
#include <fitting/Params.h>
#include "askap/StatReporter.h"
#include <casa/Logging/LogIO.h>
#include <askap/Log4cxxLogSink.h>
#include <CommandLineParser.h>
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <gridding/VisGridderFactory.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <askap/AskapUtil.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>
#include <dataaccess/MemBufferDataAccessor.h>


ASKAP_LOGGER(logger, ".tGridding");

using namespace askap;
using namespace askap::synthesis;
using namespace askap::scimath;
using namespace askap::accessors;

// Main function
int main(int argc, const char** argv)
{
    // This class must have scope outside the main try/catch block
    askap::askapparallel::AskapParallel comms(argc, argv);

    try {
        // Ensure that CASA log messages are captured
        casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
        casa::LogSink::globalSink(globalSink);

        StatReporter stats;

        // Put everything in scope to ensure that all destructors are called
        // before the final message
        {
                    cmdlineparser::Parser parser; // a command line parser
            // command line parameter
            cmdlineparser::FlaggedParameter<std::string> inputsPar("-inputs",
                    "tgridding.in");
            // this parameter is optional
            parser.add(inputsPar, cmdlineparser::Parser::return_default);

            parser.process(argc, argv);

            const std::string parsetFile = inputsPar;

            LOFAR::ParameterSet parset(parsetFile);
            LOFAR::ParameterSet subset(parset.isDefined("Cimager.gridder") ? parset.makeSubset("Cimager.") : parset);
            
            ASKAPLOG_INFO_STR(logger, "Setting up the gridder to test and the model");
            IVisGridder::ShPtr gridder = VisGridderFactory::make(subset);
            ASKAPCHECK(gridder, "Gridder is not defined");
            scimath::Params model;
            boost::shared_ptr<scimath::Params> modelPtr(&model, utility::NullDeleter());
            SynthesisParamsHelper::setUpImages(modelPtr,subset.makeSubset("Images."));
            ASKAPLOG_INFO_STR(logger, "Model contains the following elements: "<<model);
            
            const int nCycles = subset.getInt32("ncycles", 1);
            ASKAPCHECK(nCycles > 0, "Number of iterations over the dataset is supposed to be positive, you have "<<nCycles);
            const std::string dataset = subset.getString("dataset");
            ASKAPLOG_INFO_STR(logger, "Dataset "<<dataset<<" will be used");
            const int cacheSize = subset.getInt32("nUVWMachines",1);
            ASKAPCHECK(cacheSize > 0, "uvw-machine cache size should be positive");
            const double cacheTolerance = SynthesisParamsHelper::convertQuantity(subset.getString("uvwMachineDirTolerance", 
                                                   "1e-6rad"),"rad");
            const int nCopies = subset.getInt32("ncopies",1);
            ASKAPCHECK(nCopies > 0, "number of copies should be positive");            
            ASKAPLOG_INFO_STR(logger, "Will run "<<nCopies<<" copies of gridding jobs");
            
            accessors::TableDataSource ds(dataset, accessors::TableDataSource::MEMORY_BUFFERS, "DATA");
            ds.configureUVWMachineCache(size_t(cacheSize),cacheTolerance);                          
            accessors::IDataSelectorPtr sel=ds.createSelector();
            sel << subset;
            accessors::IDataConverterPtr conv=ds.createConverter();
            conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO), "Hz");
            conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
            // ensure that time is counted in seconds since 0 MJD
            conv->setEpochFrame();
            
            ASKAPLOG_INFO_STR(logger, "Instantiating and initialising gridders");
            const size_t nImages = model.names().size();
            std::vector<IVisGridder::ShPtr> gridderList(nImages * nCopies);
            for (size_t i=0; i<gridderList.size(); ++i) {
                 gridderList[i] = gridder->clone();
                 ASKAPDEBUGASSERT(gridderList[i]);
                 const size_t modelIndex = i % nImages;
                 const std::string imageName = model.names()[modelIndex];
                 const Axes axes(model.axes(imageName));
                 gridderList[i]->initialiseGrid(axes,model.value(imageName).shape(), false);
                 gridderList[i]->customiseForContext(imageName);
            }
            for (int cycle = 0; cycle < nCycles; ++cycle) {
                 ASKAPLOG_INFO_STR(logger, "-------------- 'Major cycle' number "<<(cycle + 1)<< " -----------------");
                 accessors::IDataSharedIter it=ds.createIterator(sel, conv);
                 size_t counterGrid = 0;
                 for (it.init();it.hasMore();it.next()) {
                      accessors::MemBufferDataAccessor accBuffer(*it);
                      accBuffer.rwVisibility().set(0.);
                      accBuffer.rwVisibility() -= it->visibility();
                      accBuffer.rwVisibility() *= float(-1.);
                      size_t tempCounter = 0; 
                      #ifdef _OPENMP
                      #pragma omp parallel default(shared)
                      {
                         #pragma omp for reduction(+:tempCounter)
                      #endif
                         for (size_t i = 0; i<gridderList.size(); ++i) {
                              gridderList[i]->grid(accBuffer);
                              tempCounter += accBuffer.nRow();
                         }
                      #ifdef _OPENMP
                      }
                      #endif
                      counterGrid += tempCounter;
                 }
                 ASKAPLOG_INFO_STR(logger, "Finished gridding pass, number of rows gridded is "<<counterGrid);                 
            }
            
        }
        stats.logSummary();
        ///==============================================================================
    } catch (const cmdlineparser::XParser &ex) {
        ASKAPLOG_FATAL_STR(logger, "Command line parser error, wrong arguments " << argv[0]);
        std::cerr << "Usage: " << argv[0] << " [-inputs parsetFile]"
                      << std::endl;
    } catch (const askap::AskapError& x) {
        ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
        std::cerr << "Askap error in " << argv[0] << ": " << x.what()
                      << std::endl;
        exit(1);
    } catch (const std::exception& x) {
        ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
        std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what()
                      << std::endl;
        exit(1);
    }

    return 0;
}
        