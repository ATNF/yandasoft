/// @file
///
/// @brief Generic operations-specific calibration
/// @details This application is intended for the types of calibration which 
/// cannot follow ASKAP's predict-forward approach, i.e. which require 
/// observations of various fields done in some special way. It is intended
/// for experimentation with calibration as well as some operation-specific
/// tasks like baseline and pointing refinements.  
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>
// Package level header file
#include <askap_synthesis.h>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/Application.h>
#include <askap/StatReporter.h>
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>

// local include
#include <measurementequation/ScanStats.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>



ASKAP_LOGGER(logger, ".opcal");

using namespace std;
using namespace askap;
using namespace askap::synthesis;

class OpCalApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            StatReporter stats;

            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {

                 ASKAPLOG_INFO_STR(logger, "ASKAP synthesis opcal application " << ASKAP_PACKAGE_VERSION);
               
                 if (comms.isMaster()) {
                   ASKAPLOG_INFO_STR(logger, "Parset file contents:\n" << config());
                 }
                 ASKAPCHECK(!comms.isParallel(), "This application is not intended to be used in parallel mode (at this stage)");
                 
                 ScanStats ss(3600.);
                 
                 const std::string ms = config().getString("dataset");
                 accessors::TableDataSource ds(ms, accessors::TableDataSource::MEMORY_BUFFERS/*, dataColumn()*/);
                 accessors::IDataSelectorPtr sel=ds.createSelector();
                 sel << config();
                 accessors::IDataConverterPtr conv=ds.createConverter();
                 conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO)/*getFreqRefFrame()*/, "Hz");
                 conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
                 conv->setEpochFrame(); // time in seconds since 0 MJD
                 accessors::IDataSharedIter it=ds.createIterator(sel, conv);
                 ASKAPLOG_INFO_STR(logger, "Inspecting "<<ms);
                 ss.inspect(ms, it);
                 
                 ASKAPLOG_INFO_STR(logger, "Found "<<ss.size()<<" chunks in the supplied data");                 
                 
            } catch (const askap::AskapError& e) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << e.what());
                std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
                exit(1);
            } catch (const std::exception& e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << e.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what() << std::endl;
                exit(1);
            }
            stats.logSummary();
            return 0;
        }
};

// Main function
int main(int argc, char* argv[])
{
    OpCalApp app;
    return app.main(argc, argv);
}
            