/// @file cdeconvolver.cc
///
/// @brief Image deconvolution program
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>

// Package level header file
#include <askap/askap_synthesis.h>

// System includes
#include <stdexcept>
#include <iostream>
#include <string>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".cdeconvolver");
#include <askap/askap/AskapError.h>
#include <askap/askap/Application.h>
#include <askap/askap/StatReporter.h>
#include <askap/deconvolution/DeconvolverBase.h>
#include <askap/deconvolution/DeconvolverFactory.h>
#include <askap/deconvolution/DeconvolverHelpers.h>

#include <askap/parallel/AdviseParallel.h>


using namespace askap;
using namespace askap::synthesis;

class CdeconvolverApp : public askap::Application
{
    public:
        int run(int argc, char* argv[]) final
        {
            StatReporter stats;

            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                const LOFAR::ParameterSet subset(config().makeSubset("Cdeconvolver."));

                ASKAPLOG_INFO_STR(logger, "ASKAP image deconvolver " << ASKAP_PACKAGE_VERSION);

                boost::shared_ptr<DeconvolverBase<Float, Complex> > deconvolver(DeconvolverFactory::make(subset));
                deconvolver->deconvolve();

                // Now write the model and residual to disk using the names specified in the
                // parset. We simply copy the dirty image and then write the array into
                // the resulting image.
                DeconvolverHelpers::putArrayToImage(deconvolver->model(), "model", "dirty", subset);
                DeconvolverHelpers::putArrayToImage(deconvolver->dirty(), "residual", "dirty", subset);

                Vector<Array<float> > restored(1);
                if(deconvolver->restore(restored)) {
                    DeconvolverHelpers::putArrayToImage(restored(0), "restored", "dirty", subset);
                }
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
    private:
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:" + ASKAP_PACKAGE_VERSION);
            return pkgVersion;
        }
};

// Main function
int main(int argc, char* argv[])
{
    CdeconvolverApp app;
    return app.main(argc, argv);
}
