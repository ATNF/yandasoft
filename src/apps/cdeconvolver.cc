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
#include <askap_synthesis.h>

// System includes
#include <stdexcept>
#include <iostream>
#include <string>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".cdeconvolver");
#include <askap/AskapError.h>
#include <askap/Application.h>
#include <askap/StatReporter.h>
#include <deconvolution/DeconvolverBase.h>
#include <deconvolution/DeconvolverFactory.h>
#include <deconvolution/DeconvolverHelpers.h>

using namespace askap;
using namespace askap::synthesis;

static std::string getNodeName(void)
{
    const int HOST_NAME_MAXLEN = 256;
    char name[HOST_NAME_MAXLEN];
    gethostname(name, HOST_NAME_MAXLEN);
    std::string nodename(name);

    std::string::size_type idx = nodename.find_first_of('.');
    if (idx != std::string::npos) {
        // Extract just the hostname part
        nodename = nodename.substr(0, idx);
    }

    return nodename;
}

class CdeconvolverApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            StatReporter stats;

            const LOFAR::ParameterSet subset(config().makeSubset("Cdeconvolver."));

            const std::string hostname = getNodeName();
            ASKAPLOG_REMOVECONTEXT("hostname");
            ASKAPLOG_PUTCONTEXT("hostname", hostname.c_str());

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

            stats.logSummary();
            return 0;
        }
};

// Main function
int main(int argc, char* argv[])
{
    CdeconvolverApp app;
    return app.main(argc, argv);
}
