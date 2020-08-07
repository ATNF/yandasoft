/// @file cdeconvolver.cc
///
/// @brief Image deconvolution program
///
/// Control parameters are passed in from a LOFAR ParameterSet file.
///
/// @copyright (c) 2007, 2020 CSIRO
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
/// @author Stephen Ord <stephen.ord@csiro.au>

// Package level header file
#include <askap/askap_synthesis.h>

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
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/deconvolution/DeconvolverBase.h>
#include <askap/deconvolution/DeconvolverFactory.h>
#include <askap/deconvolution/DeconvolverHelpers.h>
#include <askap/distributedimager/CubeBuilder.h>
#include <askap/askapparallel/AskapParallel.h>

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
    
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsModelCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsResidualCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsRestoredCube;
    
        
        virtual int run(int argc, char* argv[])
        {
            StatReporter stats;

            const LOFAR::ParameterSet subset(config().makeSubset("Cdeconvolver."));
            
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));
            
            ASKAPLOG_INFO_STR(logger, "ASKAP image (MPI) deconvolver " << ASKAP_PACKAGE_VERSION);
            // Need some metadata for the output cube constructions
            
            // Lets get the grid,pcf and psf cube names from the parset
            
            const std::string gridCubeName = subset.getString("grid");
            const std::string pcfCubeName = subset.getString("pcf");
            const std::string psfCubeName = subset.getString("psf");
            
            // Lets check the dimensions
            // initialise an image accessor
            // set up image handler, needed for both master and worker
            SynthesisParamsHelper::setUpImageHandler(subset);
            accessors::IImageAccess<casacore::Float>& iacc = SynthesisParamsHelper::imageHandler();
            
            // Lets load in a cube
            casacore::PagedImage<casacore::Complex> grid(gridCubeName);
            const casa::IPosition shape = grid.shape();
            
            // lets calculate the allocations ...

            casa::IPosition blc(shape.nelements(),0);
            casa::IPosition trc(shape);
            int nchanCube = trc[3];
            
            vector<IPosition> inShapeVec;
            vector<CoordinateSystem> inCoordSysVec;
            // What fraction of the full problem does this rank
            // THe units here are channels - not polarisations.
            int myFullAllocationSize = 0;
            int myFullAllocationStart = 0;
            int myFullAllocationStop = 0;
            // Where a rank is in its allocation
            int myAllocationSize = 1;
            int myAllocationStart = 0;
            int myAllocationStop = myAllocationStart + myAllocationSize;
            
            if (comms.rank() >= nchanCube) {
              ASKAPLOG_WARN_STR(logger,"Rank " << comms.rank() << " has no work to merge");
                return -1;
            }
            if (nchanCube % comms.nProcs() != 0) {
              ASKAPLOG_WARN_STR(logger,"Unbalanced allocation: num of ranks:" << comms.nProcs() << " not a factor of number of channels: "<< nchanCube);
            }
            if (comms.nProcs() >= nchanCube) {
              myFullAllocationSize = 1;
            }
            else {
              myFullAllocationSize = nchanCube/comms.nProcs();
            }
            myFullAllocationStart = comms.rank()*myFullAllocationSize;
            myFullAllocationStop = myFullAllocationStart + myFullAllocationSize;

            // unless last rank
            if (comms.rank() == comms.nProcs()-1) {
              myFullAllocationSize = trc[3] - myFullAllocationStart; // we are using End is Last
            }

            ASKAPLOG_INFO_STR(logger,"FullAllocation starts at " << myFullAllocationStart << " and is " << myFullAllocationSize << " in size");
            // So the plan is to iterate over each channel ...
            // Calculate the inShapes for each channel and file ....

            for (myAllocationStart = myFullAllocationStart; myAllocationStart < myFullAllocationStop; myAllocationStart = myAllocationStart +  myAllocationSize) {
                
                ASKAPLOG_INFO_STR(logger,"Processing Channel " << myAllocationStart);
            }

            // output Model cube
            // this->itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, freqinc, "model"));
            // output Residual cube
            // this->itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, freqinc,"residual"));
            // boost::shared_ptr<DeconvolverBase<casacore::Float, casacore::Complex> > deconvolver(DeconvolverFactory::make(subset));
            
            // deconvolver->deconvolve();

            // Now write the model and residual to disk using the names specified in the 
            // parset. We simply copy the dirty image and then write the array into 
            // the resulting image. 
            // DeconvolverHelpers::putArrayToImage(deconvolver->model(), "model", "dirty", subset);
            // DeconvolverHelpers::putArrayToImage(deconvolver->dirty(), "residual", "dirty", subset);

            // Vector<Array<float> > restored(1);
            // if(deconvolver->restore(restored)) {
            //    DeconvolverHelpers::putArrayToImage(restored(0), "restored", "dirty", subset);
            // }

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
