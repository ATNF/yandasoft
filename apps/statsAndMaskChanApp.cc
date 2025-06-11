/// @file
///
/// Demo program to mask out problematic channels of a cube based on noise statistics
///
/// @copyright (c) 2019 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 3 of the License,
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
/// @author Minh Vuong <Minh.Vuong@csiro.au>
///


// Package level header file
#include "askap/askap_synthesis.h"

// System includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdio.h>

// ASKAPsoft includes
#include <askap/askap/Application.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/StatReporter.h>
#include <askap/askapparallel/AskapParallel.h>

#include <askap/utils/StatsAndMask.h>

#include <askap/imageaccess/ImageAccessFactory.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicMath/Math.h>

ASKAP_LOGGER(logger, ".statsAndMaskChanApp");

using namespace askap;

///
/// This app does the followings :
/// (1) calculates the image's statistics
/// (2) masks the bad channels
/// (3) records the image statistics
///
class StatsAndMaskChanApp : public askap::Application
{
public:
    int run(int argc, char* argv[]) final
    {
        // This class must have scope outside the main try/catch block
        askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

        try {
            StatReporter stats;
            LOFAR::ParameterSet subset(config().makeSubset("MaskChannels."));

            std::string statsfile = subset.getString("statsFile","");
            //ASKAPLOG_INFO_STR(logger, "Reading stats file " << statsfile);

            std::string image = subset.getString("image","");
            bool editImage = subset.getBool("editImage", false);
            bool editStats = subset.getBool("editStats", false);
            std::string outputStats = subset.getString("outputStats","");
            bool useSignificance = subset.getBool("useSignificance", false);
            float threshold = subset.getFloat("significanceLevel", 10.);
            bool useNoise = subset.getBool("useNoise", true);
            float badMADFM = subset.getFloat("madfmThreshold", 1000.);
            bool maskBlank = subset.getBool("maskBlank", true);

            boost::shared_ptr<accessors::IImageAccess<casacore::Float> > iacc = accessors::imageAccessFactory(subset);
            casa::IPosition shape = iacc->shape(image);
            casa::CoordinateSystem coo = iacc->coordSys(image);
            int specAxis = coo.spectralAxisNumber();
            const auto numberOfChansInImage = shape[specAxis];
            askap::utils::StatsAndMask statisticsAndMask(comms,image,iacc);
            const int rank = comms.rank();
            const int nProcs = comms.nProcs();
            if ( rank == 0 ) {
                ASKAPLOG_INFO_STR(logger, "Shape of input image = " << shape);
                ASKAPLOG_INFO_STR(logger, "Shape size = " << shape.size());
                ASKAPLOG_INFO_STR(logger, "Spectral axis = " << specAxis);
                ASKAPLOG_INFO_STR(logger, "Size of spectral axis (i.e number of channels) = " << numberOfChansInImage);
            }
            for (unsigned int chan = 0; chan < numberOfChansInImage; chan++) {
                //if ( chan % numberOfChansInImage == rank ) {
                if ( chan % nProcs == rank ) {
                    // ASKAPLOG_INFO_STR(logger,"rank: " << rank << " running on node: " << comms.nodeName());
                    ASKAPLOG_INFO_STR(logger,"rank: " << rank << ", channel: " << chan);
                    casacore::IPosition blc(4,0,0,0,chan);
                    casacore::IPosition trc = shape - 1;
                    trc(3) = chan;
                    statisticsAndMask.calculate(blc,trc);
                }
            }
            // waits for all ranks to get here
            comms.barrier();
            //ASKAPLOG_INFO_STR(logger,"size of itsStatsPerChannelMap: " << statisticsAndMask.statsPerChannelMap().size());
            if (rank == 0 ) {
                // master collects the stats from the workers
                statisticsAndMask.receiveStats();
                ASKAPLOG_INFO_STR(logger,"size of itsStatsPerChannelMap: " << statisticsAndMask.statsPerChannelMap().size());

                // do the masking
                statisticsAndMask.maskBadChannels(image,threshold,badMADFM,maskBlank, useSignificance,
                                                  useNoise, editStats, editImage, outputStats);
                if ( editImage ) {
                    //write the stats to image table
                    statisticsAndMask.writeStatsToImageTable(image);
                }
            } else {
                // workers send the stats to the master
                statisticsAndMask.sendStats();
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
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};




int main(int argc, char *argv[])
{
    StatsAndMaskChanApp app;
    return app.main(argc, argv);
}
