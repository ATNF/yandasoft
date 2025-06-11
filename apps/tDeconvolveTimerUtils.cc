/// @file tDeconvolveTimerUtils.cc
///
/// @brief combine a number of images as a linear mosaic
/// @details This is a utility to merge images into a mosaic. Images can be set
/// explicitly or found automatically based on input tags.
///
/// @copyright (c) 2012,2014,2021 CSIRO
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

// Package level header file
#include <askap/utils/DeconvolveTimerUtils.h>
#include <askap/askapparallel/AskapParallel.h>

/// ASKAP includes
#include <askap/askap/AskapUtil.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/Application.h>


/// 3rd party
#include <Common/ParameterSet.h>

#include <map>
#include <chrono>
#include <thread>

ASKAP_LOGGER(logger, ".tdeconvolvetimerutils");

using namespace std;

#define NUM_THREADS 3

namespace askap {
class DeconvolveTimerApp : public askap::Application
{
    public:

        int run(int argc, char* argv[]) final {
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));
            try {
                askap::utils::Timer timer;
                timer.start();
                const int no_timers = 5;
                askap::utils::SectionTimer sectionTimer(no_timers);
                const int N = 5;
                for (int n = 0; n < N; n++) {
                    #pragma omp parallel num_threads(NUM_THREADS)
                    {
                        // why not use omp master and omp barrier instead of
                        // omp single ? because we want to time the master thread
                        // each time through the for loop

                        #pragma omp master
                        {
                            sectionTimer.start(0);
                            std::chrono::seconds bedTime(1);
                            std::this_thread::sleep_for(bedTime);
                            sectionTimer.stop(0);
                        }
                        #pragma omp barrier

                        double sum = 0;

                        sectionTimer.start(0);
                        #pragma omp for
                        for (int i = 0; i < 100000000; i++) {
                            // dont worry about race condition
                            //sum += i;
                            sum += sin(1.0/(i+1));
                        }
                        sectionTimer.stop(0);

                        sectionTimer.start(1);

                        #pragma omp for
                        for (int i = 0; i < 100000000; i++) {
                            // dont worry about race condition
                            sum += sin(1.0/(i+1));
                        }
                        sectionTimer.stop(1);
                        sectionTimer.start(2);

                        #pragma omp for
                        for (int i = 0; i < 100000000; i++) {
                            // dont worry about race condition
                            sum += sin(1.0/(i+1));
                        }

                        sectionTimer.stop(2);
                        #pragma omp master
                        std::cout<<"n = "<<n<<", sum = "<<sum<<std::endl;

                    }
                }
                timer.stop();
                sectionTimer.summary();
                std::cout<<"outer loop: " <<  timer.summary()<<std::endl;
                return 0;
            } catch (...) {
                return 1;
            }
        }
};

} // end namespace askap

int main(int argc, char *argv[])
{
//    MPI_Init(NULL, NULL);
    askap::DeconvolveTimerApp app;
    return app.main(argc, argv);
//    MPI_Finalize();
}

