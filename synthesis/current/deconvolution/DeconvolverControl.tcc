/// @file DeconvolverControl.tcc
/// @brief Base class for Control of Deconvolver
/// @details All the Controling is delegated to this class so that
/// more control is possible.
/// @ingroup Deconvolver
///
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
///

#include <askap_synthesis.h>

#include <casa/aips.h>
#include <askap/SignalManagerSingleton.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(decctllogger, ".deconvolution.control");

#include <deconvolution/DeconvolverState.h>
#include <deconvolution/DeconvolverControl.h>

using namespace casa;

namespace askap {

    namespace synthesis {

        template<class T>
        DeconvolverControl<T>::DeconvolverControl() :
                itsAlgorithm(""), itsTerminationCause(NOTTERMINATED), itsTargetIter(1),
                itsTargetObjectiveFunction(T(0)), itsTargetFlux(T(0.0)),
                itsGain(1.0), itsTolerance(1e-4),
                itsPSFWidth(0), itsLambda(T(100.0))
        {
            // Install a signal handler to count signals so receipt of a signal
            // can be used to terminate the minor-cycle loop
            itsOldHandler = SignalManagerSingleton::instance()->registerHandler(SIGUSR2, &itsSignalCounter);
        };

        template<class T>
        DeconvolverControl<T>::~DeconvolverControl()
        {
            itsOldHandler = SignalManagerSingleton::instance()->registerHandler(SIGUSR2, itsOldHandler);
        }

        /// Control the current state
        template<class T>
        Bool DeconvolverControl<T>::terminate(const DeconvolverState<T>& state)
        {
            // Check for convergence
            if (abs(state.objectiveFunction()) < this->itsTargetObjectiveFunction) {
                ASKAPLOG_INFO_STR(decctllogger, "Objective function " << state.objectiveFunction()
                                      << " less than target " << itsTargetObjectiveFunction);
                itsTerminationCause = CONVERGED;
                return True;
            }
            //
            if (abs(state.objectiveFunction()) < this->itsFractionalThreshold*state.initialObjectiveFunction()) {
                ASKAPLOG_INFO_STR(decctllogger, "Objective function " << state.objectiveFunction()
                                      << " less than fractional threshold " << itsFractionalThreshold
                                      << " * initialObjectiveFunction : " << state.initialObjectiveFunction());
                itsTerminationCause = CONVERGED;
                return True;
            }

            // Terminate if the target number of iterations is not set
            ASKAPCHECK(this->targetIter() > 0, "Target number of iterations not set");

            // Check for too many iterations
            if ((state.currentIter() > -1) && (this->targetIter() > 0) && (state.currentIter() >= this->targetIter())) {
                itsTerminationCause = EXCEEDEDITERATIONS;
                return True;
            }

            // Check for external signal
            if (itsSignalCounter.getCount() > 0) {
                itsTerminationCause = SIGNALED;
                itsSignalCounter.resetCount(); // This signal has been actioned, so reset
                return True;
            }
            return False;
        }

        template<class T>
        String DeconvolverControl<T>::terminationString() const
        {
            switch (itsTerminationCause) {
                case CONVERGED:
                    return String("Converged");
                    break;
                case DIVERGED:
                    return String("Diverged");
                    break;
                case EXCEEDEDITERATIONS:
                    return String("Exceeded maximum number of iterations");
                    break;
                case SIGNALED:
                    return String("Signaled to terminate");
                    break;
                case NOTTERMINATED:
                    return String("Not yet terminated");
                    break;
                case UNKNOWN:
                    return String("Termination for unknown reason");
                    break;
                default:
                    return String("Logic error in termination");
                    break;
            }
        }

        template<class T>
        void DeconvolverControl<T>::configure(const LOFAR::ParameterSet& parset)
        {
            this->setGain(parset.getFloat("gain", 0.1));
            this->setTolerance(parset.getFloat("tolerance", 1e-3));
            this->setTargetIter(parset.getInt32("niter", 100));
            this->setTargetFlux(parset.getFloat("targetflux", 0));
            this->setTargetObjectiveFunction(parset.getFloat("targetobjective", 0.0));
            this->setFractionalThreshold(parset.getFloat("fractionalthreshold", 0.0));
            this->setLambda(parset.getFloat("lambda", 0.0001));
            this->setPSFWidth(parset.getInt32("psfwidth", 0));
        }

    } // namespace synthesis

} // namespace askap
