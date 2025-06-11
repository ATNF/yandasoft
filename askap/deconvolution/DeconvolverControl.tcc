/// @file DeconvolverControl.tcc
/// @brief Base class for Control of Deconvolver
/// @details All the Controlling is delegated to this class so that
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

#include <askap/askap_synthesis.h>

#include <casacore/casa/aips.h>
#include <askap/askap/SignalManagerSingleton.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(decctllogger, ".deconvolution.control");

#include <askap/deconvolution/DeconvolverState.h>
#include <askap/deconvolution/DeconvolverControl.h>

using namespace casacore;

namespace askap {

    namespace synthesis {

        template<class T>
        DeconvolverControl<T>::DeconvolverControl() :
                itsAlgorithm(""), itsTerminationCause(NOTTERMINATED), itsTargetIter(1),
                itsTargetObjectiveFunction(T(0)),itsTargetObjectiveFunction2(T(0)),
                itsTargetFlux(T(0.0)),itsGain(1.0), itsTolerance(1e-4),
                itsFractionalThreshold(T(0.0)),itsAbsoluteThreshold(0.0),
                itsPSFWidth(0), itsDetectDivergence(False), itsDetectMildDivergence(False),
                itsDivergenceLevels({1.1f,2.0f,2.0f,0.05f}), itsDeepCleanMode(False), itsLambda(T(100.0))
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
            if (abs(state.objectiveFunction()) < targetObjectiveFunction()) {
                // Now check if we want to enter deep cleaning mode
                if (targetObjectiveFunction2()>0) {
                    if (!deepCleanMode()) {
                        ASKAPLOG_INFO_STR(decctllogger, "Starting deep cleaning phase");
                        setDeepCleanMode();
                    }
                    if (abs(state.objectiveFunction()) < targetObjectiveFunction2()) {
                        ASKAPLOG_INFO_STR(decctllogger, "Objective function " << state.objectiveFunction()
                                            << " less than 2nd target " << targetObjectiveFunction2());
                        setTerminationCause(CONVERGED);
                        return True;
                    }
                } else {
                    ASKAPLOG_INFO_STR(decctllogger, "Objective function " << state.objectiveFunction()
                                          << " less than target " << targetObjectiveFunction());
                    setTerminationCause(CONVERGED);
                    return True;
                }
            }
            //
            if (abs(state.objectiveFunction()) < fractionalThreshold()*state.initialObjectiveFunction()) {
                ASKAPLOG_INFO_STR(decctllogger, "Objective function " << state.objectiveFunction()
                                      << " less than fractional threshold " << fractionalThreshold()
                                      << " * initialObjectiveFunction : " << state.initialObjectiveFunction());
                setTerminationCause(CONVERGED);
                return True;
            }
            if (abs(state.objectiveFunction()) < absoluteThreshold()) {
                ASKAPLOG_INFO_STR(decctllogger, "Objective function " << state.objectiveFunction()
                                      << " less than absolute threshold " << absoluteThreshold());
                setTerminationCause(CONVERGED);
                return True;
            }

            // Terminate if the target number of iterations is not set
            ASKAPCHECK(targetIter() > 0, "Target number of iterations not set");

            // Check for too many iterations
            if ((state.currentIter() > -1) && (targetIter() > 0) && (state.currentIter() >= targetIter())) {
                setTerminationCause(EXCEEDEDITERATIONS);
                return True;
            }

            // Check if we have started to diverge
            if (detectDivergence()) {
                // Simplest check: next component > 2* initial residual
                if ( state.initialObjectiveFunction() > 0 &&
                     state.objectiveFunction() > itsDivergenceLevels[1] * state.initialObjectiveFunction() )
                {
                  ASKAPLOG_INFO_STR(decctllogger, "Clean diverging - Objective function " <<
                  state.objectiveFunction() << " > "<<itsDivergenceLevels[1]<<" * initialObjectiveFunction = " <<
                    itsDivergenceLevels[1] * state.initialObjectiveFunction());
                  setTerminationCause(DIVERGED);
                  return True;
                }

                // Check for major cycle divergence
                // (only need to check this once per major cycle but it fits here)
                if ( state.previousInitialObjectiveFunction() > 0 &&
                     state.initialObjectiveFunction() > itsDivergenceLevels[2] * state.previousInitialObjectiveFunction() )
                {
                  ASKAPLOG_INFO_STR(decctllogger, "Clean diverging - Initial Objective function " <<
                  state.initialObjectiveFunction() << " > "<<itsDivergenceLevels[2]<<" * previous Initial ObjectiveFunction = " <<
                    itsDivergenceLevels[2] * state.previousInitialObjectiveFunction());
                  setTerminationCause(DIVERGED);
                  return True;
                }
            }

            // Check for mild divergence - just go to next major cycle
            if (detectMildDivergence()) {
                // next component > 1.1x smallest component this cycle & have done >5% of iterations
                if ( state.objectiveFunction() > itsDivergenceLevels[0] * state.smallestObjectiveFunction() &&
                     state.currentIter() > itsDivergenceLevels[3] * targetIter())
                {
                    ASKAPLOG_INFO_STR(decctllogger, "Clean starting to diverge at "<<state.objectiveFunction()<<" > "
                    << itsDivergenceLevels[0] << "*" << state.smallestObjectiveFunction() << " - stopping minor cycles");
                    setTerminationCause(DIVERGING);
                    return True;
                }
            }

            // Check for external signal
            if (itsSignalCounter.getCount() > 0) {
                setTerminationCause(SIGNALED);
                itsSignalCounter.resetCount(); // This signal has been actioned, so reset
                return True;
            }
            return False;
        }

        template<class T>
        T DeconvolverControl<T>::level(const DeconvolverState<T>& ds, T safetyMargin) const
        {
            T level = targetObjectiveFunction();
            if (deepCleanMode()) {
                level = targetObjectiveFunction2();
            }
            T level2 = fractionalThreshold() * ds.initialObjectiveFunction();
            if (level2 > level) {
                level = level2;
            }
            ASKAPCHECK(safetyMargin > 0 && safetyMargin < 1, "safetyMargin should be between 0 and 1");
            return level * (1-safetyMargin);
        }

        template<class T>
        String DeconvolverControl<T>::terminationString() const
        {
            switch (terminationCause()) {
                case CONVERGED:
                    return String("Converged");
                    break;
                case DIVERGING:
                    return String("Starting to diverge");
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
            setGain(parset.getFloat("gain", 0.1));
            setTolerance(parset.getFloat("tolerance", 1e-3));
            setTargetIter(parset.getInt32("niter", 100));
            setTargetFlux(parset.getFloat("targetflux", 0));
            setTargetObjectiveFunction(parset.getFloat("targetobjective", 0.0));
            setTargetObjectiveFunction2(parset.getFloat("targetobjective2", 0.0));
            setFractionalThreshold(parset.getFloat("fractionalthreshold", 0.0));
            setAbsoluteThreshold(parset.getFloat("absolutethreshold", 0.0));
            setLambda(parset.getFloat("lambda", 0.0001));
            setPSFWidth(parset.getInt32("psfwidth", 0));
            setDetectDivergence(parset.getBool("detectdivergence",true));
            setDetectMildDivergence(parset.getBool("detectmilddivergence",false));
            const std::vector<casa::Float> levels = parset.getFloatVector("divergencelevels",{1.1f,2.0f,2.0f,0.05f});
            ASKAPCHECK(levels.size()==4 && levels[0]>1 && levels[1]>1 && levels[2]>1,
                "divergencelevels parameter should have 4 entries, each first three > 1.0");
            setDivergenceLevels(levels);
        }

    } // namespace synthesis

} // namespace askap
