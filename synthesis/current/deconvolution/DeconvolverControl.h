/// @file DeconvolverControl.h
/// @brief Base class for Control of Deconvolver
/// @details All the Controlling is delegated to this class. This encourages
/// a standard approach to convergence testing.
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

#ifndef ASKAP_SYNTHESIS_DECONVOLVERCONTROL_H
#define ASKAP_SYNTHESIS_DECONVOLVERCONTROL_H

#include <string>

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <boost/shared_ptr.hpp>
#include <Common/ParameterSet.h>
#include <askap/ISignalHandler.h>
#include <askap/SignalCounter.h>

#include <deconvolution/DeconvolverState.h>

namespace askap {

    namespace synthesis {

        /// @brief Base class for Control of Deconvolver
        /// @details All the controlling is delegated to this class so that
        /// standard parameters and stopping criteria are allowed.
        /// @ingroup Deconvolver
        template<typename T> class DeconvolverControl {

            public:
                typedef boost::shared_ptr<DeconvolverControl<T> > ShPtr;

                /// @brief Enumerate the possible termination causes
                enum TerminationCause {
                    CONVERGED,
                    DIVERGED,
                    EXCEEDEDITERATIONS,
                    SIGNALED,
                    NOTTERMINATED,
                    UNKNOWN
                };

                DeconvolverControl();

                virtual ~DeconvolverControl();

                /// @brief configure basic parameters of the solver
                /// @details This method encapsulates extraction of basic solver parameters
                /// from the parset.
                /// @param[in] parset parset
                virtual void configure(const LOFAR::ParameterSet &parset);

                /// @brief Terminate now?
                /// @detail The state of the deconvolver is passed via the
                /// DeconvolverState instance. Information in that is
                /// used to evaluate termination.
                casa::Bool terminate(const DeconvolverState<T>& ds);

                /// @brief Return the termination as a string
                /// @param[out] Termination cause returned as a string
                casa::String terminationString() const;

                /// @brief Set the termination cause as an enumeration
                /// @param[in] Termination cause enumeration
                void setTerminationCause(const TerminationCause cause) {
                    itsTerminationCause = cause;
                };

                /// @brief Return the termination as an enumeration
                /// @param[out] Termination cause returned as an enumeration
                TerminationCause terminationCause() const {
                    return itsTerminationCause;
                };

                /// @brief Return the desired algorithm name as a string
                /// @param[out] Name of desired algorithm e.g. MultiScale
                casa::String algorithm() const {return itsAlgorithm;};

                /// @brief Set the desired algorithm name as a string
                /// @param[in] algorithm Name of desired algorithm e.g. MultiScale
                void setAlgorithm(const casa::String algorithm) {
                    itsAlgorithm = algorithm;
                };

                /// @brief Set the desired gain
                /// @param[in] gain Desired gain
                void setGain(casa::Float gain) {
                    itsGain = gain;
                }

                casa::Float gain() const {
                    return itsGain;
                }

                /// @brief Set the desired tolerance
                /// @param[in] tolerance Desired tolerance
                void setTolerance(casa::Float tolerance) {
                    itsTolerance = tolerance;
                }

                casa::Float tolerance() const {
                    return itsTolerance;
                }

                /// @brief Set the desired target number of iterations
                /// @param[in] targetiter Desired number of iterations
                void setTargetIter(casa::Int targetiter) {
                    itsTargetIter = targetiter;
                }

                casa::Int targetIter() const {
                    return itsTargetIter;
                }

                /// @brief Set the desired Lagrange multiplier
                /// @detail This is a niche parameter! Currently it
                /// is used only in the FISTA deconvolver.
                /// @param[in] lambda The desired Lagrange multiplier
                void setLambda(T lambda) {
                    itsLambda = lambda;
                }

                T lambda() const {
                    return itsLambda;
                }

                /// @brief Set the desired target objective function
                /// @detail Algorithms can be monitored for an objective
                /// function that decreases e.g. maximum absolute residual
                /// or L1 norm. If the target objective function is set then
                /// control will stop the algorithm when the actual objective
                /// function is below the target.
                /// @param[in] objectivefunction Desired value of objective function
                void setTargetObjectiveFunction(T objectiveFunction) {
                    itsTargetObjectiveFunction = objectiveFunction;
                }

                T targetObjectiveFunction() const {
                    return itsTargetObjectiveFunction;
                }

                /// @brief Set the desired target flux
                /// @detail Some algorithms can work with a target flux (increasing
                /// with iteration). If the target flux is set then
                /// control will stop the algorithm when the actual flux is above the target
                /// flux (and the other criteria have been met.
                /// @param[in] targetFlux Desired value of flux
                void setTargetFlux(T targetFlux) {
                    itsTargetFlux = targetFlux;
                }

                T targetFlux() const {
                    return itsTargetFlux;
                }

                /// @brief Set the desired fractional threshold
                /// @detail If the fractional threshold is set then
                /// the algorithm will stop when the maximum absolute residual
                /// is less than this number times the original maximum absolute
                /// residual.
                /// @param[in] fractionalThreshold Fractional threshold (0.1 is 10%).
                void setFractionalThreshold(casa::Float fractionalThreshold) {
                    itsFractionalThreshold = fractionalThreshold;
                }

                casa::Float fractionalThreshold() {
                    return itsFractionalThreshold;
                }

                /// @brief Set the desired PSF width in pixels
                /// @detail Some algorithms can work with a fractional of the
                /// full PSF. This allows specification of the support. Note that
                /// this is the full width. The peak of the PSF will be located
                /// at psfWidth/2, psfWidth/2.
                /// @param[in] Width of the PSF in pixels (e.g. 1024).
                void setPSFWidth(const casa::Int psfWidth) {itsPSFWidth = psfWidth;}

                /// @brief Get the desired PSF width in pixels
                casa::Int psfWidth() const {return itsPSFWidth;};

            private:
                casa::String itsAlgorithm;
                TerminationCause itsTerminationCause;
                casa::Int itsTargetIter;
                T itsTargetObjectiveFunction;
                T itsTargetFlux;
                casa::Float itsFractionalThreshold;
                casa::Float itsGain;
                casa::Float itsTolerance;
                casa::Int itsPSFWidth;
                T itsLambda;
                askap::SignalCounter itsSignalCounter;
                askap::ISignalHandler* itsOldHandler;
        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/DeconvolverControl.tcc>

#endif


