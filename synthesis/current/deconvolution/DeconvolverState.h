/// @file DeconvolverState.h
/// @brief Defines the state for a deconvolve
/// @details This interface class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
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

#ifndef ASKAP_SYNTHESIS_DECONVOLVERSTATE_H
#define ASKAP_SYNTHESIS_DECONVOLVERSTATE_H

#include <string>

#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casa/Arrays/Array.h>

namespace askap {

    namespace synthesis {

        /// @brief Class to hold the state of deconvolution
        /// @details This class attempts to encapsulate the state of
        /// a deconvolver in a standard way so that termination
        /// and monitoring can be generic.
        /// @ingroup Deconvolver
        template<class T> class DeconvolverState {

            public:

                typedef boost::shared_ptr<DeconvolverState<T> > ShPtr;

                DeconvolverState();

                virtual ~DeconvolverState() {};

                void incIter() {
                    itsCurrentIter++;
                }

                casa::Int currentIter() const {
                    return itsCurrentIter;
                }

                T peakResidual() const {
                    return itsPeakResidual;
                }

                T totalFlux() const {
                    return itsTotalFlux;
                }

                casa::Int startIter() const {
                    return itsStartIter;
                }

                void setObjectiveFunction(T objectiveFunction) {
                    itsObjectiveFunction = objectiveFunction;
                    if (itsInitialObjectiveFunction <= T(0.0)) {
                        itsInitialObjectiveFunction = objectiveFunction;
                    }
                }

                T objectiveFunction() const {
                    return itsObjectiveFunction;
                }

                void setInitialObjectiveFunction(T objectiveFunction) {
                    itsInitialObjectiveFunction = objectiveFunction;
                }

                void resetInitialObjectiveFunction() {
                    itsInitialObjectiveFunction = T(0.0);
                }

                T initialObjectiveFunction() const {
                    return itsInitialObjectiveFunction;
                }

                void setCurrentIter(casa::Int currentIter) {
                    itsCurrentIter = currentIter;
                }

                void setPeakResidual(T peakResidual) {
                    itsPeakResidual = peakResidual;
                }

                void setTotalFlux(T totalFlux) {
                    itsTotalFlux = totalFlux;
                }

                void setStartIter(casa::Int startIter) {
                    itsStartIter = startIter;
                }

                /// Reset the state
                void reset();

            private:

                casa::Int itsCurrentIter;
                casa::Int itsStartIter;
                casa::Int itsEndIter;
                T itsPeakResidual;
                T itsTotalFlux;
                T itsObjectiveFunction;
                T itsInitialObjectiveFunction;
        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/DeconvolverState.tcc>

#endif
