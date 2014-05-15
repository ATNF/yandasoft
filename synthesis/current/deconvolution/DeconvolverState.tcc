/// @file DeconvolverState.tcc
/// @brief Base class for State of Deconvolver
/// @details All the Stateing is delegated to this class so that
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

#include <askap/AskapLogging.h>
#include <casa/aips.h>

#include <deconvolution/DeconvolverState.h>
#include <deconvolution/DeconvolverState.h>

using namespace casa;

namespace askap {

    namespace synthesis {

        /// The current state
        template<class T>
        DeconvolverState<T>::DeconvolverState() : itsCurrentIter(0), itsStartIter(0),
                itsEndIter(0), itsPeakResidual(T(0)),
                itsTotalFlux(T(0)), itsObjectiveFunction(T(0)),
                itsInitialObjectiveFunction(T(0))
        {};

        template<class T>
        void DeconvolverState<T>::reset()
        {
            itsCurrentIter = 0;
            itsStartIter = 0;
            itsEndIter = 0;
            itsPeakResidual = T(0);
            itsTotalFlux = T(0);
            itsObjectiveFunction = T(0);
            itsInitialObjectiveFunction = T(0);
        }

    } // namespace synthesis

} // namespace askap
