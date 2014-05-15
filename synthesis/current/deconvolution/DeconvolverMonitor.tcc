/// @file DeconvolverMonitor.tcc
/// @brief Base class for monitor of Deconvolver
/// @details All the monitoring is delegated to this class so that
/// more flexibility is possible.
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
#include <askap/AskapLogging.h>
ASKAP_LOGGER(decmonlogger, ".deconvolution.monitor");

#include <deconvolution/DeconvolverMonitor.h>

namespace askap {

    namespace synthesis {

        template<class T>
        DeconvolverMonitor<T>::DeconvolverMonitor() : itsVerbose(false),
                itsLogEvery(1)
        {
        }

        /// Monitor the current state
        template<class T>
        void DeconvolverMonitor<T>::monitor(const DeconvolverState<T>& ds)
        {
            if (itsVerbose) {
                ASKAPLOG_INFO_STR(decmonlogger, "Iteration " << ds.currentIter()
                                      << ", Peak residual " << ds.peakResidual()
                                      << ", Objective function " << ds.objectiveFunction()
                                      << ", Total flux " << ds.totalFlux());
            } else {
                if ((ds.currentIter() % itsLogEvery) == 0) {
                    ASKAPLOG_INFO_STR(decmonlogger, "Iteration " << ds.currentIter()
                                          << ", Peak residual " << ds.peakResidual()
                                          << ", Objective function " << ds.objectiveFunction()
                                          << ", Total flux " << ds.totalFlux());
                }
            }
        }

        template<class T>
        void DeconvolverMonitor<T>::configure(const LOFAR::ParameterSet& parset)
        {
            itsLogEvery = parset.getInt("logevery", 1);
            itsVerbose = parset.getBool("verbose", false);
        }

    } // namespace synthesis

} // namespace askap
