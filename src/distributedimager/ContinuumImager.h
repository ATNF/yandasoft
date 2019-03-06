/// @file ContinuumImager.h
///
/// @copyright (c) 2016 CSIRO
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
/// @author Stephen Ord <stephen.ord@csiro.au>

#ifndef ASKAP_CP_SIMAGER_CONTINUUMIMAGER_H
#define ASKAP_CP_SIMAGER_CONTINUUMIMAGER_H

// ASKAPsoft includes
#include <Common/ParameterSet.h>


// Local package includes
#include <distributedimager/CubeComms.h>

namespace askap {
namespace cp {

/// @brief Main class for the Continuum version of the distributed imager.
///        It is essentially a version of the SpectralLineImager that has a minor
///        cycle performed centrally
class ContinuumImager
{
    public:
        /// @brief Construct a Distributed Imager.
        ///
        /// @param[in]  parset  the parameter set containing
        ///                     the configuration.
        /// @param[in]  comms   an instance of IBasicComms.
        ContinuumImager(LOFAR::ParameterSet& parset,
                           CubeComms& comms);

        /// @brief Destructor.
        ~ContinuumImager();

        /// @brief Run method
        void run(void);

    private:

        // Returns true if the caller is the master process,
        // else false.
        bool isMaster(void);

        // Id of the master process
        static const int itsMaster = 0;

        // Parameter set
        LOFAR::ParameterSet& itsParset;

        // Communications class
        CubeComms& itsComms;
}; // end class

}; // end cp
}; // end askap

#endif
