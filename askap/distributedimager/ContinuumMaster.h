/// @file ContinuumMaster.h
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
/// Derived from SpectralLineMaster by @author Ben Humphreys <ben.humphreys@csiro.au>
/// @author Stephen Ord <stephen.ord@csiro.au>


#ifndef ASKAP_CP_SIMAGER_CONTINUUMMASTER_H
#define ASKAP_CP_SIMAGER_CONTINUUMMASTER_H

// System includes
#include <string>
#include <vector>

// ASKAPsoft includes
#include <boost/scoped_ptr.hpp>
#include <Common/ParameterSet.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/askap/StatReporter.h>

// Local includes

#include "askap/distributedimager/CalcCore.h"
#include "askap/distributedimager/CubeBuilder.h"
#include "askap/distributedimager/CubeComms.h"

namespace askap {
namespace cp {

class ContinuumMaster {
    public:
        ContinuumMaster(LOFAR::ParameterSet& parset,
                           CubeComms& comms, StatReporter& stats);
        ~ContinuumMaster();

        void run(void);
    private:

        struct MSInfo {
            casacore::uInt nChan;
            std::vector<casacore::Quantity> freqs;
        };

        /// @brief Utility function to get dataset names from parset.
        ///
        /// Given a ParameterSet, return a vector containing all the datasets
        /// specified. This function will look for datasets in the Cimager manner:
        /// @code
        /// Cimager.dataset        = [10uJy_stdtest_0.ms,10uJy_stdtest_1.ms]
        /// @endcode
        ///
        /// It also supports another method which is necessary for the specification
        /// of large numbers of datasets:
        /// @code
        /// Cimager.dataset0                             = 10uJy_stdtest_0.ms
        /// Cimager.dataset1                             = 10uJy_stdtest_1.ms
        /// <and so on>
        /// @endcode
        ///
        /// @return a vector containing in each element one dataset.
        std::vector<std::string> getDatasets();

        std::vector<int> getBeams(void);


        /// Parameter set
        LOFAR::ParameterSet& itsParset;

        /// Communications class
        CubeComms& itsComms;

        /// statistics
        StatReporter& itsStats;

        std::map<unsigned int, casacore::Vector<casacore::Quantum<double> > > itsBeamList;



};

};
};

#endif
