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
#include <fitting/Params.h>


// Local includes

#include "distributedimager/CalcCore.h"
#include "distributedimager/CubeBuilder.h"
#include "distributedimager/CubeComms.h"
#include "distributedimager/MSGroupInfo.h"

namespace askap {
namespace cp {

class ContinuumMaster {
    public:
        ContinuumMaster(LOFAR::ParameterSet& parset,
                           CubeComms& comms);
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
        /// @param[in] parset   the parameterset to use as input.
        /// @return a vector containing in each element one dataset.
        std::vector<std::string> getDatasets(const LOFAR::ParameterSet& itsParset);

        std::vector<int> getBeams(void);


        /// Parameter set
        LOFAR::ParameterSet& itsParset;

        /// Communications class
        CubeComms& itsComms;

        /// The rest of these are for the Spectral line case - which I haven't added yet Ord-6/2016
        MSGroupInfo isMSGroupInfo;

        boost::scoped_ptr<CubeBuilder> itsImageCube;
        boost::scoped_ptr<CubeBuilder> itsPSFCube;
        boost::scoped_ptr<CubeBuilder> itsResidualCube;
        boost::scoped_ptr<CubeBuilder> itsWeightsCube;
        boost::scoped_ptr<CubeBuilder> itsPSFimageCube;
        boost::scoped_ptr<CubeBuilder> itsRestoredCube;

        std::map<unsigned int, casacore::Vector<casacore::Quantum<double> > > itsBeamList;

      

};

};
};

#endif
