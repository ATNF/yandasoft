/// @file FlaggerFactory.h
///
/// @copyright (c) 2011 CSIRO
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>

#ifndef ASKAP_SYNTHESIS_FLAGGERFACTORY_H
#define ASKAP_SYNTHESIS_FLAGGERFACTORY_H

// ASKAPsoft includes
#include "Common/ParameterSet.h"

// Local package includes
#include "askap/flagging/IFlagger.h"

namespace askap {
namespace synthesis {

/// @brief A factory that, given a Parameter Set, will create a flagging
/// flagger instance for each flagger enabled in the parset.
class FlaggerFactory {
    public:

        /// Appends one std::vector of flaggers to another. Vector "v2" is appended
        /// to std::vector "v1"
        ///
        /// @param[in,out]  v1  the std::vector that will be appended to.
        /// @param[in]      v2  the std::vector that will be appended to v1.
        static void appendFlaggers(std::vector<std::shared_ptr<IFlagger> >& v1,
                                   const std::vector<std::shared_ptr<IFlagger> > v2);

        /// Builds flagging flagger objects based on the configuration in
        /// the parameter set.
        ///
        /// @param[in] parset   the parameter set which contains an ASCII description
        ///                     of the flagging strategies to use.
        /// @param[in] ms       name of the measurementset to flag
        /// @return a std::vector containing pointers to the flagging strategies.
        static std::vector<std::shared_ptr<IFlagger> > build(
            const LOFAR::ParameterSet& parset, const std::string &ms);
};

}
}

#endif
