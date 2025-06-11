/// @file IFlagger.h
///
/// @copyright (c) 2012-2014 CSIRO
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

#ifndef ASKAP_SYNTHESIS_IFLAGGER_H
#define ASKAP_SYNTHESIS_IFLAGGER_H

// ASKAPsoft includes
#include "askap/dataaccess/SharedIter.h"
#include "askap/dataaccess/IDataIterator.h"
// Local package includes
#include "askap/flagging/FlaggingStats.h"

namespace askap {
namespace synthesis {

// Tuple row key causes performance issues - using integers instead - this might overflow for large Nant*Nfeed
// Uncomment next line to use tuples
//#define TUPLE_INDEX
//                            fieldID    feed1             feed2         antenna1        antenna2      polarisation
#ifdef TUPLE_INDEX
using rowKey = std::tuple<casacore::Int, casacore::Int, /*casacore::Int,*/ casacore::Int, casacore::Int, casacore::Int>;
#else
using rowKey = casacore::uLong;
#endif

/// @brief An interface for classes that perform flagging on a per row basis.
class IFlagger {
    public:

        /// Destructor
        virtual ~IFlagger();

        /// Perform flagging (if necessary) for the row with index "row".
        ///
        /// @param[in,out] di  the iterator for the data
        /// @param[in] rowFlag  a vector specifying which rows are completely
        ///                     flagged already and can be skipped
        /// @param[in] pass     pass number - some flaggers collect statistics on the
        ///                     first pass (= 0) that is used in a subsequent pass
        /// @param[in] dryRun   if true the dataset will not be modified,
        ///                     however statistics will be calculated indicating
        ///                     what flagging would have been done.
        virtual void processRows(const accessors::IDataSharedIter& di,
                                 const casacore::Vector<bool>& rowFlag,
                                 const casacore::uInt pass, const bool dryRun) = 0;

        /// Returns flagging statistics
        virtual FlaggingStats stats(void) const = 0;

        /// Functions associated with multiple passees
        /// @param[in] pass     number of passes over the data already performed
        virtual casacore::Bool processingRequired(const casacore::uInt pass) const = 0;

};

}
}

#endif
