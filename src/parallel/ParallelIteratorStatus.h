/// @file 
/// @brief iteration status used in parallel iterator
/// @details This class is used in communications between master and
/// workers for the parallel iterator
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

#ifndef ASKAP_SYNTHESIS_PARALLEL_ITERATOR_STATUS_H
#define ASKAP_SYNTHESIS_PARALLEL_ITERATOR_STATUS_H

#include <casacore/casa/aips.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>

namespace askap {

namespace synthesis {

/// @brief iteration status used in parallel iterator
/// @details This class is used in communications between master and
/// workers for the parallel iterator
/// @ingroup parallel
struct ParallelIteratorStatus {

  /// @brief default constructor - end mark meaning no more data available
  ParallelIteratorStatus();

  /// @brief Output shift operator for this class
  /// @param os Output ostream
  /// @param status a reference to status class
  friend LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream& os, 
                                        const ParallelIteratorStatus &status);
  /// @brief input shift operator for this class
  /// @param[in] is Input stream
  /// @param[in] par a reference to status class to be filled
  friend LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream& is, 
                                        ParallelIteratorStatus& status); 

  // data fields

  /// @brief true if there are more data available
  bool itsHasMore;
  
  /// @brief number of rows
  casa::uInt itsNRow;
  
  /// @brief number of polarisations
  casa::uInt itsNPol;

  /// @brief operations mode
  /// @details Basically, an additional integer the meaning of which depends on the user
  int itsMode;

  /// @brief total number of channels
  /// @details Workers are in general agnostic of the full spectrum, but sometimes it is handy
  /// to pass these details through iterator code to the high-level code to avoid re-implementing
  /// the same communication pattern
  casa::uInt itsTotalNChan;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_PARALLEL_ITERATOR_STATUS_H

