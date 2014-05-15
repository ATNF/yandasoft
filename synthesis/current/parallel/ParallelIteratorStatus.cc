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

#include <parallel/ParallelIteratorStatus.h>
#include <askap/AskapError.h>

namespace askap {

namespace synthesis {

/// @brief default constructor - end mark meaning no more data available
ParallelIteratorStatus::ParallelIteratorStatus() : itsHasMore(false), itsNRow(0), itsNChan(0), itsNPol(0) {}

/// @brief Output shift operator for this class
/// @param os Output ostream
/// @param status a reference to status class
LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream& os, 
                                        const ParallelIteratorStatus &status)
{
  os.putStart("ParallelIteratorStatus",1);
  os<<status.itsHasMore<<status.itsNRow<<status.itsNChan<<status.itsNPol;
  os.putEnd();
  return os;
}                                        
                                        
/// @brief input shift operator for this class
/// @param[in] is Input stream
/// @param[in] par a reference to status class to be filled
LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream& is, 
                               ParallelIteratorStatus& status)
{
  const int version = is.getStart("ParallelIteratorStatus");
  ASKAPCHECK(version == 1, "Attempting to read a wrong version of ParallelIteratorStatus");
  is >> status.itsHasMore>>status.itsNRow>>status.itsNChan>>status.itsNPol;
  is.getEnd();
  return is;
}                               


} // namespace synthesis

} // namespace askap

