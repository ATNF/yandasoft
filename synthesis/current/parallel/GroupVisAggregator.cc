/// @file
/// 
/// @brief Object function to sum visibilities across the group of workers
/// @details If we distribute the model across multiple ranks we need to
/// sum up the results of degridding before calculation of the residual.
/// This object function can be used together with ImageFFTEquation to achieve this. 
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

#include <parallel/GroupVisAggregator.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel.groupvisaggregator");

#include <askap/AskapError.h>


namespace askap {

namespace synthesis {

/// @brief constructor, sets up communication class
/// @param[in] comms communication object
GroupVisAggregator::GroupVisAggregator(askap::askapparallel::AskapParallel& comms) : itsComms(comms), itsCommIndex(0)
{
  // we implicitly assume that casa::Complex just has two float data members and nothing else
  ASKAPDEBUGASSERT(sizeof(casa::Complex) == 2*sizeof(float));
  
  ASKAPLOG_DEBUG_STR(logger, "Set up visibility aggregator to sum degridded visibilities within each group of workers");
  const size_t group = itsComms.group();
  itsCommIndex = itsComms.interGroupCommIndex();
  ASKAPLOG_DEBUG_STR(logger, "  Worker group number "<<group<<" out of "<<itsComms.nGroups()<<
  " groups, intergroup communicator index: "<<itsCommIndex);

}  
  
/// @brief update visibility cube
/// @param[in,out] cube reference to visiblity cube to update 
void GroupVisAggregator::update(casa::Cube<casa::Complex> &cube) const
{
  ASKAPASSERT(cube.nelements() != 0); 
  ASKAPASSERT(cube.contiguousStorage());
  // not a very safe way of doing it, but this way we could benefit from MPI reduce/broadcast 
  ASKAPLOG_DEBUG_STR(logger, "about to sum over the data using comm index: "<<itsCommIndex<<" shape: "<<cube.shape()<<" (0,0,0): "<<cube(0,0,0));
  itsComms.sumAndBroadcast((float *)cube.data(), 2 * cube.nelements(), itsCommIndex);
  ASKAPLOG_DEBUG_STR(logger, "after mpi call, shape: "<<cube.shape()<<" (0,0,0): "<<cube(0,0,0));
}

/// @brief aggregate flag with the logical or operation
/// @param[in,out] flag flag to reduce
void GroupVisAggregator::aggregateFlag(bool &flag) const
{
  ASKAPLOG_DEBUG_STR(logger, "about to aggregate flag ("<<flag<<")");
  itsComms.aggregateFlag(flag, itsCommIndex);
  ASKAPLOG_DEBUG_STR(logger, "flag after aggregation is "<<flag);
}

/// @brief helper method to create an instance of this class
/// @details It checks whether the current setup has multiple groups of workers
/// and if yes, creates an instance of this class. Otherwise, an empty shared pointer
/// is returned (and therefore inter-rank communication is not done)
/// @param[in] comms communication object
/// @return shared pointer to an instance of this class  
boost::shared_ptr<GroupVisAggregator> GroupVisAggregator::create(askap::askapparallel::AskapParallel& comms)
{
  if (comms.nGroups() > 1) {
      boost::shared_ptr<GroupVisAggregator> result(new GroupVisAggregator(comms));
      return result;
  }
  ASKAPLOG_DEBUG_STR(logger, "There are no groupping of workers, inter-rank summation of degridded visibilities is not necessary");
  return boost::shared_ptr<GroupVisAggregator>();
}


} // namespace synthesis

} // namespace askap

