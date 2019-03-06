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

#ifndef GROUP_VIS_AGGREGATOR_H
#define GROUP_VIS_AGGREGATOR_H

#include <measurementequation/IVisCubeUpdate.h>
#include <askapparallel/AskapParallel.h>

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief Object function to sum visibilities across the group of workers
/// @details If we distribute the model across multiple ranks we need to
/// sum up the results of degridding before calculation of the residual.
/// This object function can be used together with ImageFFTEquation to achieve this. 
/// @ingroup parallel
class GroupVisAggregator : public IVisCubeUpdate {
public:

  /// @brief constructor, sets up communication class
  /// @param[in] comms communication object
  explicit GroupVisAggregator(askap::askapparallel::AskapParallel& comms);
  
  /// @brief update visibility cube
  /// @param[in,out] cube reference to visiblity cube to update 
  virtual void update(casa::Cube<casa::Complex> &cube) const;
  
  /// @brief aggregate flag with the logical or operation
  /// @param[in,out] flag flag to reduce
  virtual void aggregateFlag(bool &flag) const;
    
  /// @brief helper method to create an instance of this class
  /// @details It checks whether the current setup has multiple groups of workers
  /// and if yes, creates an instance of this class. Otherwise, an empty shared pointer
  /// is returned (and therefore inter-rank communication is not done)
  /// @param[in] comms communication object
  /// @return shared pointer to an instance of this class  
  static boost::shared_ptr<GroupVisAggregator> create(askap::askapparallel::AskapParallel& comms);
  
private:
  
  /// @brief class for communications
  askap::askapparallel::AskapParallel& itsComms;  
  
  /// @brief communicator index
  size_t itsCommIndex;
};

} // namespace synthesis

} // namespace askap


#endif // #ifndef GROUP_VIS_AGGREGATOR_H


