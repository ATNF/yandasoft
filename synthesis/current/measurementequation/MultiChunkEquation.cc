/// @file
/// 
/// @brief A measurement equation acting on an iterator.
/// @details This is a temporary class (I hope) to retain the existing 
/// interface for measurement equations, where these equations are applied
/// to all chunks (accessors) of the measurement set at once. It looks like 
/// in the future we need to redesign existing measurement equations to 
/// work with one iteration only (i.e. accessor instead of iterator). This 
/// class allows to simplify this transition, by factoring out the old
/// interface and implementing it via the new one.
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


// own includes
#include <measurementequation/MultiChunkEquation.h>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

/// @brief Standard constructor, which remembers data iterator.
/// @param idi data iterator
MultiChunkEquation::MultiChunkEquation(const IDataSharedIter& idi) :
                 itsSharedIterator(idi) {}

/// @brief Calculate the normal equations for the iterator
/// @details This version iterates through all chunks of data and
/// calls an abstract method declared in IMeasurementEquation for each 
//// individual accessor (each iteration of the iterator)
/// @param[in] ne Normal equations
void MultiChunkEquation::calcEquations(askap::scimath::INormalEquations& ne) const
{ 
  for (itsSharedIterator.init(); itsSharedIterator.hasMore(); 
                                 itsSharedIterator.next()) {
       calcEquations(*itsSharedIterator,ne);
  }
}

/// @brief Predict model visibility for the iterator.
/// @details This version of the predict method iterates
/// over all chunks of data and calls an abstract method declared
/// in IMeasurementEquation for each accessor. 
void MultiChunkEquation::predict() const
{
  for (itsSharedIterator.init(); itsSharedIterator.hasMore(); 
                                 itsSharedIterator.next()) {
       predict(*itsSharedIterator);
  }
}
   
/// @brief access to the iterator associated with this equation
/// @return a const reference to the iterator held by this object
const IDataSharedIter& MultiChunkEquation::iterator() const throw()  
{
  return itsSharedIterator;
}

/// @brief replace iterator to a new one
/// @details This method is probably more temporary than the whole class.
/// It is necessary for composite measurement equations to be able to 
/// substitute the iterator to something else to beat accessor-based/
/// iterator-based ME problem.
/// @param idi data iterator
void MultiChunkEquation::setIterator(const IDataSharedIter& idi)
{
  itsSharedIterator = idi;
}
