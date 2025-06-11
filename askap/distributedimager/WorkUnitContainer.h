/// @file
/// @brief Class managing a collection of work units
/// @details This class is intended to encapsulate management of a vector of work units
/// (adding elements, various preprocessing like squashing channels, iteration and breaking
/// iteration in a certain way). It allows us to separate iteration logic from the rest of
/// the algorithm and thus to have a cleaner code. 
///
/// @copyright (c) 2023 CSIRO
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

#ifndef ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CONTAINER_H
#define ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CONTAINER_H

// own includes
#include <askap/messages/ContinuumWorkUnit.h>

// std includes
#include <vector>

// boost includes
#include <boost/noncopyable.hpp>

namespace askap {

namespace synthesis {

/// @brief Class managing a collection of work units
/// @details This class is intended to encapsulate management of a vector of work units
/// (adding elements, various preprocessing like squashing channels, iteration and breaking
/// iteration in a certain way). It allows us to separate iteration logic from the rest of
/// the algorithm and thus to have a cleaner code. 
/// @ingroup distributedimager
class WorkUnitContainer : public boost::noncopyable { 
public:
   /// @brief const iterator type over all stored work units 
   typedef std::vector<cp::ContinuumWorkUnit>::const_iterator const_iterator;

   /// @brief constructor
   WorkUnitContainer();

   /// @brief stl start iterator over the whole container
   /// @return stl-compliant iterator pointing to the first element
   inline const_iterator begin() const { return itsWorkUnits.begin(); }

   /// @brief stl end iterator for iteration over the whole container
   /// @return stl-compliant iterator pointing to the end of the sequence 
   /// (i.e. an imaginary element after the last one)
   inline const_iterator end() const { return itsWorkUnits.end(); }

   // similar methods to the ones below can be implemented for other types of partitioning
   // (e.g. by direction) if needed. Perhaps, better naming would be required in this case.
   // Also it is possible to embed this functionality into the iterator instead of searching 
   // for the appropriate section on-the-fly (but the iterator type would need to be changed).
   // Such an approach would be more flexible and would make it possible to iterate over
   // non-contigous selection of work units.

   /// @brief stl start iterator over the given frequency block
   /// @details This version returns the iterator for the group of work units with unique
   /// frequency. There could be many such frequency blocks. The one desired (from 0 to N-1, 
   /// where N is the return value of numberOfFrequencyBlocks) is given as a parameter
   /// @param[in] block frequency block number
   /// @return start iterator for the section of interest
   const_iterator begin(size_t block) const;

   /// @brief stl end iterator for the given frequency block
   /// @details This version returns the end iterator for the group of work units with unique
   /// frequency. There could be many such frequency blocks. The one desired (from 0 to N-1, 
   /// where N is the return value of numberOfFrequencyBlocks) is given as a parameter
   /// @param[in] block frequency block number
   /// @return end iterator for the section of interest
   const_iterator end(size_t block) const;

   /// @brief add a work unit to the container
   /// @details The new unit is prepended to the existing vector of work units
   /// (to match the code behaviour prior to refactoring)
   /// @param[in] wu work unit to add
   void add(const cp::ContinuumWorkUnit &wu);

   /// @brief squash work units with adjacent channels into one work unit
   /// @details It modifies the container in situ by merging work units corresponding to 
   /// adjacent channels. This allows us to save on processing in the continuum case.
   void mergeAdjacentChannels(); 

   /// @brief obtain the number of work units in the container
   /// @return the current number of stored workunits
   inline size_t size() const { return itsWorkUnits.size(); }

   /// @brief return the number of unique frequencies
   /// @details The work units are groupped by frequency channels, but may contain different beams, epochs.
   /// This method returns the number of unique frequencies which can be used together with frequency-specific
   /// begin and end methods which require the zero-based sequence number of such frequency block.
   /// @return the number of unique frequencies across all stored work units. Zero is returned for an empty container.
   size_t numberOfFrequencyBlocks() const;

   /// @brief access work unit by index
   /// @details This is hopefully a temporary method to enable early integration with ContinuumWorker. It
   /// mimics the direct access as for the std::vector. In the future, it is expected that iterator interfaces
   /// would be sufficient.
   /// @param[in] index index to work with (should be less than size)
   /// @return const reference to the work unit stored at the given index
   const cp::ContinuumWorkUnit& operator[](size_t index) const;

private:
   /// @brief helper method to populate itsFreqBoundaries if it needs an update
   /// @details It goes over all stored workunits and appends a pointer (in the form of iterator)
   /// if the frequency changes. Adding new elements invalidates itsFreqBoundaries, so we don't
   /// need to keep track the validity of individual iterators stored in the vector.
   /// @note The method has been declared const because it works with mutable fields (caching scenario)
   void updateFreqBoundariesIfNecessary() const;

   /// @brief vector of work units 
   /// @note It may be more appropriate to use some other container type here based on the usage, but at this stage 
   /// mimic the original code prior to refactoring to minimise work
   std::vector<cp::ContinuumWorkUnit> itsWorkUnits;

   /// @brief true if itsFreqBoundaries is valid
   mutable bool itsFreqBoundariesValid;

   /// @brief iterators to the of the first element with each new frequency
   /// @details The number of sections with unique frequency is one more than
   /// the size of this vector (i.e. the very first is pointed to by the iterator 
   /// returned by begin() method), so if all work units correspond to the same frequency
   /// this vector will be empty. It is only valid if itsFreqBoundariesValid is true
   mutable std::vector<const_iterator> itsFreqBoundaries;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CONTAINER_H

