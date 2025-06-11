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

// own includes
#include "askap/distributedimager/WorkUnitContainer.h"
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>

ASKAP_LOGGER(logger, ".WorkUnitContainer");

namespace askap {

namespace synthesis {

/// @brief constructor
WorkUnitContainer::WorkUnitContainer() : itsFreqBoundariesValid(false) {}

/// @brief helper method to populate itsFreqBoundaries if it needs an update
/// @details It goes over all stored workunits and appends a pointer (in the form of iterator)
/// if the frequency changes. Adding new elements invalidates itsFreqBoundaries, so we don't
/// need to keep track the validity of individual iterators stored in the vector.
/// @note The method has been declared const because it works with mutable fields (caching scenario)
void WorkUnitContainer::updateFreqBoundariesIfNecessary() const
{
   if (!itsFreqBoundariesValid) {
       itsFreqBoundariesValid = true;
       itsFreqBoundaries.clear();
       if (itsWorkUnits.size() == 0u) {
           return;
       }
       itsFreqBoundaries.reserve(itsWorkUnits.size());
       const_iterator ci = itsWorkUnits.begin();
       for (double currentFreq = (ci++)->get_channelFrequency(); ci != itsWorkUnits.end(); ++ci) {
            // although it is usually bad to rely on floating point exact comparison, it was like that prior to refactoring
            // (and is ok here because this frequency is assigned to the appropriate field without further math)
            if (currentFreq != ci->get_channelFrequency()) {
                itsFreqBoundaries.push_back(ci);
                currentFreq = ci->get_channelFrequency();
            }
       }
   }
}

/// @brief return the number of unique frequencies
/// @details The work units are groupped by frequency channels, but may contain different beams, epochs.
/// This method returns the number of unique frequencies which can be used together with frequency-specific
/// begin and end methods which require the zero-based sequence number of such frequency block.
/// @return the number of unique frequencies across all stored work units. Zero is returned for an empty container.
size_t WorkUnitContainer::numberOfFrequencyBlocks() const
{
   if (itsWorkUnits.size() == 0u) {
       return 0u;
   }
   updateFreqBoundariesIfNecessary();
   return itsFreqBoundaries.size() + 1u;
}

/// @brief stl start iterator over the given frequency block
/// @details This version returns the iterator for the group of work units with unique
/// frequency. There could be many such frequency blocks. The one desired (from 0 to N-1, 
/// where N is the return value of numberOfFrequencyBlocks) is given as a parameter
/// @param[in] block frequency block number
/// @return start iterator for the section of interest
WorkUnitContainer::const_iterator WorkUnitContainer::begin(size_t block) const
{
   // numberOfFrequencyBlocks() will call updateFreqBoundariesIfNecessary
   ASKAPCHECK(block < numberOfFrequencyBlocks(), "Requested frequency block "<<block<<" exceeds the number available");
   if (block == 0u) {
       return itsWorkUnits.begin();
   }
   return itsFreqBoundaries[block - 1];
}

/// @brief stl end iterator for the given frequency block
/// @details This version returns the end iterator for the group of work units with unique
/// frequency. There could be many such frequency blocks. The one desired (from 0 to N-1, 
/// where N is the return value of numberOfFrequencyBlocks) is given as a parameter
/// @param[in] block frequency block number
/// @return end iterator for the section of interest
WorkUnitContainer::const_iterator WorkUnitContainer::end(size_t block) const
{
   // numberOfFrequencyBlocks() will call updateFreqBoundariesIfNecessary
   ASKAPCHECK(block < numberOfFrequencyBlocks(), "Requested frequency block "<<block<<" exceeds the number available");
   if (block == itsFreqBoundaries.size()) {
       return itsWorkUnits.end();
   }
   return itsFreqBoundaries[block];
}

/// @brief add a work unit to the container
/// @details The new unit is prepended to the existing vector of work units
/// (to match the code behaviour prior to refactoring)
/// @param[in] wu work unit to add
void WorkUnitContainer::add(const cp::ContinuumWorkUnit &wu) {
   itsWorkUnits.insert(itsWorkUnits.begin(),wu); 
   itsFreqBoundariesValid = false;
}

/// @brief squash work units with adjacent channels into one work unit
/// @details It modifies the container in situ by merging work units corresponding to 
/// adjacent channels. This allows us to save on processing in the continuum case.
void WorkUnitContainer::mergeAdjacentChannels() 
{
   // current working element, everything before that one has already been processed
   std::vector<cp::ContinuumWorkUnit>::iterator cursorIt = itsWorkUnits.begin();
   for (; cursorIt != itsWorkUnits.end(); ++cursorIt) {
        size_t contiguousCount = 1u;
        int sign = 1;
        std::vector<cp::ContinuumWorkUnit>::const_iterator testIt = cursorIt;
        const unsigned int cursorChan = cursorIt->get_localChannel();
        for (++testIt; testIt != itsWorkUnits.end(); ++testIt, ++contiguousCount) {
             // break if beam or dataset changes (can add more conditions like that later on)
             if ((cursorIt->get_dataset() != testIt->get_dataset()) || (cursorIt->get_beam() != testIt->get_beam())) {
                 break;
             }
             const unsigned int testChan = testIt->get_localChannel();
             ASKAPDEBUGASSERT(cursorChan != testChan);
             if (contiguousCount == 1u) {
                 sign = testChan > cursorChan ? 1 : -1;
             }
             // gap condition
             if (testChan != cursorChan + sign * contiguousCount) {
                 break;
             }
        }
        if (contiguousCount > 1u) {
            // can merge all work units up to (but not including) testIt into the one pointed by cursorIt
            // For now leave frequencies, width, etc untouched and just update the number of channels. This matches the 
            // old behaviour prior to refactoring
            cursorIt->set_nchan(contiguousCount);
            if (sign < 0) {
                // always keep the lowest channel number for the group. By default the unit pointed to by cursorIt will have the first
                ASKAPDEBUGASSERT(cursorChan + 1 >= contiguousCount);
                cursorIt->set_localChannel(cursorChan + 1 - contiguousCount);
            }
            contiguousCount = 1u;
            std::vector<cp::ContinuumWorkUnit>::const_iterator obsoletePartIt = cursorIt;
            // we keep the element pointed by cursorIt, but delete the ones after that
            ++obsoletePartIt;
            ASKAPDEBUGASSERT(obsoletePartIt != itsWorkUnits.end());
            // the following invalidates iterators from obsoletePartIt onwards, including itsWorkUnits.end(), but we don't 
            // cache it (so the call as part of the loop condition is ok). The cursorIt iterator remains valid.
            itsWorkUnits.erase(obsoletePartIt, testIt);
        }
   }
  
   // the original code prior to refactoring updated parset as well from inside this procedure (to select a group of channels rather than one)
   // This needs to be performed elsewhere. I leave this original code here commented out as a reminder. 
   // string ChannelParam = "["+toString(contiguousCount)+","+
   //       toString(compressedWorkUnit.get_localChannel())+"]";
   //       ASKAPLOG_DEBUG_STR(logger, "compressWorkUnit: ChannelParam = "<<ChannelParam);
   //       itsParset.replace("Channels",ChannelParam);
}

/// @brief access work unit by index
/// @details This is hopefully a temporary method to enable early integration with ContinuumWorker. It
/// mimics the direct access as for the std::vector. In the future, it is expected that iterator interfaces
/// would be sufficient.
/// @param[in] index index to work with (should be less than size)
/// @return const reference to the work unit stored at the given index
const cp::ContinuumWorkUnit& WorkUnitContainer::operator[](size_t index) const
{
   ASKAPASSERT(index < itsWorkUnits.size());
   return itsWorkUnits[index];
}

} // namespace synthesis

} // namespace askap
