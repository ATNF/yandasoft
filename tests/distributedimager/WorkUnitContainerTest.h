/// @file
///
/// Unit test for WorkUnitContainer class (used in distributed imager)
///
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

#include <askap/distributedimager/WorkUnitContainer.h>

#include <cppunit/extensions/HelperMacros.h>

namespace askap {

namespace synthesis {

class WorkUnitContainerTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(WorkUnitContainerTest);
   CPPUNIT_TEST(testAdd);
   CPPUNIT_TEST(testChannelMerge);
   CPPUNIT_TEST(testChannelMergeWithInversion);
   CPPUNIT_TEST(testChannelMergeGap);
   CPPUNIT_TEST(testChannelMergeDifferentDatasets);
   CPPUNIT_TEST(testChannelMergeDifferentBeams);
   CPPUNIT_TEST(testChannelMergeNotAdjacent);
   CPPUNIT_TEST(testChannelPartition);
   CPPUNIT_TEST_SUITE_END();
public:
   
   void testAdd() {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_dataset("test1.ms");
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(0u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() == wuc.end());
      wuc.add(wu);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
      addWorkUnits(wuc, wu, 5u, 1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(6u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
   }

   void testChannelMerge() {
      testSimpleChannelMergeImpl(1);
   }

   void testChannelMergeWithInversion() {
      testSimpleChannelMergeImpl(-1);
   }

   void testChannelMergeGap() {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_localChannel(1u);
      wu.set_dataset("test1.ms");
      wu.set_beam(0u);
      addWorkUnits(wuc, wu, 5u, 1);
      wu.set_localChannel(10u);
      addWorkUnits(wuc, wu, 5u, 1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(10u), wuc.size());
      wuc.mergeAdjacentChannels();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(2u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
      WorkUnitContainer::const_iterator ci = wuc.begin();
      CPPUNIT_ASSERT(ci != wuc.end());
      const cp::ContinuumWorkUnit firstWu = *ci;
      CPPUNIT_ASSERT_EQUAL(wu.get_dataset(), firstWu.get_dataset());
      CPPUNIT_ASSERT_EQUAL(wu.get_beam(), firstWu.get_beam());
      CPPUNIT_ASSERT_EQUAL(5u, firstWu.get_nchan());
      // local channel should point to the base of the group
      // because wuc.add prepends work units the order is backwards
      CPPUNIT_ASSERT_EQUAL(10u, firstWu.get_localChannel());
      ++ci;
      CPPUNIT_ASSERT(ci != wuc.end());
      const cp::ContinuumWorkUnit secondWu = *ci;
      CPPUNIT_ASSERT_EQUAL(wu.get_dataset(), secondWu.get_dataset());
      CPPUNIT_ASSERT_EQUAL(wu.get_beam(), secondWu.get_beam());
      CPPUNIT_ASSERT_EQUAL(5u, secondWu.get_nchan());
      // local channel should point to the base of the group
      // because wuc.add prepends work units the order is backwards
      CPPUNIT_ASSERT_EQUAL(1u, secondWu.get_localChannel());
      ++ci;
      CPPUNIT_ASSERT(ci == wuc.end());
   }

   void testChannelMergeDifferentDatasets() {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_localChannel(1u);
      wu.set_dataset("test1.ms");
      wu.set_beam(0u);
      addWorkUnits(wuc, wu, 5u, 1);
      wu.set_localChannel(6u);
      wu.set_dataset("test2.ms");
      addWorkUnits(wuc, wu, 5u, 1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(10u), wuc.size());
      wuc.mergeAdjacentChannels();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(2u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
      WorkUnitContainer::const_iterator ci = wuc.begin();
      CPPUNIT_ASSERT(ci != wuc.end());
      // note, the order of workunits is backwards because wuc.add prepands new items
      const cp::ContinuumWorkUnit firstWu = *ci;
      CPPUNIT_ASSERT_EQUAL(wu.get_dataset(), firstWu.get_dataset());
      CPPUNIT_ASSERT_EQUAL(wu.get_beam(), firstWu.get_beam());
      CPPUNIT_ASSERT_EQUAL(5u, firstWu.get_nchan());
      // local channel should point to the base of the group
      CPPUNIT_ASSERT_EQUAL(6u, firstWu.get_localChannel());
      ++ci;
      CPPUNIT_ASSERT(ci != wuc.end());
      const cp::ContinuumWorkUnit secondWu = *ci;
      CPPUNIT_ASSERT_EQUAL(std::string("test1.ms"), secondWu.get_dataset());
      CPPUNIT_ASSERT_EQUAL(wu.get_beam(), secondWu.get_beam());
      CPPUNIT_ASSERT_EQUAL(5u, secondWu.get_nchan());
      // local channel should point to the base of the group
      CPPUNIT_ASSERT_EQUAL(1u, secondWu.get_localChannel());
      ++ci;
      CPPUNIT_ASSERT(ci == wuc.end());
   }

   void testChannelMergeDifferentBeams() {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_localChannel(1u);
      wu.set_dataset("test1.ms");
      wu.set_beam(0u);
      addWorkUnits(wuc, wu, 5u, 1);
      wu.set_localChannel(6u);
      wu.set_beam(1u);
      addWorkUnits(wuc, wu, 5u, 1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(10u), wuc.size());
      wuc.mergeAdjacentChannels();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(2u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
      WorkUnitContainer::const_iterator ci = wuc.begin();
      CPPUNIT_ASSERT(ci != wuc.end());
      // note, the order of workunits is backwards because wuc.add prepands new items
      const cp::ContinuumWorkUnit firstWu = *ci;
      CPPUNIT_ASSERT_EQUAL(wu.get_dataset(), firstWu.get_dataset());
      CPPUNIT_ASSERT_EQUAL(wu.get_beam(), firstWu.get_beam());
      CPPUNIT_ASSERT_EQUAL(5u, firstWu.get_nchan());
      // local channel should point to the base of the group
      CPPUNIT_ASSERT_EQUAL(6u, firstWu.get_localChannel());
      ++ci;
      CPPUNIT_ASSERT(ci != wuc.end());
      const cp::ContinuumWorkUnit secondWu = *ci;
      CPPUNIT_ASSERT_EQUAL(wu.get_dataset(), secondWu.get_dataset());
      CPPUNIT_ASSERT_EQUAL(0u, secondWu.get_beam());
      CPPUNIT_ASSERT_EQUAL(5u, secondWu.get_nchan());
      // local channel should point to the base of the group
      CPPUNIT_ASSERT_EQUAL(1u, secondWu.get_localChannel());
      ++ci;
      CPPUNIT_ASSERT(ci == wuc.end());
   }

   void testChannelMergeNotAdjacent() {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_localChannel(1u);
      wu.set_dataset("test1.ms");
      wu.set_beam(0u);
      addWorkUnits(wuc, wu, 5u, 2);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(5u), wuc.size());
      wuc.mergeAdjacentChannels();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(5u), wuc.size());
   }

   void testChannelPartition() {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_localChannel(10u);
      wu.set_dataset("test1.ms");
      wu.set_beam(0u);
      wu.set_channelFrequency(1.421e9);
      addWorkUnits(wuc, wu, 5u, -1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wuc.numberOfFrequencyBlocks());
      CPPUNIT_ASSERT_EQUAL(5u, verifyFrequencyInWorkUnits(wuc.begin(0u), wuc.end(0u), 1.421e9));
      wu.set_channelFrequency(1.420e9);
      addWorkUnits(wuc, wu, 5u, -1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(2u), wuc.numberOfFrequencyBlocks());
      // we're prepending new data, so the index will also change
      CPPUNIT_ASSERT_EQUAL(5u, verifyFrequencyInWorkUnits(wuc.begin(1u), wuc.end(1u), 1.421e9));
      CPPUNIT_ASSERT_EQUAL(5u, verifyFrequencyInWorkUnits(wuc.begin(0u), wuc.end(0u), 1.420e9));
      wu.set_channelFrequency(1.419e9);
      addWorkUnits(wuc, wu, 5u, -1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(3u), wuc.numberOfFrequencyBlocks());
      CPPUNIT_ASSERT_EQUAL(5u, verifyFrequencyInWorkUnits(wuc.begin(2u), wuc.end(2u), 1.421e9));
      CPPUNIT_ASSERT_EQUAL(5u, verifyFrequencyInWorkUnits(wuc.begin(1u), wuc.end(1u), 1.420e9));
      CPPUNIT_ASSERT_EQUAL(5u, verifyFrequencyInWorkUnits(wuc.begin(0u), wuc.end(0u), 1.419e9));
   }

protected:
   /// @brief helper method testing common frequency in the set of work units
   /// @param[in] begin start iterator
   /// @param[in] end end iterator
   /// @param[in] freq expected frequency
   /// @return the number of elements checked
   /// @note the frequency is checked with 1e-13 tolerance, although exact comparison would be acceptable
   /// in the current circumstances because no math is done with the numbers.
   template<typename It>
   static unsigned int verifyFrequencyInWorkUnits(It begin, It end, double freq) {
       unsigned int count = 0u;
       for (It it = begin; it != end; ++it, ++count) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(freq, it->get_channelFrequency(), 1e-13);
       }
       return count;
   }

   /// @brief helper method doing the test of work units merge
   /// @param[in] increment channel increment to use (+1 or -1)
   void testSimpleChannelMergeImpl(int increment) {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_localChannel(increment > 0 ? 1u : 5u);
      wu.set_dataset("test1.ms");
      wu.set_beam(0u);
      CPPUNIT_ASSERT_EQUAL(1u, wu.get_nchan());
      addWorkUnits(wuc, wu, 5u, increment);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(5u), wuc.size());
      wuc.mergeAdjacentChannels();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
      const cp::ContinuumWorkUnit newWu = *wuc.begin();
      CPPUNIT_ASSERT_EQUAL(wu.get_dataset(), newWu.get_dataset());
      CPPUNIT_ASSERT_EQUAL(wu.get_beam(), newWu.get_beam());
      CPPUNIT_ASSERT_EQUAL(5u, newWu.get_nchan());
   }
 
   /// @brief helper method to add work units for a range of local channels
   /// @details 
   /// @param[in] wuc container to work with
   /// @param[in] wuTemplate  template work unit for the units to be added (only local channel 
   ///                        is altered by adding increments)
   /// @param[in] nchan number of channels 
   /// @param[in] increment local channel increment for the added work units
   static void addWorkUnits(WorkUnitContainer& wuc, const cp::ContinuumWorkUnit& wuTemplate, unsigned int nchan, int increment) {
      for (unsigned int chan = 0; chan < nchan; ++chan) {
           // add backwards to account for the way wuc.add works
           const int adjustment = increment * chan;
           const unsigned int startChan = wuTemplate.get_localChannel();
           if (adjustment < 0) {
               CPPUNIT_ASSERT(startChan >= -adjustment);
           }
           cp::ContinuumWorkUnit newWu(wuTemplate);
           newWu.set_localChannel(startChan + adjustment);
           wuc.add(newWu);
      }
   }
   

};
    
} // namespace synthesis

} // namespace askap

