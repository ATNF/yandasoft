/// @file
///
/// Unit test for the frequency axis mapping class
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

#include <gridding/FrequencyMapper.h>
#include <cppunit/extensions/HelperMacros.h>

#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/IPosition.h>
#include <casa/BasicSL/Complex.h>

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

class FrequencyMapperTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(FrequencyMapperTest);
   CPPUNIT_TEST(testFreqMapping);
   CPPUNIT_TEST(testSetupImage);
   CPPUNIT_TEST(testApproxFreqMapping);
   CPPUNIT_TEST(testMFSFreqMapping);
   CPPUNIT_TEST_SUITE_END();
public:

   void setUp() {
       scimath::Axes axis;
       axis.add("FREQUENCY",1.308e9,1.42e9);
       itsFreqMapper.reset(new FrequencyMapper(axis, 8));
   }
   
   void testSetupImage() {
       casa::Vector<casa::Double> freqs(10);
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
            freqs[i] = 1.308e9+1.6e7*(double(i)-1.);
       } 
       CPPUNIT_ASSERT(itsFreqMapper);
       // MFS mode first
       itsFreqMapper->setupSinglePlaneGridding();
       itsFreqMapper->setupMapping(freqs);
       
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
                CPPUNIT_ASSERT(itsFreqMapper->isMapped(i));                
                CPPUNIT_ASSERT_EQUAL(0u,(*itsFreqMapper)(i));       
       }
       // now we do new setup image call, which should revert the
       // frequency mapper back to the spectral line mode
       scimath::Axes axis;
       axis.add("FREQUENCY",1.308e9,1.42e9);
       itsFreqMapper->setupImage(axis,8);
       itsFreqMapper->setupMapping(freqs);
       // test the mapping
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
            if (i==0 || i+1 == freqs.nelements()) {
                CPPUNIT_ASSERT(!itsFreqMapper->isMapped(i));                
            } else {
                CPPUNIT_ASSERT(itsFreqMapper->isMapped(i));                
                CPPUNIT_ASSERT_EQUAL(i,(*itsFreqMapper)(i) + 1);                
            }
       }
   }

   void testApproxFreqMapping() {
       casa::Vector<casa::Double> freqs(10);
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
            // this time with a negative increment and small wobble in frequency
            freqs[i] = 1.436e9-1.6e7*double(i) + (i%2 == 0? +1. : -1.)*1e6;
       } 
       CPPUNIT_ASSERT(itsFreqMapper);
       itsFreqMapper->setupMapping(freqs);
       // test the mapping
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
            if (i==0 || i+1 == freqs.nelements()) {
                CPPUNIT_ASSERT(!itsFreqMapper->isMapped(i));                
            } else {
                CPPUNIT_ASSERT(itsFreqMapper->isMapped(i));                
                CPPUNIT_ASSERT_EQUAL(casa::uInt(freqs.nelements())-i-1u,(*itsFreqMapper)(i) + 1);                
            }
       }
   }

   void testMFSFreqMapping() {
       casa::Vector<casa::Double> freqs(20);
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
            //small increment, so two measured channels map into one image channel
            freqs[i] = 1.304e9+8e6*(double(i)-2.);
       } 
       CPPUNIT_ASSERT(itsFreqMapper);
       itsFreqMapper->setupMapping(freqs);
       // test the mapping
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
            if (i<2 || i+2 >= freqs.nelements()) {
                CPPUNIT_ASSERT(!itsFreqMapper->isMapped(i));                
            } else {
                CPPUNIT_ASSERT(itsFreqMapper->isMapped(i));                
                CPPUNIT_ASSERT_EQUAL(i/2,(*itsFreqMapper)(i) + 1);                
            }
       }
   }

   void testFreqMapping() {
       casa::Vector<casa::Double> freqs(10);
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
            freqs[i] = 1.308e9+1.6e7*(double(i)-1.);
       } 
       CPPUNIT_ASSERT(itsFreqMapper);
       itsFreqMapper->setupMapping(freqs);
       // test the mapping
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
            if (i==0 || i+1 == freqs.nelements()) {
                CPPUNIT_ASSERT(!itsFreqMapper->isMapped(i));                
            } else {
                CPPUNIT_ASSERT(itsFreqMapper->isMapped(i));                
                CPPUNIT_ASSERT_EQUAL(i,(*itsFreqMapper)(i) + 1);                
            }
       }
       
       // test MFS mode
       itsFreqMapper->setupSinglePlaneGridding();
       itsFreqMapper->setupMapping(freqs);
       
       for (casa::uInt i=0; i<freqs.nelements(); ++i) {
                CPPUNIT_ASSERT(itsFreqMapper->isMapped(i));                
                CPPUNIT_ASSERT_EQUAL(0u,(*itsFreqMapper)(i));       
       }
   }
   
   
private:
   /// @brief frequency mapping class
   boost::shared_ptr<FrequencyMapper> itsFreqMapper;         
};
    
} // namespace synthesis

} // namespace askap

