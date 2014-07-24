/// @file
/// 
/// @brief Unit tests for ObservationDescription class.
/// @details ObservationDescription encapsulates parameters of a scan,
/// i.e. a collection of data with the same pointing and beam and
/// contiguous in time.
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

#ifndef SYNTHESIS_OBSERVATION_DESCRIPTION_TEST_H
#define SYNTHESIS_OBSERVATION_DESCRIPTION_TEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <measurementequation/ObservationDescription.h>

namespace askap
{
  namespace synthesis
  {
    
    class ObservationDescriptionTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(ObservationDescriptionTest);
      CPPUNIT_TEST(testCreate);
      CPPUNIT_TEST_EXCEPTION(testAccessUninitialised, askap::CheckError);
      CPPUNIT_TEST_EXCEPTION(testIllegalUpdate, askap::CheckError);
      CPPUNIT_TEST_SUITE_END();
    public:
       
      void testCreate() {
         ObservationDescription obs;
         CPPUNIT_ASSERT(!obs.isValid());
         const casa::MVDirection dir(casa::Quantity(1.01,"rad"), casa::Quantity(-1.1,"rad"));
         ObservationDescription obs1("test",5,3600.5, 3, dir, 1e9);
         CPPUNIT_ASSERT(obs1.name() == "test");
         CPPUNIT_ASSERT(obs1.isValid());
         obs.set("another",3, 1800.1, 0, dir, 8.5e8);
         CPPUNIT_ASSERT(obs.isValid());
         CPPUNIT_ASSERT(obs.name() == "another");
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3600.5, obs1.startTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3600.5, obs1.endTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., obs1.direction().separation(dir), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1e9, obs1.frequency(), 1);
         CPPUNIT_ASSERT_EQUAL(5u, obs1.startCycle());
         CPPUNIT_ASSERT_EQUAL(5u, obs1.endCycle());         
         CPPUNIT_ASSERT_EQUAL(3u, obs1.beam());
         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1800.1, obs.startTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1800.1, obs.endTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., obs.direction().separation(dir), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(8.5e8, obs.frequency(), 1);
         CPPUNIT_ASSERT_EQUAL(3u, obs.startCycle());
         CPPUNIT_ASSERT_EQUAL(3u, obs.endCycle());
         CPPUNIT_ASSERT_EQUAL(0u, obs.beam());
         
         obs.update(6,1900.1);
         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1800.1, obs.startTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1900.1, obs.endTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., obs.direction().separation(dir), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(8.5e8, obs.frequency(), 1);         
         CPPUNIT_ASSERT_EQUAL(3u, obs.startCycle());
         CPPUNIT_ASSERT_EQUAL(6u, obs.endCycle());
         CPPUNIT_ASSERT_EQUAL(0u, obs.beam());
         
         // assignment operator
         obs1 = obs;
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1800.1, obs1.startTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1900.1, obs1.endTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., obs1.direction().separation(dir), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(8.5e8, obs1.frequency(), 1);         
         CPPUNIT_ASSERT_EQUAL(3u, obs1.startCycle());
         CPPUNIT_ASSERT_EQUAL(6u, obs1.endCycle());
         CPPUNIT_ASSERT_EQUAL(0u, obs1.beam());
         CPPUNIT_ASSERT(obs1.name() == "another");
         
         
         // copy constructor
         ObservationDescription obs2(obs);
         CPPUNIT_ASSERT(obs2.isValid());         
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1800.1, obs2.startTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1900.1, obs2.endTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., obs2.direction().separation(dir), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(8.5e8, obs2.frequency(), 1);                  
         CPPUNIT_ASSERT_EQUAL(3u, obs2.startCycle());
         CPPUNIT_ASSERT_EQUAL(6u, obs2.endCycle());
         CPPUNIT_ASSERT_EQUAL(0u, obs2.beam());
         CPPUNIT_ASSERT(obs2.name() == "another");                  
      }
      
      void testAccessUninitialised() {
         ObservationDescription obs;
         obs.startCycle();
      }
      
      void testIllegalUpdate() {
         const casa::MVDirection dir(casa::Quantity(1.01,"rad"), casa::Quantity(-1.1,"rad"));
         ObservationDescription obs("test",5,3600.5, 3, dir, 1e9);
         CPPUNIT_ASSERT(obs.isValid());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3600.5, obs.startTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3600.5, obs.endTime(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., obs.direction().separation(dir), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1e9, obs.frequency(), 1);         
         CPPUNIT_ASSERT_EQUAL(5u, obs.startCycle());
         CPPUNIT_ASSERT_EQUAL(5u, obs.endCycle());         
         CPPUNIT_ASSERT_EQUAL(3u, obs.beam());
         obs.update(6, 3000.);
      }
    
    };

  } // namespace synthesis
} // namespace askap

#endif // #ifndef SYNTHESIS_OBSERVATION_DESCRIPTION_TEST_H

