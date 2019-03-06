/// @file
/// 
/// @brief Unit tests for RestoringBeamHelper.
/// @details RestoringBeamHelper encapsulates operations with the beam
/// parameters to allow us either fit them or set up explicitly.
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

#ifndef RESTORING_BEAM_HELPER_TEST_H
#define RESTORING_BEAM_HELPER_TEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <measurementequation/RestoringBeamHelper.h>

namespace askap
{
  namespace synthesis
  {
    
    class RestoringBeamHelperTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(RestoringBeamHelperTest);
      CPPUNIT_TEST(testCreateExplicit);
      CPPUNIT_TEST(testCreateForDelayedFit);
      CPPUNIT_TEST_SUITE_END();
    public:
       
      void testCreateExplicit() {
         RestoringBeamHelper rbh;
         CPPUNIT_ASSERT(!rbh.valid());
         casa::Vector<casa::Quantum<double> > beam(3);
         beam[0] = casa::Quantity(0.001,"rad");
         beam[1] = casa::Quantity(0.002,"rad");
         beam[2] = casa::Quantity(-1.1,"rad");         
         rbh.assign(beam);
         CPPUNIT_ASSERT(rbh.valid());
         CPPUNIT_ASSERT(!rbh.fitRequired());

         RestoringBeamHelper rbh2(beam);
         CPPUNIT_ASSERT(rbh2.valid());
         CPPUNIT_ASSERT(!rbh2.fitRequired());
         beam.resize(0);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1e-3, rbh.value()[0].getValue(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2e-3, rbh.value()[1].getValue(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.1, rbh.value()[2].getValue(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1e-3, rbh2.value()[0].getValue(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2e-3, rbh2.value()[1].getValue(), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.1, rbh2.value()[2].getValue(), 1e-6);
      }
    
      void testCreateForDelayedFit() {
         RestoringBeamHelper rbh;
         CPPUNIT_ASSERT(!rbh.valid());
         rbh.configureFit(0.05);
         CPPUNIT_ASSERT(rbh.valid());
         CPPUNIT_ASSERT(rbh.fitRequired());
      }
    };

  } // namespace synthesis
} // namespace askap

#endif // #ifndef RESTORING_BEAM_HELPER_TEST_H

