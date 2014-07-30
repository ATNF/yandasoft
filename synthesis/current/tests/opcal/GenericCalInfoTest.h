/// @file
/// 
/// @brief Unit tests for GenericCalInfo class.
/// @details GenericCalInfo encapsulates results of running calibration 
/// procedures on a set of data
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

#ifndef SYNTHESIS_GENERIC_CAL_INFO_TEST_H
#define SYNTHESIS_GENERIC_CAL_INFO_TEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <opcal/GenericCalInfo.h>
#include <askap/AskapError.h>

namespace askap
{
  namespace synthesis
  {
    
    class GenericCalInfoTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(GenericCalInfoTest);
      CPPUNIT_TEST(testCreate);
      CPPUNIT_TEST_EXCEPTION(testAccessUninitialised1, askap::CheckError);
      CPPUNIT_TEST_EXCEPTION(testAccessUninitialised2, askap::CheckError);
      CPPUNIT_TEST_SUITE_END();
    public:
       
      void testCreate() {
         GenericCalInfo info1;
         CPPUNIT_ASSERT(!info1.gainDefined());
         CPPUNIT_ASSERT(!info1.delayDefined());
         
         GenericCalInfo info2(casa::Complex(-1.,2.),3.);
         CPPUNIT_ASSERT(info2.gainDefined());
         CPPUNIT_ASSERT(info2.delayDefined());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-1., static_cast<double>(real(info2.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2., static_cast<double>(imag(info2.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3., info2.delay(), 1e-6);
         
         GenericCalInfo info3(info2);
         CPPUNIT_ASSERT(info3.gainDefined());
         CPPUNIT_ASSERT(info3.delayDefined());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-1., static_cast<double>(real(info3.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2., static_cast<double>(imag(info3.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3., info3.delay(), 1e-6);
         
         info3 = info1;
         CPPUNIT_ASSERT(!info3.gainDefined());
         CPPUNIT_ASSERT(!info3.delayDefined());

         info3.setGain(casa::Complex(2., -1.));
         CPPUNIT_ASSERT(info3.gainDefined());
         CPPUNIT_ASSERT(!info3.delayDefined());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2., static_cast<double>(real(info3.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-1., static_cast<double>(imag(info3.gain())), 1e-6);
         info3.setDelay(-5.);
         CPPUNIT_ASSERT(info3.delayDefined());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2., static_cast<double>(real(info3.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-1., static_cast<double>(imag(info3.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-5., info3.delay(), 1e-6);
                  
         info3 = info2;
         CPPUNIT_ASSERT(info3.gainDefined());
         CPPUNIT_ASSERT(info3.delayDefined());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-1., static_cast<double>(real(info3.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2., static_cast<double>(imag(info3.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3., info3.delay(), 1e-6); 
         
         info3.invalidate();
         CPPUNIT_ASSERT(!info3.gainDefined());
         CPPUNIT_ASSERT(!info3.delayDefined());         
      }
      
      void testAccessUninitialised1() {
         GenericCalInfo info1;
         info1.setGain(casa::Complex(-1.,2.));
         CPPUNIT_ASSERT(info1.gainDefined());
         CPPUNIT_ASSERT(!info1.delayDefined());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-1., static_cast<double>(real(info1.gain())), 1e-6);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2., static_cast<double>(imag(info1.gain())), 1e-6);
         // the following should trigger an exception
         info1.delay();         
      }

      void testAccessUninitialised2() {
         GenericCalInfo info1;
         info1.setDelay(3.);
         CPPUNIT_ASSERT(!info1.gainDefined());
         CPPUNIT_ASSERT(info1.delayDefined());
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3., info1.delay(), 1e-6);         
         // the following should trigger an exception
         info1.gain();         
      }          
    };

  } // namespace synthesis
} // namespace askap

#endif // #ifndef SYNTHESIS_GENERIC_CAL_INFO_TEST_H


