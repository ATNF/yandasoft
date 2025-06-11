/// @file
///
/// Unit test for the deconvolution base class
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>

#include <askap/deconvolution/DeconvolverBase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <casacore/casa/BasicSL/Complex.h>

#include <boost/shared_ptr.hpp>

using namespace casa;

namespace askap {

namespace synthesis {

class DeconvolverBaseTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(DeconvolverBaseTest);
  CPPUNIT_TEST(testCasacoreAssumptions);
  CPPUNIT_TEST(testCasacoreAssumptions2);
  CPPUNIT_TEST(testCreate);
  CPPUNIT_TEST_EXCEPTION(testWrongShape, AskapError);
  CPPUNIT_TEST_SUITE_END();
public:

  void setUp() {
    IPosition dimensions(2,100,100);
    itsDirty.reset(new Array<Float>(dimensions));
    itsDirty->set(0.0);
    itsPsf.reset(new Array<Float>(dimensions));
    itsPsf->set(0.0);
    (*itsPsf)(IPosition(2,50,50))=1.0;
    itsDB = DeconvolverBase<Float,Complex>::ShPtr(new DeconvolverBase<Float, Complex>(*itsDirty, *itsPsf));
    CPPUNIT_ASSERT(itsDB);
    CPPUNIT_ASSERT(itsDB->control());
    CPPUNIT_ASSERT(itsDB->monitor());
    CPPUNIT_ASSERT(itsDB->state());
    itsWeight.reset(new Array<Float>(dimensions));
    itsWeight->set(10.0);
    itsDB->setWeight(*itsWeight);
  }

  void tearDown() {
      // Ensure arrays are destroyed last
      itsDB.reset();
      itsWeight.reset();
      itsPsf.reset();
      itsDirty.reset();
  }

  void testCasacoreAssumptions() {
     // MV: this test is not about DeconvolverBase, but rather about assumptions
     // about casacore methods used throughout the deconvolve package which were not documented well
     IPosition shape(2,2,3);
     Array<Float> realArr(shape, 1.f);
     Array<Complex> complexArr(shape, Complex(2.,3.));
     CPPUNIT_ASSERT(realArr.contiguousStorage());
     CPPUNIT_ASSERT(complexArr.contiguousStorage());
     for (Array<Float>::const_iterator ci = realArr.begin(); ci != realArr.end(); ++ci) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1., *ci, 1e-10);
     }
     for (Array<Complex>::const_iterator ci = complexArr.begin(); ci != complexArr.end(); ++ci) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(2., real(*ci), 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(3., imag(*ci), 1e-10);
     }
     // this method overwrites just real part leaving imaginary part intact
     setReal(complexArr, realArr);
     for (Array<Complex>::const_iterator ci = complexArr.begin(); ci != complexArr.end(); ++ci) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1., real(*ci), 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(3., imag(*ci), 1e-10);
     }
     // this method just sets all values to the constant
     realArr.set(5.);
     for (Array<Float>::const_iterator ci = realArr.begin(); ci != realArr.end(); ++ci) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(5., *ci, 1e-10);
     }
     // and it shouldn't propagate to the complex array (i.e. no reference semantics here)
     for (Array<Complex>::const_iterator ci = complexArr.begin(); ci != complexArr.end(); ++ci) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1., real(*ci), 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(3., imag(*ci), 1e-10);
     }

     // special array with twice as many elements on the first axis
     Array<Float> realArr2(IPosition(2,shape(0)*2, shape(1)), 5.f);
     int index = 0;
     for (Array<Float>::iterator it = realArr2.begin(); it != realArr2.end(); ++it, ++index) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(5., *it, 1e-10);
          if (index % 2 == 0) {
              *it = 10.;
          }
     }

     // this method packs real array of double length into a complex array
     // we also effectively test that
     RealToComplex(complexArr, realArr2);
     for (Array<Complex>::const_iterator ci = complexArr.begin(); ci != complexArr.end(); ++ci) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(10., real(*ci), 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(5., imag(*ci), 1e-10);
     }

  }

  void testCasacoreAssumptions2() {
     // MV: this test is not about DeconvolverBase, but rather about assumptions
     // about casacore methods used throughout the deconvolve package, in particular
     // about underlying memory layout
     const uInt nplanes = 13;
     const uInt size = 100;
     Cube<Float> tmp(size,size,nplanes,-1.);
     for (uInt plane=0; plane < nplanes; ++plane) {
          CPPUNIT_ASSERT(tmp.xyPlane(plane).contiguousStorage());
          // slice is the reference to a particular zplane
          Matrix<Float> slice = tmp.xyPlane(plane);
          CPPUNIT_ASSERT(slice.contiguousStorage());
          CPPUNIT_ASSERT(slice.column(0).contiguousStorage());
          CPPUNIT_ASSERT(!slice.row(0).contiguousStorage());
          slice.set(static_cast<Float>(plane));
          slice.column(0).set(static_cast<Float>(plane)*plane);
     }
     for (uInt row=0; row < size; ++row) {
          for (uInt column=0; column < size; ++column) {
               for (uInt plane=0; plane < nplanes; ++plane) {
                    const Float expected = column == 0 ? static_cast<Float>(plane) * plane : static_cast<Float>(plane);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, tmp(row, column, plane), 1e-6);
               }
          }
     }
  }

  void testCreate() {
    Array<Float> newDirty(IPosition(4,100,100,1,1));
    itsDB->updateDirty(newDirty);
  }
  void testWrongShape() {
    Array<Float> newDirty(IPosition(4,200,200,1,1));
    itsDB->updateDirty(newDirty);
  }
private:

private:
  boost::shared_ptr< Array<Float> > itsDirty;
  boost::shared_ptr< Array<Float> > itsPsf;
  boost::shared_ptr< Array<Float> > itsWeight;

   /// @brief DeconvolutionBase class
  boost::shared_ptr<DeconvolverBase<Float, Complex> > itsDB;
};

} // namespace synthesis

} // namespace askap
