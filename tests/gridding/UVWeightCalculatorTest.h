/// @file
///
/// Unit test for various procedures computing uv-weight or doing other operations with it
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

#include <askap/gridding/ConjugatesAdderFFT.h>
#include <askap/gridding/RobustUVWeightCalculator.h>
#include <askap/gridding/CompositeUVWeightCalculator.h>
#include <askap/gridding/ReciprocalUVWeightCalculator.h>
#include <cppunit/extensions/HelperMacros.h>

// for saveAsCasaImage (it was needed for debugging)
//#include <askap/scimath/utils/ImageUtils.h>

namespace askap {

namespace synthesis {

class UVWeightCalculatorTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(UVWeightCalculatorTest);
   CPPUNIT_TEST(testConjugatesAdderFFT);
   CPPUNIT_TEST(testRobustWeights);
   CPPUNIT_TEST(testCompositeUVWeightCalculator);
   CPPUNIT_TEST(testReciprocalUVWeightCalculator);
   CPPUNIT_TEST_SUITE_END();
public:

   void testConjugatesAdderFFT() {
        casacore::Int pos[4][2] = {{200, 233}, {242, 324}, {270, 300}, {300, 222}};
        casacore::Matrix<float> buf(512, 512, 0.f);
        for (casacore::uInt pt = 0; pt < 4; ++pt) {
             buf(pos[pt][0], pos[pt][1]) = 2.f * pt / 6;
        }
 
        ConjugatesAdderFFT calc;
        calc.process(buf);
        
        for (casacore::uInt pt = 0; pt < 4; ++pt) {
             // factor of 2 still remains because we added conjugates as opposed to just taking the real part, so this way of doing it preserves the values
             CPPUNIT_ASSERT_DOUBLES_EQUAL(2.f*pt/6, buf(pos[pt][0], pos[pt][1]), 1e-6);
             // check conjugates as well (note, we didn't simulate any point on the 0th row or column which do not have conjugates with even-sized grids)
             CPPUNIT_ASSERT_DOUBLES_EQUAL(2.f*pt/6, buf(512 - pos[pt][0], 512 - pos[pt][1]), 1e-6);
        }
   }

   void testRobustWeights() {
      testRobustWeights(-2.);
      testRobustWeights(0.);
      testRobustWeights(+2.);
   }

   void testRobustWeights(float robustness) {
        casacore::Matrix<float> buf(512, 512, 0.f);

        for (casacore::uInt row = 0; row < buf.nrow(); ++row) {
             const float x = (static_cast<float>(row) - buf.nrow() / 2) / 2. / 50.;
             for (casacore::uInt col = 0; col < buf.ncolumn(); ++col) {
                  const float y = (static_cast<float>(col) - buf.ncolumn() / 2) / 2. / 75.;
                  buf(row, col) = exp(-x*x - y*y);
             }
        }
        const float avgWt = casacore::sum(buf*buf) / casacore::sum(buf);       
        
        RobustUVWeightCalculator calc(robustness);
        calc.process(buf);

        // need to check the result
        // first undo the inversion (we put it there because the weights are applied via multiplication)
        // by nature of the Wiener filter all values should be non-zero. But we check just in case.
        for (casacore::uInt row = 0; row < buf.nrow(); ++row) {
             for (casacore::uInt col = 0; col < buf.ncolumn(); ++col) {
                  const float val = buf(row, col);
                  CPPUNIT_ASSERT(val > 0.f);
                  buf(row,col) = 1.f / val;
             }
        }
        // now unroll the formula for the filter
        buf -= 1.f;
        // normalise to simplify the comparison
        const float peak = casacore::max(buf);       
        CPPUNIT_ASSERT(peak > 0.);
        buf /= peak;
        for (casacore::uInt row = 0; row < buf.nrow(); ++row) {
             const float x = (static_cast<float>(row) - buf.nrow() / 2) / 2. / 50.;
             for (casacore::uInt col = 0; col < buf.ncolumn(); ++col) {
                  const float y = (static_cast<float>(col) - buf.ncolumn() / 2) / 2. / 75.;
                  const float expected = exp(-x*x - y*y);
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, buf(row,col), 2e-5);
             }
        }
        // check normalisation factor to match formula for Robust weighting with the given robustness
        // generous tolerance accounts for the fact that the math is done with single precision for the test
        CPPUNIT_ASSERT_DOUBLES_EQUAL(25.f / avgWt * pow(10., -2.*robustness), peak, 1e-4 * peak);
        
        //scimath::saveAsCasaImage("tst.img", buf);
   }

   void testReciprocalUVWeightCalculator() {
        casacore::Matrix<float> buf(128, 64, 1.f);
        for (casacore::uInt row = 0; row < buf.nrow(); ++row) {
             const float x = (static_cast<float>(row) - buf.nrow() / 2);
             for (casacore::uInt col = 0; col < buf.ncolumn(); ++col) {
                  const float y = (static_cast<float>(col) - buf.ncolumn() / 2);
                  buf(row, col) += x*x + y*y;
             }
        }

        ReciprocalUVWeightCalculator calc(3.f);
        calc.process(buf);

        // checking the result
        casacore::uInt zeroCount = 0u;
        for (casacore::uInt row = 0; row < buf.nrow(); ++row) {
             const float x = (static_cast<float>(row) - buf.nrow() / 2);
             for (casacore::uInt col = 0; col < buf.ncolumn(); ++col) {
                  const float y = (static_cast<float>(col) - buf.ncolumn() / 2);
                  const float val = buf(row, col);
                  if (val > 0.f) {
                      const float expected = (1.f + x*x + y*y);
                      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, 1.f / val, 2e-6 * expected);
                  } else {
                      ++zeroCount;
                  }
             }
        }
        CPPUNIT_ASSERT_EQUAL(9u, zeroCount);
   }

   void testCompositeUVWeightCalculator() {
        casacore::Matrix<float> buf(128, 128, 0.f);

        CompositeUVWeightCalculator calc;
        boost::shared_ptr<DummyUVWeightCalculator> dummy1(new DummyUVWeightCalculator(1., 2.));
        boost::shared_ptr<DummyUVWeightCalculator> dummy2(new DummyUVWeightCalculator(2., 3.));
        boost::shared_ptr<DummyUVWeightCalculator> dummy3(new DummyUVWeightCalculator(-1., 0.5));
        calc.add(dummy1);
        calc.add(dummy2);
        calc.add(dummy3);
        calc.process(buf);
        for (casacore::uInt row = 0; row < buf.nrow(); ++row) {
             for (casacore::uInt col = 0; col < buf.ncolumn(); ++col) {
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.5f, buf(row,col), 1e-6);
             }
        }
   }

private:
   struct DummyUVWeightCalculator : virtual public IUVWeightCalculator {
  
      /// @bief constructor to set up the linear transformation (add a number, multiply by another)
      /// @details It is handy to do it this way (as opposed to simple multiplication and addition), 
      /// so we can test the order of operations as well because in
      /// general such operations will not commute
      /// @param[in] add the number of to add
      /// @param[in] mul the number to multiply by
      explicit DummyUVWeightCalculator(float add, float mul) : itsNumber2Add(add), itsNumber2Mul(mul) {} 

      /// @brief dummy processing, just adds the number it was setup with to all elements and multiplies by another one     
      /// @param[in] wt weight to work with (it is modified in situ).
      /// @note The shape is supposed to stay intact.
      virtual void process(casacore::Matrix<float> &wt) const override {
         wt += itsNumber2Add;
         wt *= itsNumber2Mul;
      }
   private:
      /// @brief dummy weight number to add
      const float itsNumber2Add;

      /// @brief dummy number to multiply by
      const float itsNumber2Mul;
   };

};
    
} // namespace synthesis

} // namespace askap

