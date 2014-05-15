/// @file
///
/// Unit test for non-lionear sampling of the w-space
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

#ifndef NON_LINEAR_W_SAMPLING_TEST_H
#define NON_LINEAR_W_SAMPLING_TEST_H

#include <gridding/GaussianWSampling.h>
#include <gridding/PowerWSampling.h>
#include <cppunit/extensions/HelperMacros.h>

namespace askap {

namespace synthesis {

class NonLinearWSamplingTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(NonLinearWSamplingTest);
   CPPUNIT_TEST(testGaussian);
   CPPUNIT_TEST(testPowerLaw);
   CPPUNIT_TEST_SUITE_END();
public:
   void setUp() {}
   
   void testPowerLaw() {
      const int nTrials = 10;
      for (int i = 0; i<nTrials; ++i) {
           if (i*2 == nTrials) {
               break;
           }
           const double exponent = 2.*double(i - nTrials/2)/double(nTrials);
           PowerWSampling smp(exponent);
           doConversionTest(smp);
      }
   }
   
   void testGaussian()  {
      const int nTrials = 10;
      // starting from 3 because too narrow gaussian would genuinly fail the test due to numerical dynamic range of 
      // double (calculating gaussian and inverting it)
      for (int i = 3; i<=nTrials; ++i) {
           const double wplanes50 = double(i)/double(nTrials+1)/sqrt(2); // in the allowed range: 0<wplanes50<1/sqrt(2)
           GaussianWSampling smp(wplanes50);
           doConversionTest(smp);
           // check that wplanes50 parameter indeed means the number of planes in the (-0.5,0.5) range
           CPPUNIT_ASSERT(fabs(wplanes50 - itsMidRange)<0.02);
      }   
   }
   
   void doConversionTest(const IWSampling &smp) {
       const int nTrials = 100;
       itsMidRange = 0.;
       for (int i = 0; i<=nTrials; ++i) {
            const double plane = 2.*double(i - nTrials/2)/double(nTrials);
            const double mapped = smp.map(plane);
            if ((mapped >= -0.5) && (mapped <= 0.5)) {
                 itsMidRange += 1./double(nTrials+1);
            }
            const double result = smp.index(mapped);
            CPPUNIT_ASSERT(fabs(plane - result)<1e-6);
       }
       CPPUNIT_ASSERT(fabs(smp.index(0.))<1e-6);
       CPPUNIT_ASSERT(fabs(smp.index(-1.)+1.)<1e-6);
       CPPUNIT_ASSERT(fabs(smp.index(1.)-1.)<1e-6);
   }
private:
   // counter of the mapped results in the (-0.5,0.5) range
   double itsMidRange;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef NON_LINEAR_W_SAMPLING_TEST_H


