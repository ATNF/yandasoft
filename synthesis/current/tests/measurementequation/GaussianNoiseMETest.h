/// @file
/// 
/// @brief Unit tests for GaussianNoiseME
/// @details GaussianNoiseME is the measurement equation designed to simulate gaussian 
/// noise in visibilities. This class checks that the statistics are as expected.
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

#ifndef GAUSSIAN_NOISE_ME_TEST_H
#define GAUSSIAN_NOISE_ME_TEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <dataaccess/DataAccessorStub.h>
#include <measurementequation/GaussianNoiseME.h>

#include <casa/BasicSL/Complex.h>
#include <casa/Arrays/Vector.h>


namespace askap
{
namespace synthesis
{

  class GaussianNoiseMETest : public CppUnit::TestFixture
  {
     CPPUNIT_TEST_SUITE(GaussianNoiseMETest);
     CPPUNIT_TEST(testNoiseStatistics);
     CPPUNIT_TEST_SUITE_END();
  public:
     void testNoiseStatistics() 
     { 
       accessors::DataAccessorStub acc(true);
       const double variance = 0.1;
       const casa::uInt nRuns = 100;
       GaussianNoiseME me(variance);
       double s = 0., s2 = 0., riCovar = 0.;
       casa::uInt cnt = 0;
       for (casa::uInt run=0; run<nRuns; ++run) {
            me.predict(acc);
            const casa::Vector<casa::Complex> visVec = 
                  acc.rwVisibility().reform(casa::IPosition(1,acc.visibility().nelements()));
            for (casa::uInt elem = 0; elem<visVec.nelements(); ++elem) {
                 const casa::Complex cVal = visVec[elem];
                 s += casa::real(cVal)+casa::imag(cVal);
                 s2 += casa::norm(cVal);
                 riCovar += casa::real(cVal)*casa::imag(cVal);
                 ++cnt;
            }
       }
       // we don't really care about the number of elements in the accessor
       // however, it affects how close the numbers are to the expectations
       // If the default setting for the accessor stub is changed for some reason,
       // the expected values may need adjustment (hence, the following assert to simplify debugging)
       CPPUNIT_ASSERT(cnt == 348000);
       // real and imaginary are two separate random numbers
       s /= 2.*double(cnt);
       s2 /= 2.*double(cnt);
       riCovar /= double(cnt);
       // real and imaginary parts don't correlate
       CPPUNIT_ASSERT(fabs(riCovar)<1e-3);
       // mean is zero
       CPPUNIT_ASSERT(fabs(s)<1e-2);
       // variance specified is sigma squared
       CPPUNIT_ASSERT(fabs((s2-s*s)*double(2*cnt)/double(2*cnt-1)-variance)<1e-2);
     }
  };

} // namespace synthesis

} // namespace askap

#endif // #ifndef GAUSSIAN_NOISE_ME_TEST_H

