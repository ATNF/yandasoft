/// @file
///
/// Unit test for BandpassDelayHelper class which is intended to solve for and modify delays in bandpass solutions
///
///
/// @copyright (c) 2025 CSIRO
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

#include <askap/utils/BandpassDelayHelper.h>
#include <askap/calibaccess/CachedCalSolutionAccessor.h>
#include <askap/scimath/fitting/CalParamNameHelper.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/askap/AskapUtil.h>

#include <cppunit/extensions/HelperMacros.h>

namespace askap {

namespace utils {

class BandpassDelayHelperTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(BandpassDelayHelperTest);
   CPPUNIT_TEST(testLoad);
   CPPUNIT_TEST(testStore);
   CPPUNIT_TEST(testIdeal);
   CPPUNIT_TEST(testDelayCalcAndApply);
   CPPUNIT_TEST(testDelayCalc);
   CPPUNIT_TEST(testAddDelay);
   CPPUNIT_TEST_SUITE_END();
public:

   void testLoad() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 9u;
      const casacore::uInt nChan = 16u;
      BandpassDelayHelper bdh(nAnt, nBeam, nChan);

      // instead of table-based solution accessor we use CachedCalSolutionAccessor adapter which stores info in Params class
      // first fill it with data and then try to load it into the helper class.
      scimath::Params buffer;
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                // we'll have only one polarisation defined for testing depending on whether ant == beam or otherwise
                const std::string parName = scimath::CalParamNameHelper::paramName(ant, beam, ant == beam ? casacore::Stokes::YY : casacore::Stokes::XX);
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     // ensure we have different parameters for each antenna / beam / channel combination
                     buffer.add(scimath::CalParamNameHelper::addChannelInfo(parName, chan), casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1));
                }
           }
      }

      // setup params-based solution accessor
      boost::shared_ptr<scimath::Params> bufferPtr(&buffer, utility::NullDeleter());
      accessors::CachedCalSolutionAccessor acc(bufferPtr);

      // load bandpass into the helper class
      bdh.loadBandpass(acc, 1e6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1e6, bdh.itsResolution, 1.);

      // now check the content of the helper class to verify that bandpass has been loaded correctly along with the validity flags
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     const casacore::IPosition validPoint(4, chan, beam == ant ? 1 : 0, beam, ant);
                     const casacore::IPosition invalidPoint(4, chan, beam == ant ? 0 : 1, beam, ant);
                     CPPUNIT_ASSERT(bdh.itsBandpassValid(validPoint));
                     CPPUNIT_ASSERT(!bdh.itsBandpassValid(invalidPoint));
                     const casacore::Complex expectedVal = casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::real(expectedVal), casacore::real(bdh.itsBandpass(validPoint)), 1e-6);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::imag(expectedVal), casacore::imag(bdh.itsBandpass(validPoint)), 1e-6);
                }
           }
      }
   }

   void testStore() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 9u;
      const casacore::uInt nChan = 16u;
      BandpassDelayHelper bdh(nAnt, nBeam, nChan);

      // fill the values in the buffer directly
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     const casacore::IPosition validPoint(4, chan, beam == ant ? 1 : 0, beam, ant);
                     const casacore::IPosition invalidPoint(4, chan, beam == ant ? 0 : 1, beam, ant);
                     bdh.itsBandpassValid(validPoint) = true;
                     // check that invalid point is invalid by default instead of assigning
                     CPPUNIT_ASSERT(!bdh.itsBandpassValid(invalidPoint));
                     bdh.itsBandpass(validPoint) = casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1);
                }
           }
      }
      // params-based accessor to store the bandpass, Params object is held internally
      accessors::CachedCalSolutionAccessor acc;
      bdh.storeBandpass(acc);

      // check the results
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                // we'll have only one polarisation defined for testing depending on whether ant == beam or otherwise
                const std::string parName = scimath::CalParamNameHelper::paramName(ant, beam, ant == beam ? casacore::Stokes::YY : casacore::Stokes::XX);
                const std::string otherPolParName = scimath::CalParamNameHelper::paramName(ant, beam, ant == beam ? casacore::Stokes::XX : casacore::Stokes::YY);
                // the actual parameter names have channel number encoded, so these ones should be missing
                CPPUNIT_ASSERT(!acc.cache().has(otherPolParName));
                CPPUNIT_ASSERT(!acc.cache().has(parName));
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     const std::string fullParName = scimath::CalParamNameHelper::addChannelInfo(parName, chan);
                     CPPUNIT_ASSERT(acc.cache().has(fullParName));
                     CPPUNIT_ASSERT(!acc.cache().has(scimath::CalParamNameHelper::addChannelInfo(otherPolParName, chan)));
                     const casacore::Complex expected = casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1);
                     const casacore::Complex obtained = acc.cache().complexValue(scimath::CalParamNameHelper::addChannelInfo(parName, chan));
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::real(expected), casacore::real(obtained), 1e-6);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::imag(expected), casacore::imag(obtained), 1e-6);
                }
           }
      }
   }

   void testIdeal() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 3u;
      const casacore::uInt nChan = 10u;
      BandpassDelayHelper bdh(nAnt, nBeam, nChan);

      bdh.setIdealBandpass(1e6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1e6, bdh.itsResolution, 1.);

      // check the results
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                for (casacore::uInt pol = 0; pol < 2u; ++pol) {
                     for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                          const casacore::IPosition index(4, chan, pol, beam, ant);
                          CPPUNIT_ASSERT(bdh.itsBandpassValid(index));
                          CPPUNIT_ASSERT_DOUBLES_EQUAL(1., casacore::real(bdh.itsBandpass(index)), 1e-6);
                          CPPUNIT_ASSERT_DOUBLES_EQUAL(0., casacore::imag(bdh.itsBandpass(index)), 1e-6);
                     }
                }
           }
      }
   }

   void testDelayCalcAndApply() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 9u;
      const casacore::uInt nChan = 64u;
      BandpassDelayHelper bdh(nAnt, nBeam, nChan);
      CPPUNIT_ASSERT(bdh.itsDelay.shape() == casacore::IPosition(3, 2, nBeam, nAnt));

      bdh.setIdealBandpass(1e6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1e6, bdh.itsResolution, 1.);

      // fill in delays only with some validity flags
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                const casacore::IPosition validPoint(3, beam == ant ? 1 : 0, beam, ant);
                const casacore::IPosition invalidPoint(3, beam == ant ? 0 : 1, beam, ant);
                // check that invalid point is invalid by default instead of assigning it
                CPPUNIT_ASSERT(!bdh.itsDelayValid(invalidPoint));
                bdh.itsDelayValid(validPoint) = true;
                // some unique number of nanoseconds for each antenna and beam combination
                bdh.itsDelay(validPoint) = (ant * nBeam + beam) * 1e-9;
           }
      } 

      // now apply the delays set above to the ideal bandpass
      bdh.applyDelays();

      // reset delays and validity flags
      bdh.itsDelay.set(0.f);
      bdh.itsDelayValid.set(false);

      // compute delays and we should get figures close to what we put in originally
      bdh.calcDelays();

      // check the results
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                const casacore::IPosition validPoint(3, beam == ant ? 1 : 0, beam, ant);
                const casacore::IPosition invalidPoint(3, beam == ant ? 0 : 1, beam, ant);
                CPPUNIT_ASSERT(!bdh.itsDelayValid(invalidPoint));
                CPPUNIT_ASSERT(bdh.itsDelayValid(validPoint));
                const float expected = (ant * nBeam + beam) * 1e-9;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, bdh.itsDelay(validPoint), 1e-6);
           }
      } 
      
   }

   void testDelayCalc() {
      // this test ensures correctness of the actual delay calculation
      // only have 1 antenna in the bandpass table for simplicity, although this is impossible in practice
      // Moreover, leave the second polarisation undefined, so it can be tested too that this flag propagates to delays
      const casacore::uInt nChan = 100u;
      BandpassDelayHelper bdh(1u, 1u, nChan);
      // 1 Hz spectral resolution and 100 channels, one wrap across the band is equivalent to 10ms delay
      // (note, in practice delays will be in nanoseconds - microseconds range, but for a test it is ok)
      bdh.setIdealBandpass(1.);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1., bdh.itsResolution, 1e-6);
      // generate phase slope from first principles 
      casacore::Vector<casacore::Complex> bp = casacore::Matrix<casacore::Complex>(bdh.itsBandpass.reform(casacore::IPosition(2,nChan, 2))).column(0);
      casacore::Matrix<bool> bpValidMatr = bdh.itsBandpassValid.reform(casacore::IPosition(2,nChan, 2));
      for (casacore::uInt chan = 0; chan < nChan; ++chan) {
           // after setIdealBandpass all points are valid
           CPPUNIT_ASSERT(bpValidMatr(chan,0));
           CPPUNIT_ASSERT(bpValidMatr(chan,1));
           // only make the second polarisation invalid
           bpValidMatr(chan,1) = false;
           // have the amplitude ramping up too, it should be irrelevant for the delay
           bp[chan] = casacore::polar(static_cast<float>(1+chan), static_cast<float>(casacore::C::_2pi * 0.01 * chan));
           // few flagged channels across the band shouldn't cause a problem (but we need a moderate tolerance when we compare the result as the
           // match wouldn't be perfect)
           if (chan % 10u == 1u) {
               bpValidMatr(chan,0) = false;
           }
      }
      // unflag one channel for the otherwise completely flagged polarisation to check that we require a minimum of 2 channels for a valid delay
      // solution (although of course it will still be dodgy, but we leave it up to the user; stricter limits could be implemented if we want)
      bpValidMatr(nChan / 2,1) = true;

      // compute delays
      bdh.calcDelays();
      
      // check the result
      CPPUNIT_ASSERT(bdh.itsDelayValid.shape() == casacore::IPosition(3, 2, 1, 1));
      CPPUNIT_ASSERT(bdh.itsDelay.shape() == casacore::IPosition(3, 2, 1, 1));
      CPPUNIT_ASSERT(bdh.itsDelayValid(0, 0, 0));
      CPPUNIT_ASSERT(!bdh.itsDelayValid(1, 0, 0));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.01f, bdh.itsDelay(0,0,0), 1e-5);

      // check operations of negateDelays as well
      bdh.negateDelays();

      // check the result
      CPPUNIT_ASSERT(bdh.itsDelayValid.shape() == casacore::IPosition(3, 2, 1, 1));
      CPPUNIT_ASSERT(bdh.itsDelay.shape() == casacore::IPosition(3, 2, 1, 1));
      CPPUNIT_ASSERT(bdh.itsDelayValid(0, 0, 0));
      CPPUNIT_ASSERT(!bdh.itsDelayValid(1, 0, 0));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.01f, bdh.itsDelay(0,0,0), 1e-5);

      // check delay zeroing preserving the flags (while we have a suitably setup helper class)
      bdh.zeroDelays();
      CPPUNIT_ASSERT(bdh.itsDelayValid.shape() == casacore::IPosition(3, 2, 1, 1));
      CPPUNIT_ASSERT(bdh.itsDelay.shape() == casacore::IPosition(3, 2, 1, 1));
      CPPUNIT_ASSERT(bdh.itsDelayValid(0, 0, 0));
      CPPUNIT_ASSERT(!bdh.itsDelayValid(1, 0, 0));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, bdh.itsDelay(0,0,0), 1e-5);
   }

   void testAddDelay() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 9u;
      const casacore::uInt nChan = 64u;
      BandpassDelayHelper bdh(nAnt+1u, nBeam+1u, nChan*2u);
      CPPUNIT_ASSERT(bdh.itsDelay.shape() == casacore::IPosition(3, 2, nBeam+1u, nAnt+1u));
      BandpassDelayHelper bdh2(nAnt, nBeam, nChan);
      CPPUNIT_ASSERT(bdh2.itsDelay.shape() == casacore::IPosition(3, 2, nBeam, nAnt));

      bdh.setIdealBandpass(5e5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(5e5, bdh.itsResolution, 1.);
      bdh2.setIdealBandpass(1e6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1e6, bdh2.itsResolution, 1.);

      // fill in delays only with some validity flags
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                const casacore::IPosition validPoint(3, beam == ant ? 1 : 0, beam, ant);
                const casacore::IPosition invalidPoint(3, beam == ant ? 0 : 1, beam, ant);
                // check that invalid point is invalid by default instead of assigning it
                CPPUNIT_ASSERT(!bdh2.itsDelayValid(invalidPoint));
                CPPUNIT_ASSERT(!bdh.itsDelayValid(invalidPoint));
                // flag one beam completely only in one of the instances
                bdh2.itsDelayValid(validPoint) = (beam != nBeam / 2);
                const float value = static_cast<float>(ant + nAnt * beam) * (beam % 2 ? 1.f: -1.f) * 1e-9;
                bdh2.itsDelay(validPoint) = value;
                bdh.itsDelayValid(validPoint) = true;
                bdh.itsDelay(validPoint) = value;
           }
           // the first helper class has one more antenna and beam - define some of these delays
           bdh.itsDelay(0, nBeam, ant) = 1e-9;
           bdh.itsDelayValid(1, nBeam, ant) = true;
      }
      bdh.itsDelayValid(1, nBeam, nAnt) = true;
      bdh.itsDelay(0, nBeam, nAnt) = -1e-9;

      // apply delays, although strictly speaking this is irrelevant for this particular test which only verifies delay manipulation
      bdh.applyDelays();
      bdh2.applyDelays();

      bdh2.negateDelays();
      // the following should make all valid delays in bdh2 zero, extra beam and antenna shouldn't matter
      bdh2.addDelays(bdh);
      // the following should invalidate delays related to extra beam and antenna, other values should remain the same
      bdh.addDelays(bdh2);

      // now check
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                const casacore::IPosition validPoint(3, beam == ant ? 1 : 0, beam, ant);
                const casacore::IPosition invalidPoint(3, beam == ant ? 0 : 1, beam, ant);
                CPPUNIT_ASSERT(!bdh2.itsDelayValid(invalidPoint));
                CPPUNIT_ASSERT(!bdh.itsDelayValid(invalidPoint));
                // account for one flagged beam
                if (beam != nBeam / 2) {
                    CPPUNIT_ASSERT(bdh2.itsDelayValid(validPoint));
                    CPPUNIT_ASSERT(bdh.itsDelayValid(validPoint));
                    const float expected = static_cast<float>(ant + nAnt * beam) * (beam % 2 ? 1.f: -1.f) * 1e-9;
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.f, bdh2.itsDelay(validPoint), 1e-14);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, bdh.itsDelay(validPoint), 1e-14);
                } else {
                    CPPUNIT_ASSERT(!bdh2.itsDelayValid(validPoint));
                    CPPUNIT_ASSERT(!bdh.itsDelayValid(validPoint));
                }
           }
           CPPUNIT_ASSERT(!bdh.itsDelayValid(0, nBeam, ant));
           CPPUNIT_ASSERT(!bdh.itsDelayValid(1, nBeam, ant));
      }
      CPPUNIT_ASSERT(!bdh.itsDelayValid(0, nBeam, nAnt));
      CPPUNIT_ASSERT(!bdh.itsDelayValid(1, nBeam, nAnt));
   }
};
    
} // namespace utils

} // namespace askap

