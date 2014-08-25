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

#include <measurementequation/ComponentEquation.h>
#include <measurementequation/Calibrator1934.h>
#include <dataaccess/DataIteratorStub.h>

#include <cppunit/extensions/HelperMacros.h>

#include <askap/AskapError.h>

#include<algorithm>

namespace askap {

namespace synthesis {

class Calibrator1934Test : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(Calibrator1934Test);
   CPPUNIT_TEST(testFluxModel);
   CPPUNIT_TEST(testVisModel);
   CPPUNIT_TEST(testFullStokesModel);   
   CPPUNIT_TEST_SUITE_END();
public:
   void testFluxModel() {
      // a number of points in the 500 MHz to 2 GHz frequency range obtained with the miriad's calplot task
      CPPUNIT_ASSERT_DOUBLES_EQUAL(15.138, Calibrator1934::fluxDensity(1201.4),1e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(13.181, Calibrator1934::fluxDensity(769.7),1e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(10.631, Calibrator1934::fluxDensity(598.9),1e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(14.389, Calibrator1934::fluxDensity(1595.1),1e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(13.678, Calibrator1934::fluxDensity(1803.9),2e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(12.939, Calibrator1934::fluxDensity(2003.1),2e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(14.926, Calibrator1934::fluxDensity(1391.1),1e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(14.768, Calibrator1934::fluxDensity(990.3),1e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(14.292, Calibrator1934::fluxDensity(897.8),1e-3);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(15.055, Calibrator1934::fluxDensity(1094.6),1e-3);      
   };
   
   void testVisModel() {
      accessors::IDataSharedIter idi = accessors::IDataSharedIter(new accessors::DataIteratorStub(1));
      CPPUNIT_ASSERT(idi);
      CPPUNIT_ASSERT_EQUAL(8u, idi->nChannel());
      accessors::DataAccessorStub &acc = dynamic_cast<accessors::DataAccessorStub&>(*idi);
      // fill frequency axis with some points for which we know the flux of 1934-638 from miriad's calplot task
      const double freqs[] = {598.9e6,769.7e6,897.8e6,990.3e6,1094.6e6,1201.4e6,1391.1e6,1595.1e6};
      std::copy(freqs, freqs + 8, acc.itsFrequency.begin());
      askap::scimath::Params ip;
      ip.add("calibrator.1934-638");      
      ComponentEquation ce(ip,idi);
      ce.predict();
      const double fluxes[] = {10.631, 13.181, 14.292, 14.768, 15.055, 15.138, 14.926, 14.389};
      for (casa::uInt ch = 0; ch<acc.nChannel(); ++ch) {
           const casa::Matrix<casa::Complex> vis = acc.itsVisibility.xzPlane(ch);
           const double expectedFlux = fluxes[ch];
           CPPUNIT_ASSERT_EQUAL(acc.nRow(), vis.nrow());
           CPPUNIT_ASSERT_EQUAL(acc.nPol(), vis.ncolumn());
           
           for (casa::uInt row = 0; row<vis.nrow(); ++row) {
                for (casa::uInt col = 0; col<vis.ncolumn(); ++col) {
                     const casa::Complex simFlux = vis(row,col); 
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedFlux, static_cast<double>(real(simFlux)), 1e-3);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(0., static_cast<double>(imag(simFlux)), 1e-3);                     
                }
           }
      }      
   }

   void testFullStokesModel() {
      accessors::IDataSharedIter idi = accessors::IDataSharedIter(new accessors::DataIteratorStub(1));
      CPPUNIT_ASSERT(idi);
      accessors::DataAccessorStub &acc = dynamic_cast<accessors::DataAccessorStub&>(*idi);
      ASKAPASSERT(acc.itsStokes.nelements() == 1);
          
      casa::Vector<casa::Stokes::StokesTypes> stokes(4);
      stokes[0] = casa::Stokes::XX;
      stokes[1] = casa::Stokes::XY;
      stokes[2] = casa::Stokes::YX;
      stokes[3] = casa::Stokes::YY;         
          
      acc.itsStokes.assign(stokes.copy());
      acc.itsVisibility.resize(acc.nRow(), acc.nChannel() ,4);
      acc.itsVisibility.set(casa::Complex(-10.,15.));
      acc.itsNoise.resize(acc.nRow(),acc.nChannel(),acc.nPol());
      acc.itsNoise.set(1.);
      acc.itsFlag.resize(acc.nRow(),acc.nChannel(),acc.nPol());
      acc.itsFlag.set(casa::False);

      // fill frequency axis with some points for which we know the flux of 1934-638 from miriad's calplot task
      const double freqs[] = {598.9e6,769.7e6,897.8e6,990.3e6,1094.6e6,1201.4e6,1391.1e6,1595.1e6};
      CPPUNIT_ASSERT_EQUAL(8u, idi->nChannel());      
      std::copy(freqs, freqs + 8, acc.itsFrequency.begin());
      askap::scimath::Params ip;
      ip.add("calibrator.1934-638");      
      ComponentEquation ce(ip,idi);
      ce.predict();
      CPPUNIT_ASSERT_EQUAL(4u, acc.nPol());
      const double fluxes[] = {10.631, 13.181, 14.292, 14.768, 15.055, 15.138, 14.926, 14.389};
      for (casa::uInt ch = 0; ch<acc.nChannel(); ++ch) {
           const casa::Matrix<casa::Complex> vis = acc.itsVisibility.xzPlane(ch);
           CPPUNIT_ASSERT_EQUAL(acc.nRow(), vis.nrow());
           CPPUNIT_ASSERT_EQUAL(acc.nPol(), vis.ncolumn());
           for (casa::uInt pol = 0; pol<acc.nPol(); ++pol) {
                // we simulate XX and YY; x-pols should all be zeros 
                const double expectedFlux = (pol % 3 == 0) ? fluxes[ch] * 0.5 : 0.;
           
                for (casa::uInt row = 0; row<vis.nrow(); ++row) {
                     const casa::Complex simFlux = vis(row,pol); 
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedFlux, static_cast<double>(real(simFlux)), 1e-3);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(0., static_cast<double>(imag(simFlux)), 1e-3);                     
                }
           }
      }      
   }
};

} // namespace synthesis

} // namespace askap

