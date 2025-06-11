/// @file
///
/// @brief Unit tests for gain calibration interpolation.
/// @details The tests gathered in this file interpolate gains and
/// apply them to test data
///
///
/// @copyright (c) 2022 CSIRO
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
/// @author Mark Wieringa <mark.wieringa@csiro.au>

#ifndef GAIN_INTERPOLATE_TEST_H
#define GAIN_INTERPOLATE_TEST_H

#include <askap/measurementequation/ComponentEquation.h>
#include <askap/measurementequation/CalibrationME.h>
#include <askap/measurementequation/PreAvgCalMEBase.h>
#include <askap/measurementequation/LeakageTerm.h>
#include <askap/measurementequation/NoXPolGain.h>
#include <askap/measurementequation/Product.h>
#include <askap/scimath/fitting/ComplexDiffMatrix.h>
#include <askap/scimath/fitting/Params.h>

#include <askap/scimath/fitting/LinearSolver.h>
#include <askap/dataaccess/DataIteratorStub.h>
#include <askap/measurementequation/CalibrationApplicatorME.h>
#include <askap/calibaccess/TableCalSolutionSource.h>
#include <askap/calibaccess/CalSolutionSourceStub.h>
#include <askap/measurementequation/CalibParamsMEAdapter.h>
#include <cppunit/extensions/HelperMacros.h>

#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>

#include <boost/shared_ptr.hpp>


namespace askap
{
  namespace synthesis
  {
    using namespace accessors;

    class GainInterpolateTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(GainInterpolateTest);
      CPPUNIT_TEST(testInterpolatedApplication);
      CPPUNIT_TEST_SUITE_END();

      static boost::shared_ptr<ICalSolutionSource> rwSource(bool doRemove, uint nAnt, uint nChan) {
          const std::string fname("calibdata.tab");
          if (doRemove) {
              TableCalSolutionSource::removeOldTable(fname);
          }
          boost::shared_ptr<TableCalSolutionSource> css(new TableCalSolutionSource(fname, nAnt, 1, nChan));
          CPPUNIT_ASSERT(css);
          return css;
       }

     public:

     void fillGains(uint nAnt, uint nChan) {
         // set up gains for 2 solutions at t=0 and t=60
         boost::shared_ptr<ICalSolutionSource> css = rwSource(true, nAnt, nChan);
         long newID = css->newSolutionID(0.);
         CPPUNIT_ASSERT_EQUAL(0l, newID);
         boost::shared_ptr<ICalSolutionAccessor> acc = css->rwSolution(newID);
         const float pi4 = 3.14159265358/4;
         for (uint i = 0; i < nAnt; i++) {
             acc->setGain(JonesIndex(i,0u),JonesJTerm(1.6f*exp(casacore::Complex(0.,pi4*(-1.+2*(i%2)))),true,1.6f*exp(casacore::Complex(0.,pi4/2)),true));
         }
         // reuse existing table
         acc.reset();
         css.reset();
         css = rwSource(false, nAnt, nChan);
         newID = css->newSolutionID(60.);
         CPPUNIT_ASSERT_EQUAL(1l, newID);
         acc = css->rwSolution(newID);
         for (uint i = 0; i < nAnt; i++) {
             acc->setGain(JonesIndex(i,0u),JonesJTerm(0.4f*exp(casacore::Complex(0,pi4*(1-2*(i%2)))),true,0.4f*exp(casacore::Complex(0.,-pi4/2)),true));
         }
     }
     void setUp() {
         itsIter = boost::shared_ptr<accessors::DataIteratorStub>(new accessors::DataIteratorStub(1));
         accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*itsIter);
         ASKAPASSERT(da.itsStokes.nelements() == 1);

         casacore::Vector<casacore::Stokes::StokesTypes> stokes(4);
         stokes[0] = casacore::Stokes::XX;
         stokes[1] = casacore::Stokes::XY;
         stokes[2] = casacore::Stokes::YX;
         stokes[3] = casacore::Stokes::YY;

         da.itsStokes.assign(stokes.copy());
         da.itsVisibility.resize(4,8,da.nRow());
         da.itsVisibility.set(casacore::Complex(1.0,-1.0));
         da.itsNoise.resize(da.nPol(),da.nChannel(),da.nRow());
         da.itsNoise.set(1.);
         da.itsFlag.resize(da.nPol(),da.nChannel(),da.nRow());
         da.itsFlag.set(casacore::False);
         da.itsFrequency.resize(da.nChannel());
         for (casacore::uInt ch = 0; ch < da.nChannel(); ++ch) {
              da.itsFrequency[ch] = 1.4e9 + 20e6*double(ch);
         }
         da.itsTime = 30.0;
     }

     casacore::Cube<casacore::Complex> interpolatedApplication(bool interpolate = true, double time = 30.0) {
         setUp();
         // check that everything is set up for full stokes
         CPPUNIT_ASSERT(itsIter);
         accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*itsIter);
         CPPUNIT_ASSERT(da.itsStokes.nelements() == 4);
         da.itsTime = time;

         uint nAnt = 30;
         // Stub iterator is set up with 30 antennas
         CPPUNIT_ASSERT((nAnt*(nAnt-1))/2 == da.nRow());
         uint nChan = 8;
         fillGains(nAnt, nChan);

         boost::shared_ptr<ICalSolutionSource> css = rwSource(false, nAnt, nChan);
         CalibrationApplicatorME calME(css);
         calME.interpolateTime(interpolate);
         calME.correct(da);
         return da.visibility();
     }

     void testInterpolatedApplication() {
         const casacore::Cube<casacore::Complex>& vis1 = interpolatedApplication(false);
         const casacore::Cube<casacore::Complex>& vis2 = interpolatedApplication(true, 0.01);
         const casacore::Cube<casacore::Complex>& vis3 = interpolatedApplication(true, 30.0);
         const casacore::Cube<casacore::Complex>& vis4 = interpolatedApplication(true, 59.99);
         const casacore::Cube<casacore::Complex>& vis5 = interpolatedApplication(true, 60.0);
         // difference between no interpolation (using previous valid solution) and tiny amount of interpolation
         //std::cout<< "Vis1 -Vis2 <" <<max(abs(vis1 - vis2))<< std::endl;
         CPPUNIT_ASSERT(max(abs(vis1 - vis2)) < 0.001);
         CPPUNIT_ASSERT(max(abs(vis4 - vis5)) < 0.02);
         // difference between no interpolation (using next valid solution) and tiny amount of interpolation
         //std::cout<< "Vis4 -Vis5 <" <<max(abs(vis4 - vis5))<< std::endl;
         // values of midpoint interpolation for row 0 and 1
         //std::cout<< "Vmid(0) = "<< vis3(0,0,0) <<" "<< vis3(0,0,1)<<" "<<vis3(0,0,2)<<" "<<vis3(0,0,3)<<std::endl;
         //std::cout<< "Vmid(1) = "<< vis3(1,0,0) <<" "<< vis3(1,0,1)<<" "<<vis3(1,0,2)<<" "<<vis3(1,0,3)<<std::endl;
         CPPUNIT_ASSERT_DOUBLES_EQUAL(real(vis3(0,0,0)), 1.18151, 0.001);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(real(vis3(0,0,2)), 1.18151, 0.001);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(vis3(0,0,0)), -0.7772, 0.001);
         CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(vis3(0,0,2)), -0.7772, 0.001);
         for (uint pol = 1; pol < 3; pol+=2) {
             CPPUNIT_ASSERT_DOUBLES_EQUAL(real(vis3(0,0,pol)),1.0, 0.001);
             CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(vis3(0,0,pol)),-1.0, 0.001);
         }
         for (uint pol = 0; pol < 3; pol++) {
             CPPUNIT_ASSERT_DOUBLES_EQUAL(real(vis3(1,0,pol)),1.0, 0.001);
             CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(vis3(1,0,pol)),-1.0, 0.001);
         }
     }

     private:
       accessors::SharedIter<accessors::DataIteratorStub> itsIter;

   };
  } // namespace synthesis

} // namespace askap


#endif // #ifndef GAIN_INTERPOLATE_TEST_H
