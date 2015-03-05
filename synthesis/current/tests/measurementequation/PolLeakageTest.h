/// @file
/// 
/// @brief Unit tests for polarisation leakage calibration.
/// @details The tests gathered in this file predict visibility data
/// with some calibration errors and then solve for them.
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

#ifndef POL_LEAKAGE_TEST_H
#define POL_LEAKAGE_TEST_H

#include <measurementequation/ComponentEquation.h>
#include <measurementequation/CalibrationME.h>
#include <measurementequation/PreAvgCalMEBase.h>
#include <measurementequation/LeakageTerm.h>
#include <measurementequation/NoXPolGain.h>
#include <measurementequation/Product.h>
#include <fitting/ComplexDiffMatrix.h>
#include <fitting/Params.h>

#include <fitting/LinearSolver.h>
#include <dataaccess/DataIteratorStub.h>
#include <measurementequation/CalibrationApplicatorME.h>
#include <calibaccess/CachedCalSolutionAccessor.h>
#include <calibaccess/CalSolutionSourceStub.h>
#include <measurementequation/CalibParamsMEAdapter.h>
#include <cppunit/extensions/HelperMacros.h>

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>

#include <boost/shared_ptr.hpp>


namespace askap
{
  namespace synthesis
  {

    class PolLeakageTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(PolLeakageTest);
      CPPUNIT_TEST(testBuildCDM);
      CPPUNIT_TEST(testBETAMuellerMatrix);      
      CPPUNIT_TEST(testSolve);
      CPPUNIT_TEST(testSolvePreAvg);
      CPPUNIT_TEST(testApplication);
      CPPUNIT_TEST(testSimulation);            
      CPPUNIT_TEST_SUITE_END();
     
     public:
      void setUp() {
          itsIter = boost::shared_ptr<accessors::DataIteratorStub>(new accessors::DataIteratorStub(1));
          accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*itsIter);
          ASKAPASSERT(da.itsStokes.nelements() == 1);
          
          casa::Vector<casa::Stokes::StokesTypes> stokes(4);
          stokes[0] = casa::Stokes::XX;
          stokes[1] = casa::Stokes::XY;
          stokes[2] = casa::Stokes::YX;
          stokes[3] = casa::Stokes::YY;         
          
          da.itsStokes.assign(stokes.copy());
          da.itsVisibility.resize(da.nRow(), 2 ,4);
          da.itsVisibility.set(casa::Complex(-10.,15.));
          da.itsNoise.resize(da.nRow(),da.nChannel(),da.nPol());
          da.itsNoise.set(1.);
          da.itsFlag.resize(da.nRow(),da.nChannel(),da.nPol());
          da.itsFlag.set(casa::False);
          da.itsFrequency.resize(da.nChannel());
          for (casa::uInt ch = 0; ch < da.nChannel(); ++ch) {
               da.itsFrequency[ch] = 1.4e9 + 20e6*double(ch);
          }
                    
          const casa::uInt nAnt = 30;
          
          // leakages, assume g12=g21
          const double realD[nAnt] = {0.1, -0.1, 0.05, -0.13, 0.333,
                                          0.1, 0.0, 0.0, -0.2, 0.03, 
                                         -0.05, 0.1, -0.1, -0.02, 0.03,
                                         -0.03, -0.1, 0.1, 0.1, 0.05,
                                          0.0, -0.03, 0.1, 0.03, 0.08,
                                          0.05, -0.07, 0.054, 0.0, 0.1}; 
          const double imagD[nAnt] = {0.0, 0., -0.05, 0.0587, 0.,
                                          0., -0.1, 0.02, -0.1, 0.84, 
                                          0.086, 0.1, 0.1, 0., 0.03,
                                         -0.084, 0., 0., -0.1, -0.05,
                                          0.02, 0.09, 0.1, 0.03, -0.1,
                                         -0.09, 0.072, -0.04, 0.05, -0.1}; 
      
          itsParams1.reset(new scimath::Params);
          itsParams1->add("flux.i.cena", 100.);
          itsParams1->add("direction.ra.cena", 0.5*casa::C::arcsec);
          itsParams1->add("direction.dec.cena", -0.3*casa::C::arcsec);
          itsParams1->add("shape.bmaj.cena", 3.0e-3*casa::C::arcsec);
          itsParams1->add("shape.bmin.cena", 2.0e-3*casa::C::arcsec);
          itsParams1->add("shape.bpa.cena", -55*casa::C::degree);
          for (casa::uInt ant=0; ant<nAnt; ++ant) {
               itsParams1->add("leakage.d12."+toString(ant)+".0",
                            casa::Complex(realD[ant],imagD[ant]));
               itsParams1->add("leakage.d21."+toString(ant)+".0",
                            casa::Complex(realD[ant],imagD[ant]));
          }
          //

          itsCE1.reset(new ComponentEquation(*itsParams1, itsIter));
          itsEq1.reset(new METype(*itsParams1,itsIter,itsCE1));
      
          
          itsParams2.reset(new scimath::Params);
          itsParams2->add("flux.i.cena", 100.);
          itsParams2->add("direction.ra.cena", 0.50000*casa::C::arcsec);
          itsParams2->add("direction.dec.cena", -0.30000*casa::C::arcsec);
          itsParams2->add("shape.bmaj.cena", 3.0e-3*casa::C::arcsec);
          itsParams2->add("shape.bmin.cena", 2.0e-3*casa::C::arcsec);
          itsParams2->add("shape.bpa.cena", -55*casa::C::degree);
          for (casa::uInt ant=0; ant<nAnt; ++ant) {
               itsParams2->add("leakage.d12."+toString(ant)+".0",casa::Complex(0.));
               itsParams2->add("leakage.d21."+toString(ant)+".0",casa::Complex(0.));
          }
       
          itsCE2.reset(new ComponentEquation(*itsParams2, itsIter));
      
      }
      /// @brief creates dummy parameters
      /// @details This method resets and fills itsParams1 with some dummy values
      /// to be used in individual unit tests
      void fillGainsAndLeakages() {
          const casa::uInt nAnt = 30;
          // use the following values to form both gains and leakages
          const double realGains[nAnt] = {1.1, 0.9, 1.05, 0.87, 1.333,
                                          1.1, 1.0, 1.0, -1.0, 0.3, 
                                         -0.5, 1.1, 0.9, 0.98, 1.03,
                                         -0.3, -1.1, 0.9, 1.1, 1.05,
                                          1.0, -0.3, 1.1, 0.3, 1.8,
                                          0.5, -0.7, 1.054, 1.0, 1.1}; 
          const double imagGains[nAnt] = {0.0, 0., -0.05, 0.587, 0.,
                                          0., -0.1, 0.02, -0.1, 0.84, 
                                          0.86, 0.1, 0.1, 0., 0.03,
                                         -0.84, 0., 0., -0.1, -0.05,
                                          0.2, 0.9, 1.1, 0.3, -0.1,
                                         -0.9, 0.72, -0.04, 0.05, -0.1}; 
          itsParams1.reset(new scimath::Params);
          itsParams1->add("flux.i.cena", 1.);
          itsParams1->add("direction.ra.cena", 0.*casa::C::arcsec);
          itsParams1->add("direction.dec.cena", 0.*casa::C::arcsec);
          for (casa::uInt ant=0; ant<nAnt; ++ant) {
               itsParams1->add(accessors::CalParamNameHelper::paramName(ant,0,casa::Stokes::XX),
                            casa::Complex(realGains[ant],imagGains[ant]));
               itsParams1->add(accessors::CalParamNameHelper::paramName(ant,0,casa::Stokes::YY),
                            casa::Complex(realGains[nAnt - 1 - ant],imagGains[nAnt - 1 - ant]));
                            
               itsParams1->add(accessors::CalParamNameHelper::paramName(ant,0,casa::Stokes::XY),
                            casa::Complex(realGains[ant] - 1.,imagGains[ant])/casa::Complex(10.,0.));
               itsParams1->add(accessors::CalParamNameHelper::paramName(ant,0,casa::Stokes::YX),
                            casa::Complex(realGains[ant] - 1.,imagGains[ant])/casa::Complex(10.,0.));                                           
                            
          }
      }
      
      void testBETAMuellerMatrix() {
          itsIter = boost::shared_ptr<accessors::DataIteratorStub>(new accessors::DataIteratorStub(1));
          CPPUNIT_ASSERT(itsIter);          
          accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*itsIter);
          CPPUNIT_ASSERT_EQUAL(size_t(1u), da.itsStokes.nelements());
          
          casa::Vector<casa::Stokes::StokesTypes> stokes(4);
          stokes[0] = casa::Stokes::XX;
          stokes[1] = casa::Stokes::XY;
          stokes[2] = casa::Stokes::YX;
          stokes[3] = casa::Stokes::YY;         
          
          da.itsStokes.assign(stokes.copy());
          const casa::uInt nAnt = 6;
          const casa::uInt nBaselines = nAnt * (nAnt - 1) / 2;
          da.itsVisibility.resize(nBaselines, 1 , stokes.nelements());
          da.itsVisibility.set(casa::Complex(-10.,15.));
          CPPUNIT_ASSERT_EQUAL(nBaselines, da.nRow());
          CPPUNIT_ASSERT_EQUAL(1u, da.nChannel());
          CPPUNIT_ASSERT_EQUAL(4u, da.nPol());
          da.itsNoise.resize(da.nRow(),da.nChannel(),da.nPol());
          da.itsNoise.set(1.);
          da.itsFlag.resize(da.nRow(),da.nChannel(),da.nPol());
          da.itsFlag.set(casa::False);
          da.itsFrequency.resize(da.nChannel());
          da.itsFrequency.set(1.4e9);
          da.itsAntenna1.resize(da.nRow());
          da.itsAntenna2.resize(da.nRow());
          da.itsFeed1.resize(da.nRow());
          da.itsFeed2.resize(da.nRow());
          da.itsFeed1.set(0u);
          da.itsFeed2.set(0u);
          for (casa::uInt ant1=0, row=0; ant1<nAnt; ++ant1) {
               for (casa::uInt ant2=0; ant2<ant1; ++ant2,++row) {
                    CPPUNIT_ASSERT(row < da.nRow());
                    da.itsAntenna1[row] = ant1;
                    da.itsAntenna2[row] = ant2;
               }
          }     
          // other fields are not used in the code we're testing, so can leave them uninitialised
          
          // results obtained from BETA SB 619-621, see ASKAPSDP-1633 (rough, for a representative channel)
          // order ak06,ak01,ak03,ak15,ak08,ak09 as in the real system
          const float d12[nAnt] = {0.025, 0.021, 0.08, 0.01, 0.045, 0.035};
          const float d21[nAnt] = {0.018, 0.015, 0.09, 0.01, 0.035, 0.02};
          const float xyPhase[nAnt] = {80., 140., -120., 140., 50., -175.};

          // fill parameters
          itsParams1.reset(new scimath::Params);
          ASKAPDEBUGASSERT(itsParams1);
          for (casa::uInt ant=0; ant<nAnt; ++ant) {
               itsParams1->add(accessors::CalParamNameHelper::paramName(ant,0,casa::Stokes::XX),
                            casa::Complex(1.,0.));
               itsParams1->add(accessors::CalParamNameHelper::paramName(ant,0,casa::Stokes::YY),
                            casa::Complex(1.,0.));
               
               const casa::Complex d12Term = casa::polar(d12[ant], float(xyPhase[ant]/180.*casa::C::pi));              
               itsParams1->add(accessors::CalParamNameHelper::paramName(ant,0,casa::Stokes::XY), d12Term);
               const casa::Complex d21Term = casa::polar(d21[ant], -float(xyPhase[ant]/180.*casa::C::pi));              
               itsParams1->add(accessors::CalParamNameHelper::paramName(ant,0,casa::Stokes::YX), d21Term);                            
          }
          
          typedef Product<NoXPolGain,LeakageTerm> EffectType; 
          const EffectType effect(itsParams1);
          for (casa::uInt testRow = 0; testRow < itsIter->nRow(); ++testRow) {                    
               const scimath::ComplexDiffMatrix cdm = effect.get(*itsIter,testRow);
               const casa::uInt testAnt1 = itsIter->antenna1()[testRow];
               const casa::uInt testAnt2 = itsIter->antenna2()[testRow];
               casa::Matrix<casa::Complex> mueller(cdm.nRow(),cdm.nColumn());
               ASKAPDEBUGASSERT(mueller.nrow() == mueller.ncolumn());
               casa::Vector<float> buf(mueller.nrow() * mueller.ncolumn() - mueller.nrow());
               casa::uInt count = 0;
               for (casa::uInt row=0; row<cdm.nRow(); ++row) {
                    for (casa::uInt col=0; col<cdm.nColumn(); ++col) {
                         mueller(row,col) = cdm(row,col).value();
                         if (row != col) {
                             CPPUNIT_ASSERT(count < buf.nelements());
                             buf(count++) = casa::abs(mueller(row,col));
                         }
                    }
               }
               CPPUNIT_ASSERT_EQUAL(buf.nelements(), size_t(count));
               const float expectedMax = ((testAnt1 == 2) || (testAnt2 == 2)) ? 0.1 : 0.05;
               CPPUNIT_ASSERT(casa::max(buf) < expectedMax); 
               //std::cout<<testAnt1<<" "<<testAnt2<<" "<<casa::min(buf)<<" "<<casa::max(buf)<<" "<<casa::median(buf)<<std::endl;
          }
      }
      
      void testBuildCDM() {
          fillGainsAndLeakages();
          CPPUNIT_ASSERT(itsParams1);
          typedef Product<NoXPolGain,LeakageTerm> EffectType; 
          
          const EffectType effect(itsParams1);
          CPPUNIT_ASSERT(itsIter);
          const casa::uInt nPol = 4;
          CPPUNIT_ASSERT(itsIter->stokes().nelements() == nPol);
          accessors::CachedCalSolutionAccessor acc(itsParams1);                    
          CPPUNIT_ASSERT(itsIter->nRow() > 0);
          CPPUNIT_ASSERT(itsIter->nPol() == nPol);
          for (casa::uInt testRow = 0; testRow < itsIter->nRow(); ++testRow) {
          
               const scimath::ComplexDiffMatrix cdm = effect.get(*itsIter,testRow);

               const casa::uInt testAnt1 = itsIter->antenna1()[testRow];
               const casa::uInt testAnt2 = itsIter->antenna2()[testRow];
               const casa::uInt testBeam1 = itsIter->feed1()[testRow];
               const casa::uInt testBeam2 = itsIter->feed2()[testRow];

               const casa::SquareMatrix<casa::Complex, 2> jones1 = acc.jones(testAnt1,testBeam1, 0);
               const casa::SquareMatrix<casa::Complex, 2> jones2 = acc.jones(testAnt2,testBeam2, 0);
               
               for (casa::uInt i = 0; i < nPol; ++i) {
                    for (casa::uInt j = 0; j < nPol; ++j) {
                         // element of the Mueller matrix obtained from first principles (outer product)
                         const casa::Complex expected = jones1(i / 2, j / 2) * conj(jones2(i % 2, j % 2));
                         CPPUNIT_ASSERT(i < cdm.nRow());
                         CPPUNIT_ASSERT(j < cdm.nColumn());                         
                         const casa::Complex obtained = cdm(i,j).value();
                         CPPUNIT_ASSERT_DOUBLES_EQUAL(real(expected),real(obtained),1e-6);
                         CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(expected),imag(obtained),1e-6);                         
                    }
               }
          }                    
      }
     
      void testSolve() {
          // Predict with the "perfect" parameters"
          CPPUNIT_ASSERT(itsEq1);
          itsEq1->predict();
          std::vector<std::string> freeNames = itsParams2->freeNames();
          for (std::vector<std::string>::const_iterator it = freeNames.begin();
               it!=freeNames.end();++it) {
               if (it->find("leakage") == std::string::npos) {
                   itsParams2->fix(*it);
               }
          }
          
          for (size_t iter=0; iter<10; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne; //(*params2);
            
               itsEq2.reset(new METype(*itsParams2,itsIter,itsCE2));
            
               itsEq2->calcEquations(ne);
               Quality q;
               LinearSolver solver1;
               solver1.addNormalEquations(ne);
               solver1.setAlgorithm("SVD");
               solver1.solveNormalEquations(*itsParams2,q);  
               //std::cout<<q<<std::endl;               
          }
          
          freeNames = itsParams2->freeNames();
          for (std::vector<std::string>::const_iterator it = freeNames.begin();
               it!=freeNames.end();++it) {
               CPPUNIT_ASSERT(itsParams2->has(*it));
               CPPUNIT_ASSERT(itsParams1->has(*it));       
               CPPUNIT_ASSERT_DOUBLES_EQUAL(casa::real(itsParams2->complexValue(*it) -
                                            itsParams1->complexValue(*it)), 0., 5e-3);
               CPPUNIT_ASSERT_DOUBLES_EQUAL(casa::imag(itsParams2->complexValue(*it) -
                                            itsParams1->complexValue(*it)), 0., 5e-3);
          }
      } 
      
      void testSolvePreAvg() {
          // Predict with the "perfect" parameters"
          CPPUNIT_ASSERT(itsEq1);
          itsEq1->predict();
          std::vector<std::string> freeNames = itsParams2->freeNames();
          for (std::vector<std::string>::const_iterator it = freeNames.begin();
               it!=freeNames.end();++it) {
               if (it->find("leakage") == std::string::npos) {
                   itsParams2->fix(*it);
               }
          }
    
          typedef CalibrationME<LeakageTerm,PreAvgCalMEBase> PreAvgMEType;
          boost::shared_ptr<PreAvgMEType> preAvgEq(new PreAvgMEType());
          CPPUNIT_ASSERT(preAvgEq);
          itsIter.init();          
          preAvgEq->accumulate(itsIter,itsCE2);

          for (size_t iter=0; iter<10; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne; 
            
               preAvgEq->setParameters(*itsParams2);            
               preAvgEq->calcEquations(ne);
               
               Quality q;
               LinearSolver solver1;
               solver1.addNormalEquations(ne);
               solver1.setAlgorithm("SVD");
               solver1.solveNormalEquations(*itsParams2,q);  
          }

          freeNames = itsParams2->freeNames();
          for (std::vector<std::string>::const_iterator it = freeNames.begin();
               it!=freeNames.end();++it) {
               CPPUNIT_ASSERT(itsParams2->has(*it));
               CPPUNIT_ASSERT(itsParams1->has(*it));                              
               CPPUNIT_ASSERT_DOUBLES_EQUAL(casa::real(itsParams2->complexValue(*it) -
                                            itsParams1->complexValue(*it)), 0., 5e-3);
               CPPUNIT_ASSERT_DOUBLES_EQUAL(casa::imag(itsParams2->complexValue(*it) -
                                            itsParams1->complexValue(*it)), 0., 5e-3);
          }                   
      }
      
      
      void testApplication() {        
          // check that everything is set up for full stokes
          CPPUNIT_ASSERT(itsIter);
          accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*itsIter);          
          CPPUNIT_ASSERT(da.itsStokes.nelements() == 4);
          da.rwVisibility().set(0.);          
          
          fillGainsAndLeakages();
          CPPUNIT_ASSERT(itsParams1);
          
          itsCE1.reset(new ComponentEquation(*itsParams1, itsIter));
          typedef CalibrationME<Product<NoXPolGain,LeakageTerm> > METype2;
          
          boost::shared_ptr<METype2> eq1(new METype2(*itsParams1,itsIter,itsCE1));
          eq1->predict();
          
          accessors::CachedCalSolutionAccessor acc(itsParams1);                    
          accessors::CalSolutionSourceStub src(boost::shared_ptr<accessors::CachedCalSolutionAccessor>(&acc,utility::NullDeleter()));
          CalibrationApplicatorME calME(boost::shared_ptr<accessors::CalSolutionSourceStub>(&src,utility::NullDeleter()));
          calME.correct(da);

          // check visibilities after calibration application
          const casa::Cube<casa::Complex>& vis = da.visibility();
          for (casa::uInt row = 0; row < da.nRow(); ++row) {
               for (casa::uInt chan = 0; chan < da.nChannel(); ++chan) {
                    for (casa::uInt pol = 0; pol < da.nPol(); ++pol) {
                         CPPUNIT_ASSERT_DOUBLES_EQUAL(pol % 3 == 0 ? 0.5 : 0., real(vis(row,chan,pol)),1e-6);
                         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., imag(vis(row,chan,pol)),1e-6);                         
                    }
               }
          }
        }
        
        void checkTwoParamsClasses(const scimath::Params &param1, const scimath::Params &param2) {
            const std::vector<string> names = param1.names();
            CPPUNIT_ASSERT_EQUAL(names.size(), param2.names().size());
            for (std::vector<string>::const_iterator ci = names.begin(); ci!=names.end(); ++ci) {
                 CPPUNIT_ASSERT(param1.has(*ci) && param2.has(*ci));
                 CPPUNIT_ASSERT_DOUBLES_EQUAL(real(param1.complexValue(*ci)),real(param2.complexValue(*ci)), 1e-6);
                 CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(param1.complexValue(*ci)),imag(param2.complexValue(*ci)), 1e-6);                 
            }
        }
        
        void testSimulation() {
          // this test is similar to testApplication, but simulation is done using parameters obtained
          // via calibration solution interface
          // check that everything is set up for full stokes
          CPPUNIT_ASSERT(itsIter);
          accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*itsIter);          
          CPPUNIT_ASSERT(da.itsStokes.nelements() == 4);
          da.rwVisibility().set(0.);          
          
          fillGainsAndLeakages();
          CPPUNIT_ASSERT(itsParams1);
          typedef CalibrationME<Product<NoXPolGain,LeakageTerm> > METype;
          itsCE2.reset(new ComponentEquation(*itsParams1, itsIter));
          boost::shared_ptr<METype> firstPrinciplesEqn(new METype(*itsParams1,itsIter,itsCE2));
          
          // itsParams1 has been copied inside firstPrinciplesEqn and itsCE2, can 
          // leave only gains/leakages in there
          itsParams1->remove("flux.i.cena");
          itsParams1->remove("direction.ra.cena");
          itsParams1->remove("direction.dec.cena");          
          scimath::Params tmpParams1(*itsParams1);          
          
          accessors::CachedCalSolutionAccessor acc(itsParams1);                    
          accessors::CalSolutionSourceStub css(boost::shared_ptr<accessors::CachedCalSolutionAccessor>(&acc,utility::NullDeleter()));
          
          itsParams2.reset(new scimath::Params);
          itsParams2->add("flux.i.cena", 1.);
          itsParams2->add("direction.ra.cena", 0.*casa::C::arcsec);
          itsParams2->add("direction.dec.cena", 0.*casa::C::arcsec);
          
          itsCE1.reset(new ComponentEquation(*itsParams2, itsIter));
          
          boost::shared_ptr<METype> eqn(new METype(*itsParams2,itsIter,itsCE1));
          boost::shared_ptr<accessors::CalSolutionSourceStub> cssPtr(&css, utility::NullDeleter());
          scimath::Params tmpParams2(*itsParams2);
          
          CalibParamsMEAdapter adapter(eqn,cssPtr,itsIter);                    
          // this should simulate 1 Jy point source in the phase centre with gains and leakages applied
          adapter.predict();
          //eqn->predict();
          
          // check that itsParams1 and 2 are intact
          checkTwoParamsClasses(tmpParams2, *itsParams2);
          checkTwoParamsClasses(tmpParams1, *itsParams1);
          const casa::Cube<casa::Complex> corruptedVis(da.visibility().copy());
          
          // now correct using the same solution source          
          CalibrationApplicatorME calME(cssPtr);
          calME.correct(da);

          
          // check visibilities after calibration application
          const casa::Cube<casa::Complex>& vis = da.visibility();
          for (casa::uInt row = 0; row < da.nRow(); ++row) {
               for (casa::uInt chan = 0; chan < da.nChannel(); ++chan) {
                    for (casa::uInt pol = 0; pol < da.nPol(); ++pol) {
                         CPPUNIT_ASSERT_DOUBLES_EQUAL(pol % 3 == 0 ? 0.5 : 0., real(vis(row,chan,pol)),1e-6);
                         CPPUNIT_ASSERT_DOUBLES_EQUAL(0., imag(vis(row,chan,pol)),1e-6);                         
                    }
               }
          }
          // simulate corrupted visibilities again from "first principles", i.e. using explicitly
          // defined gains and leakages in the parameters of the measurement equation
          CPPUNIT_ASSERT(firstPrinciplesEqn);
          firstPrinciplesEqn->predict();
          // check that the result is the same as with the ME adapter
          CPPUNIT_ASSERT_EQUAL(corruptedVis.shape(),vis.shape());
          for (casa::uInt row = 0; row < da.nRow(); ++row) {
               for (casa::uInt chan = 0; chan < da.nChannel(); ++chan) {
                    for (casa::uInt pol = 0; pol < da.nPol(); ++pol) {
                         CPPUNIT_ASSERT_DOUBLES_EQUAL(real(vis(row,chan,pol)),real(corruptedVis(row,chan,pol)),1e-6);
                         CPPUNIT_ASSERT_DOUBLES_EQUAL(imag(vis(row,chan,pol)),imag(corruptedVis(row,chan,pol)),1e-6);                         
                    }
               }
          }           
        }
      
     private:
      typedef CalibrationME<LeakageTerm> METype;
      boost::shared_ptr<ComponentEquation> itsCE1, itsCE2;
      boost::shared_ptr<METype> itsEq1,itsEq2;
      boost::shared_ptr<scimath::Params> itsParams1, itsParams2;
      accessors::SharedIter<accessors::DataIteratorStub> itsIter;      
    };
  } // namespace synthesis

} // namespace askap


#endif // #ifndef POL_LEAKAGE_TEST_H



