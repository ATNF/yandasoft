/// @file
/// 
/// @brief Unit tests for PreAvgCalBuffer.
/// @details PreAvgCalBuffer accululates partial sums for a number of
/// visibility groups (indexed by baseline and beam), which are then used
/// in the least square problem avoiding the iteration over the original dataset.
/// This file contains unit tests of this class
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

#ifndef PRE_AVG_CAL_BUFFER_TEST_H
#define PRE_AVG_CAL_BUFFER_TEST_H

#include <dataaccess/DataIteratorStub.h>
#include <cppunit/extensions/HelperMacros.h>
#include <measurementequation/PreAvgCalBuffer.h>
#include <measurementequation/ComponentEquation.h>
#include <fitting/PolXProducts.h>

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>

#include <boost/shared_ptr.hpp>


namespace askap {

namespace synthesis   {

/// @brief unit tests of PreAvgCalBuffer
class PreAvgCalBufferTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(PreAvgCalBufferTest);
  CPPUNIT_TEST(testInitByAccessor);
  CPPUNIT_TEST(testInitExplicit);
  CPPUNIT_TEST(testBeamIndependent);
  CPPUNIT_TEST(testPolIndex);
  CPPUNIT_TEST(testAccumulate);
  CPPUNIT_TEST(testFDPAccumulate);
  CPPUNIT_TEST(testFDPInitExplicit);
  CPPUNIT_TEST(testAccumulateXPol);
  CPPUNIT_TEST_SUITE_END();
      
  private:
     boost::shared_ptr<ComponentEquation> itsME;
     boost::shared_ptr<Params> itsParams;
     accessors::SharedIter<accessors::DataIteratorStub> itsIter;
     
  public:
     void setUp() {
         itsParams.reset(new Params);
         itsParams->add("flux.i.src", 100.);
         itsParams->add("direction.ra.src", 0.5*casa::C::arcsec);
         itsParams->add("direction.dec.src", -0.3*casa::C::arcsec);
         itsParams->add("shape.bmaj.src", 3.0e-3*casa::C::arcsec);
         itsParams->add("shape.bmin.src", 2.0e-3*casa::C::arcsec);
         itsParams->add("shape.bpa.src", -55*casa::C::degree);
         
         itsIter = accessors::SharedIter<accessors::DataIteratorStub>(new accessors::DataIteratorStub(1));
         accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*itsIter);
         ASKAPASSERT(da.itsStokes.nelements() == 1);
         da.itsStokes[0] = casa::Stokes::I;   
         da.itsNoise.set(1.0);
         
         itsME.reset(new ComponentEquation(*itsParams, itsIter));               
     }
     void testInitByAccessor() {
         PreAvgCalBuffer pacBuf;
         pacBuf.initialise(*itsIter);
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         
         CPPUNIT_ASSERT_EQUAL(itsIter->nRow(),pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(1u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(itsIter->nPol(),pacBuf.nPol());
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nRow()),size_t(pacBuf.flag().nrow()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nPol()),size_t(pacBuf.flag().nplane()));
         CPPUNIT_ASSERT_EQUAL(1u,casa::uInt(pacBuf.flag().ncolumn()));
         CPPUNIT_ASSERT_EQUAL(1u,casa::uInt(pacBuf.stokes().nelements()));
         CPPUNIT_ASSERT_EQUAL(casa::Stokes::I, pacBuf.stokes()[0]);
                  
         for (casa::uInt row=0; row < pacBuf.nRow(); ++row)  {
              CPPUNIT_ASSERT_EQUAL(itsIter->antenna1()[row],pacBuf.antenna1()[row]);
              CPPUNIT_ASSERT_EQUAL(itsIter->antenna2()[row],pacBuf.antenna2()[row]);
              CPPUNIT_ASSERT_EQUAL(itsIter->feed1()[row],pacBuf.feed1()[row]);
              CPPUNIT_ASSERT_EQUAL(itsIter->feed2()[row],pacBuf.feed2()[row]);
              for (casa::uInt pol=0; pol<pacBuf.nPol(); ++pol) {
                   CPPUNIT_ASSERT_EQUAL(true,pacBuf.flag()(row,0,pol));
              }              
         }
     }
     
     /// @brief helper method to change the beam index for all accessor stub
     /// @param[in] beam new beam index
     void setBeamIndex(const casa::uInt beam) {
         boost::shared_ptr<accessors::IDataAccessor> origDA(itsIter.operator->(),utility::NullDeleter());
         CPPUNIT_ASSERT(origDA);
         const boost::shared_ptr<accessors::DataAccessorStub> da = 
               boost::dynamic_pointer_cast<accessors::DataAccessorStub>(origDA);
         CPPUNIT_ASSERT(da);
         da->itsFeed1.set(beam);
         da->itsFeed2.set(beam);
     }
     
     void testInitExplicit() {
         // 20 antennas instead of 30 available, 2 beams instead of 1 available in the stubbed
         // accessor
         PreAvgCalBuffer pacBuf(20,2);
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         
         // 20 antennas and 2 beams give 380 rows; 4 polarisation by default
         CPPUNIT_ASSERT_EQUAL(380u,pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(1u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(4u,pacBuf.nPol());
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nRow()),size_t(pacBuf.flag().nrow()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nPol()),size_t(pacBuf.flag().nplane()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nChannel()),size_t(pacBuf.flag().ncolumn()));
         CPPUNIT_ASSERT_EQUAL(4u,casa::uInt(pacBuf.stokes().nelements()));
         CPPUNIT_ASSERT(scimath::PolConverter::isLinear(pacBuf.stokes()));
         for (casa::uInt pol = 0; pol<4; ++pol) {
              CPPUNIT_ASSERT_EQUAL(pol,scimath::PolConverter::getIndex(pacBuf.stokes()[pol]));
         }
         
         CPPUNIT_ASSERT(itsME);
         CPPUNIT_ASSERT(itsIter);
         
         // simulate visibilities
         itsME->predict(*itsIter);
         
         pacBuf.accumulate(*itsIter, itsME);
         
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         // (435 - 190) * 8 = 1960 samples unaccounted for 
         // (accessor has 1 polarisation)
         CPPUNIT_ASSERT_EQUAL(1960u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         

         // change the beam index and accumulate again
         setBeamIndex(2);
         pacBuf.accumulate(*itsIter, itsME);

         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         // (2*435 - 190) * 8 = 5440 samples unaccounted for 
         // (the second beam is not mapped)
         CPPUNIT_ASSERT_EQUAL(5440u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         

         const PolXProducts& pxp = pacBuf.polXProducts();
         for (casa::uInt row=0; row<pacBuf.nRow(); ++row) {
              CPPUNIT_ASSERT_EQUAL(pacBuf.feed1()[row], pacBuf.feed2()[row]);
              for (casa::uInt pol=0; pol<pacBuf.nPol(); ++pol) {
                   if ((pol == 0) && (pacBuf.feed1()[row] == 0)) {
                       CPPUNIT_ASSERT_EQUAL(false, pacBuf.flag()(row,0,pol));                       
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(double(real(pxp.getModelProduct(row,0,pol,pol))),
                                              double(real(pxp.getModelMeasProduct(row,0,pol,pol))),1e-2);
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelMeasProduct(row,0,pol,pol))),1e-5);
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelProduct(row,0,pol,pol))),1e-5);
                       // 8 channels and 100 Jy source give sums of 80000 per accessor summed in
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(80000., double(real(pxp.getModelProduct(row,0,pol,pol))),1e-2);
                   } else {
                       // nothing should be found in the accessor, so the appropriate samples should be flagged
                       CPPUNIT_ASSERT_EQUAL(true, pacBuf.flag()(row,0,pol));
                   }
              }
         }
         
     }

     void testBeamIndependent() {
         // 20 antennas instead of 30 available, 2 beams instead of 1 available in the stubbed
         // accessor
         PreAvgCalBuffer pacBuf(20);
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         
         // 20 antennas and 1 aggregate beam give 190 rows; 4 polarisation by default
         CPPUNIT_ASSERT_EQUAL(190u,pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(1u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(4u,pacBuf.nPol());
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nRow()),size_t(pacBuf.flag().nrow()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nPol()),size_t(pacBuf.flag().nplane()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nChannel()),size_t(pacBuf.flag().ncolumn()));
         CPPUNIT_ASSERT_EQUAL(4u,casa::uInt(pacBuf.stokes().nelements()));
         CPPUNIT_ASSERT(scimath::PolConverter::isLinear(pacBuf.stokes()));
         for (casa::uInt pol = 0; pol<4; ++pol) {
              CPPUNIT_ASSERT_EQUAL(pol,scimath::PolConverter::getIndex(pacBuf.stokes()[pol]));
         }
         
         CPPUNIT_ASSERT(itsME);
         CPPUNIT_ASSERT(itsIter);
         
         // simulate visibilities
         itsME->predict(*itsIter);
                  
         pacBuf.accumulate(*itsIter, itsME);
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         // (435 - 190) * 8 = 1960 samples unaccounted for 
         // (accessor has 1 polarisation)
         CPPUNIT_ASSERT_EQUAL(1960u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         


         // change the beam index and accumulate again
         setBeamIndex(2);
         pacBuf.accumulate(*itsIter, itsME);
         
         
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         // (435 - 190) * 8 * 2 = 3920 samples are now unaccounted for 
         CPPUNIT_ASSERT_EQUAL(3920u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         

         const PolXProducts& pxp = pacBuf.polXProducts();
         for (casa::uInt row=0; row<pacBuf.nRow(); ++row) {
              CPPUNIT_ASSERT_EQUAL(pacBuf.feed1()[row], pacBuf.feed2()[row]);
              for (casa::uInt pol=0; pol<pacBuf.nPol(); ++pol) {
                   if ((pol == 0) && (pacBuf.feed1()[row] == 0)) {
                       CPPUNIT_ASSERT_EQUAL(false, pacBuf.flag()(row,0,pol));                       
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(double(real(pxp.getModelProduct(row,0,pol,pol))),
                                              double(real(pxp.getModelMeasProduct(row,0,pol,pol))),1e-2);
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelMeasProduct(row,0,pol,pol))),1e-5);
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelProduct(row,0,pol,pol))),1e-5);
                       // 8 channels, 2 beams and 100 Jy source give sums of 160000 per accessor summed in
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(160000., double(real(pxp.getModelProduct(row,0,pol,pol))),1e-2);
                   } else {
                       // nothing should be found in the accessor, so the appropriate samples should be flagged
                       CPPUNIT_ASSERT_EQUAL(true, pacBuf.flag()(row,0,pol));
                   }
              }
         }
         
     }
     
     void testPolIndex() {
         // 20 antennas, 1 beam + 4 polarisations by default
         PreAvgCalBuffer pacBuf(20,1);
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         
         // 20 antennas and 1 beam give 190 rows; 4 polarisation by default
         CPPUNIT_ASSERT_EQUAL(190u,pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(1u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(4u,pacBuf.nPol());
         const scimath::PolXProducts &pxp = pacBuf.polXProducts();
         for (casa::uInt pol1 = 0; pol1 < pacBuf.nPol(); ++pol1) {
              for (casa::uInt pol2 = 0; pol2 <= pol1; ++pol2) {
                   // now check polarisation indexing inside pacBuf (now handled by PolXProducts class)
                   pxp.getModelProduct(0,0,pol1,pol2);
              }
         }
     }
     
     void testResults(const PreAvgCalBuffer &pacBuf, const int run = 1) {
         for (casa::uInt row=0; row<pacBuf.nRow(); ++row) {
              CPPUNIT_ASSERT_EQUAL(pacBuf.feed1()[row], pacBuf.feed2()[row]);
              CPPUNIT_ASSERT(pacBuf.nPol() > 0);
              
              const scimath::PolXProducts &pxp = pacBuf.polXProducts();              
              for (casa::uInt pol=0; pol<pacBuf.nPol(); ++pol) {
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(double(real(pxp.getModelProduct(row,0,pol,pol))),
                                                double(real(pxp.getModelMeasProduct(row,0,pol,pol))),1e-2*run);
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelMeasProduct(row,0,pol,pol))),1e-5);
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelProduct(row,0,pol,pol))),1e-5);
                   // 8 channels and 100 Jy source give sums of 80000 per accessor summed in
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(pol%3 == 0 ? 80000.*run : 0., double(real(pxp.getModelProduct(row,0,pol,pol))),1e-2*run);
                   CPPUNIT_ASSERT_EQUAL(false, pacBuf.flag()(row,0,pol));                                     
    
                  // checking cross-pol terms, if any
                  for (casa::uInt pol2 = 0; pol2<pol; ++pol2) {
                       const casa::Complex expected = ((pol == 3) && (pol2 == 0)) ? casa::Complex(80000.*run,0.) : casa::Complex(0.,0.);
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,casa::abs(expected - pxp.getModelProduct(row,0,pol,pol2)),1e-2*run);
                       CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,casa::abs(expected - pxp.getModelMeasProduct(row,0,pol,pol2)),1e-2*run);
                  }
              }
         }
     }

     void testFrequencyDependentResults(const PreAvgCalBuffer &pacBuf, const int run = 1) {
         for (casa::uInt row=0; row<pacBuf.nRow(); ++row) {
              CPPUNIT_ASSERT_EQUAL(pacBuf.feed1()[row], pacBuf.feed2()[row]);
              CPPUNIT_ASSERT(pacBuf.nPol() > 0);
              
              const scimath::PolXProducts &pxp = pacBuf.polXProducts();              
              for (casa::uInt pol=0; pol<pacBuf.nPol(); ++pol) {
                   for (casa::uInt chan=0; chan<pacBuf.nChannel(); ++chan) {
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(double(real(pxp.getModelProduct(row,chan,pol,pol))),
                                                     double(real(pxp.getModelMeasProduct(row,chan,pol,pol))),1e-2*run);
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelMeasProduct(row,chan,pol,pol))),1e-5);
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelProduct(row,chan,pol,pol))),1e-5);
                        // 1 channel and 100 Jy source give sums of 10000 per accessor summed in
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(pol%3 == 0 ? 10000.*run : 0., double(real(pxp.getModelProduct(row,chan,pol,pol))),1e-2*run);
                        CPPUNIT_ASSERT_EQUAL(false, pacBuf.flag()(row,chan,pol));                                     
    
                        // checking cross-pol terms, if any
                        for (casa::uInt pol2 = 0; pol2<pol; ++pol2) {
                             const casa::Complex expected = ((pol == 3) && (pol2 == 0)) ? casa::Complex(10000.*run,0.) : casa::Complex(0.,0.);
                             CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,casa::abs(expected - pxp.getModelProduct(row,chan,pol,pol2)),1e-2*run);
                             CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,casa::abs(expected - pxp.getModelMeasProduct(row,chan,pol,pol2)),1e-2*run);
                        }
                        
                   }
              }
         }
     }
     
     void testAccumulate() {
         PreAvgCalBuffer pacBuf;
         CPPUNIT_ASSERT(itsME);
         CPPUNIT_ASSERT(itsIter);
         
         // simulate visibilities
         itsME->predict(*itsIter);
         
         // buffer should be initialised by the first encountered accessor
         pacBuf.accumulate(*itsIter, itsME);

         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());
         CPPUNIT_ASSERT_EQUAL(itsIter->nRow(),pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(1u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(itsIter->nPol(),pacBuf.nPol());
         testResults(pacBuf,1);                  

         // add up another accessor
         pacBuf.accumulate(*itsIter, itsME);         
         testResults(pacBuf,2);                  
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         
     }

     void testFDPAccumulate() {
         PreAvgCalBuffer pacBuf;
         CPPUNIT_ASSERT(itsME);
         CPPUNIT_ASSERT(itsIter);
         
         // simulate visibilities
         itsME->predict(*itsIter);
         
         // buffer should be initialised by the first encountered accessor,
         // frequency-dependent flag is on. We do accummulation number of
         // channel times to be able to reuse the same testing code
         pacBuf.accumulate(*itsIter, itsME, true);

         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());
         CPPUNIT_ASSERT_EQUAL(itsIter->nRow(),pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(8u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(/*itsIter->nPol()*/1u,pacBuf.nPol());
         testFrequencyDependentResults(pacBuf,1);                  

         // add up another accessor
         pacBuf.accumulate(*itsIter, itsME, true);
         testFrequencyDependentResults(pacBuf,2);                  
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());
         CPPUNIT_ASSERT_EQUAL(8u,pacBuf.nChannel());                  
     }
     
     void testFDPInitExplicit() {
         // 20 antennas instead of 30 available, 2 beams instead of 1 available in the stubbed, 8 channels
         // accessor
         PreAvgCalBuffer pacBuf(20,2,8);
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         
         // 20 antennas and 2 beams give 380 rows; 4 polarisation by default
         CPPUNIT_ASSERT_EQUAL(380u,pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(8u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(4u,pacBuf.nPol());
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nRow()),size_t(pacBuf.flag().nrow()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nPol()),size_t(pacBuf.flag().nplane()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nChannel()),size_t(pacBuf.flag().ncolumn()));     
         CPPUNIT_ASSERT(itsME);
         CPPUNIT_ASSERT(itsIter);
         
         // simulate visibilities
         itsME->predict(*itsIter);
         
         pacBuf.accumulate(*itsIter, itsME, true);
         
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         // (435 - 190) * 8 = 1960 samples unaccounted for 
         // (accessor has 1 polarisation)
         CPPUNIT_ASSERT_EQUAL(1960u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());         
         CPPUNIT_ASSERT_EQUAL(380u,pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(8u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(4u,pacBuf.nPol());
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nRow()),size_t(pacBuf.flag().nrow()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nPol()),size_t(pacBuf.flag().nplane()));
         CPPUNIT_ASSERT_EQUAL(size_t(pacBuf.nChannel()),size_t(pacBuf.flag().ncolumn()));

         //
         const scimath::PolXProducts &pxp = pacBuf.polXProducts();                            
         for (casa::uInt row=0; row<pacBuf.nRow(); ++row) {
              CPPUNIT_ASSERT_EQUAL(pacBuf.feed1()[row], pacBuf.feed2()[row]);
              const bool dataExpected = (pacBuf.feed1()[row] == 0);
              
              // test pattern is different from testFrequencyDependentResults as only the first 
              // polarisation product as the dummy accessor has data
              for (casa::uInt pol=0; pol<pacBuf.nPol(); ++pol) {
                   for (casa::uInt chan=0; chan<pacBuf.nChannel(); ++chan) {
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(double(real(pxp.getModelProduct(row,chan,pol,pol))),
                                                     double(real(pxp.getModelMeasProduct(row,chan,pol,pol))),1e-2);
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelMeasProduct(row,chan,pol,pol))),1e-5);
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(0,double(imag(pxp.getModelProduct(row,chan,pol,pol))),1e-5);
                        // 1 channel and 100 Jy source give sums of 10000 per accessor summed in
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(dataExpected && (pol == 0) ? 10000. : 0., double(real(pxp.getModelProduct(row,chan,pol,pol))),1e-2);
                        CPPUNIT_ASSERT_EQUAL(!dataExpected || (pol > 0), pacBuf.flag()(row,chan,pol));                                     
    
                        // checking cross-pol terms, if any
                        for (casa::uInt pol2 = 0; pol2<pol; ++pol2) {
                             //const casa::Complex expected = ((pol == 3) && (pol2 == 0)) ? casa::Complex(10000.*run,0.) : casa::Complex(0.,0.);
                             const casa::Complex expected(0.,0.);
                             CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,casa::abs(expected - pxp.getModelProduct(row,chan,pol,pol2)),1e-2);
                             CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,casa::abs(expected - pxp.getModelMeasProduct(row,chan,pol,pol2)),1e-2);
                        }
                        
                   }
              }
         }
         
     }
     
     void testAccumulateXPol() {
         casa::Vector<casa::Stokes::StokesTypes> stokes(4);
         stokes[0] = casa::Stokes::XX;
         stokes[1] = casa::Stokes::XY;
         stokes[2] = casa::Stokes::YX;
         stokes[3] = casa::Stokes::YY;         

         CPPUNIT_ASSERT(itsIter);       
         accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*itsIter);          
         da.itsStokes.assign(stokes.copy());
         da.itsVisibility.resize(da.nRow(), da.nChannel() ,4);
         da.itsVisibility.set(casa::Complex(-10.,15.));
         da.itsNoise.resize(da.nRow(),da.nChannel(),da.nPol());
         da.itsNoise.set(1.);
         da.itsFlag.resize(da.nRow(),da.nChannel(),da.nPol());
         da.itsFlag.set(casa::False);

         CPPUNIT_ASSERT(itsParams);
         // increase Stokes I flux, so individual polarisation products would get the same
         // value as the Stokes I in other tests (so test methods could be reused). We need to
         // regenerate the measurement equation here because the paramter values are copied at the construction
         itsParams->update("flux.i.src",200.);
         itsME.reset(new ComponentEquation(*itsParams, itsIter));
         CPPUNIT_ASSERT(itsME);
         // simulate visibilities
         itsME->predict(*itsIter);
         PreAvgCalBuffer pacBuf;

         // buffer should be initialised by the first encountered accessor
         pacBuf.accumulate(*itsIter, itsME);
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());
         CPPUNIT_ASSERT_EQUAL(itsIter->nRow(),pacBuf.nRow());
         CPPUNIT_ASSERT_EQUAL(1u,pacBuf.nChannel());
         CPPUNIT_ASSERT_EQUAL(itsIter->nPol(),pacBuf.nPol());
         testResults(pacBuf,1);                  

         // add up another accessor
         pacBuf.accumulate(*itsIter, itsME);         
         testResults(pacBuf,2);                  
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToType());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredNoMatch());
         CPPUNIT_ASSERT_EQUAL(0u,pacBuf.ignoredDueToFlags());              
     }
};
} // namespace synthesis

} // namespace askap

#endif // #ifndef PRE_AVG_CAL_BUFFER_TEST_H

