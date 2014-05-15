/// @file
/// 
/// @brief Unit tests for GainCalibrationEquation.
/// @details GainCalibrationEquation just multiplies by a gain matrix
/// visibilities produced by another measurement equation. It also generates
/// normal equations, which allow to solve for unknowns in the gain matrix.
/// The tests given in this file attempt to predict a visibility data set
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

#ifndef CALIBRATION_ME_TEST_H
#define CALIBRATION_ME_TEST_H

#include <measurementequation/ComponentEquation.h>
#include <measurementequation/CalibrationME.h>
#include <measurementequation/PreAvgCalMEBase.h>
#include <measurementequation/NoXPolGain.h>
#include <measurementequation/NoXPolFreqDependentGain.h>
#include <measurementequation/IdentityComponent.h>
#include <measurementequation/Product.h>
#include <measurementequation/Sum.h>
#include <measurementequation/ZeroComponent.h>
#include <fitting/LinearSolver.h>
#include <dataaccess/DataIteratorStub.h>
#include <calibaccess/CalParamNameHelper.h>

#include <cppunit/extensions/HelperMacros.h>

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>

#include <boost/shared_ptr.hpp>


namespace askap
{
  namespace synthesis
  {
    using utility::toString;
    
    class CalibrationMETest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(CalibrationMETest);
      /*
      CPPUNIT_TEST(testSolveNoPreAvg);
      CPPUNIT_TEST(testSolveBPNoPreAvg);
      CPPUNIT_TEST(testSolvePreAvg);      
      CPPUNIT_TEST(testSolvePreAvg2);
      */
      CPPUNIT_TEST(testSolveBPPreAvg);      
      CPPUNIT_TEST_SUITE_END();
      
      private:
        typedef CalibrationME<Sum<Product<NoXPolGain, IdentityComponent,
                  IdentityComponent>, ZeroComponent> > METype;
        typedef CalibrationME<NoXPolFreqDependentGain> BPMEType;
        boost::shared_ptr<ComponentEquation> p1, p2;
        boost::shared_ptr<METype> eq1,eq2;
        boost::shared_ptr<Params> params1, params2;
        accessors::SharedIter<accessors::DataIteratorStub> idi;

      protected:        
        
        /// @brief helper method to take care of absolute phase uncertainty
        /// @param[in] shared pointer to parameters to update
        /// @param[in] chan if non-negative, rotation will be done for this channel and bandpass 
        /// calibration is assumed. Otherwise, normal gain calibration is implied
        static void rotatePhase(const boost::shared_ptr<Params> &params, const int chan = -1) {
          CPPUNIT_ASSERT(params);
          const std::string baseName = (chan >= 0 ? accessors::CalParamNameHelper::bpPrefix() : std::string()) + "gain";
          // taking care of the absolute phase uncertainty
          const casa::uInt refAnt = 0;
          const std::string refParamName = baseName + ".g11."+toString(refAnt)+".0" + (chan < 0 ? std::string() :
                      std::string(".") + toString(chan));
          const casa::Complex refPhaseTerm = casa::polar(1.f, 
                  -arg(params->complexValue(refParamName)));
                       
          std::vector<std::string> freeNames(params->freeNames());
          for (std::vector<std::string>::const_iterator it=freeNames.begin();
                                              it!=freeNames.end();++it)  {
               const std::string parname = *it;
               if (parname.find("gain") == 0) {
                   CPPUNIT_ASSERT(params->has(parname));                    
                   params->update(parname,
                        params->complexValue(parname)*refPhaseTerm);                                 
               } else if (parname.find("bp.gain") == 0) {
                   CPPUNIT_ASSERT(chan >= 0);
                   if (accessors::CalParamNameHelper::extractChannelInfo(parname).first == casa::uInt(chan)) {
                       params->update(parname,
                           params->complexValue(parname)*refPhaseTerm);                                 
                   }
               }
          }
        }
        
        /// @brief check gain parameters
        /// @details This method checks that all gain parameters in params2 are
        /// equal to the values from params1.
        /// @param[in] isBP if true, the bandpass solution is implied
        void checkSolution(const bool isBP = false) {
          CPPUNIT_ASSERT(params1);
          CPPUNIT_ASSERT(params2);
          const std::string baseName = (isBP ? accessors::CalParamNameHelper::bpPrefix() : std::string()) + "gain";
          
          // checking that solved gains should be close to 1 for g11 
          // and to 0.9 for g22 (we don't have data to solve for the second
          // polarisation, so it should be left unchanged)
          std::vector<std::string> completions(params2->completions(baseName));
          for (std::vector<std::string>::const_iterator it=completions.begin();
                                                it!=completions.end();++it)  {
               const std::string parname = baseName+*it;                                 
                              
               if (it->find(".g22") == 0) {
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,params2->scalarValue(parname),3e-7);
               } else if (it->find(".g11") == 0) {
                   const casa::Complex diff = params2->complexValue(parname)- 
                          (isBP ? params1->complexValue(accessors::CalParamNameHelper::extractChannelInfo("gain"+*it).second) : 
                          params1->complexValue(parname));
                   //std::cout<<parname<<" "<<diff<<" "<<abs(diff)<<std::endl;        
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,abs(diff),3e-7);
               } else {
                 ASKAPTHROW(AskapError, "an invalid gain parameter "<<parname<<" has been detected");
               }
          }
        }
        
        /// @brief prepare data and parameters
        /// @details This method predicts the data using "perfect" gains and fixes parameters which
        /// are not to be solved for. Technically, these operations can be put into setUp, but
        /// doing them in each unit test may be better because some exceptions may be thrown         
        void initDataAndParameters() {
           // Predict with the "perfect" parameters"
           eq1->predict();
           std::vector<std::string> freeNames = params2->freeNames();
           for (std::vector<std::string>::const_iterator it = freeNames.begin();
                it!=freeNames.end();++it) {
                if (it->find("gain") != 0) {
                    params2->fix(*it);
                }
           }         
        }
        
      public:        
        void setUp()
        {
          idi = boost::shared_ptr<accessors::DataIteratorStub>(new accessors::DataIteratorStub(1));
          accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*idi);
          ASKAPASSERT(da.itsStokes.nelements() == 1);
          da.itsStokes[0] = casa::Stokes::XX;
          
          const casa::uInt nAnt = 30;
          //const casa::uInt nAnt1 = 6;
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
          
          params1.reset(new Params);
          params1->add("flux.i.cena", 100.);
          params1->add("direction.ra.cena", 0.5*casa::C::arcsec);
          params1->add("direction.dec.cena", -0.3*casa::C::arcsec);
          params1->add("shape.bmaj.cena", 3.0e-3*casa::C::arcsec);
          params1->add("shape.bmin.cena", 2.0e-3*casa::C::arcsec);
          params1->add("shape.bpa.cena", -55*casa::C::degree);
          for (casa::uInt ant=0; ant<nAnt; ++ant) {
               params1->add("gain.g11."+toString(ant)+".0",
                            casa::Complex(realGains[ant],imagGains[ant]));
               params1->add("gain.g22."+toString(ant)+".0",1.);
          }

          p1.reset(new ComponentEquation(*params1, idi));
          eq1.reset(new METype(*params1,idi,p1));

          params2.reset(new Params);
          params2->add("flux.i.cena", 100.);
          params2->add("direction.ra.cena", 0.50000*casa::C::arcsec);
          params2->add("direction.dec.cena", -0.30000*casa::C::arcsec);
          params2->add("shape.bmaj.cena", 3.0e-3*casa::C::arcsec);
          params2->add("shape.bmin.cena", 2.0e-3*casa::C::arcsec);
          params2->add("shape.bpa.cena", -55*casa::C::degree);
          for (casa::uInt ant=0; ant<nAnt; ++ant) {
               params2->add("gain.g11."+toString(ant)+".0",casa::Complex(1.0,0.0));
               params2->add("gain.g22."+toString(ant)+".0",1.0);
               params2->fix("gain.g22."+toString(ant)+".0");
          }
       
          p2.reset(new ComponentEquation(*params2, idi));
          //eq2.reset(new METype(*params2,idi,p2));

        }
        
        void testSolvePreAvg() 
        {
          initDataAndParameters();
          typedef CalibrationME<NoXPolGain, PreAvgCalMEBase> PreAvgMEType;
          // preaverage and iterate over the data
          boost::shared_ptr<PreAvgMEType> preAvgEq(new PreAvgMEType(*params2));
          CPPUNIT_ASSERT(preAvgEq);
          idi.init();
          preAvgEq->accumulate(idi,p2);
          // major cycles detached from iteration over data
          for (size_t iter=0; iter<5; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne;
               preAvgEq->calcEquations(ne);
               Quality q;
               LinearSolver solver; 
               solver.addNormalEquations(ne);
               solver.setAlgorithm("SVD");
               const boost::shared_ptr<Params> params = preAvgEq->rwParameters();
               CPPUNIT_ASSERT(params);
               solver.solveNormalEquations(*params,q);  
               rotatePhase(params);
          }
          params2 = preAvgEq->parameters().clone();
          checkSolution();
        }

        /// @brief another way to handle parameters in major cycle
        /// @note Actual calculations are supposed to be identical to testSolvePreAvg2
        void testSolvePreAvg2() 
        {
          initDataAndParameters();
          typedef CalibrationME<NoXPolGain, PreAvgCalMEBase> PreAvgMEType;
          // preaverage and iterate over the data
          boost::shared_ptr<PreAvgMEType> preAvgEq(new PreAvgMEType());
          CPPUNIT_ASSERT(preAvgEq);
          idi.init();
          preAvgEq->accumulate(idi,p2);
          // major cycles detached from iteration over data
          for (size_t iter=0; iter<5; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne;
               preAvgEq->setParameters(*params2);
               preAvgEq->calcEquations(ne);
               Quality q;
               LinearSolver solver; 
               solver.addNormalEquations(ne);
               solver.setAlgorithm("SVD");
               solver.solveNormalEquations(*params2,q);  
               rotatePhase(params2);
               //std::cout<<iter<<" "<<params2->complexValue("gain.g11.0.0")<<std::endl;
          }
          checkSolution();
        }

        void testSolveNoPreAvg()
        {
          initDataAndParameters();
          for (size_t iter=0; iter<5; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne;
            
               eq2.reset(new METype(*params2,idi,p2));
            
               eq2->calcEquations(ne);
               Quality q;
               LinearSolver solver1;
               solver1.addNormalEquations(ne);
               solver1.setAlgorithm("SVD");
               solver1.solveNormalEquations(*params2,q);  
               //std::cout<<q<<std::endl;               
                              
               // taking care of the absolute phase uncertainty
               rotatePhase(params2);
          //std::cout<<*params2<<std::endl;
          }
          checkSolution();        
        }

        void testSolveBPPreAvg()
        {
          initDataAndParameters();

          std::vector<std::string> freeNames = params2->freeNames();
          for (std::vector<std::string>::const_iterator it = freeNames.begin();
               it!=freeNames.end();++it) {
               if (it->find("gain") == 0) {
                   for (casa::uInt chan = 0; chan < idi->nChannel(); ++chan) {
                        params2->add(accessors::CalParamNameHelper::addChannelInfo(accessors::CalParamNameHelper::bpPrefix()+*it,chan),
                                  params2->complexValue(*it));
                        params2->fix(*it);
                   }
               }
          }
          
          // with pre-averaging we have just one iteration over data, i.e. outside the loop
          CalibrationME<NoXPolFreqDependentGain, PreAvgCalMEBase> bpEq;
          idi.init();
          bpEq.accumulate(idi,p2);
            
          for (size_t iter=0; iter<5; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne;
                           
               bpEq.setParameters(*params2);
               bpEq.calcEquations(ne);
               Quality q;
               LinearSolver solver1;
               solver1.addNormalEquations(ne);
               solver1.setAlgorithm("SVD");
               solver1.solveNormalEquations(*params2,q);  
               //std::cout<<q<<std::endl;               
                              
               // taking care of the absolute phase uncertainty
               for (casa::uInt chan = 0; chan<idi->nChannel(); ++chan) {               
                    rotatePhase(params2,int(chan));
               }
          }
          checkSolution(true);                           
        }
        
        void testSolveBPNoPreAvg()
        {
          initDataAndParameters();

          std::vector<std::string> freeNames = params2->freeNames();
          for (std::vector<std::string>::const_iterator it = freeNames.begin();
               it!=freeNames.end();++it) {
               if (it->find("gain") == 0) {
                   for (casa::uInt chan = 0; chan < idi->nChannel(); ++chan) {
                        params2->add(accessors::CalParamNameHelper::addChannelInfo(accessors::CalParamNameHelper::bpPrefix()+*it,chan),
                                  params2->complexValue(*it));
                        params2->fix(*it);
                   }
               }
          }
          for (size_t iter=0; iter<5; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne;
            
               BPMEType bpEq(*params2,idi,p2);
            
               bpEq.calcEquations(ne);
               Quality q;
               LinearSolver solver1;
               solver1.addNormalEquations(ne);
               solver1.setAlgorithm("SVD");
               solver1.solveNormalEquations(*params2,q);  
               //std::cout<<q<<std::endl;               
                              
               // taking care of the absolute phase uncertainty
               for (casa::uInt chan = 0; chan<idi->nChannel(); ++chan) {               
                    rotatePhase(params2,int(chan));
               }
          //std::cout<<*params2<<std::endl;
          }
          checkSolution(true);                           
        }
        
   };
    
  } // namespace synthesis
} // namespace askap

#endif // #ifndef CALIBRATION_ME_TEST_H

