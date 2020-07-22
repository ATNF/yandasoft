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

#ifndef CALIBRATION_DD_TEST_H
#define CALIBRATION_DD_TEST_H

#include <askap/measurementequation/ComponentEquation.h>
#include <askap/measurementequation/CalibrationME.h>
#include <askap/measurementequation/PreAvgDDCalMEBase.h>
#include <askap/measurementequation/NoXPolGain.h>
#include <askap/measurementequation/NoXPolFreqDependentGain.h>
#include <askap/measurementequation/IdentityComponent.h>
#include <askap/measurementequation/Product.h>
#include <askap/measurementequation/Sum.h>
#include <askap/measurementequation/ZeroComponent.h>
#include <askap/scimath/fitting/LinearSolver.h>
#include <askap/scimath/fitting/JonesIndex.h>
#include <askap/dataaccess/DataIteratorStub.h>
#include <askap/dataaccess/IConstDataAccessor.h>
#include <askap/dataaccess/IDataAccessor.h>
#include <askap/dataaccess/DDCalBufferDataAccessor.h>
#include <askap/calibaccess/CalParamNameHelper.h>

#include <cppunit/extensions/HelperMacros.h>

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>

#include <boost/shared_ptr.hpp>


namespace askap
{
  namespace synthesis
  {
    using utility::toString;
    
    class CalibrationDDTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(CalibrationDDTest);
      CPPUNIT_TEST(testSolveAmp);      
      CPPUNIT_TEST(testSolvePhase);
      CPPUNIT_TEST_SUITE_END();
      
      private:
        typedef CalibrationME<Sum<Product<NoXPolGain, IdentityComponent,
                  IdentityComponent>, ZeroComponent>> METype;
        boost::shared_ptr<ComponentEquation> p1, p2, p3;
        boost::shared_ptr<METype> eq1;
        boost::shared_ptr<Params> params1, params2, params3;
        accessors::SharedIter<accessors::DataIteratorStub> idi;
        casacore::Vector<casa::Float> ampError, lOffset, mOffset;
        const casacore::uInt nAnt = 30;
        const casacore::uInt nDir = 2;

      protected:        
        
        /// @brief helper method to take care of absolute phase uncertainty
        /// @param[in] shared pointer to parameters to update
        /// @param[in] chan if non-negative, rotation will be done for this channel and bandpass 
        /// calibration is assumed. Otherwise, normal gain calibration is implied
        static void rotatePhase(const boost::shared_ptr<Params> &params, const int chan = -1) {
          CPPUNIT_ASSERT(params);
          const std::string baseName = (chan >= 0 ? accessors::CalParamNameHelper::bpPrefix() : std::string()) + "gain";
          // taking care of the absolute phase uncertainty
          const casacore::uInt refAnt = 0;
          const std::string refParamName = baseName + ".g11."+toString(refAnt)+".0" + (chan < 0 ? std::string() :
                      std::string(".") + toString(chan));
          const casacore::Complex refPhaseTerm = casacore::polar(1.f, 
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
                   if (accessors::CalParamNameHelper::extractChannelInfo(parname).first == casacore::uInt(chan)) {
                       params->update(parname,
                           params->complexValue(parname)*refPhaseTerm);                                 
                   }
               }
          }
        }
          
        /// @brief reset amplitudes
        /// @details This method resets the amplitude of params3 gains to the params1 values
        /// This is a simple hack to force phase-only solutions.
        void resetAmplitudes(const boost::shared_ptr<Params> &params) {
          CPPUNIT_ASSERT(params1);
          CPPUNIT_ASSERT(params);
          const std::string baseName = "gain";
          
          std::vector<std::string> completions(params->completions(baseName));
          for (std::vector<std::string>::const_iterator it=completions.begin();
                                                it!=completions.end();++it)  {
               const std::string parname = baseName+*it;                                 
               const JonesIndex Jindex = accessors::CalParamNameHelper::parseParam(parname).first;

               if (it->find(".g22") == 0) {
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,params->scalarValue(parname),3e-7);
               } else if (it->find(".g11") == 0) {
                   const std::string refname = "gain.g11."+toString(Jindex.antenna())+".0";
                   params->update(parname, params->complexValue(parname) *
                        abs(params1->complexValue(refname)) / abs(params->complexValue(parname)));
               } else {
                 ASKAPTHROW(AskapError, "an invalid gain parameter "<<parname<<" has been detected");
               }
          }
        }
      
        /// @brief check amplitudes
        /// @details This method checks that the amplitude of gain parameters in params2
        /// relative to params1 is consistent with the multiplicative model error.
        void checkAmplitudes() {
          CPPUNIT_ASSERT(params1);
          CPPUNIT_ASSERT(params2);
          const std::string baseName = "gain";
          
          std::vector<std::string> completions(params2->completions(baseName));
          for (std::vector<std::string>::const_iterator it=completions.begin();
                                                it!=completions.end();++it)  {
               const std::string parname = baseName+*it;                                 
               const JonesIndex Jindex = accessors::CalParamNameHelper::parseParam(parname).first;

               if (it->find(".g22") == 0) {
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,params2->scalarValue(parname),3e-7);
               } else if (it->find(".g11") == 0) {
                   const std::string refname = "gain.g11."+toString(Jindex.antenna())+".0";
                   const casacore::Float gainRatio =
                        abs(params1->complexValue(refname)) / abs(params2->complexValue(parname));
                   //std::cout<<parname<<" (abs("<<params1->complexValue(refname)<<") / "<<
                   //                      "abs("<<params2->complexValue(parname)<<"))^2 = "<<
                   //                      gainRatio*gainRatio<<" => "<<
                   //                      gainRatio*gainRatio - (1.+ampError[Jindex.beam()])<<std::endl;
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,abs(gainRatio*gainRatio - (1.+ampError[Jindex.beam()])),1e-5);
               } else {
                 ASKAPTHROW(AskapError, "an invalid gain parameter "<<parname<<" has been detected");
               }
          }
        }
       
        /// @brief check phases
        /// @details This method checks that the phase of gain parameters in params3
        /// relative to params1 is consistent with the position error in the model.
        void checkPhases() {
          CPPUNIT_ASSERT(params1);
          CPPUNIT_ASSERT(params3);
          const std::string baseName = "gain";
          
          // 
          std::vector<std::string> completions(params3->completions(baseName));
          for (std::vector<std::string>::const_iterator it=completions.begin();
                                                it!=completions.end();++it)  {
               const std::string parname = baseName+*it;                                 
               const JonesIndex Jindex = accessors::CalParamNameHelper::parseParam(parname).first;

               if (it->find(".g22") == 0) {
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,params3->scalarValue(parname),3e-7);
               } else if (it->find(".g11") == 0) {
                   const std::string refname = "gain.g11."+toString(Jindex.antenna())+".0";
                   const casacore::Float gainRatio =
                        abs(params1->complexValue(refname)) / abs(params3->complexValue(parname));
                   //std::cout<<parname<<" abs("<<params1->complexValue(refname)<<") / "<<
                   //                     "abs("<<params3->complexValue(parname)<<") = "<< gainRatio<<std::endl;
                   // ensure that the amplitudes have been correctly normalised
                   CPPUNIT_ASSERT_DOUBLES_EQUAL(gainRatio,1.,3e-7);
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
                    params3->fix(*it);
                }
           }         
        }
        
      public:        
        void setUp()
        {
          // iterator over data.
          idi = boost::shared_ptr<accessors::DataIteratorStub>(new accessors::DataIteratorStub(1));
          // data accessor associated with the iterator
          accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*idi);
          ASKAPASSERT(da.itsStokes.nelements() == 1);
          da.itsStokes[0] = casacore::Stokes::XX;
          
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
          params1->add("flux.i.src1", 100.);
          params1->add("direction.ra.src1", 5.0*casacore::C::arcmin);
          params1->add("direction.dec.src1", -13.0*casacore::C::arcmin);
          params1->add("shape.bmaj.src1", 3.0e-3*casacore::C::arcsec);
          params1->add("shape.bmin.src1", 2.0e-3*casacore::C::arcsec);
          params1->add("shape.bpa.src1", -55*casacore::C::degree);
          params1->add("flux.i.src2", 100.);
          params1->add("direction.ra.src2", -126.0*casacore::C::arcmin);
          params1->add("direction.dec.src2", 32.0*casacore::C::arcmin);
          for (casacore::uInt ant=0; ant<nAnt; ++ant) {
               //params1->add("gain.g11."+toString(ant)+".0",casacore::Complex(1.0,0.0));
               params1->add("gain.g11."+toString(ant)+".0",
                            casacore::Complex(realGains[ant],imagGains[ant]));
               params1->add("gain.g22."+toString(ant)+".0",1.);
          }

          p1.reset(new ComponentEquation(*params1, idi));
          eq1.reset(new METype(*params1,idi,p1));

          // multiplicative gain errors, e.g. from a primary beam model error
          ampError.resize(nDir);
          ampError[0] =  0.1;
          ampError[1] = -0.1;

          params2.reset(new Params);
          params2->add("flux.i.src1", 100.*(1.0 + ampError[0]));
          params2->add("direction.ra.src1", 5.0*casacore::C::arcmin);
          params2->add("direction.dec.src1", -13.0*casacore::C::arcmin);
          params2->add("shape.bmaj.src1", 3.0e-3*casacore::C::arcsec);
          params2->add("shape.bmin.src1", 2.0e-3*casacore::C::arcsec);
          params2->add("shape.bpa.src1", -55*casacore::C::degree);
          params2->add("flux.i.src2", 100.*(1.0 + ampError[1]));
          params2->add("direction.ra.src2", -126.0*casacore::C::arcmin);
          params2->add("direction.dec.src2", 32.0*casacore::C::arcmin);
          for (casacore::uInt dir=0; dir<nDir; ++dir) {
              for (casacore::uInt ant=0; ant<nAnt; ++ant) {
                  //params2->add("gain.g11."+toString(ant)+"."+toString(dir),casacore::Complex(1.0,0.0));
                  params2->add("gain.g11."+toString(ant)+"."+toString(dir),
                               casacore::Complex(realGains[ant],imagGains[ant]));
                  params2->add("gain.g22."+toString(ant)+"."+toString(dir),1.0);
                  params2->fix("gain.g22."+toString(ant)+"."+toString(dir));
              }
          }
       
          p2.reset(new ComponentEquation(*params2, idi));

          // additive position errors, e.g. an ionospheric shift at a single frequency
          // synthesised beam is about 25", so set offset based on this
          lOffset.resize(nDir);
          mOffset.resize(nDir);
          lOffset[0] =  0.0*casacore::C::arcsec;
          mOffset[0] =  0.0*casacore::C::arcsec;
          lOffset[1] =  2.0*casacore::C::arcsec;
          mOffset[1] = -5.0*casacore::C::arcsec;

          params3.reset(new Params);
          params3->add("flux.i.src1", 100.);
          params3->add("direction.ra.src1", 5.0*casacore::C::arcmin + lOffset[0]);
          params3->add("direction.dec.src1", -13.0*casacore::C::arcmin + mOffset[0]);
          params3->add("shape.bmaj.src1", 3.0e-3*casacore::C::arcsec);
          params3->add("shape.bmin.src1", 2.0e-3*casacore::C::arcsec);
          params3->add("shape.bpa.src1", -55*casacore::C::degree);
          params3->add("flux.i.src2", 100.);
          params3->add("direction.ra.src2", -126.0*casacore::C::arcmin + lOffset[1]);
          params3->add("direction.dec.src2", 32.0*casacore::C::arcmin + mOffset[1]);
          for (casacore::uInt dir=0; dir<nDir; ++dir) {
              for (casacore::uInt ant=0; ant<nAnt; ++ant) {
                  //params3->add("gain.g11."+toString(ant)+"."+toString(dir),casacore::Complex(1.0,0.0));
                  params3->add("gain.g11."+toString(ant)+"."+toString(dir),
                               casacore::Complex(realGains[ant],imagGains[ant]));
                  params3->add("gain.g22."+toString(ant)+"."+toString(dir),1.0);
                  params3->fix("gain.g22."+toString(ant)+"."+toString(dir));
              }
          }
       
          p3.reset(new ComponentEquation(*params3, idi));

        }
        
        void testSolveAmp() 
        {

          initDataAndParameters();
          boost::shared_ptr<PreAvgDDCalMEBase> preAvgEq;
          preAvgEq.reset(new CalibrationME<NoXPolGain, PreAvgDDCalMEBase>(*params2));
          preAvgEq->initialise(nAnt, nDir, 1);
          CPPUNIT_ASSERT(preAvgEq);
          idi.init();
          p2->setNDir(nDir);
          preAvgEq->accumulate(idi,p2);
          // major cycles detached from iteration over data
          for (size_t iter=0; iter<10; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne;
               preAvgEq->calcEquations(ne);
               Quality q;
               LinearSolver solver; 
               solver.addNormalEquations(ne);
               solver.setAlgorithm("LSQR");
               const boost::shared_ptr<Params> params = preAvgEq->rwParameters();
               CPPUNIT_ASSERT(params);
               solver.solveNormalEquations(*params,q);  
               rotatePhase(params);
          }
          params2 = preAvgEq->parameters().clone();
          checkAmplitudes();

        }

        /// @brief another way to handle parameters in major cycle
        /// @note Actual calculations are supposed to be identical to testSolvePhase
        void testSolvePhase() 
        {

          initDataAndParameters();
          boost::shared_ptr<PreAvgDDCalMEBase> preAvgEq;
          preAvgEq.reset(new CalibrationME<NoXPolGain, PreAvgDDCalMEBase>(*params3));
          preAvgEq->initialise(nAnt, nDir, 1);
          CPPUNIT_ASSERT(preAvgEq);
          idi.init();
          p3->setNDir(nDir);
          preAvgEq->accumulate(idi,p3);
          // major cycles detached from iteration over data
          for (size_t iter=0; iter<10; ++iter) {
               // Calculate gradients using "imperfect" parameters"
               GenericNormalEquations ne;
               preAvgEq->calcEquations(ne);
               Quality q;
               LinearSolver solver; 
               solver.addNormalEquations(ne);
               solver.setAlgorithm("LSQR");
               const boost::shared_ptr<Params> params = preAvgEq->rwParameters();
               CPPUNIT_ASSERT(params);
               solver.solveNormalEquations(*params,q);  
               //rotatePhase(params);
               resetAmplitudes(params);
          }
          params3 = preAvgEq->parameters().clone();
          checkPhases();

        }

      private:

   };
    
  } // namespace synthesis
} // namespace askap

#endif // #ifndef CALIBRATION_DD_TEST_H

