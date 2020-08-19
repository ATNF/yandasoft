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

#ifndef CALIBRATION_IONO_TEST_H
#define CALIBRATION_IONO_TEST_H

#include <askap/measurementequation/ComponentEquation.h>
#include <askap/measurementequation/CalibrationME.h>
#include <askap/measurementequation/PreAvgDDCalMEBase.h>
#include <askap/measurementequation/IonosphericTerm.h>
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
#include <complex>

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>

#include <boost/shared_ptr.hpp>


namespace askap
{
  namespace synthesis
  {
    using utility::toString;
    
    class CalibrationIonoTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(CalibrationIonoTest);
      //CPPUNIT_TEST(testSolveAmp);      
      CPPUNIT_TEST(testSolveOffset);
      CPPUNIT_TEST_SUITE_END();
      
      private:
        typedef CalibrationME<Sum<Product<IonosphericTerm, IdentityComponent,
                  IdentityComponent>, ZeroComponent>> METype;
        boost::shared_ptr<ComponentEquation> p1, p2;
        boost::shared_ptr<METype> eq1;
        boost::shared_ptr<Params> params1, params2;
        accessors::SharedIter<accessors::DataIteratorStub> idi;
        const casacore::uInt nAnt = 30;
        const casacore::uInt nDir = 2;

      protected:        
          
        /// @brief check ionospheric offsets
        /// @details This method checks that the ionospheric offset parameters
        /// in params2 match those of params1
        void checkOffsets() {
            CPPUNIT_ASSERT(params1);
            CPPUNIT_ASSERT(params2);

            const std::string baseName = "ionosphere";
            std::vector<std::string> completions(params2->completions(baseName));

            for (std::vector<std::string>::const_iterator it=completions.begin(); it!=completions.end();++it)  {
                 const std::string parname = baseName+*it;                                 
           
                 std::cout <<parname<<" "<<params1->complexValue(parname)/float(casacore::C::arcmin)<<" - "<<
                                           params2->complexValue(parname)/float(casacore::C::arcmin)<<" = "<<
                                           params1->complexValue(parname)/float(casacore::C::arcmin) - 
                                           params2->complexValue(parname)/float(casacore::C::arcmin) << std::endl;
           
                 // check that the real parts are equal
                 CPPUNIT_ASSERT_DOUBLES_EQUAL(params2->complexValue(parname).real(),
                                              params2->complexValue(parname).real(),1e-9);
                 // check that the imaginary parts are equal (should be zero)
                 CPPUNIT_ASSERT_DOUBLES_EQUAL(params2->complexValue(parname).imag(),
                                              params2->complexValue(parname).imag(),1e-9);
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
                if (it->find("ionosphere") != 0) {
                    params2->fix(*it);
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

          boost::shared_ptr<IonosphericTerm> iono;
          casacore::Vector<casacore::RigidVector<double,3> > xyz = iono->getAntennaPositions(da);
          CPPUNIT_ASSERT(xyz.shape()==nAnt);

          params1.reset(new Params);
          params1->add("source.src1", 0);
          params1->add("flux.i.src1", 100.);
          params1->add("direction.ra.src1", 0.0*casacore::C::arcmin);
          params1->add("direction.dec.src1", -90.0*casacore::C::arcmin);
          params1->add("shape.bmaj.src1", 3.0e-3*casacore::C::arcsec);
          params1->add("shape.bmin.src1", 2.0e-3*casacore::C::arcsec);
          params1->add("shape.bpa.src1", -55*casacore::C::degree);
          params1->add("source.src2", 1);
          params1->add("flux.i.src2", 100.);
          params1->add("direction.ra.src2", -126.0*casacore::C::arcmin);
          params1->add("direction.dec.src2", 32.0*casacore::C::arcmin);
          // antenna-dependent fixed parameters
          for (casacore::uInt ant=0; ant<nAnt; ++ant) {
               params1->add("antenna.position.x."+toString(ant), xyz(ant)(0));
               params1->add("antenna.position.y."+toString(ant), xyz(ant)(1));
               params1->add("antenna.position.z."+toString(ant), xyz(ant)(2));
               params1->add("gain.g11."+toString(ant)+".0",
                            casacore::Complex(realGains[ant],imagGains[ant]));
               params1->add("gain.g22."+toString(ant)+".0",1.);
          }
          // free parameters. work in lambda^2 once the basic setup is established
          for (casacore::uInt dir=0; dir<nDir; ++dir) {
               params1->add("ionosphere.coeff.l."+toString(dir), 1.0*casacore::C::arcmin);
               params1->add("ionosphere.coeff.m."+toString(dir), 2.0*casacore::C::arcmin);
          }

// uvw info is not always available in IonosphericTerm, and eventually xyz staton position data will be needed.
// Add that now as parameters?
// uij = Xi - Xj
// one option:
//  - take the first N-1 visibilities
//  - call station 0 the reference and then add: X0' = 0; Xj' = u0j
//  - then go through and subtract the average
//  - however, that will all be relative to the phase centre at the first time step
//  - just do it like this for now and fix later
// could also require that array lat,long are supplied and calculate rotatedUVW for that plane
//  - 

          p1.reset(new ComponentEquation(*params1, idi));
          eq1.reset(new METype(*params1,idi,p1));

          params2.reset(new Params);
          params2->add("source.src1", 0);
          params2->add("flux.i.src1", 100.);
          params2->add("direction.ra.src1", 0.0*casacore::C::arcmin);
          params2->add("direction.dec.src1", -90.0*casacore::C::arcmin);
          params2->add("shape.bmaj.src1", 3.0e-3*casacore::C::arcsec);
          params2->add("shape.bmin.src1", 2.0e-3*casacore::C::arcsec);
          params2->add("shape.bpa.src1", -55*casacore::C::degree);
          params2->add("source.src2", 1);
          params2->add("flux.i.src2", 100.);
          params2->add("direction.ra.src2", -126.0*casacore::C::arcmin);
          params2->add("direction.dec.src2", 32.0*casacore::C::arcmin);
          // antenna-dependent fixed parameters
          for (casacore::uInt ant=0; ant<nAnt; ++ant) {
              for (casacore::uInt dir=0; dir<nDir; ++dir) {
                  params2->add("gain.g11."+toString(ant)+"."+toString(dir),
                               casacore::Complex(realGains[ant],imagGains[ant]));
                  params2->add("gain.g22."+toString(ant)+"."+toString(dir),1.0);
              }
              params2->add("antenna.position.x."+toString(ant), xyz(ant)(0));
              params2->add("antenna.position.y."+toString(ant), xyz(ant)(1));
              params2->add("antenna.position.z."+toString(ant), xyz(ant)(2));
          }
          // free parameters. work in lambda^2 once the basic setup is established
          for (casacore::uInt dir=0; dir<nDir; ++dir) {
              params2->add("ionosphere.coeff.l."+toString(dir), 0.0*casacore::C::arcmin);
              params2->add("ionosphere.coeff.m."+toString(dir), 0.0*casacore::C::arcmin);
          }

          p2.reset(new ComponentEquation(*params2, idi));

        }
        
        /// @brief another way to handle parameters in major cycle
        /// @note Actual calculations are supposed to be identical to testSolveOffset
        void testSolveOffset() 
        {

          initDataAndParameters();
          boost::shared_ptr<PreAvgDDCalMEBase> preAvgEq;
          preAvgEq.reset(new CalibrationME<IonosphericTerm, PreAvgDDCalMEBase>(*params2));
          preAvgEq->initialise(nAnt, nDir, 1);
          CPPUNIT_ASSERT(preAvgEq);
          idi.init();
          p2->setNDir(nDir);
          preAvgEq->accumulate(idi,p2);
          for (size_t iter=0; iter<5; ++iter) {
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
          }
          params2 = preAvgEq->parameters().clone();
          checkOffsets();

        }

      private:

   };
    
  } // namespace synthesis
} // namespace askap

#endif // #ifndef CALIBRATION_IONO_TEST_H

