/// @file
///
/// @brief Tests of the functionality provided by VectorOperations
/// @details This file contains appropriate unit tests
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
/// @author Max Voronkov <maxim.voronkov@csiro.au
///



// own includes
#include <measurementequation/VectorOperations.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <scimath/Mathematics/AutoDiff.h>
#include <casa/BasicSL/Complex.h>
#include <cppunit/extensions/HelperMacros.h>
#include <askap/AskapError.h>
#include <fitting/ComplexDiff.h>


// std includes
#include <cmath>
#include <vector>

#ifndef VECTOR_OPERATIONS_TEST_H
#define VECTOR_OPERATIONS_TEST_H

using std::abs;

#include <boost/shared_ptr.hpp>

using namespace askap;
using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class VectorOperationsTest : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(VectorOperationsTest);
      CPPUNIT_TEST(testCopy);
      CPPUNIT_TEST(testSubtract);
      CPPUNIT_TEST(testAdd);
      CPPUNIT_TEST(testAddScaled);
      CPPUNIT_TEST(test1); // temporary
      CPPUNIT_TEST_SUITE_END();

      
    public:
        void testCopy()
        {
          casa::Matrix<casa::Double> in(2,2,1.);
          vector<double> vec(2,3.);
          vec[0]=-3.;
          copyVector(vec,in.row(0));
          CPPUNIT_ASSERT(fabs(in(0,0)+3.)<1e-10 && fabs(in(0,1)-3.)<1e-10);
          casa::Vector<casa::Complex> complexVec(1,casa::Complex(-1.,-2.));
          copyVector(complexVec,in.row(1));
          CPPUNIT_ASSERT(fabs(in(1,0)+1.)<1e-10 && fabs(in(1,1)+2.)<1e-10);
          vector<casa::AutoDiff<double> > autoDiffVec(2, 
                                          casa::AutoDiff<double>(0.,1));
          autoDiffVec[0]=sin(casa::AutoDiff<double>(0.,1,0));
          autoDiffVec[1]=1.+cos(casa::AutoDiff<double>(casa::C::_2pi/4.,1,0));
          copyVector(autoDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0])<1e-10 && fabs(vec[1]-1)<1e-10);
          copyDerivativeVector(0,autoDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0]-1)<1e-10 && fabs(vec[1]+1)<1e-10);
          vector<scimath::ComplexDiff> complexDiffVec(1,
                        scimath::ComplexDiff("par1",casa::Complex(0.,-1.)));
          complexDiffVec[0] *= scimath::ComplexDiff("par2", casa::Complex(2.,0.));
          copyVector(complexDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0])<1e-10 && fabs(vec[1]+2.)<1e-10);
          copyReDerivativeVector("par1",complexDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0]-2)<1e-10 && fabs(vec[1])<1e-10);
          copyReDerivativeVector("par2",complexDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0])<1e-10 && fabs(vec[1]+1)<1e-10);
          copyImDerivativeVector("par1",complexDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0])<1e-10 && fabs(vec[1]-2)<1e-10);
          copyImDerivativeVector("par2",complexDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0]-1)<1e-10 && fabs(vec[1])<1e-10);
        }
        
        void testSubtract()
        {
          casa::Matrix<casa::Double> in(2,2,1.);
          vector<double> vec(2,3.);
          vec[0]=-3.;
          subtractVector(vec,in.row(1));
          CPPUNIT_ASSERT(fabs(in(0,0)-1.)<1e-10 && fabs(in(0,1)-1.)<1e-10);
          CPPUNIT_ASSERT(fabs(in(1,0)-4.)<1e-10 && fabs(in(1,1)+2.)<1e-10);
          casa::Vector<casa::Complex> complexVec(1,casa::Complex(-1.,-2.));
          subtractVector(complexVec,in.row(1));
          CPPUNIT_ASSERT(fabs(in(1,0)-5.)<1e-10 && fabs(in(1,1))<1e-10);
          vector<casa::AutoDiff<double> > autoDiffVec(2, 
                                          casa::AutoDiff<double>(0.,1));
          autoDiffVec[0]=sin(casa::AutoDiff<double>(0.,1,0));
          autoDiffVec[1]=1.+cos(casa::AutoDiff<double>(casa::C::_2pi/4.,1,0));
          subtractVector(autoDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0]+3.)<1e-10 && fabs(vec[1]-2.)<1e-10);
          subtractVector(complexVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0]+2.)<1e-10 && fabs(vec[1]-4.)<1e-10);
        }
        
        void testAdd()
        {
          casa::Matrix<casa::Double> in(2,2,1.);
          vector<double> vec(2,3.);
          vec[0]=-3.;
          addVector(vec,in.row(1));
          CPPUNIT_ASSERT(fabs(in(0,0)-1.)<1e-10 && fabs(in(0,1)-1.)<1e-10);
          CPPUNIT_ASSERT(fabs(in(1,0)+2.)<1e-10 && fabs(in(1,1)-4.)<1e-10);
          casa::Vector<casa::Complex> complexVec(1,casa::Complex(-1.,-2.));
          addVector(complexVec,in.row(1));
          CPPUNIT_ASSERT(fabs(in(1,0)+3.)<1e-10 && fabs(in(1,1)-2.)<1e-10);
          vector<casa::AutoDiff<double> > autoDiffVec(2, 
                                          casa::AutoDiff<double>(0.,1));
          autoDiffVec[0]=sin(casa::AutoDiff<double>(0.,1,0));
          autoDiffVec[1]=1.+cos(casa::AutoDiff<double>(casa::C::_2pi/4.,1,0));
          addVector(autoDiffVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0]+3.)<1e-10 && fabs(vec[1]-4.)<1e-10);
          addVector(complexVec,vec);
          CPPUNIT_ASSERT(fabs(vec[0]+4.)<1e-10 && fabs(vec[1]-2.)<1e-10);
        }

        void testAddScaled()
        {
          casa::Matrix<casa::Double> in(2,2,1.);
          vector<double> vec(2,3.);
          vec[0]=-3.;
          addScaledVector(vec,in.row(1),0.5);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1., in(0,0), 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1., in(0,1), 1e-10);                    
          CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, in(1,0), 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5, in(1,1), 1e-10);          
          casa::Vector<casa::Complex> complexVec(1,casa::Complex(-1.,-2.));
          addScaledVector(complexVec,in.row(1),0.5);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(-1, in(1,0), 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5, in(1,1), 1e-10);          
          vector<casa::AutoDiff<double> > autoDiffVec(2, 
                                          casa::AutoDiff<double>(0.,1));
          autoDiffVec[0]=sin(casa::AutoDiff<double>(0.,1,0));
          autoDiffVec[1]=1.+cos(casa::AutoDiff<double>(casa::C::_2pi/4.,1,0));
          addScaledVector(autoDiffVec,vec,0.5);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(-3, vec[0], 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(3.5, vec[1], 1e-10);          
          
          addScaledVector(complexVec,vec,0.5);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.5, vec[0], 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5, vec[1], 1e-10);          
          addScaledVector(vec, complexVec, casa::Complex(-1., 2.));
          CPPUNIT_ASSERT_DOUBLES_EQUAL(-2.5, real(complexVec[0]), 1e-10);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(-11.5, imag(complexVec[0]), 1e-10);
        }
         
        // temporary method for autodiff experiments
        void test1()
        {
          vector<casa::AutoDiff<casa::Complex> > autoDiffVec(2, 
                                          casa::AutoDiff<casa::Complex>(casa::Complex(0.,0.),1));
          autoDiffVec[0]=sin(casa::AutoDiff<casa::Complex>(0.,1,0))+casa::Complex(0,1.)*cos(casa::AutoDiff<casa::Complex>(0.,1,0));
          //autoDiffVec[1]=casa::Complex(1.)+cos(casa::AutoDiff<casa::Complex>(casa::C::_2pi/4.,1,0));
          autoDiffVec[1]=exp(casa::AutoDiff<casa::Complex>(casa::C::_2pi/4.,1,0)*casa::Complex(0,1.));
          //std::cout<<autoDiffVec[0].value()<<" "<<autoDiffVec[0].derivative(0)<<std::endl;
          //std::cout<<autoDiffVec[1].value()<<" "<<autoDiffVec[1].derivative(0)<<std::endl;
          
        }

    };

  } // namespace synthesis
} // namespace askap

#endif // #ifndef VECTOR_OPERATIONS_TEST_H
