/// @file
///
/// Unit test for the basis function
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>

#include <deconvolution/BasisFunction.h>
#include <deconvolution/PointBasisFunction.h>
#include <deconvolution/MultiScaleBasisFunction.h>
#include <cppunit/extensions/HelperMacros.h>

#include <casa/BasicSL/Complex.h>

#include <boost/shared_ptr.hpp>

using namespace casa;

namespace askap {
  
  namespace synthesis {
    
    class BasisFunctionTest : public CppUnit::TestFixture
    {
      CPPUNIT_TEST_SUITE(BasisFunctionTest);
      CPPUNIT_TEST(testSetup);
      CPPUNIT_TEST(testPoint);
      CPPUNIT_TEST(testMultiScale);
      CPPUNIT_TEST(testMultiScaleSetShape);
      CPPUNIT_TEST_SUITE_END();
    public:
      
      void setUp() {
        itsBasisFunctionShape=IPosition(3,50,50,1);
	Vector<Float> scales(3);
	scales[0]=0.0;
	scales[1]=3.0;
	scales[2]=6.0;
	itsBasisFunction=boost::shared_ptr<BasisFunction<Float> >(new MultiScaleBasisFunction<Float>(scales));
	itsBasisFunction->initialise(itsBasisFunctionShape);
      }
      void testSetup() {
        {
          CPPUNIT_ASSERT(itsBasisFunction->basisFunction().shape()==IPosition(3,50,50,3));
          CPPUNIT_ASSERT(itsBasisFunction->numberBases()==3);
        }
      }
      void testPoint() {
        {
          itsBasisFunction.reset();
          itsBasisFunction=boost::shared_ptr<BasisFunction<Float> >(new PointBasisFunction<Float>());
	  itsBasisFunction->initialise(itsBasisFunctionShape);
          CPPUNIT_ASSERT(itsBasisFunction->basisFunction().shape()==IPosition(3,50,50,1));
          CPPUNIT_ASSERT(itsBasisFunction->numberBases()==1);
          IPosition centre(3,25,25,0);
          CPPUNIT_ASSERT(abs(itsBasisFunction->basisFunction()(centre)-1.0)<1e-6);
        }
      }
      void testMultiScale() {
        {
          itsBasisFunction.reset();
          Vector<Float> scales(3);
          scales[0]=0.0;
          scales[1]=3.0;
          scales[2]=6.0;
          itsBasisFunction=boost::shared_ptr<BasisFunction<Float> >(new MultiScaleBasisFunction<Float>(scales));
	  itsBasisFunction->initialise(itsBasisFunctionShape);
          CPPUNIT_ASSERT(itsBasisFunction->basisFunction().shape()==IPosition(3,50,50,3));
          CPPUNIT_ASSERT(itsBasisFunction->numberBases()==3);
          IPosition centre(3,25,25,0);
          centre(2)=0;
          CPPUNIT_ASSERT(abs(itsBasisFunction->basisFunction()(centre)-1.0)<1e-5);
          centre(2)=1;
          CPPUNIT_ASSERT(abs(itsBasisFunction->basisFunction()(centre)-0.192449)<1e-05);
          centre(2)=2;
          CPPUNIT_ASSERT(abs(itsBasisFunction->basisFunction()(centre)-0.048122)<1e-5);
        }
      }
      void testMultiScaleSetShape() {
        {
          itsBasisFunction.reset();
          Vector<Float> scales(3);
          scales[0]=0.0;
          scales[1]=3.0;
          scales[2]=6.0;
          itsBasisFunction=boost::shared_ptr<BasisFunction<Float> >(new MultiScaleBasisFunction<Float>(scales));
	  itsBasisFunction->initialise(IPosition(2,20,20));
          CPPUNIT_ASSERT(itsBasisFunction->basisFunction().shape()==IPosition(3,20,20,3));
          CPPUNIT_ASSERT(itsBasisFunction->numberBases()==3);
          IPosition centre(3,10,10,0);
          centre(2)=0;
          CPPUNIT_ASSERT(abs(itsBasisFunction->basisFunction()(centre)-1.0)<1e-5);
          centre(2)=1;
          CPPUNIT_ASSERT(abs(itsBasisFunction->basisFunction()(centre)-0.192449)<1e-05);
          centre(2)=2;
          CPPUNIT_ASSERT(abs(itsBasisFunction->basisFunction()(centre)-0.048122)<1e-5);
        }
      }
      void tearDown() {
        itsBasisFunction.reset();
      }
      
    private:
      boost::shared_ptr<BasisFunction<Float> > itsBasisFunction;
      
      IPosition itsBasisFunctionShape;

    };
    
  } // namespace synthesis
  
} // namespace askap

