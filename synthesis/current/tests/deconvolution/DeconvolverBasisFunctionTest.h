/// @file
///
/// Unit test for the deconvolution base class
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

#include <deconvolution/DeconvolverBasisFunction.h>
#include <deconvolution/MultiScaleBasisFunction.h>
#include <deconvolution/PointBasisFunction.h>
#include <cppunit/extensions/HelperMacros.h>

#include <casa/BasicSL/Complex.h>

#include <boost/shared_ptr.hpp>

using namespace casa;

namespace askap {

namespace synthesis {

class DeconvolverBasisFunctionTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(DeconvolverBasisFunctionTest);
  CPPUNIT_TEST(testCreate);
  CPPUNIT_TEST(testDeconvolveCenter);
  CPPUNIT_TEST_EXCEPTION(testWrongShape, casa::ArrayShapeError);
  CPPUNIT_TEST_EXCEPTION(testDeconvolveOffsetPSF, AskapError);
  CPPUNIT_TEST_SUITE_END();
public:
   
  void setUp() {
    IPosition dimensions(4,100,100,1,1);
    itsDirty.reset(new Array<Float>(dimensions));
    itsDirty->set(0.0);
    itsPsf.reset(new Array<Float>(dimensions));
    itsPsf->set(0.0);
    (*itsPsf)(IPosition(4,50,50,0,0))=1.0;
    itsDB = DeconvolverBasisFunction<Float,Complex>::ShPtr(new DeconvolverBasisFunction<Float, Complex>(*itsDirty, *itsPsf));
    Vector<Float> scales(3);
    scales[0]=0.0;
    scales[1]=3.0;
    scales[2]=6.0;
    itsBasisFunction=boost::shared_ptr<BasisFunction<Float> >(new MultiScaleBasisFunction<Float>::MultiScaleBasisFunction(IPosition(4,100,100,1,1), scales));
    itsDB->setBasisFunction(itsBasisFunction);
    CPPUNIT_ASSERT(itsDB);
    CPPUNIT_ASSERT(itsDB->control());
    CPPUNIT_ASSERT(itsDB->monitor());
    CPPUNIT_ASSERT(itsDB->state());
    CPPUNIT_ASSERT(itsDB->basisFunction());
    boost::shared_ptr<DeconvolverControl<Float> > DC(new DeconvolverControl<Float>::DeconvolverControl());
    CPPUNIT_ASSERT(itsDB->setControl(DC));
    boost::shared_ptr<DeconvolverMonitor<Float> > DM(new DeconvolverMonitor<Float>::DeconvolverMonitor());
    CPPUNIT_ASSERT(itsDB->setMonitor(DM));
    boost::shared_ptr<DeconvolverState<Float> > DS(new DeconvolverState<Float>::DeconvolverState());
    CPPUNIT_ASSERT(itsDB->setControl(DC));
    itsWeight.reset(new Array<Float>(dimensions));
    itsWeight->set(10.0);
    itsDB->setWeight(*itsWeight);
    itsDB->state()->setCurrentIter(0);
    itsDB->control()->setTargetIter(10);
    itsDB->control()->setGain(1.0);
    itsDB->control()->setTargetObjectiveFunction(0.01); 
    itsDB->dirty().set(0.0);
  }

  void tearDown() {
      // Ensure arrays are destroyed last
      itsDB.reset();
      itsWeight.reset();
      itsPsf.reset();
      itsDirty.reset();
      itsBasisFunction.reset();
  }

  void testCreate() {
    Array<Float> newDirty(IPosition(4,100,100,1,1));
    itsDB->updateDirty(newDirty); 
  }
  void testWrongShape() {
    Array<Float> newDirty(IPosition(4,200,200,1,1));
    itsDB->updateDirty(newDirty);
  }
  void testDeconvolveOffsetPSF() {
    itsPsf.reset(new Array<Float>(dimensions));
    itsPsf->set(0.0);
    (*itsPsf)(IPosition(4,70,70,0,0))=1.0;
    CPPUNIT_ASSERT(itsDB->deconvolve());
    CPPUNIT_ASSERT(itsDB->control()->terminationCause()==DeconvolverControl<Float>::CONVERGED);
  }
   
  void testDeconvolveCenter() {
    itsDB->dirty()(IPosition(4,50,50,0,0))=1.0;
    CPPUNIT_ASSERT(itsDB->deconvolve());
    CPPUNIT_ASSERT(itsDB->control()->terminationCause()==DeconvolverControl<Float>::CONVERGED);
  }
   
private:

  boost::shared_ptr< Array<Float> > itsDirty;
  boost::shared_ptr< Array<Float> > itsPsf;
  boost::shared_ptr< Array<Float> > itsWeight;

   /// @brief DeconvolutionBasisFunction class
  boost::shared_ptr<DeconvolverBasisFunction<Float, Complex> > itsDB;

  boost::shared_ptr<BasisFunction<Float> > itsBasisFunction;
};
    
} // namespace synthesis

} // namespace askap

