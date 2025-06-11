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

#include <askap/deconvolution/DeconvolverMultiTermBasisFunction.h>
#include <askap/deconvolution/MultiScaleBasisFunction.h>
#include <askap/deconvolution/PointBasisFunction.h>
#include <cppunit/extensions/HelperMacros.h>

#include <casacore/casa/BasicSL/Complex.h>

#include <boost/shared_ptr.hpp>

using namespace casa;

namespace askap {

namespace synthesis {

class DeconvolverMultiTermBasisFunctionTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(DeconvolverMultiTermBasisFunctionTest);
  CPPUNIT_TEST(testCreate);
  CPPUNIT_TEST(testDeconvolveCenter1);
  CPPUNIT_TEST(testDeconvolveCenter2);
  CPPUNIT_TEST_EXCEPTION(testWrongShape, AskapError);
  CPPUNIT_TEST_EXCEPTION(testDeconvolveOffsetPSF, AskapError);
  CPPUNIT_TEST(testDeconvolvePSFSubset);
  CPPUNIT_TEST(testScaleMask);
  CPPUNIT_TEST_SUITE_END();
public:

  void setUp() {
    ASKAPLOG_INFO_STR(logger,"test setUp");
    IPosition dimensions(2,100,100);
    itsDirty.reset(new Array<Float>(dimensions));
    itsDirty->set(0.0);
    itsPsf.reset(new Array<Float>(dimensions));
    itsPsf->set(0.0);
    (*itsPsf)(IPosition(2,50,50))=1.0;
    itsDB = DeconvolverMultiTermBasisFunction<Float,Complex>::ShPtr(new DeconvolverMultiTermBasisFunction<Float, Complex>(*itsDirty, *itsPsf));
    Vector<Float> scales(3);
    scales[0]=0.0;
    scales[1]=3.0;
    scales[2]=6.0;
    itsBasisFunction=boost::shared_ptr<BasisFunction<Float> >(new MultiScaleBasisFunction<Float>(IPosition(4,100,100,1,1), scales));
    itsDB->setBasisFunction(itsBasisFunction);
    CPPUNIT_ASSERT(itsDB);
    CPPUNIT_ASSERT(itsDB->control());
    CPPUNIT_ASSERT(itsDB->monitor());
    CPPUNIT_ASSERT(itsDB->state());
    CPPUNIT_ASSERT(itsDB->basisFunction());
    boost::shared_ptr<DeconvolverControl<Float> > DC(new DeconvolverControl<Float>());
    CPPUNIT_ASSERT(itsDB->setControl(DC));
    boost::shared_ptr<DeconvolverMonitor<Float> > DM(new DeconvolverMonitor<Float>());
    CPPUNIT_ASSERT(itsDB->setMonitor(DM));
    boost::shared_ptr<DeconvolverState<Float> > DS(new DeconvolverState<Float>());
    CPPUNIT_ASSERT(itsDB->setControl(DC));
    itsWeight.reset(new Array<Float>(dimensions));
    itsWeight->set(10.0);
    itsDB->setWeight(*itsWeight);
    itsDB->state()->setCurrentIter(0);
    itsDB->control()->setTargetIter(10);
    itsDB->control()->setGain(1.0);
    itsDB->control()->setTargetObjectiveFunction(0.001);
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
    ASKAPLOG_INFO_STR(logger,"test Create");
    itsDirty.reset(new Array<Float>(IPosition(2,100,100)));
    itsDirty->set(1.0);
    itsDB->updateDirty(*itsDirty);
  }
  void testWrongShape() {
    ASKAPLOG_INFO_STR(logger,"test WrongShape");
    itsDirty.reset(new Array<Float>(IPosition(2,200,200)));
    itsDirty->set(0.0);
    itsDB->updateDirty(*itsDirty);
  }
  void testDeconvolveOffsetPSF() {
    ASKAPLOG_INFO_STR(logger,"test DeconvolveOffset");
    itsDB->dirty()(IPosition(2,30,20))=1.0;
    itsDB->psf().set(0.0);
    itsDB->psf()(IPosition(2,70,70))=1.0;
    CPPUNIT_ASSERT(itsDB->deconvolve());
    CPPUNIT_ASSERT(itsDB->control()->terminationCause()==DeconvolverControl<Float>::CONVERGED);
  }


  void testDeconvolveCenter1() {
      testDeconvolveCenter("MAXBASE");
  }
  void testDeconvolveCenter2() {
      testDeconvolveCenter("MAXCHISQ");
  }

  void testDeconvolveCenter(std::string solutionType) {
    ASKAPLOG_INFO_STR(logger,"test DeconvolveCentre - "+solutionType);
    itsDB->setSolutionType(solutionType);
    itsDB->dirty().set(0.0);
    itsDB->dirty()(IPosition(2,50,50))=1.0;
    CPPUNIT_ASSERT(itsDB->deconvolve());
    CPPUNIT_ASSERT(itsDB->control()->terminationCause()==DeconvolverControl<Float>::CONVERGED);
  }

  void testDeconvolvePSFSubset() {
    // this test was added in response to a bug where reuse of the deconvolver triggered shape mismatch when
    // a subset of the PSF was used for deconvolution
    ASKAPLOG_INFO_STR(logger,"test DeconvolvePSFSubset");
    itsDB->dirty().set(0.0);
    itsDB->dirty()(IPosition(2,50,50))=1.0;
    itsDB->control()->setPSFWidth(50);
    CPPUNIT_ASSERT(itsDB->deconvolve());
    CPPUNIT_ASSERT(itsDB->control()->terminationCause()==DeconvolverControl<Float>::CONVERGED);
    // second "major cycle"
    itsDB->dirty().set(0.0);
    itsDB->dirty()(IPosition(2,49,49))=1.0;
    itsDB->state()->setCurrentIter(0);
    CPPUNIT_ASSERT(itsDB->deconvolve());
    CPPUNIT_ASSERT(itsDB->control()->terminationCause()==DeconvolverControl<Float>::CONVERGED);
  }

  void testScaleMask() {
      ASKAPLOG_INFO_STR(logger,"test ScaleMask");
      itsDB->dirty()(IPosition(2,45,45),IPosition(2,55,55))=1.0;
      itsDB->control()->setTargetObjectiveFunction2(0.1);
      itsDB->control()->setTargetIter(15);
      CPPUNIT_ASSERT(itsDB->dirty()(IPosition(2,50,50))==1.0);
      CPPUNIT_ASSERT(itsDB->deconvolve());
      CPPUNIT_ASSERT(itsDB->control()->terminationCause()==DeconvolverControl<Float>::EXCEEDEDITERATIONS);
      Array<float> mask = itsDB->scaleMask();
      itsDB->setScaleMask(mask);
      Array<float> mask2 = itsDB->scaleMask();
      for( uint j=0; j<mask.nelements(); j++) {
        CPPUNIT_ASSERT_EQUAL(mask.data()[j],mask2.data()[j]);
      }
  }

private:

  boost::shared_ptr< Array<Float> > itsDirty;
  boost::shared_ptr< Array<Float> > itsPsf;
  boost::shared_ptr< Array<Float> > itsWeight;

   /// @brief DeconvolutionBasisFunction class
  boost::shared_ptr<DeconvolverMultiTermBasisFunction<Float, Complex> > itsDB;

  boost::shared_ptr<BasisFunction<Float> > itsBasisFunction;
};

} // namespace synthesis

} // namespace askap
