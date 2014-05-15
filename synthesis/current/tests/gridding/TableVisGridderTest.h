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

//#include <gridding/BoxVisGridder.h>
#include <gridding/SphFuncVisGridder.h>
#include <gridding/AWProjectVisGridder.h>
#include <gridding/AProjectWStackVisGridder.h>
#include <gridding/WStackVisGridder.h>
#include <gridding/AWProjectVisGridder.h>
#include <gridding/WProjectVisGridder.h>
#include <gridding/ATCAIllumination.h>
#include <gridding/DiskIllumination.h>
#include <fitting/Params.h>
#include <measurementequation/ComponentEquation.h>
#include <dataaccess/DataIteratorStub.h>
#include <casa/aips.h>
#include <casa/Arrays/Matrix.h>
#include <measures/Measures/MPosition.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/MVPosition.h>
#include <casa/BasicSL/Constants.h>
#include <askap/AskapError.h>
#include <gridding/VisGridderFactory.h>

#include <cppunit/extensions/HelperMacros.h>

#include <stdexcept>
#include <boost/shared_ptr.hpp>

using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class TableVisGridderTest : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(TableVisGridderTest);
      //      CPPUNIT_TEST(testForwardBox);
      //      CPPUNIT_TEST(testReverseBox);
      CPPUNIT_TEST_EXCEPTION(testCreateAbstract,AskapError);      
      CPPUNIT_TEST_EXCEPTION(testUnknownGridder,AskapError);      
      CPPUNIT_TEST(testForwardSph);
      CPPUNIT_TEST(testReverseSph);
      CPPUNIT_TEST(testForwardAWProject);
      CPPUNIT_TEST(testReverseAWProject);
      CPPUNIT_TEST(testForwardWProject);
      CPPUNIT_TEST(testReverseWProject);
      CPPUNIT_TEST(testForwardWStack);
      CPPUNIT_TEST(testReverseWStack);
      CPPUNIT_TEST(testForwardAWProject);
      CPPUNIT_TEST(testReverseAWProject);
      CPPUNIT_TEST(testForwardAProjectWStack);
      CPPUNIT_TEST(testReverseAProjectWStack);
      CPPUNIT_TEST(testForwardATCAIllumination);
      CPPUNIT_TEST(testReverseATCAIllumination);
      CPPUNIT_TEST_SUITE_END();

  private:
      //      boost::shared_ptr<BoxVisGridder> itsBox;
      boost::shared_ptr<SphFuncVisGridder> itsSphFunc;
      boost::shared_ptr<AWProjectVisGridder> itsAWProject;
      boost::shared_ptr<WProjectVisGridder> itsWProject;
      boost::shared_ptr<WStackVisGridder> itsWStack;
      boost::shared_ptr<AProjectWStackVisGridder> itsAProjectWStack;

      accessors::IDataSharedIter idi;
      boost::shared_ptr<Axes> itsAxes;
      boost::shared_ptr<casa::Array<double> > itsModel;
      boost::shared_ptr<casa::Array<double> > itsModelPSF;
      boost::shared_ptr<casa::Array<double> > itsModelWeights;

  public:
      void setUp()
      {
        idi = accessors::IDataSharedIter(new accessors::DataIteratorStub(1));

        Params ip;
        ip.add("flux.i.cena", 100.0);
        ip.add("direction.ra.cena", 0.5);
        ip.add("direction.dec.cena", -0.3);
        ip.add("shape.bmaj.cena", 0.0);
        ip.add("shape.bmin.cena", 0.0);
        ip.add("shape.bpa.cena", 0.0);

        ComponentEquation ce(ip, idi);
        ce.predict();

	//        itsBox.reset(new BoxVisGridder());
        itsSphFunc.reset(new SphFuncVisGridder());
        boost::shared_ptr<IBasicIllumination> illum(new DiskIllumination(120.0, 10.0));
        // illumination models can easily be shared between gridders if parameters,
        // like dish size and blockage are the same
        itsAWProject.reset(new AWProjectVisGridder(illum, 10000.0, 9, 1e-3, 1, 128, 1));
        itsAProjectWStack.reset(new AProjectWStackVisGridder(illum, 10000.0, 9, 0., 1, 128, 1));
        itsWProject.reset(new WProjectVisGridder(10000.0, 9, 1e-3, 1, 128, 0, ""));
        itsWStack.reset(new WStackVisGridder(10000.0, 9));

        double cellSize=10*casa::C::arcsec;

        casa::Matrix<double> xform(2,2,0.);
        xform.diagonal().set(1.);
               
        itsAxes.reset(new Axes());
        itsAxes->addDirectionAxis(casa::DirectionCoordinate(casa::MDirection::J2000, 
                     casa::Projection(casa::Projection::SIN), 0.,0.,cellSize,cellSize,xform,256.,256.));
        
        itsModel.reset(new casa::Array<double>(casa::IPosition(4,512,512,1,1)));
        itsModel->set(0.0);
        itsModelPSF.reset(new casa::Array<double>(casa::IPosition(4,512,512,1,1)));
        itsModelPSF->set(0.0);
        itsModelWeights.reset(new casa::Array<double>(casa::IPosition(4,512,512,1,1)));
        itsModelWeights->set(0.0);
      }

      void tearDown()
      {
      }

      //      void testReverseBox()
      //      {
      //        itsBox->initialiseGrid(*itsAxes, itsModel->shape(), true);
      //        itsBox->grid(idi);
      //        itsBox->finaliseGrid(*itsModel);
      //        itsBox->finalisePSF(*itsModelPSF);
      //        itsBox->finaliseWeights(*itsModelWeights);
      //      }
      //      void testForwardBox()
      //      {
      //        itsBox->initialiseDegrid(*itsAxes, *itsModel);
      //        itsBox->degrid(idi);
      //      }
      
      void testCreateAbstract()
      {
         // TableVisGridder is still an abstract class
         // calling createGridder static method should raise an
         // exception
         LOFAR::ParameterSet parset;
         TableVisGridder::createGridder(parset);
      }
      
      void testUnknownGridder()
      {
         // here we attempt to create a gridder, which is not recognized by a factory
         // The code should try to load a properly named dynamic library and because
         // it is not there an exception will be thrown 
         LOFAR::ParameterSet parset;
         parset.add("gridder","AGridderWithSuchNameShouldDefinitelyBeUnknownToTheSystsem");
         VisGridderFactory::make(parset);
      }
      
      void testReverseSph()
      {
        itsSphFunc->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsSphFunc->grid(*idi);
        itsSphFunc->finaliseGrid(*itsModelPSF);
        itsSphFunc->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsSphFunc->grid(*idi);
        itsSphFunc->finaliseGrid(*itsModel);
        itsSphFunc->finaliseWeights(*itsModelWeights);
      }
      void testForwardSph()
      {
        itsSphFunc->initialiseDegrid(*itsAxes, *itsModel);
        itsSphFunc->degrid(*idi);
      }
      void testReverseAWProject()
      {
        itsAWProject->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsAWProject->grid(*idi);
        itsAWProject->finaliseGrid(*itsModel);
        itsAWProject->finaliseWeights(*itsModelWeights);
        itsAWProject->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsAWProject->grid(*idi);
        itsAWProject->finaliseGrid(*itsModelPSF);
      }
      void testForwardAWProject()
      {
        itsAWProject->initialiseDegrid(*itsAxes, *itsModel);
        itsAWProject->degrid(*idi);
      }
      void testReverseWProject()
      {
        itsWProject->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsWProject->grid(*idi);
        itsWProject->finaliseGrid(*itsModel);
        itsWProject->finaliseWeights(*itsModelWeights);
        itsWProject->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsWProject->grid(*idi);
        itsWProject->finaliseGrid(*itsModelPSF);
      }
      void testForwardWProject()
      {
        itsWProject->initialiseDegrid(*itsAxes, *itsModel);
        itsWProject->degrid(*idi);
      }
      void testReverseWStack()
      {
        itsWStack->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsWStack->grid(*idi);
        itsWStack->finaliseGrid(*itsModel);
        itsWStack->finaliseWeights(*itsModelWeights);
        itsWStack->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsWStack->grid(*idi);
        itsWStack->finaliseGrid(*itsModelPSF);
      }
      void testForwardWStack()
      {
        itsWStack->initialiseDegrid(*itsAxes, *itsModel);
        itsWStack->degrid(*idi);
      }
      void testReverseAProjectWStack()
      {
        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsAProjectWStack->grid(*idi);
        itsAProjectWStack->finaliseGrid(*itsModel);
        itsAProjectWStack->finaliseWeights(*itsModelWeights);
        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsAProjectWStack->grid(*idi);
        itsAProjectWStack->finaliseGrid(*itsModelPSF);
      }
      void testForwardAProjectWStack()
      {
        itsAProjectWStack->initialiseDegrid(*itsAxes, *itsModel);
        itsAProjectWStack->degrid(*idi);
      }
      void testReverseATCAIllumination()
      {
        boost::shared_ptr<ATCAIllumination> illum(new ATCAIllumination(120.0, 10.0));
        illum->simulateTapering(1.0);
        itsAProjectWStack.reset(new AProjectWStackVisGridder(illum, 10000.0, 9, 0., 1, 128, 1));
        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsAProjectWStack->grid(*idi);
        itsAProjectWStack->finaliseGrid(*itsModel);
        itsAProjectWStack->finaliseWeights(*itsModelWeights);
        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsAProjectWStack->grid(*idi);
        itsAProjectWStack->finaliseGrid(*itsModelPSF);
      }
      void testForwardATCAIllumination()
      {
        boost::shared_ptr<ATCAIllumination> illum(new ATCAIllumination(12.0, 2.0));
        illum->simulateTapering(1.0);
        itsAProjectWStack.reset(new AProjectWStackVisGridder(illum, 10000.0, 9, 0., 1, 128, 1));
        itsAProjectWStack->initialiseDegrid(*itsAxes, *itsModel);
        itsAProjectWStack->degrid(*idi);
      }
    };

  }
}
