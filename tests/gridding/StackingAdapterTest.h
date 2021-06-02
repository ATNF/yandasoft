/// @copyright (c) 2021 CSIRO
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

#include <askap/gridding/WProjectVisGridder.h>
#include <askap/gridding/StackingGridderAdapter.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/measurementequation/ComponentEquation.h>
#include <askap/dataaccess/DataIteratorStub.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVPosition.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <askap/askap/AskapError.h>
#include <askap/gridding/VisGridderFactory.h>

#include <cppunit/extensions/HelperMacros.h>

#include <stdexcept>
#include <boost/shared_ptr.hpp>

using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class StackingAdapterTest : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(StackingAdapterTest);
      //      CPPUNIT_TEST(testForwardBox);
      //      CPPUNIT_TEST(testReverseBox);
      CPPUNIT_TEST(testForwardStackingGridderAdapter);
      CPPUNIT_TEST(testReverseStackingGridderAdapter);
      CPPUNIT_TEST_SUITE_END();

  private:
      //      boost::shared_ptr<BoxVisGridder> itsBox;
      boost::shared_ptr<StackingGridderAdapter> itsStacker;
      boost::shared_ptr<WProjectVisGridder> itsWProject;

      accessors::IDataSharedIter idi;
      boost::shared_ptr<Axes> itsAxes;
      boost::shared_ptr<casa::Array<imtype> > itsModel;
      boost::shared_ptr<casa::Array<imtype> > itsModelPSF;
      boost::shared_ptr<casa::Array<imtype> > itsModelWeights;

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

        itsWProject.reset(new WProjectVisGridder(10000.0, 9, 1e-3, 1, 128, 0, ""));
        itsStacker.reset(new StackingGridderAdapter(itsWProject));

        double cellSize=10*casa::C::arcsec;

        casa::Matrix<double> xform(2,2,0.);
        xform.diagonal().set(1.);

        itsAxes.reset(new Axes());
        itsAxes->addDirectionAxis(casa::DirectionCoordinate(casa::MDirection::J2000,
                     casa::Projection(casa::Projection::SIN), 0.,0.,cellSize,cellSize,xform,256.,256.));

        itsModel.reset(new casa::Array<imtype>(casa::IPosition(4,512,512,1,1)));
        itsModel->set(0.0);
        itsModelPSF.reset(new casa::Array<imtype>(casa::IPosition(4,512,512,1,1)));
        itsModelPSF->set(0.0);
        itsModelWeights.reset(new casa::Array<imtype>(casa::IPosition(4,512,512,1,1)));
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

        itsSphFunc->finaliseWeights(*itsModelWeights);
      }
      void testForwardStackingGridderAdapter()
      {
        itsStacker->initialiseDegrid(*itsAxes, *itsModel);
        itsStacker->degrid(*idi);
      }      
      void testReverseStackingGridderAdapter()
      {
        itsStacker->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsStacker->grid(*idi);
        itsStacker->finaliseGrid(*itsModel);
        itsStacker->finaliseWeights(*itsModelWeights);
        itsStacker->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsStacker->grid(*idi);
        itsStacker->finaliseGrid(*itsModelPSF);
      }

    };

  }
}
