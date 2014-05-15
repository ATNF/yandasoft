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

#include <measurementequation/ImageDFTEquation.h>
#include <fitting/LinearSolver.h>
#include <dataaccess/DataIteratorStub.h>
#include <casa/aips.h>
#include <casa/Arrays/Matrix.h>
#include <measures/Measures/MPosition.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/MVPosition.h>
#include <casa/BasicSL/Constants.h>

#include <cppunit/extensions/HelperMacros.h>

#include <askap/AskapError.h>

#include <cmath>

using std::abs;

#include <boost/shared_ptr.hpp>

using namespace askap;
using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class ImageDFTEquationTest : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(ImageDFTEquationTest);
      CPPUNIT_TEST(testPredict);
      CPPUNIT_TEST(testSVD);
      CPPUNIT_TEST_EXCEPTION(testFixed, AskapError);
      CPPUNIT_TEST_SUITE_END();

      private:
        ImageDFTEquation *p1, *p2;
        Params *params1, *params2;
        accessors::IDataSharedIter idi;

      public:
        void setUp()
        {
          idi = accessors::IDataSharedIter(new accessors::DataIteratorStub(1));

          uint npix=16;
          Axes imageAxes;
          double arcsec=casa::C::pi/(3600.0*180.0);
          casa::Matrix<double> xform(2,2,0.);
          xform.diagonal().set(1.);
               
          imageAxes.addDirectionAxis(casa::DirectionCoordinate(casa::MDirection::J2000, 
                     casa::Projection(casa::Projection::SIN), 0.,0.,15.*arcsec,15.*arcsec,xform,npix/2.,npix/2.));
          
          params1 = new Params;
          casa::Array<double> imagePixels1(casa::IPosition(2, npix, npix));
          imagePixels1.set(0.0);
          imagePixels1(casa::IPosition(2, npix/2, npix/2))=1.0;
          imagePixels1(casa::IPosition(2, 12, 3))=0.7;
          params1->add("image.i.cena", imagePixels1, imageAxes);

          p1 = new ImageDFTEquation(*params1, idi);

          params2 = new Params;
          casa::Array<double> imagePixels2(casa::IPosition(2, npix, npix));
          imagePixels2.set(0.0);
          imagePixels2(casa::IPosition(2, npix/2, npix/2))=0.9;
          imagePixels2(casa::IPosition(2, 12, 3))=0.75;
          params2->add("image.i.cena", imagePixels2, imageAxes);

          p2 = new ImageDFTEquation(*params2, idi);

        }

        void tearDown()
        {
          delete p1;
          delete p2;
        }

        void testPredict()
        {
          p1->predict();
        }

        void testSVD()
        {
// Calculate gradients using "imperfect" parameters"
          p1->predict();
// Predict with the "perfect" parameters"
          GenericNormalEquations ne; //(*params2);
          p2->calcEquations(ne);
          {
            LinearSolver solver1;
            solver1.addNormalEquations(ne);
            Quality q;
            solver1.setAlgorithm("SVD");
            solver1.solveNormalEquations(*params2,q);
            casa::Array<double> improved = params2->value("image.i.cena");
            uint npix=16;
            std::cout << q << std::endl;
            CPPUNIT_ASSERT(std::abs(q.cond()/1115634013709.060-1.0)<0.0001);
            CPPUNIT_ASSERT(std::abs(improved(casa::IPosition(2, npix/2, npix/2))-1.0)<0.003);
            std::cout << improved(casa::IPosition(2, npix/2, npix/2))-1.0 << std::endl;
            CPPUNIT_ASSERT(std::abs(improved(casa::IPosition(2, 12, 3))-0.700)<0.003);
            std::cout << improved(casa::IPosition(2, 12, 3))-0.700 << std::endl;
          }
        }

        void testFixed()
        {
          p1->predict();
          GenericNormalEquations ne; //(*params2);
          p2->calcEquations(ne);
          Quality q;
          params2->fix("image.i.cena");
          LinearSolver solver1;
          solver1.addNormalEquations(ne);
          solver1.solveNormalEquations(*params2,q);
        }
    };

  }
}
