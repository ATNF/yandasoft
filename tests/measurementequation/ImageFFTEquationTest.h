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

#include <askap/measurementequation/ImageFFTEquation.h>
//#include <askap/gridding/BoxVisGridder.h>
#include <askap/gridding/SphFuncVisGridder.h>
#include <askap/gridding/AWProjectVisGridder.h>
#include <askap/gridding/DiskIllumination.h>
#include <askap/measurementequation/ImageSolver.h>
#include <askap/measurementequation/ImageMultiScaleSolver.h>
#include <askap/dataaccess/DataIteratorStub.h>
#include <askap/scimath/fitting/ParamsCasaTable.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVPosition.h>
#include <casacore/casa/BasicSL/Constants.h>

#include <askap/askap/AskapError.h>

#include <cppunit/extensions/HelperMacros.h>

//#include <askap/measurementequation/SynthesisParamsHelper.h>
//#include <casacore/casa/Arrays/ArrayMath.h>

#include <stdexcept>

#include <boost/shared_ptr.hpp>

using namespace askap;
using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class ImageFFTEquationTest : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(ImageFFTEquationTest);
      CPPUNIT_TEST(testCopyAndReference);
      CPPUNIT_TEST(testPredict);
      CPPUNIT_TEST(testSolveSphFun);
      CPPUNIT_TEST(testSolveAntIllum);
      CPPUNIT_TEST_EXCEPTION(testFixed, CheckError);
      CPPUNIT_TEST(testFullPol);
      CPPUNIT_TEST_SUITE_END();

  private:
      boost::shared_ptr<ImageFFTEquation> p1, p2;
      boost::shared_ptr<Params> params1, params2;
      accessors::IDataSharedIter idi;
      uint npix;

  public:
      void setUp()
      {
        idi = accessors::IDataSharedIter(new accessors::DataIteratorStub(1));

        npix=1024;
        Axes imageAxes;
        double arcsec=casacore::C::pi/(3600.0*180.0);
        double cell=8.0*arcsec;
        casacore::Matrix<double> xform(2,2,0.);
        xform.diagonal().set(1.);

        imageAxes.addDirectionAxis(casacore::DirectionCoordinate(casacore::MDirection::J2000,
                     casacore::Projection(casacore::Projection::SIN), 0.,0.,cell,cell,xform,npix/2.,npix/2.));

        imageAxes.addStokesAxis(casacore::Vector<casacore::Stokes::StokesTypes>(1,casacore::Stokes::I));
        imageAxes.add("FREQUENCY",1.4e9,1.4e9);

        params1.reset(new Params(true));
        casacore::Array<double> imagePixels1(casacore::IPosition(4, npix, npix, 1, 1));
        imagePixels1.set(0.0);
        imagePixels1(casacore::IPosition(4, npix/2, npix/2, 0, 0))=1.0;
        imagePixels1(casacore::IPosition(4, 3*npix/8, 7*npix/16, 0, 0))=0.7;
        params1->add("image.i.cena", imagePixels1, imageAxes);

        p1.reset(new ImageFFTEquation(*params1, idi));

        params2.reset(new Params(true));
        casacore::Array<double> imagePixels2(casacore::IPosition(4, npix, npix, 1, 1));
        imagePixels2.set(0.0);
        imagePixels2(casacore::IPosition(4, npix/2, npix/2, 0, 0))=0.9;
        imagePixels2(casacore::IPosition(4, 3*npix/8, 7*npix/16, 0, 0))=0.75;
        params2->add("image.i.cena", imagePixels2, imageAxes);
        p2.reset(new ImageFFTEquation(*params2, idi));

      }

      void testCopyAndReference()
      {
          boost::shared_ptr<ImageFFTEquation> p3, p4;
          boost::shared_ptr<Params> params3, params4;
          params3.reset(new Params(*params1));
          params4.reset(new Params(*params1));

          IVisGridder::ShPtr gridder=IVisGridder::ShPtr(new SphFuncVisGridder());

          // create using copied params
          p3.reset(new ImageFFTEquation(*params3,idi, gridder));
          *params3 = *params2;
          casacore::IPosition pos(4, npix/2, npix/2, 0, 0);
          // check values are not affected
          CPPUNIT_ASSERT(abs(p3->parameters().valueT("image.i.cena")(pos)-1.0)<0.001);
          // create using referenced params
          p4.reset(new ImageFFTEquation(params4,idi, gridder));
          *params4 = *params2;
          // check values are affected (reference)
          CPPUNIT_ASSERT(abs(p4->parameters().valueT("image.i.cena")(pos)-0.9)<0.001);
      }

      void testPredict()
      {
        //        {
        //          ParamsCasaTable pt("ImageFFTEquationTest_original.tab", false);
        //          pt.setParameters(*params1);
        //        }
        p1->predict();
      }

      void testFullPol() {
         // tests of the full stokes simulation
         CPPUNIT_ASSERT(params1);
         CPPUNIT_ASSERT(params2);

         casacore::Array<imtype> pix = params1->valueT("image.i.cena").copy().reform(casacore::IPosition(4,npix/2,npix/2,4,1));
         pix.set(0.);
         pix(casacore::IPosition(4, npix/4, npix/4, 0, 0))=1.0;
         pix(casacore::IPosition(4, npix/4, npix/4, 1, 0))=0.01;
         pix(casacore::IPosition(4, npix/4, npix/4, 2, 0))=-0.01;
         pix(casacore::IPosition(4, npix/4, npix/4, 3, 0))=0.9;
         //pix(casacore::IPosition(4, 3*npix/16, 7*npix/32, 0, 0))=0.7;
         casacore::Vector<casacore::Stokes::StokesTypes> stokes(4);
         stokes[0] = casacore::Stokes::XX;
         stokes[1] = casacore::Stokes::XY;
         stokes[2] = casacore::Stokes::YX;
         stokes[3] = casacore::Stokes::YY;
         params1->axes("image.i.cena").addStokesAxis(stokes);
         // overwrite direction axis because we now have a smaller image
         const double arcsec=casacore::C::pi/(3600.0*180.0);
         const double cell=8.0*arcsec;
         casacore::Matrix<double> xform(2,2,0.);
         xform.diagonal().set(1.);
         params1->axes("image.i.cena").addDirectionAxis(casacore::DirectionCoordinate(casacore::MDirection::J2000,
                     casacore::Projection(casacore::Projection::SIN), 0.,0.,cell,cell,xform,npix/4.,npix/4.));

         params1->valueT("image.i.cena").assign(pix);
         p1.reset(new ImageFFTEquation(*params1, idi));

         accessors::DataAccessorStub &da = dynamic_cast<accessors::DataAccessorStub&>(*idi);
         da.itsStokes.assign(stokes.copy());
         //da.itsVisibility.resize(da.nRow(), 2 ,4);
         da.itsVisibility.resize(4,2,da.nRow());
         da.itsVisibility.set(casacore::Complex(-10.,15.));
         //da.itsNoise.resize(da.nRow(),da.nChannel(),da.nPol());
         da.itsNoise.resize(da.nPol(),da.nChannel(),da.nRow());
         da.itsNoise.set(1.);
         //da.itsFlag.resize(da.nRow(),da.nChannel(),da.nPol());
         da.itsFlag.resize(da.nPol(),da.nChannel(),da.nRow());
         da.itsFlag.set(casacore::False);


         p1->predict();
         CPPUNIT_ASSERT(da.nPol() == 4);
         CPPUNIT_ASSERT(da.nChannel() == 2);

         for (casacore::uInt row=0; row<da.nRow(); ++row) {
              for (casacore::uInt ch=0; ch<da.nChannel(); ++ch) {
                   //CPPUNIT_ASSERT(casacore::abs(casacore::DComplex(1.,0.) - casacore::DComplex(da.visibility()(row,ch,0)))<1e-5);
                   //CPPUNIT_ASSERT(casacore::abs(casacore::DComplex(0.01,0.) - casacore::DComplex(da.visibility()(row,ch,1)))<1e-5);
                   //CPPUNIT_ASSERT(casacore::abs(casacore::DComplex(-0.01,0.) - casacore::DComplex(da.visibility()(row,ch,2)))<1e-5);
                   //CPPUNIT_ASSERT(casacore::abs(casacore::DComplex(0.9,0.) - casacore::DComplex(da.visibility()(row,ch,3)))<1e-5);

                   CPPUNIT_ASSERT(casacore::abs(casacore::DComplex(1.,0.) - casacore::DComplex(da.visibility()(0,ch,row)))<1e-5);
                   CPPUNIT_ASSERT(casacore::abs(casacore::DComplex(0.01,0.) - casacore::DComplex(da.visibility()(1,ch,row)))<1e-5);
                   CPPUNIT_ASSERT(casacore::abs(casacore::DComplex(-0.01,0.) - casacore::DComplex(da.visibility()(2,ch,row)))<1e-5);
                   CPPUNIT_ASSERT(casacore::abs(casacore::DComplex(0.9,0.) - casacore::DComplex(da.visibility()(3,ch,row)))<1e-5);
              }
         }
      }

      void testSolveSphFun()
      {
        // Predict with the "perfect" parameters"
        p1->predict();
        // perform a given number of major cycles
        const size_t nMajCycles = 1;
        casacore::Array<imtype> improved; // buffer for the result
        for (size_t cycle = 0; cycle<nMajCycles; ++cycle) {
             // Calculate gradients using "imperfect" parameters"
             ImagingNormalEquations ne(*params2);
             p2->calcEquations(ne);
             Quality q;
             ImageMultiScaleSolver solver1;
             solver1.setAlgorithm("Hogbom");
             solver1.addNormalEquations(ne);
             solver1.solveNormalEquations(*params2,q);
             improved = params2->valueT("image.i.cena");
        }
        /*
        casacore::Array<float> dbg(improved.shape());
        casacore::convertArray<float,double>(dbg,improved);
        SynthesisParamsHelper::saveAsCasaImage("dbg.img",dbg);
        */
        // This only works for the pixels with emission but it's a good test nevertheless
        CPPUNIT_ASSERT(abs(improved(casacore::IPosition(4, npix/2, npix/2, 0, 0))
            -1.0)<0.003);
        CPPUNIT_ASSERT(abs(improved(casacore::IPosition(4, 3*npix/8, 7*npix/16, 0,
            0))-0.700)<0.003);
      }

      void testSolveAntIllum()
      {
        // Predict with the "perfect" parameters"
        boost::shared_ptr<IBasicIllumination> illum(new DiskIllumination(12.0, 1.0));

        IVisGridder::ShPtr gridder=IVisGridder::ShPtr(new AWProjectVisGridder(illum,
                                          8000, 9, 1e-3, 8, 512, 0));
        p1.reset(new ImageFFTEquation(*params1, idi, gridder));
        p2.reset(new ImageFFTEquation(*params2, idi, gridder));
        p1->predict();
        // Calculate gradients using "imperfect" parameters"
        ImagingNormalEquations ne(*params2);
        p2->calcEquations(ne);
        Quality q;
        ImageSolver solver1;
        solver1.addNormalEquations(ne);
        solver1.solveNormalEquations(*params2,q);
        const casacore::Array<imtype> improved = params2->valueT("image.i.cena");
        // This only works for the pixels with emission but it's a good test nevertheless
        CPPUNIT_ASSERT(abs(improved(casacore::IPosition(4, npix/2, npix/2, 0, 0))
            -1.0)<0.005);
        CPPUNIT_ASSERT(abs(improved(casacore::IPosition(4, 3*npix/8, 7*npix/16, 0,
            0))-0.700)<0.005);
      }

      void testFixed()
      {
        ImagingNormalEquations ne(*params1);
        p1->predict();
        p2->calcEquations(ne);
        Quality q;
        params2->fix("image.i.cena");
        ImageSolver solver1;
        solver1.addNormalEquations(ne);
        solver1.solveNormalEquations(*params2,q);
      }
    };

  }
}
