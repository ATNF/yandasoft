/// @file
///
/// @brief Tests of the preconditioners
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
#include <askap/measurementequation/GaussianTaperPreconditioner.h>
#include <askap/measurementequation/GaussianTaperCache.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <askap/askap/AskapError.h>
#include <casacore/lattices/LatticeMath/Fit2D.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#include <casacore/lattices/LEL/LatticeExpr.h>


#include <cppunit/extensions/HelperMacros.h>


#ifndef PRECONDITIONER_TESTS_H
#define PRECONDITIONER_TESTS_H

using std::abs;

#include <boost/shared_ptr.hpp>

using namespace askap;
using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class PreconditionerTests : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(PreconditionerTests);
      CPPUNIT_TEST(testGaussianTaper);
      CPPUNIT_TEST(testGaussianTaperCache);
      CPPUNIT_TEST_SUITE_END();


    public:
        void testGaussianTaperCache()
        {
          GaussianTaperCache gtc(25.,15.,-M_PI/18.);
          casacore::Array<float> taper = gtc.taper(casacore::IPosition(4,128,128,1,1));

          std::vector<std::string> direction(3);
          direction[0]="12h30m00.0";
          direction[1]="-15.00.00.00";
          direction[2]="J2000";

          std::vector<int> shape(2,128);
          std::vector<std::string> cellsize(2,"1arcsec");
          casacore::Vector<casacore::Stokes::StokesTypes> stokes(1, casacore::Stokes::I);

          scimath::Params params(true);
          SynthesisParamsHelper::add(params,"psf.test",direction,cellsize,shape,false,1.4e9,
                              1.4e9,1,stokes);

          #ifdef ASKAP_FLOAT_IMAGE_PARAMS
          params.update("psf.test",taper);
          #else
          casacore::Array<double> temp(taper.shape());
          casacore::convertArray<double,float>(temp, taper);
          params.update("psf.test",temp);
          #endif

          casacore::Vector<casacore::Quantum<double> > fit = SynthesisParamsHelper::fitBeam(params,0.05,101,"psf.test");
          CPPUNIT_ASSERT(fit.nelements() == 3);
          // the cell size is 1 arcsec, so the tolerance of 0.1 arcsec seems good enough
          CPPUNIT_ASSERT(fabs(fit[0].getValue("arcsec")-25.)<0.1);
          CPPUNIT_ASSERT(fabs(fit[1].getValue("arcsec")-15.)<0.1);
          CPPUNIT_ASSERT(fabs(fit[2].getValue("rad") + M_PI/18.)<0.1);


        }
        void testGaussianTaper()
        {
          GaussianTaperPreconditioner gtp(25.,15.,M_PI/18.);
          casacore::IPosition shape(2,128,128);
          casacore::Array<float> psf(shape), dirty(shape), pcf;
          psf.set(0.); dirty.set(0.);
          dirty(casacore::IPosition(2,64,64)) = 1.;

          casacore::IPosition index(2);
          const double fwhm2sigma = sqrt(8.*log(2.));
          for (index[0] = 0; index[0]<128; ++index[0]) {
               for (index[1] = 0; index[1]<128; ++index[1]) {
                    const double xOffset = (double(index[0])-64.);
                    const double yOffset = (double(index[1])-64.);
                    const double expFactor = exp(-casacore::square(xOffset/2*fwhm2sigma)/2.-
                            casacore::square(yOffset/1.3*fwhm2sigma)/2.);
                    psf(index) = expFactor;
               }
          }

          gtp.doPreconditioning(psf,dirty,pcf);

          casacore::ArrayLattice<float> psfLattice(psf);
          casacore::ArrayLattice<casacore::Complex> scratch(psf.shape());
          scratch.copyData(casacore::LatticeExpr<casacore::Complex>(toComplex(psfLattice)));
          casacore::LatticeFFT::cfft2d(scratch, true);
          psfLattice.copyData(casacore::LatticeExpr<float> ( real(scratch) ));


          casacore::LogIO logger;
          casacore::Fit2D fitter(logger);
          casacore::Vector<casacore::Double> param = fitter.estimate(casacore::Fit2D::GAUSSIAN, psf);
          fitter.addModel(casacore::Fit2D::GAUSSIAN, param);
          casacore::Array<float> sigma(psf.shape());
          sigma.set(1.);
          CPPUNIT_ASSERT(fitter.fit(psf,sigma) == casacore::Fit2D::OK);
          param = fitter.availableSolution();

          /// @todo We need to revisit normalisation factors at some stage,
          /// I am still not happy.
          //std::cout<<param<<std::endl;
          CPPUNIT_ASSERT(param.size() == 6);
          CPPUNIT_ASSERT(std::abs(param[1]-64.)<1e-5);
          CPPUNIT_ASSERT(std::abs(param[2]-64.)<1e-5);
          CPPUNIT_ASSERT(std::abs(param[3]-25.)<1);
          CPPUNIT_ASSERT(std::abs(param[4]-15.)<1);
          CPPUNIT_ASSERT(std::abs(param[5]/M_PI*180.-100.)<1);
        }
    };

  } // namespace synthesis
} // namespace askap

#endif // #ifndef PRECONDITIONER_TESTS_H
