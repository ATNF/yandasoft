/// @file
///
/// @brief Tests of routines needed for Nyquist gridding
/// @details This file contains appropriate unit tests
///
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
/// @author Daniel Mitchell <daniel.mitchell@csiro.au
///

// own includes
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


#ifndef NYQUISTGRIDDING_TESTS_H
#define NYQUISTGRIDDING_TESTS_H

#include <boost/shared_ptr.hpp>

using namespace askap;
using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class NyquistGriddingTests : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(NyquistGriddingTests);
      CPPUNIT_TEST(testReversibility);      // test reversibility for perfectly band-limited uv plane
      //CPPUNIT_TEST(testCoordinates);        // test coordinate reversibility? Not sure they're handled at this level
      //CPPUNIT_TEST(testPreconditioning);    // test that preconditioning is unaffected by downsampling
      //CPPUNIT_TEST(testHogbom);             // test that Hogbom cleaning is unaffected by downsampling
      //CPPUNIT_TEST(testCube);               // test that resolution changes work for cubes?
      CPPUNIT_TEST_SUITE_END();

    public:
        const bool printDetails = false;

    public:
        void testReversibility()
        {

          // Generate an empty uv plane
          const int n_over = 128;
          casacore::IPosition shape(2,n_over,n_over);
          casacore::Array<casacore::Complex> scratch(shape);
          casacore::ArrayLattice<casacore::Complex> scratchLattice(scratch);
          scratchLattice.set(0.);

          // Fill central region of the uv plane
          //  - don't go below row/col 32 or above row/col 96
          scratchLattice.putAt(1., casacore::IPosition(2,54,66));
          scratchLattice.putAt(1., casacore::IPosition(2,44,75));
          scratchLattice.putAt(1., casacore::IPosition(2,40,45));
          scratchLattice.putAt(1., casacore::IPosition(2,55,37));
          scratchLattice.putAt(1., casacore::IPosition(2,86,45));
          scratchLattice.putAt(1., casacore::IPosition(2,74,85));
          scratchLattice.putAt(1., casacore::IPosition(2,74,85));
          const float normFactor = shape.product() / casacore::sum(real(scratch));

          // iFFT to the image domain
          casacore::LatticeFFT::cfft2d(scratchLattice, false);

          // set image
          casacore::Array<float> dirty_orig(shape);
          casacore::ArrayLattice<float> dirtyLattice(dirty_orig);
          dirtyLattice.copyData(casacore::LatticeExpr<float> ( real(scratchLattice) ));
          dirty_orig *= normFactor;
          if (printDetails) {
              std::cout << std::endl;
              std::cout << "Original shape: "<<dirty_orig.shape() << std::endl;
          }

          for (int n_down = n_over/2; n_down < n_over/2+8; ++n_down) {

              // set the desired oversampling factor (ratio between downsampled and final output/cleaning resolution)
              const double extraOversampleFactor = double(n_over) / double(n_down);

              // save a copy
              casacore::Array<float> dirty = dirty_orig;

              // downsample
              SynthesisParamsHelper::downsample(dirty,extraOversampleFactor);
              if (printDetails) {
                  std::cout << "Nyquist shape:  "<<dirty.shape() << std::endl;
                  std::cout << " - oversampling factor = "<<extraOversampleFactor << std::endl;
              }

              // oversample
              SynthesisParamsHelper::oversample(dirty,extraOversampleFactor,false);
              if (printDetails) {
                  std::cout << " - final shape = "<<dirty.shape() << std::endl;
              }

              // compare images
              CPPUNIT_ASSERT(dirty.shape().isEqual(shape));
              if (printDetails) {
                  std::cout << " - error(n/2,n/2) = "<<std::abs(dirty(casacore::IPosition(2,64,64))-1.) << std::endl;
                  std::cout << " - max(error)     = "<<casacore::max(casacore::abs(dirty-dirty_orig)) << std::endl;
              }
              CPPUNIT_ASSERT(std::abs(dirty(casacore::IPosition(2,64,64))-1.)<1e-6);
              CPPUNIT_ASSERT(casacore::max(casacore::abs(dirty-dirty_orig))<1e-6);

          }

        }
/*
        void testPreconditioning()
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
*/
    };

  } // namespace synthesis
} // namespace askap

#endif // #ifndef NYQUISTGRIDDING_TESTS_H
