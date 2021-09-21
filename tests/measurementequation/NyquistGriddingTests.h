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
#include <askap/measurementequation/WienerPreconditioner.h>
#include <Common/ParameterSet.h>
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
      CPPUNIT_TEST(testPreconditioning);    // test that preconditioning is unaffected by downsampling
      //CPPUNIT_TEST(testIO);                 // test that imported sky model matches original?
      //CPPUNIT_TEST(testCube);               // test that resolution changes work for cubes?
      CPPUNIT_TEST_SUITE_END();

    private:
        const bool printDetails = false;
        const int n = 128;
        casacore::Array<float> dirty_orig;
        casacore::Array<float> psf_orig;
        casacore::Array<float> pcf_orig;

    public:
        void setUp()
        {
            // Generate an empty uv plane
            const casacore::IPosition shape(2,n,n);
            casacore::Array<casacore::Complex> scratch(shape);
            casacore::ArrayLattice<casacore::Complex> scratchLattice(scratch);
           
            // Fill central region of the uv plane
            //  - don't go below row/col n/4 or above row/col 3*n/4
            const std::vector<int> u = {54,44,44,40,55,86,74,74,74,74,74,74,74,74,74};
            const std::vector<int> v = {66,75,75,45,37,45,85,85,85,85,85,85,85,85,85};
            const float normFactor = shape.product() / float(u.size());
           
            // Form the dirty image, starting with gridded visibilities
            scratchLattice.set(0.);
            for (int k=0; k<u.size(); ++k) {
                scratch(casacore::IPosition(2,u[k],v[k])) += 1.; // could include a phase term
            }
            casacore::LatticeFFT::cfft2d(scratchLattice, false); // iFFT to the image domain
            dirty_orig = casacore::Array<float>(shape);
            casacore::ArrayLattice<float> dirtyLattice(dirty_orig);
            dirtyLattice.copyData(casacore::LatticeExpr<float> ( real(scratchLattice) ));
            dirty_orig *= normFactor;
           
            // Form the PSF image, starting with gridded visibilities
            scratchLattice.set(0.);
            for (int k=0; k<u.size(); ++k) {
                scratch(casacore::IPosition(2,u[k],v[k])) += 1.;
            }
            casacore::LatticeFFT::cfft2d(scratchLattice, false); // iFFT to the image domain
            psf_orig = casacore::Array<float>(shape);
            casacore::ArrayLattice<float> psfLattice(psf_orig);
            psfLattice.copyData(casacore::LatticeExpr<float> ( real(scratchLattice) ));
            psf_orig *= normFactor;
           
            // Form the PCF image, starting with gridded visibilities
            scratchLattice.set(0.);
            for (int k=0; k<u.size(); ++k) {
                // weighted w kernel size is accumulated in the imaginary part, but needs conjugate symmetry
                const float w = v > n/2 ? 1. : -1.; // more complicated if v==n/2, but that isn't the case
                scratch(casacore::IPosition(2,u[k],v[k])) += casacore::Complex(1.,w);
            }
            casacore::LatticeFFT::cfft2d(scratchLattice, false); // iFFT to the image domain
            pcf_orig = casacore::Array<float>(shape);
            casacore::ArrayLattice<float> pcfLattice(pcf_orig);
            pcfLattice.copyData(casacore::LatticeExpr<float> ( real(scratchLattice) ));
            pcf_orig *= normFactor;
           
            if (printDetails) {
                std::cout << std::endl;
                std::cout << "Original shape: "<<dirty_orig.shape() << std::endl;
            }

        }

        void testReversibility()
        {

            // step through several different pixel ratios to make sure there are no issues commuting
            for (int n_down = n/2; n_down < n/2 + 8; ++n_down) {
                // set the oversampling factor
                const double extraOversampleFactor = double(n) / double(n_down);
                // make a copy of the original
                casacore::Array<float> dirty = dirty_orig.copy();
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
                CPPUNIT_ASSERT(dirty.shape().isEqual(dirty_orig.shape()));
                if (printDetails) {
                    std::cout << " - error_peak = "<<std::abs(dirty(casacore::IPosition(2,n/2,n/2))-1.) << std::endl;
                    std::cout << " - error_max  = "<<casacore::max(casacore::abs(dirty-dirty_orig)) << std::endl;
                }
                CPPUNIT_ASSERT(std::abs(dirty(casacore::IPosition(2,n/2,n/2))-1.)<1e-6);
                CPPUNIT_ASSERT(casacore::max(casacore::abs(dirty-dirty_orig))<1e-6);
            }

        }
        void testPreconditioning()
        {

            const float robustness = -2.;
            const int n_down = n/2 + 2;
            const double extraOversampleFactor = double(n) / double(n_down);
           
            casacore::Array<float> dirty1 = dirty_orig.copy();
            casacore::Array<float> psf1 = psf_orig.copy();
            casacore::Array<float> pcf1 = pcf_orig.copy();
                     
            casacore::LogIO logger;

            WienerPreconditioner wp1(robustness);
            wp1.doPreconditioning(psf1,dirty1,pcf1);
            if (printDetails) {
                std::cout << "After preconditioning (full resolution)" << std::endl;
                std::cout << " - shape = "<<dirty1.shape() << std::endl;
                std::cout << " - mid = "<<dirty1(casacore::IPosition(2,n/2,n/2)) << std::endl;
                std::cout << " - uni - nat max error = "<<casacore::max(casacore::abs(dirty1-dirty_orig)) << std::endl;
            }
           
            casacore::Array<float> dirty2 = dirty_orig.copy();
            casacore::Array<float> psf2 = psf_orig.copy();
            casacore::Array<float> pcf2 = pcf_orig.copy();
            SynthesisParamsHelper::downsample(dirty2,extraOversampleFactor);
            SynthesisParamsHelper::downsample(psf2,extraOversampleFactor);
            SynthesisParamsHelper::downsample(pcf2,extraOversampleFactor);
            if (printDetails) {
                std::cout << "Nyquist preconditioning:" << std::endl;
                std::cout << " - shape = "<<dirty2.shape() << std::endl;
            }

            // use the createPreconditioner function to set itsUseCachedPcf = false (since size will have changed)
            LOFAR::ParameterSet parset;
            parset.add("robustness", utility::toString(robustness));
            boost::shared_ptr<WienerPreconditioner> wp2 = WienerPreconditioner::createPreconditioner(parset, false);
            wp2->doPreconditioning(psf2,dirty2,pcf2);

            SynthesisParamsHelper::oversample(dirty2,extraOversampleFactor,false);
            SynthesisParamsHelper::oversample(psf2,extraOversampleFactor,false);
            if (printDetails) {
                std::cout << "After preconditioning and oversampling" << std::endl;
                std::cout << " - shape = "<<dirty2.shape() << std::endl;
                std::cout << " - mid = "<<dirty2(casacore::IPosition(2,n/2,n/2)) << std::endl;
                std::cout << " - uni - nat max error = "<<casacore::max(casacore::abs(dirty2-dirty_orig)) << std::endl;
                std::cout << " - nyq - uni max error = "<<casacore::max(casacore::abs(dirty2-dirty1)) << std::endl;
            }
            CPPUNIT_ASSERT(casacore::max(casacore::abs(dirty2-dirty1))<1e-6);
            CPPUNIT_ASSERT(casacore::max(casacore::abs(psf2-psf1))<1e-6);

        }
    };

  } // namespace synthesis
} // namespace askap

#endif // #ifndef NYQUISTGRIDDING_TESTS_H
