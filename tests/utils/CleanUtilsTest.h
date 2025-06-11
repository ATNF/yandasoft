/// @file CleanUtilsTest.cc
///
/// @copyright (c) 2024 CSIRO
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
/// @author Mark Wieringa <mark.wieringa@csiro.au>

// CPPUnit includes
#include <cppunit/extensions/HelperMacros.h>

// Support classes
#include "askap/askap/AskapError.h"
#include <askap/scimath/fitting/Params.h>
#include <boost/optional.hpp>

// Classes to test
#include "askap/utils/CleanUtils.h"

using namespace casacore;
using namespace askap::scimath;

namespace askap {
namespace utils {

class CleanUtilsTest : public CppUnit::TestFixture {
        CPPUNIT_TEST_SUITE(CleanUtilsTest);
        CPPUNIT_TEST(testOverlapMaskEmpty);
        CPPUNIT_TEST(testOverlapMask);
        CPPUNIT_TEST_SUITE_END();

    public:
        void setUp() {
        };

        void tearDown() {
        }

        /// Tests building all flaggers
        void testOverlapMaskEmpty() {
            // test with empty inputs
            scimath::Params ip(true);
            std::map<std::string, int> taylorMap;
            boost::optional<float> extraOversamplingFactor;
            CPPUNIT_ASSERT_EQUAL(0ul, overlapMask(ip,taylorMap,extraOversamplingFactor).size());

            // with just taylorMap defined - expect exception from Params
            taylorMap["image.first"] = 1;
            CPPUNIT_ASSERT_THROW(overlapMask(ip,taylorMap,extraOversamplingFactor), askap::AskapError);
        }

        void testOverlapMask() {
            scimath::Params ip(true);
            std::map<std::string, int> taylorMap;
            boost::optional<float> extraOversamplingFactor;
            // make large and small image with varying offset & overlap
            int npix1 = 1024, npix2 = 256;

            int offset = 1;
            imtype result = npix1 * npix1 - npix2 * npix2;
            makeImages(ip,taylorMap,offset,npix1,npix2);
            CPPUNIT_ASSERT_EQUAL(result, sum(overlapMask(ip,taylorMap,extraOversamplingFactor)));

            offset = 129;
            result = npix1 * npix1 - npix2 * npix2;
            makeImages(ip,taylorMap,offset,npix1,npix2);
            CPPUNIT_ASSERT_EQUAL(result, sum(overlapMask(ip,taylorMap,extraOversamplingFactor)));

            offset = 385;
            result = npix1 * npix1 - 65280; // due to geometry
            makeImages(ip,taylorMap,offset,npix1,npix2);
            CPPUNIT_ASSERT_EQUAL(result, sum(overlapMask(ip,taylorMap,extraOversamplingFactor)));

            offset = 513;
            result = npix1 * npix1 - 32512; // due to only half overlapping
            makeImages(ip,taylorMap,offset,npix1,npix2);
            CPPUNIT_ASSERT_EQUAL(result, sum(overlapMask(ip,taylorMap,extraOversamplingFactor)));

            offset = 641;
            result = 0; // due to no overlap
            makeImages(ip,taylorMap,offset,npix1,npix2);
            CPPUNIT_ASSERT_EQUAL(result, sum(overlapMask(ip,taylorMap,extraOversamplingFactor)));

        }

        // make two images with given offset
        void makeImages(scimath::Params& ip, std::map<std::string, int>& taylorMap, int offset,
            int npix1, int npix2)
        {
            ip.reset();
            taylorMap.clear();

            // make some images
            Axes imageAxes;
            const double cell = 8.0*C::arcsec;
            casacore::Matrix<double> xform(2,2,0.);
            xform.diagonal().set(1.);

            imageAxes.addDirectionAxis(casacore::DirectionCoordinate(casacore::MDirection::J2000,
                         casacore::Projection(casacore::Projection::SIN), 0.,0.,cell,cell,xform,npix1/2.,npix1/2.));
            imageAxes.addStokesAxis(casacore::Vector<casacore::Stokes::StokesTypes>(1,casacore::Stokes::I));
            imageAxes.add("FREQUENCY",1.4e9,1.4e9);

            casacore::Array<float> imagePixels1(casacore::IPosition(4, npix1, npix1, 1, 1),0.f);
            ip.add("image.i.main", imagePixels1, imageAxes);

            double off = offset*cell;
            Axes imageAxes2;

            imageAxes2.addDirectionAxis(casacore::DirectionCoordinate(casacore::MDirection::J2000,
                         casacore::Projection(casacore::Projection::SIN), 0.,off,cell,cell,xform,npix2/2.,npix2/2.));
            imageAxes2.addStokesAxis(casacore::Vector<casacore::Stokes::StokesTypes>(1,casacore::Stokes::I));
            imageAxes2.add("FREQUENCY",1.4e9,1.4e9);

            casacore::Array<float> imagePixels2(casacore::IPosition(4, npix2, npix2, 1, 1),0.f);
            ip.add("image.i.offset", imagePixels2, imageAxes2);

            // fill map
            taylorMap["image.i.main"] = 1;
            taylorMap["image.i.offset"] = 1;
        }

    private:
};

}   // End namespace synthesis
}   // End namespace askap
