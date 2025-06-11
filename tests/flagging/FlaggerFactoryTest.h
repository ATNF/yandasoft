/// @file FlaggerFactoryTest.cc
///
/// @copyright (c) 2013 CSIRO
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>

// CPPUnit includes
#include <cppunit/extensions/HelperMacros.h>

// Support classes
#include "casacore/ms/MeasurementSets/MeasurementSet.h"
#include "askap/askap/AskapError.h"

// Classes to test
#include "askap/flagging/FlaggerFactory.h"

using namespace casacore;

namespace askap {
namespace synthesis {

class FlaggerFactoryTest : public CppUnit::TestFixture {
        CPPUNIT_TEST_SUITE(FlaggerFactoryTest);
        CPPUNIT_TEST(testBuildAll);
        CPPUNIT_TEST(testBuildStokesV);
        CPPUNIT_TEST(testBuildElevation);
        CPPUNIT_TEST(testBuildAmplitude);
        CPPUNIT_TEST(testBuildSelection);
        CPPUNIT_TEST_SUITE_END();

    public:
        void setUp() {
        };

        void tearDown() {
        }

        /// Tests building all flaggers
        void testBuildAll() {
            LOFAR::ParameterSet parset;
            parset.add("dataset", "target.ms");

            parset.add("stokesv_flagger.enable", "true");
            parset.add("stokesv_flagger.threshold", "5.0");

            parset.add("elevation_flagger.enable", "true");
            parset.add("elevation_flagger.high", "60.0");
            parset.add("elevation_flagger.low", "20.0");

            parset.add("amplitude_flagger.enable", "true");
            parset.add("amplitude_flagger.high", "1500.0");
            parset.add("amplitude_flagger.low", "1e-15");
            parset.add("amplitude_flagger.stokes", "[XX, YY]");

            parset.add("selection_flagger.rules", "[rule1, rule2]");
            parset.add("selection_flagger.rule1.spw", "35:0~15;288~303");
            parset.add("selection_flagger.rule2.antenna", "Pad01");

            // NOTE: If the FlaggerFactory ever needs a "real" measurement set, then this
            // is going to break. The below constructor does not actually initialise a valid
            // measurement set.
            std::string ms;
            CPPUNIT_ASSERT_EQUAL(5ul, FlaggerFactory::build(parset, ms).size());
        }

        void testBuildStokesV() {
            // The Stokes-V flagger will work with defaults
            LOFAR::ParameterSet parset;
            parset.add("stokesv_flagger.enable", "true");
            std::string ms;
            CPPUNIT_ASSERT_EQUAL(1ul, FlaggerFactory::build(parset, ms).size());
        }

        void testBuildElevation() {
            // The elevation flagger will work with defaults
            LOFAR::ParameterSet parset;
            parset.add("elevation_flagger.enable", "true");
            std::string ms;
            CPPUNIT_ASSERT_EQUAL(1ul, FlaggerFactory::build(parset, ms).size());
        }

        void testBuildAmplitude() {
            // The amplitude based flagger doesn't have defaults, so throws
            // an exception if enabled without being configured.
            LOFAR::ParameterSet parset;
            parset.add("amplitude_flagger.enable", "true");
            std::string ms;
            CPPUNIT_ASSERT_THROW(FlaggerFactory::build(parset, ms), askap::AskapError);

            // Either high or low are necessary, not both
            LOFAR::ParameterSet parsetHigh = parset;
            parsetHigh.add("amplitude_flagger.high", "1500.0");
            CPPUNIT_ASSERT_EQUAL(1ul, FlaggerFactory::build(parsetHigh, ms).size());

            LOFAR::ParameterSet parsetLow = parset;
            parsetLow.add("amplitude_flagger.low", "1e-15");
            CPPUNIT_ASSERT_EQUAL(1ul, FlaggerFactory::build(parsetLow, ms).size());
        }

        void testBuildSelection() {
            // If a rule is listed, it is expected some criteria will also be specified
            LOFAR::ParameterSet parset;
            parset.add("selection_flagger.rules", "[rule1, rule2]");
            std::string ms;
            CPPUNIT_ASSERT_THROW(FlaggerFactory::build(parset, ms), askap::AskapError);

            // Configure one rule (of two specified), still expect an exception
            parset.add("selection_flagger.rule1.spw", "35:0~15;288~303");
            CPPUNIT_ASSERT_THROW(FlaggerFactory::build(parset, ms), askap::AskapError);
        }

    private:
};

}   // End namespace synthesis
}   // End namespace askap
