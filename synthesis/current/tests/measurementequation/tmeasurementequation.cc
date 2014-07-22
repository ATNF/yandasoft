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

// ASKAPsoft includes
#include <AskapTestRunner.h>

// just to avoid template compilation which will not work without logging
#define A_PROJECT_GRIDDER_BASE_TCC

// Test includes
#include <ComponentEquationTest.h>
#include <Calibrator1934Test.h>
#include <VectorOperationsTest.h>
#include <ImageDFTEquationTest.h>
#include <ImageFFTEquationTest.h>
#include <CalibrationMETest.h>
#include <PreconditionerTests.h>
#include <SynthesisParamsHelperTest.h>
#include <ImageParamsHelperTest.h>
#include <GaussianNoiseMETest.h>
#include <PolLeakageTest.h>
#include <PreAvgCalBufferTest.h>
#include <RestoringBeamHelperTest.h>
#include <VisMetaDataStatsTest.h>
#include <ObservationDescriptionTest.h>

int main( int argc, char **argv)
{
    askapdev::testutils::AskapTestRunner runner(argv[0]);
    
    runner.addTest(askap::synthesis::VectorOperationsTest::suite());
    runner.addTest(askap::synthesis::ComponentEquationTest::suite()); 
    runner.addTest(askap::synthesis::Calibrator1934Test::suite()); 
    runner.addTest(askap::synthesis::PreAvgCalBufferTest::suite()); 
    runner.addTest(askap::synthesis::CalibrationMETest::suite()); 
    //runner.addTest(askap::synthesis::ImageDFTEquationTest::suite());
    runner.addTest(askap::synthesis::ImageFFTEquationTest::suite());
    runner.addTest(askap::synthesis::PreconditionerTests::suite());
    runner.addTest(askap::synthesis::SynthesisParamsHelperTest::suite());
    runner.addTest(askap::synthesis::ImageParamsHelperTest::suite());
    runner.addTest(askap::synthesis::GaussianNoiseMETest::suite());
    runner.addTest(askap::synthesis::PolLeakageTest::suite()); 
    runner.addTest(askap::synthesis::RestoringBeamHelperTest::suite());
    runner.addTest(askap::synthesis::VisMetaDataStatsTest::suite());
    runner.addTest(askap::synthesis::ObservationDescriptionTest::suite());
    
    const bool wasSucessful = runner.run();

    return wasSucessful ? 0 : 1;
}

