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
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>

ASKAP_LOGGER(logger, ".deconvolution");

// Test includes
#include <EntropyTest.h>
#include <BasisFunctionTest.h>
#include <DeconvolverBaseTest.h>
#include <DeconvolverFistaTest.h>
#include <DeconvolverHogbomTest.h>
#include <DeconvolverMultiTermBasisFunctionTest.h>
#include <DeconvolverControlTest.h>
#include <DeconvolverMonitorTest.h>
#include <DeconvolverStateTest.h>

int main(int argc, char *argv[])
{
    askapdev::testutils::AskapTestRunner runner(argv[0]);

    runner.addTest( askap::synthesis::DeconvolverBaseTest::suite());
    runner.addTest( askap::synthesis::DeconvolverFistaTest::suite());
    runner.addTest( askap::synthesis::DeconvolverHogbomTest::suite());
    runner.addTest( askap::synthesis::DeconvolverMultiTermBasisFunctionTest::suite());
    runner.addTest( askap::synthesis::DeconvolverControlTest::suite());
    runner.addTest( askap::synthesis::DeconvolverMonitorTest::suite());
    runner.addTest( askap::synthesis::DeconvolverStateTest::suite());
    runner.addTest( askap::synthesis::EntropyTest::suite());
    runner.addTest( askap::synthesis::BasisFunctionTest::suite());
    bool wasSuccessful = runner.run();

    return wasSuccessful ? 0 : 1;
}
