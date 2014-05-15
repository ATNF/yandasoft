/// @file
///
/// Unit test for the deconvolution control class
///
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>

#include <deconvolution/DeconvolverState.h>
#include <deconvolution/DeconvolverState.h>
#include <cppunit/extensions/HelperMacros.h>

#include <casa/BasicSL/Complex.h>

#include <boost/shared_ptr.hpp>

using namespace casa;

namespace askap {

namespace synthesis {

class DeconvolverStateTest : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE(DeconvolverStateTest);
   CPPUNIT_TEST(testSetGet);
   CPPUNIT_TEST(testReset);
   CPPUNIT_TEST_SUITE_END();
public:
   
  void setUp() {
    itsDS=boost::shared_ptr<DeconvolverState<Float> >(new DeconvolverState<Float>());
  }
   void testSetGet() {
     {
       itsDS->setCurrentIter(100);
       CPPUNIT_ASSERT(itsDS->currentIter()==100);
       itsDS->setPeakResidual(1.4);
       CPPUNIT_ASSERT(itsDS->peakResidual()==Float(1.4));
       itsDS->setTotalFlux(0.3);
       CPPUNIT_ASSERT(itsDS->totalFlux()==Float(0.3));
       itsDS->setObjectiveFunction(10.3);
       CPPUNIT_ASSERT(itsDS->objectiveFunction()==Float(10.3));
       itsDS->setCurrentIter(0);
       CPPUNIT_ASSERT(itsDS->currentIter()==0);
       itsDS->incIter();
       CPPUNIT_ASSERT(itsDS->currentIter()==1);
     }
   }
   void testReset() {
     {
       itsDS->setCurrentIter(100);
       CPPUNIT_ASSERT(itsDS->currentIter()==100);
       itsDS->reset();
       CPPUNIT_ASSERT(itsDS->currentIter()==0);
     }
   }
  void tearDown() {
    itsDS->reset();
  }
   
private:
   /// @brief DeconvolutionState class
  boost::shared_ptr<DeconvolverState<Float> > itsDS;
};
    
} // namespace synthesis

} // namespace askap

