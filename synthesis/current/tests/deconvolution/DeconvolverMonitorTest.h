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
#include <deconvolution/DeconvolverMonitor.h>

#include <cppunit/extensions/HelperMacros.h>

#include <boost/shared_ptr.hpp>

using namespace casa;

namespace askap {

namespace synthesis {

class DeconvolverMonitorTest : public CppUnit::TestFixture
{
   CPPUNIT_TEST_SUITE(DeconvolverMonitorTest);
   CPPUNIT_TEST(testMonitor);
   CPPUNIT_TEST_SUITE_END();
public:
   
  void setUp() {
    itsDM=boost::shared_ptr<DeconvolverMonitor<Float> >(new DeconvolverMonitor<Float>());
  }
   void testMonitor() {
     {
       DeconvolverState<Float> ds;
       for (int iter=0;iter<10;iter++) {
         ds.setCurrentIter(iter);
         ds.setObjectiveFunction(1.0/Float(iter+1));
         itsDM->monitor(ds);
       }
     }
   }
  void tearDown() {
    itsDM.reset();
  }
   
private:
   /// @brief DeconvolutionMonitor class
  boost::shared_ptr<DeconvolverMonitor<Float> > itsDM;
};
    
} // namespace synthesis

} // namespace askap

