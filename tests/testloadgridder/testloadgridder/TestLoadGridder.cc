/// @copyright (c) 2009 CSIRO
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

#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".testloadgridder");

#include <testloadgridder/TestLoadGridder.h>
#include <gridding/VisGridderFactory.h>

namespace askap
{
  namespace synthesis
  {

    TestLoadGridder::TestLoadGridder()
    {
    }

    TestLoadGridder::~TestLoadGridder()
    {
    }

    /// Clone a copy of this Gridder
    IVisGridder::ShPtr TestLoadGridder::clone() 
    {
      return IVisGridder::ShPtr(new TestLoadGridder(*this));
    }

    IVisGridder::ShPtr TestLoadGridder::makeGridder
    (const LOFAR::ParameterSet&)
    {
      std::cout << "TestLoadGridder::makeGridder" << std::endl;
      ASKAPLOG_INFO_STR (logger, "TestLoadGridder::makeGridder");
      return IVisGridder::ShPtr(new TestLoadGridder());
    }

    const std::string& TestLoadGridder::gridderName()
    {
      static std::string name("TestLoadGridder");
      return name;
    }

    void TestLoadGridder::registerGridder()
    {
      VisGridderFactory::registerGridder (gridderName(), &makeGridder);
    }

    void TestLoadGridder::initIndices(const accessors::IConstDataAccessor&) 
    {
    }

    void TestLoadGridder::initConvolutionFunction(const accessors::IConstDataAccessor&)
    {
      itsSupport=0;
      itsOverSample=1;
      const int cSize=2*(itsSupport+1)*itsOverSample+1; // 3
      const int cCenter=(cSize-1)/2; // 1
      itsConvFunc.resize(1);
      itsConvFunc[0].resize(cSize, cSize); // 3, 3, 1
      itsConvFunc[0].set(0.0);
      itsConvFunc[0](cCenter,cCenter)=1.0; // 1,1,0 = 1
    }
    
    void TestLoadGridder::correctConvolution(casa::Array<double>& image)
    {
    }

  }
}
