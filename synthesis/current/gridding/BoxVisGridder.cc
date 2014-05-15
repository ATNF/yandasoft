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

#include <gridding/BoxVisGridder.h>

namespace askap
{
  namespace synthesis
  {

    BoxVisGridder::BoxVisGridder()
    {
    }

    BoxVisGridder::~BoxVisGridder()
    {
    }

		/// Clone a copy of this Gridder
		IVisGridder::ShPtr BoxVisGridder::clone() 
		{
			return IVisGridder::ShPtr(new BoxVisGridder(*this));
		}
		
    /// @brief static method to create gridder
    /// @details Each gridder should have a static factory method, which is
    /// able to create a particular type of the gridder and initialise it with
    /// the parameters taken form the given parset. It is assumed that the 
  	/// method receives a subset of parameters where the gridder name is already
    /// taken out. 
    /// @return a shared pointer to the gridder instance					 
    IVisGridder::ShPtr BoxVisGridder::createGridder(const LOFAR::ParameterSet&)
    {
      return IVisGridder::ShPtr(new BoxVisGridder());  
    }
		

		void BoxVisGridder::initIndices(const accessors::IConstDataAccessor&) 
		{
		}

    void BoxVisGridder::initConvolutionFunction(const accessors::IConstDataAccessor&)
    {
      itsSupport=1;
      itsOverSample=1;
      const int cSize=2*itsSupport+1; // 3
      const int cCenter=(cSize-1)/2; // 1
      itsConvFunc.resize(1);
      itsConvFunc[0].resize(cSize, cSize); // 3, 3
      itsConvFunc[0].set(0.0);
      itsConvFunc[0](cCenter,cCenter)=1.0; // 1,1 = 1
    }
    
		void BoxVisGridder::correctConvolution(casa::Array<double>& /*image*/)
		{
		}

  }
}
