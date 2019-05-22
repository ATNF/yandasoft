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

#include "IVisGridder.h"
#include <askap/AskapError.h>

namespace askap
{
  namespace synthesis
  {

    IVisGridder::IVisGridder()
    {
    }

    IVisGridder::~IVisGridder()
    {
    }
    
    /// @brief static method to create gridder
	/// @details Each gridder should have a static factory method, which is
	/// able to create a particular type of the gridder and initialise it with
	/// the parameters taken form the given parset. It is assumed that the 
	/// method receives a subset of parameters where the gridder name is already
	/// taken out. 
	/// @note This method just throws an exception in this basic interface. It is 
	/// added to ensure that all derived classes have this method defined. We have 
	/// to use a static method as opposed to pure virtual function because we plan to use it
	/// to create a brand new instance of the gridder (and hence no object would
	/// exist at that stage)			 
	IVisGridder::ShPtr IVisGridder::createGridder(const LOFAR::ParameterSet&) 
    {
       ASKAPTHROW(AskapError, "createGridder is supposed to be defined for every derived gridder, "
                              "IVisGridder::createGridder should never be called");           
       return IVisGridder::ShPtr();                                            
    }
  }
}
