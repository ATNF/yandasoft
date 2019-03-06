/// @file
///
/// @brief Abstract component
/// @details
/// IComponent is a base class for components working with ComponentEquation
/// examples of components include, e.g. Gaussian or point sources.
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>
/// 

#include <measurementequation/IComponent.h>

/// virtual destructor to keep the compiler happy
askap::synthesis::IComponent::~IComponent()  {} 

/// @brief convert StokesTypes into an index 0..3
/// @details It is decided that all components have to be defined in
/// terms of IQUV stokes parameters. It is not prohibited that the 
/// constructors of actual components accept other stokes parameters like
/// XX, etc. However, in the latter case, these parameters should be converted
/// to IQUV at the time of the object construction. Most likely actual 
/// components will hold an array of fluxes for each stokes parameter. Therefore,
/// it is necessary to convert quickly from StokesTypes to the index.
/// This method gives a mapping of I to 0, Q to 1, U to 2 and V to 3. For
/// other values an exception is thrown.
/// @param[in] pol required polarization
/// @return an index (I: 0, Q: 1, U: 2 and V: 3)
size_t askap::synthesis::IComponent::stokesIndex(casa::Stokes::StokesTypes pol) 
                    throw(AskapError)
{
  ASKAPDEBUGASSERT(pol!=casa::Stokes::Undefined && pol<=casa::Stokes::V);
  return static_cast<size_t>(pol-casa::Stokes::I);
}
