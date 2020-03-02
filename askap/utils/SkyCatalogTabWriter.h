/// @file SkyCatalogTabWriter.h
///
/// @copyright (c) 2014 CSIRO
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

#include <casacore/casa/aips.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/BasicSL/String.h>

/// An auxiliary class to write annotation table from the supplied data
class SkyCatalogTabWriter {
  casa::uInt ncomp;  // current number of components
                     // (to be able to add a correct name)
  casa::Table tab;
public:
     SkyCatalogTabWriter(const casa::String &fname) throw(casa::AipsError);
     // lng, lat - longitude and latitude in degrees,
     // flux - flux density in Jy
     // type - e.g. J2000
     // label - annotation (default is a sequence number)          
     void addComponent(casa::Double lng, casa::Double lat,
                       casa::Double flux = 1., const casa::String &type = "J2000",
		       const casa::String &label = "") throw(casa::AipsError);
};
