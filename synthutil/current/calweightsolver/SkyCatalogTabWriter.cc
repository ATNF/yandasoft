/// @file SkyCatalogTabWriter.cc
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

#include "SkyCatalogTabWriter.h"
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/casa/Containers/RecordFieldId.h>
#include <casacore/casa/Containers/RecordInterface.h>


using namespace casa;
#include <sstream>
#include <string>
using namespace std;

//  An auxiliary class to write annotation table from the supplied data
SkyCatalogTabWriter::SkyCatalogTabWriter(const casa::String &fname)
            throw(casa::AipsError) : ncomp(0)
{
  TableDesc td("Skycatalog",TableDesc::Scratch);
  td.addColumn(ScalarColumnDesc<String>("Type"));
  td.addColumn(ScalarColumnDesc<Double>("Long"));
  td.addColumn(ScalarColumnDesc<Double>("Lat"));
  td.addColumn(ScalarColumnDesc<String>("Annotation"));
  td.addColumn(ScalarColumnDesc<Float>("Flux"));  
  td.rwColumnDesc("Long").rwKeywordSet().define("UNIT","deg");
  td.rwColumnDesc("Lat").rwKeywordSet().define("UNIT","deg");
  td.rwColumnDesc("Flux").rwKeywordSet().define("UNIT","Jy");
  SetupNewTable newtab(fname,td,Table::New);
  tab=Table(newtab,Table::Plain,0);
  tab.tableInfo().setType("Skycatalog");
}

void SkyCatalogTabWriter::addComponent(casa::Double lng, casa::Double lat,
                       casa::Double flux, const casa::String &type,
		       const casa::String &label) throw(casa::AipsError)
{
  tab.addRow();
  ScalarColumn<String>(tab,"Type").put(ncomp,type);
  ScalarColumn<Double>(tab,"Long").put(ncomp,lng);
  ScalarColumn<Double>(tab,"Lat").put(ncomp,lat);
  String annotation(label);
  if (annotation=="") {          
     ostringstream os;
     os<<ncomp+1;
     if (!os) throw AipsError("SkyCatalogTabWriter - Unable to convert sequence number");          
     annotation=os.str();     
  }
  ScalarColumn<String>(tab,"Annotation").put(ncomp,annotation);
  ScalarColumn<Float>(tab,"Flux").put(ncomp,flux);
  ++ncomp;
  //if (ncomp==1) addComponent(lng,lat,0.9*flux,type,label);
}

//
///////////////////////////////////////////////////////////////////////////////
