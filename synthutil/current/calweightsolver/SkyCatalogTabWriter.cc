///////////////////////////////////////////////////////////////////////////////
//
//  An auxiliary class to write annotation table from the supplied data
//

#include "SkyCatalogTabWriter.h"
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/ScalarColumn.h>
#include <casa/Containers/RecordFieldId.h>
#include <casa/Containers/RecordInterface.h>


using namespace casa;
#include <sstream>
#include <string>
using namespace std;

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
