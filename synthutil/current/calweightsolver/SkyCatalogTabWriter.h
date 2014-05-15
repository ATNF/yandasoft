///////////////////////////////////////////////////////////////////////////////
//
//  An auxiliary class to write annotation table from the supplied data
//

#include <casa/aips.h>
#include <casa/Exceptions/Error.h>
#include <tables/Tables/Table.h>
#include <casa/BasicSL/String.h>

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

//
///////////////////////////////////////////////////////////////////////////////
