/// @file 
///
/// @brief Experiments with the measurement set
/// @details This is not a general purpose program. 
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

/// casa includes
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/TableColumn.h>
#include <casa/Arrays/ArrayMath.h>

// own includes
#include <askap/AskapError.h>

// std
#include <stdexcept>
#include <list>
#include <string>
#include <algorithm>

using namespace askap;

template<typename T>
void copyArrayColumn(casa::Table &tab, const std::string &name, const casa::uInt nOrigRows)
{
  const casa::TableDesc td = tab.actualTableDesc();
  ASKAPASSERT(tab.nrow() == 2*nOrigRows);
  const casa::ColumnDesc cd = td[name];
  ASKAPASSERT(cd.isArray());
  casa::ArrayColumn<T> col(tab, name);
  casa::Array<T> buf;
  for (casa::uInt row=0; row<nOrigRows; ++row) {       
       col.get(row,buf,casa::True);
       col.put(nOrigRows + row, buf);
  }  
}

bool has(const std::list<std::string> & lst, const std::string &name) {
  return std::find(lst.begin(), lst.end(), name) != lst.end();
}

int main(int argc, const char **argv) {
try {
  if (argc!=2) {
      throw std::runtime_error("usage mstabtest ms");
  }
  casa::Table mytab(argv[1],casa::Table::Update);
  casa::TableDesc td = mytab.actualTableDesc();
  const casa::uInt nOrigRows = mytab.nrow();
  std::cout<<"Table has "<<nOrigRows<<" rows and "<<td.ncolumn()<<" columns"<<std::endl;
  casa::Vector<casa::String> colList = td.columnNames();
  mytab.addRow(nOrigRows);
  {
    std::list<std::string> floatColumns;
    floatColumns.push_back("WEIGHT");
    floatColumns.push_back("SIGMA");

    std::list<std::string> boolColumns;
    boolColumns.push_back("FLAG");
    
    std::list<std::string> skipColumns;
    skipColumns.push_back("UVW");
    skipColumns.push_back("DATA");
    skipColumns.push_back("FLAG_CATEGORY");
    
  
    casa::TableColumn outcol;
    casa::ROTableColumn incol;
    for (casa::uInt col = 0; col<colList.nelements(); ++col) {
         const std::string colName = colList[col];
         std::cout<<"col = "<<col<<": "<<colName<<std::endl;       
         if (has(skipColumns, colName)) {
             std::cout<<"skipping..."<<std::endl;
         } else if (has(floatColumns,colName)) {
             std::cout<<"copy as array of floats..."<<std::endl;             
             copyArrayColumn<casa::Float>(mytab,colName, nOrigRows);
         } else if (has(boolColumns,colName)) {
             std::cout<<"copy as array of bool..."<<std::endl;             
             copyArrayColumn<casa::Bool>(mytab,colName, nOrigRows);
         } else {
            std::cout<<"default processing..."<<std::endl;
            incol.attach(mytab,col);
            outcol.attach(mytab,col);
            for (casa::uInt row = 0; row<nOrigRows; ++row) {
                outcol.put(nOrigRows + row, incol, row);
            }       
         }
    }
  }
  std::cout<<"processing UVW and DATA columns"<<std::endl;
  casa::ArrayColumn<double> uvwCol(mytab,"UVW");
  casa::ArrayColumn<casa::Complex> dataCol(mytab,"DATA");
  casa::Array<double> dBuf;
  casa::Array<casa::Complex> cBuf;
  for (casa::uInt row=0; row<nOrigRows; ++row) {       
       uvwCol.get(row,dBuf,casa::True);
       uvwCol.put(nOrigRows + row, -1.*dBuf);
       dataCol.get(row,cBuf,casa::True);
       dataCol.put(nOrigRows + row, casa::conj(cBuf));   
  }
}  
catch (const std::exception &ex) {
  std::cerr<<ex.what()<<std::endl;
  return -1;
}

return 0;
}

