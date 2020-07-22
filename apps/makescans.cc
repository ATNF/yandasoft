/// @file makescans.cc
///
/// test code to help with various non-standard experiments where synchronisation
/// between some parameter change and scans is not possible by other means than time
/// This tool reads an external 3 column file with scan number start and stop MJDs and
/// rewrites the MS's scan column to match the time intervals. All data outside the intervals
/// are flagged. It is checked that the input MS has scan = 0 for all rows
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


#include <askap/AskapLogging.h>
#include <askapparallel/AskapParallel.h>

ASKAP_LOGGER(logger, ".makescans");

#include <askap/AskapError.h>
#include <dataaccess/ParsetInterface.h>
#include <askap/AskapUtil.h>


// casa
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/ArrayMath.h>

// LOFAR
#include <Common/ParameterSet.h>

// std
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <string>

// boost
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>

using namespace askap;

void loadScans(std::vector<boost::tuple<casa::Int, double, double> > &buf, const std::string &fname)
{
   ASKAPDEBUGASSERT(buf.size() == 0);
   std::ifstream is(fname.c_str());
   while (is) {
       int scan = -1; 
       double startMJD = 0., endMJD = 0.;
       is >> scan >> startMJD >> endMJD;
       if (is) {
           ASKAPCHECK(scan >= 0, "Expect non-negative scan number in "<<fname);
           buf.push_back(boost::tuple<casa::Int, double, double>(scan, startMJD, endMJD));
       }
   }
   std::cerr<<"Found "<<buf.size()<<" records in "<<fname<<std::endl;
}

int findScan(const std::vector<boost::tuple<casa::Int, double, double> > &buf, double mjd) 
{
   // we could've optimise it assuming the list is always sorted
   for (std::vector<boost::tuple<casa::Int, double, double> >::const_iterator ci = buf.begin(); ci != buf.end(); ++ci) {
        if (mjd >= ci->get<1>() && mjd <= ci->get<2>()) { 
            return ci->get<0>();
        }
   }
   return -1;
}

void processOneMS(const std::vector<boost::tuple<casa::Int, double, double> > &buf, const std::string &fname)
{
  std::cerr<<"Updating scan information in "<<fname<<std::endl;
  casa::Table ms(fname, casa::Table::Update);
  casa::ArrayColumn<casa::Bool> flagCol(ms, "FLAG");
  casa::ScalarColumn<casa::Int> scanCol(ms, "SCAN_NUMBER");
  casa::ScalarColumn<casa::Double> timeCol(ms, "TIME_CENTROID");

  std::cerr<<"Total number of rows in the measurement set: "<<ms.nrow()<<std::endl;

  casa::uInt rowsOutside = 0;
  std::set<casa::Int> scansWithData;
  casa::Int prevScan = 0;
 
  for (casa::uInt row = 0; row<ms.nrow(); ++row) {
       //casa::Array<casa::Bool> flagBuf;
       //flagCol.get(row,flagBuf);
       casa::Int scan = scanCol.get(row);
       ASKAPCHECK(scan == 0, "Expect scan=0 in the input MS before modification");
       // quick and dirty way to extract mjd, we could've done the same more neatly through the measures column
       scan = findScan(buf, timeCol.get(row)/86400.);
       if (scan < 0) {
           ++rowsOutside;
           // flag it
           casa::Array<casa::Bool> flagBuf;
           flagCol.get(row,flagBuf);
           flagBuf.set(casa::True);
           flagCol.put(row, flagBuf);
           scanCol.put(row, prevScan);
       } else {
          scansWithData.insert(scan);
          scanCol.put(row, scan);
          prevScan = scan;
          //std::cout<<"row = "<<row<<" scan="<<scan<<" "<<std::setprecision(15)<<timeCol.get(row)/86400.<<std::endl;
       }
  }
  std::cerr<<"Rows outside defined scans: "<<rowsOutside<<" out of "<<ms.nrow()<<" total number of rows"<<std::endl;
  std::cerr<<"Matched "<<scansWithData.size()<<" scans in the MS out of "<<buf.size()<<" scans defined"<<std::endl;
}

int main(int argc, const char **argv) {
  // This class must have scope outside the main try/catch block
  ///askapparallel::AskapParallel comms(argc, argv);

  try {
     if (argc!=3) {
         std::cerr<<"Usage: "<<argv[0]<<" MS scans.txt"<<std::endl;
	 return -2;
     }
     

     casa::Timer timer;

     timer.mark();
     std::cerr<<"Initialization: "<<timer.real()<<std::endl;
     timer.mark();
     std::vector<boost::tuple<casa::Int, double, double> > scans;
     loadScans(scans, argv[2]);
     processOneMS(scans, argv[1]);
     std::cerr<<"Job: "<<timer.real()<<std::endl;
     
  }
  catch(const AskapError &ce) {
     std::cerr<<"AskapError has been caught. "<<ce.what()<<std::endl;
     return -1;
  }
  catch(const std::exception &ex) {
     std::cerr<<"std::exception has been caught. "<<ex.what()<<std::endl;
     return -1;
  }
  catch(...) {
     std::cerr<<"An unexpected exception has been caught"<<std::endl;
     return -1;
  }
  return 0;
}
