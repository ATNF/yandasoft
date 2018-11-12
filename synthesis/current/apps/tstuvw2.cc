/// @file tstuvw2.cc
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

#include <dataaccess/TableDataSource.h>
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, "");

#include <askap/AskapError.h>
#include <dataaccess/SharedIter.h>
#include <dataaccess/ParsetInterface.h>

#include <dataaccess/TableManager.h>
#include <dataaccess/IDataConverterImpl.h>
#include <utils/EigenDecompose.h>
#include <askap/AskapUtil.h>


// casa
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <measurementequation/SynthesisParamsHelper.h>


// std
#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

using namespace askap;
using namespace synthesis;
using namespace askap::accessors;

void doReadOnlyTest(const IConstDataSource &ds) {
  IDataSelectorPtr sel=ds.createSelector();
  sel->chooseCrossCorrelations();
  //sel->chooseFeed(1);
  //sel<<LOFAR::ParameterSet("test.in").makeSubset("TestSelection.");
  IDataConverterPtr conv=ds.createConverter();  
  conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"Hz");
  const double refMJD = 58287.5;
  conv->setEpochFrame(casa::MEpoch(casa::Quantity(refMJD,"d"),
                      casa::MEpoch::Ref(casa::MEpoch::UTC)),"s");
  conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));                    
  /*
  const casa::MVDirection tangent(SynthesisParamsHelper::convertQuantity("12h30m00.000","rad"),
                                  SynthesisParamsHelper::convertQuantity("-45.00.00.000","rad"));
  const casa::MDirection tangentDir(tangent, casa::MDirection::J2000);
  */
  double firstTimeStamp = -1;
  size_t counter = 0;
  for (IConstDataSharedIter it=ds.createConstIterator(sel,conv);it!=it.end();++it,++counter) {  
       
       const IConstDataAccessor &acc = *it;
       if (firstTimeStamp < 0) {
           firstTimeStamp = acc.time();
       }
       casa::MEpoch epoch(casa::Quantity(refMJD,"d")+casa::Quantity(acc.time(),"s"), casa::MEpoch::Ref(casa::MEpoch::UTC));
       std::cout<<"time: "<<epoch<<" or "<<acc.time() - firstTimeStamp<<" seconds since start, cycle "<<counter + 1<<std::endl;
       for (casa::uInt row = 0; row < acc.nRow(); ++row) {
            std::cout<<row<<" "<<acc.antenna1()[row]<<" "<<acc.antenna2()[row]<<" "<<acc.feed1()[row]<<" "<<
                  acc.uvw()[row](0)<<" "<<acc.uvw()[row](1)<<" "<<acc.uvw()[row](2)<<std::endl;
            acc.pointingDir1();
       }
  }
}

int main(int argc, char **argv) {
  try {
     if (argc!=2) {
         cerr<<"Usage "<<argv[0]<<" measurement_set"<<endl;
	 return -2;
     }

     casa::Timer timer;

     timer.mark();
     //TableDataSource ds(argv[1],TableDataSource::REMOVE_BUFFERS |
     //                           TableDataSource::MEMORY_BUFFERS);     
     //TableDataSource ds(argv[1],TableDataSource::MEMORY_BUFFERS | TableDataSource::WRITE_PERMITTED);     
     TableDataSource ds(argv[1],TableDataSource::MEMORY_BUFFERS);     
     std::cerr<<"Initialization: "<<timer.real()<<std::endl;
     timer.mark();
     doReadOnlyTest(ds);
     std::cerr<<"Job: "<<timer.real()<<std::endl;
     
  }
  catch(const AskapError &ce) {
     cerr<<"AskapError has been caught. "<<ce.what()<<endl;
     return -1;
  }
  catch(const std::exception &ex) {
     cerr<<"std::exception has been caught. "<<ex.what()<<endl;
     return -1;
  }
  catch(...) {
     cerr<<"An unexpected exception has been caught"<<endl;
     return -1;
  }
  return 0;
}
