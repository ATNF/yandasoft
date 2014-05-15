/// @file 
/// @brief test of parallel iterator
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


#include <dataaccess/TableDataSource.h>
#include <askap_accessors.h>
#include <askap/AskapLogging.h>
#include <casa/Logging/LogIO.h>
#include <askap/Log4cxxLogSink.h>
ASKAP_LOGGER(logger, "");

#include <askap/AskapError.h>
#include <dataaccess/SharedIter.h>
#include <dataaccess/IDataConverterImpl.h>
#include <parallel/ParallelWriteIterator.h>
#include <askapparallel/AskapParallel.h>

// casa
#include <measures/Measures/MFrequency.h>
#include <tables/Tables/Table.h>
#include <casa/OS/Timer.h>
#include <casa/Arrays/ArrayMath.h>


// std
#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

using namespace askap;
using namespace askap::accessors;
using namespace askap::synthesis;



void doReadWriteTest(askap::askapparallel::AskapParallel &comms, const std::string &name) {
  if (comms.isMaster()) {
     TableDataSource ds(name,TableDataSource::MEMORY_BUFFERS | TableDataSource::WRITE_PERMITTED);     
     IDataSelectorPtr sel=ds.createSelector();
     //sel->chooseFeed(1);  
     IDataConverterPtr conv=ds.createConverter();
     conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"MHz");
     conv->setEpochFrame(casa::MEpoch(casa::Quantity(53635.5,"d"),
                         casa::MEpoch::Ref(casa::MEpoch::UTC)),"s");
     IDataSharedIter it=ds.createIterator(sel,conv);
     ParallelWriteIterator::masterIteration(comms,it);
     ASKAPLOG_INFO_STR(logger, "Master has finished its job");
  }
  if (comms.isWorker()) {
      ParallelWriteIterator it(comms);
      size_t cnt=0;
      for (it.init();it.hasMore();it.next(),++cnt) {
           it->frequency();
           it->pointingDir1();
           it->time();
           it->antenna1();
           it->feed1();
           it->uvw();
           for (casa::uInt chan=0; chan<it->nChannel(); ++chan) {
                it->rwVisibility().xzPlane(chan).set(casa::Complex(it->frequency()[chan]-1420.,0.));
           }
           const double l=0., m=0.003975472185;
           for (casa::uInt row = 0; row<it->nRow(); ++row) {
                for (casa::uInt chan=0; chan<it->nChannel(); ++chan) {
                     const double phase = 2.*casa::C::pi*(it->uvw()(row)(0)*l+it->uvw()(row)(1)*m)/casa::C::c*it->frequency()[chan]*1e6;
                     const casa::Complex phasor(cos(phase),sin(phase));
                     casa::Array<casa::Complex> tmp = it->rwVisibility().yzPlane(row).row(chan);
                     tmp *= phasor;
                }
           }
      }
      ASKAPLOG_INFO_STR(logger, "Worker at rank ="<<comms.rank()<<" completed "<<cnt<<" iterations");
  }
}

int main(int argc, const char **argv) {
  // This class must have scope outside the main try/catch block
  askap::askapparallel::AskapParallel comms(argc, argv);
  
  try {
     // Ensure that CASA log messages are captured
     casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
     casa::LogSink::globalSink(globalSink);
    
     casa::Timer timer;
     timer.mark();
    
     if (argc!=2) {
         ASKAPLOG_FATAL_STR(logger, "Total times - user:   " << timer.user() << " system: " << timer.system()
                             << " real:   " << timer.real());     
         return -2;
     }
     ASKAPCHECK(comms.isParallel(), "This test could only be run as a parallel MPI job");
     
     ASKAPLOG_INFO_STR(logger, "Initialization: "<<timer.real());
     timer.mark();
     doReadWriteTest(comms,argv[1]);    
     ASKAPLOG_INFO_STR(logger, "Job: "<<timer.real());
     
  }
  catch(const AskapError &ce) {
     ASKAPLOG_FATAL_STR(logger, "AskapError has been caught. "<<ce.what());
     comms.abort();
  }
  catch(const std::exception &ex) {
     ASKAPLOG_FATAL_STR(logger, "std::exception has been caught. "<<ex.what());
     comms.abort();
  }
  catch(...) {
     ASKAPLOG_FATAL_STR(logger, "an unexpected exception has been caught");
     comms.abort();
  }
  return 0;
}
