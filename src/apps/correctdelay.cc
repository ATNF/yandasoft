///
/// @file correctdelay.cc : Small utility to apply delay to the existing measurement set
/// @details This utility is largely intended for MRO experiments. It may evolve into 
/// something bigger over time. But ideally we need to implement the operation via 
/// the measurement equation.
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
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, "");

#include <askap/AskapError.h>
#include <dataaccess/SharedIter.h>
#include <dataaccess/ParsetInterface.h>

#include <dataaccess/TableManager.h>
#include <dataaccess/IDataConverterImpl.h>

// casa
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>

#include <CommandLineParser.h>

#include <measurementequation/IMeasurementEquation.h>
//#include <dataaccess/MemBufferDataAccessor.h>
#include <askapparallel/MPIComms.h>

// std
#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

using namespace askap;
using namespace synthesis;
using namespace askap::accessors;

/// @brief apply delays to the given data source
/// @details ds data source to work with
/// @details delays vector with delays (one per antenna) (in ns)
void correctDelays(const IDataSource &ds, const casa::Vector<double> &delays) {
  IDataSelectorPtr sel=ds.createSelector();
  //sel->chooseFeed(1);  
  IDataConverterPtr conv=ds.createConverter();
  conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"Hz");
  conv->setEpochFrame(casa::MEpoch(casa::Quantity(53635.5,"d"),
                      casa::MEpoch::Ref(casa::MEpoch::UTC)),"s");
  IDataSharedIter it=ds.createIterator(sel,conv);
  for (it.init();it!=it.end();it.next()) {
       const casa::Vector<double> freqs = it->frequency();
       ASKAPDEBUGASSERT(freqs.nelements() == it->nChannel());
       for (casa::uInt row=0; row< it->nRow(); ++row) {
            casa::Matrix<casa::Complex> thisRow = it->rwVisibility().yzPlane(row);
            const casa::uInt ant1 = it->antenna1()[row];
            const casa::uInt ant2 = it->antenna2()[row];
            ASKAPCHECK((ant1 < delays.nelements()) && (ant2 < delays.nelements()), "Encountered antenna with underfined delay: baseline "<<
                        ant1<<" - "<<ant2<<" at row="<<row<<", delays defined for "<<delays.nelements()<<" antennas");
            // minus sign because we're correcting
            const double delayBy2pi = -2.*casa::C::pi*(delays[ant2] - delays[ant1])*1e-9; // delays are in ns
            for (casa::uInt ch=0; ch<freqs.nelements(); ++ch) {
                 ASKAPDEBUGASSERT(ch < thisRow.nrow());
                 casa::Vector<casa::Complex> allPols = thisRow.row(ch);
                 const float phase = static_cast<float> (delayBy2pi * freqs[ch]);
                 const casa::Complex phasor(cos(phase),sin(phase));
                 allPols *= phasor;
            }
       }
  }
}

int main(int argc, char **argv) {
  try {
     casa::Timer timer;

     timer.mark();
     
     cmdlineparser::Parser parser; // a command line parser
     
     // command line parameters
     cmdlineparser::GenericParameter<std::string> msName;
     parser.add(msName, cmdlineparser::Parser::throw_exception);
     parser.process(argc, argv);
     
     // Initialize MPI (also succeeds if no MPI available).
     askap::askapparallel::MPIComms comms(argc, argv);

     TableDataSource ds(msName.getValue(),TableDataSource::MEMORY_BUFFERS | TableDataSource::WRITE_PERMITTED);     
     
     // delays are hardcoded for now
     casa::Vector<double> delays(3,0.);
     //delays[0] = 0.;
     //delays[1] = 0.;
     //delays[2] = 0.;

     std::cerr<<"Initialization: "<<timer.real()<<std::endl;

     timer.mark();
     correctDelays(ds,delays);
     std::cerr<<"Job: "<<timer.real()<<std::endl;
     
  }
  catch(const cmdlineparser::XParser &) {
     std::cerr<<"Usage "<<argv[0]<<" measurement_set"<<std::endl;
     std::cerr<<"Note, measurement set is replaced in situ!"<<std::endl;
	 return -2;    
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
