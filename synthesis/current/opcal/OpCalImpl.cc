/// @file
///
/// @brief Generic operations-specific calibration
/// @details This class is intended for the types of calibration which 
/// cannot follow ASKAP's predict-forward approach, i.e. which require 
/// observations of various fields done in some special way. It is intended
/// for experimentation with calibration as well as some operation-specific
/// tasks like baseline and pointing refinements.  
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>

ASKAP_LOGGER(logger, ".OpCalImpl");


#include <opcal/OpCalImpl.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>


namespace askap {

namespace synthesis {

/// @brief Constructor from ParameterSet
/// @details The parset is used to construct the internal state. We could
/// also support construction from a python dictionary (for example).
/// @param comms communication object 
/// @param parset ParameterSet for inputs
/// @note We don't use parallel aspect at this stage, the code expects a single rank if compiled with MPI.
OpCalImpl::OpCalImpl(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset) : 
   itsConfig(parset), itsScanStats(parset.getFloat("maxtime",-1.), parset.getInt("maxcycles",-1))
{
   ASKAPCHECK(!comms.isParallel(), "This application is not intended to be used in parallel mode (at this stage)");  
}


/// @brief main entry point
void OpCalImpl::run()
{
  inspectData();
  ASKAPLOG_INFO_STR(logger, "Found "<<itsScanStats.size()<<" chunks in the supplied data");                   
}
   
   
/// @brief gather scan statistics
/// @details This method iterates over data for all supplied MSs and fills itsScanStats at the end.
/// Optional parameters describing how to break long observations are taken from the parset.
void OpCalImpl::inspectData()
{
   const std::string ms = config().getString("dataset");
   accessors::TableDataSource ds(ms, accessors::TableDataSource::MEMORY_BUFFERS/*, dataColumn()*/);
   accessors::IDataSelectorPtr sel=ds.createSelector();
   sel << config();
   accessors::IDataConverterPtr conv=ds.createConverter();
   conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO)/*getFreqRefFrame()*/, "Hz");
   conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
   conv->setEpochFrame(); // time in seconds since 0 MJD
   accessors::IDataSharedIter it=ds.createIterator(sel, conv);
   ASKAPLOG_INFO_STR(logger, "Inspecting "<<ms);
   itsScanStats.inspect(ms, it);                 
}


} // namespace synthesis

} // namespace askap

