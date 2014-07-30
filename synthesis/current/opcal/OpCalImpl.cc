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


// ASKAPsoft includes

#include <opcal/OpCalImpl.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>

// std includes

#include <string>
#include <vector>
#include <set>

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
   if (itsScanStats.timeLimit() > 0.) {
       ASKAPLOG_INFO_STR(logger, "Chunks will be limited to "<<itsScanStats.timeLimit()<<" seconds"); 
   }  
   if (itsScanStats.cycleLimit() > 0) {
       ASKAPLOG_INFO_STR(logger, "Chunks will be limited to "<<itsScanStats.cycleLimit()<<" correlator cycles in length"); 
   }     
}


/// @brief main entry point
void OpCalImpl::run()
{
  inspectData();
  ASKAPLOG_INFO_STR(logger, "Found "<<itsScanStats.size()<<" chunks in the supplied data");
  runCalibration();                   
}
   
   
/// @brief gather scan statistics
/// @details This method iterates over data for all supplied MSs and fills itsScanStats at the end.
/// Optional parameters describing how to break long observations are taken from the parset.
void OpCalImpl::inspectData()
{
   const std::vector<std::string> msList = config().getStringVector("dataset");
   for (std::vector<std::string>::const_iterator ci = msList.begin(); ci!=msList.end(); ++ci) {
        
        const size_t sizeBefore = itsScanStats.size();
        accessors::TableDataSource ds(*ci, accessors::TableDataSource::MEMORY_BUFFERS);
        accessors::IDataSelectorPtr sel=ds.createSelector();
        sel << config();
        accessors::IDataConverterPtr conv=ds.createConverter();
        conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO), "Hz");
        conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
        conv->setEpochFrame(); // time in seconds since 0 MJD
        accessors::IDataSharedIter it=ds.createIterator(sel, conv);
        ASKAPLOG_INFO_STR(logger, "Inspecting "<<*ci);
        itsScanStats.inspect(*ci, it);
        ASKAPLOG_INFO_STR(logger, "   - found "<<itsScanStats.size() - sizeBefore<<" chunks");
   }                 
}

/// @brief perform calibration for every scan
/// @details This method runs calibration procedure for each scan in itsScanStats, initialises and
/// fills itsCalData
void OpCalImpl::runCalibration()
{   
   const casa::uInt nAnt = config().getUint("nant",6);
   ASKAPCHECK(nAnt > 0, "Expect a positive number of antennas");
   itsCalData.resize(itsScanStats.size(), nAnt);
   // ensure each element is undefined although this is not necessary, strictly speaking
   for (casa::uInt scan=0; scan < itsCalData.nrow(); ++scan) {
        for (casa::uInt ant=0; ant < itsCalData.ncolumn(); ++ant) {
             itsCalData(scan,ant).invalidate();
        }
   }
   // now obtain calibration information for every "scan"
   // first, build a set of all names to have more structured access to the data (in the hope of having a speed up)
   std::set<std::string> names;
   for (ScanStats::const_iterator ci = itsScanStats.begin(); ci != itsScanStats.end(); ++ci) {
        names.insert(ci->name()); 
   }
   for (std::set<std::string>::const_iterator ci = names.begin(); ci != names.end(); ++ci) {
        ASKAPLOG_INFO_STR(logger, "Performing calibration for "<<*ci);
        // 
   }
} 



} // namespace synthesis

} // namespace askap

