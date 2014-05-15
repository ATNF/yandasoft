/// @file
///
/// MEParallelApp: Support for parallel applications using the measurement 
/// equation classes. This code implements common behavior for imaging, calibration and
/// continuum subtraction. Unlike MEParallel it has some application-specific code in
/// addition to parallelism.
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
///

// logging stuff
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

// own includes
#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <parallel/MEParallelApp.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <gridding/VisGridderFactory.h>
#include <gridding/TableVisGridder.h>


using namespace askap;
using namespace askap::synthesis;

/// @brief constructor 
/// @details sets communication object and parameter set
/// @param[in] comms communication object
/// @param[in] parset parameter set
MEParallelApp::MEParallelApp(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset) : 
   MEParallel(comms,parset),   
   itsUVWMachineCacheSize(1), itsUVWMachineCacheTolerance(1e-6)   
{
   // set up image handler, needed for both master and worker
   SynthesisParamsHelper::setUpImageHandler(parset);

   // set up default reference frame
   SynthesisParamsHelper::setDefaultFreqFrame(getFreqRefFrame());

   if (itsComms.isWorker()) {
       /// Get the list of measurement sets and the column to use.
       itsDataColName = parset.getString("datacolumn", "DATA");
       itsMs = parset.getStringVector("dataset");
        
       ASKAPCHECK(itsMs.size()>0, "Need dataset specification");
       const int nProcs = itsComms.nProcs();
       
       if (itsMs.size() == 1) {
           const string tmpl=itsMs[0];
           if (nProcs>2) {
               itsMs.resize(nProcs-1);
           }
           for (int i=0; i<nProcs-1; ++i) {
                itsMs[i] = substitute(tmpl);
                if ((itsComms.rank() - 1) == i) {
                    ASKAPLOG_INFO_STR(logger, "Measurement set "<<tmpl<<" for rank "<<i+1<<" is substituted by "<<itsMs[i]);
                }
           }
       } else {
          ASKAPLOG_INFO_STR(logger, "Skip measurment set substitution, names are given explicitly: "<<itsMs);
       }
       if (nProcs>1) {
           ASKAPCHECK(int(itsMs.size()) == (nProcs-1),
              "When running in parallel, need one data set per node");
       } 
       
       // configure uvw-machine cache parameters (to be set up via Data Source)
       const int cacheSize = parset.getInt32("nUVWMachines",1);
       ASKAPCHECK(cacheSize > 0 ,"Cache size is supposed to be a positive number, you have "<<cacheSize);
       itsUVWMachineCacheSize = size_t(cacheSize);
       itsUVWMachineCacheTolerance = 
           SynthesisParamsHelper::convertQuantity(parset.getString("uvwMachineDirTolerance", 
                                                   "1e-6rad"),"rad");
         
       ASKAPLOG_DEBUG_STR(logger, "UVWMachine cache will store "<<itsUVWMachineCacheSize<<" machines");
       ASKAPLOG_DEBUG_STR(logger, "Tolerance on the directions is "<<itsUVWMachineCacheTolerance/casa::C::pi*180.*3600.<<" arcsec");
        
       // Create the gridder using a factory acting on a parameterset
       itsGridder = createGridder(comms, parset);
       ASKAPCHECK(itsGridder, "Gridder is not defined correctly");              
   }
}

