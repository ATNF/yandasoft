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
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

// own includes
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/parallel/MEParallelApp.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/gridding/VisGridderFactory.h>
#include <askap/gridding/TableVisGridder.h>


using namespace askap;
using namespace askap::synthesis;

/// @brief constructor
/// @details sets communication object and parameter set
/// @param[in] comms communication object
/// @param[in] parset parameter set
MEParallelApp::MEParallelApp(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset, bool useFloat) :
   MEParallel(comms,parset,useFloat), itsDataColName(parset.getString("datacolumn", "DATA")),
   itsUVWMachineCacheSize(1), itsUVWMachineCacheTolerance(1e-6), itsMasterDoesWork(parset.getBool("masterDoesWork",false))
{
   // set up image handler, needed for both master and worker
   SynthesisParamsHelper::setUpImageHandler(parset);

   // set up default reference frame
   SynthesisParamsHelper::setDefaultFreqFrame(getFreqRefFrame());

   // MV: we used to have column selection inside the following if-statement as
   // only workers access measurement sets in proper master-worker design.
   // New imager violates it, hence technically it should not have been derived from
   // mw-framework classes! Some technical debt here. Moving its initialisation to the
   // initalisation at construction solves the immediate problem, however ms substitution
   // still will not work correctly (and never was).
   if (doWork()) {
       const int nProcs = nWorkers();
       const int rank = workerRank();
       /// Get the list of measurement sets
       itsMs = parset.getStringVector("dataset");
//       ASKAPCHECK(itsMs.size()>0, "Need dataset specification");

       if (itsMs.size() == 0) {
           ASKAPLOG_WARN_STR(logger,"dataset not present or empty");
       }
       else {

           if (itsMs.size() == 1) {
               const string tmpl=itsMs[0];
               if (nProcs > 1) {
                   itsMs.resize(nProcs);
               }
               for (int i=0; i<nProcs; ++i) {
                   itsMs[i] = substitute(tmpl);
                   if (rank == i) {
                       ASKAPLOG_DEBUG_STR(logger, "Measurement set "<<tmpl<<
                               " for rank "<<itsComms.rank()<<" is substituted by "<<itsMs[i]);
                   }
               }
           } else {
               ASKAPLOG_DEBUG_STR(logger,
                       "Skip measurement set substitution, names are given explicitly: "<<itsMs);
           }
           // commented out warning: joint imaging many MSs is now a standard processing mode
           // if (nProcs>1) {
           //     if (int(itsMs.size()) != (nProcs)) {
           //         ASKAPLOG_WARN_STR(logger,"Running in parallel, dataset per rank usually required");
           //     }
           // }
       }
       // configure uvw-machine cache parameters (to be set up via Data Source)
       const int cacheSize = parset.getInt32("nUVWMachines",1);
       ASKAPCHECK(cacheSize > 0 ,
           "Cache size is supposed to be a positive number, you have "<<cacheSize);
       itsUVWMachineCacheSize = size_t(cacheSize);
       itsUVWMachineCacheTolerance =
           SynthesisParamsHelper::convertQuantity(parset.getString("uvwMachineDirTolerance",
                                                   "1e-6rad"),"rad");

       ASKAPLOG_DEBUG_STR(logger, "UVWMachine cache will store "<<
           itsUVWMachineCacheSize<<" machines");
       ASKAPLOG_DEBUG_STR(logger, "Tolerance on the directions is "<<
           itsUVWMachineCacheTolerance/casacore::C::pi*180.*3600.<<" arcsec");

       // Create the gridder using a factory acting on a parameterset
       if (parset.isDefined("gridder")) {
           itsGridder = createGridder(comms, parset);
           ASKAPCHECK(itsGridder, "Gridder is not defined correctly");
       }
   }
}
