/// @file
/// MV: Note, this file is not currently used (and excluded from cmake). It is left checked in just in case we need it in the future.
/// @brief Class managing caching of a work unit
/// @details This class is created as part of ContinuumWorker refactoring. It provides caching logic where a part of
/// the dataset can be stored on a disk (e.g. a ramdisk) and accessed from there. This class manages corresponding
/// change of the file name, local channel, etc hiding the details from the end user. The code was moved from
/// the original ContinuumWorker pretty much as it was, this class only makes the code more structured. 
///
/// @copyright (c) 2023 CSIRO
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
/// refactoring by Max Voronkov <maxim.voronkov@csiro.au>, not sure who wrote the original
/// code which was moved essentially unchanged


// own includes
#include <askap/distributedimager/WorkUnitCache.h>
#include <askap/distributedimager/MSSplitter.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>

// std includes
#include <string>
#include <sstream>
#include <sys/stat.h>

// boost includes
#include "boost/filesystem.hpp"

ASKAP_LOGGER(logger, ".WorkUnitCache");

namespace askap {

namespace synthesis {

/// @brief constructor
/// @param[in] doCache if true caching will happen, otherwise operator() passes the reference unchanged
/// @param[in] cachePath path to the directory where the caching is done (e.g. a RAM disk), unused if doCache is false
/// @param[in] comms reference to communication object. It is required to provide synchronisation in the case 
/// when caching is done. Note, in this case all ranks should execute this code and have consistent doCache.
/// It might be possible to abstract out this dependency somehow, but for now the same approach is used as it was
/// prior to the refactoring
/// @param[in] bucketSize bucket size for the cached MS (only used if doCache is true)
/// @param[in] tileNCorr number of correlations in the tile for the cached MS (only used if doCache is true)
/// @param[in] tileNChan number of channels in the tile for the cached MS (only used if doCache is true)
/// @note There are additional parameters which influence the behaviour of MSSplitter indirectly via the parset.
/// It is not clear to me (MV), however, whether this is essential for current modes of operation. But it could 
/// lead to unexpected corner cases.
WorkUnitCache::WorkUnitCache(bool doCache, const std::string& cachePath, cp::CubeComms& comms,
                 casacore::uInt bucketSize, casacore::uInt tileNCorr, casacore::uInt tileNChan) :
       itsCacheEnabled(doCache), itsCachePath(cachePath), itsComms(comms), itsBucketSize(bucketSize),
       itsTileNCorr(tileNCorr), itsTileNChan(tileNChan) {}     

/// @brief destructor
WorkUnitCache::~WorkUnitCache()
{
   reset();
}

/// @brief revert to the state immediately after construction
/// @details This method clears all the caches the same way the destructor does and can be used to execute clean up
/// at a particular time (rather than when the object goes out of scope).
void WorkUnitCache::reset()
{
  // MV: the original code prior to refactoring didn't have inter-rank coordination for deleting files. I think this is wrong.
  // We need at least the same kind of approach (which is not ideal either) as during caching when these files are written
  if (itsComms.inGroup(0)) {
      for (std::string ms : itsFilesForCleanup) {
           struct stat buffer;

           if (stat(ms.c_str(), &buffer) == 0) {
               ASKAPLOG_DEBUG_STR(logger, "Split file " << ms << " exists - deleting");
               boost::filesystem::remove_all(ms.c_str());
           } else {
               ASKAPLOG_DEBUG_STR(logger, "Split file " << ms << " does not exist - nothing to do");
           }
      }
  }
  itsFilesForCleanup.clear();
}

/// @brief main entry - cache given work units if necessary
/// @param[in] wu const reference to a work unit to process
/// @return const reference to either given work unit or the cached one
const cp::ContinuumWorkUnit& WorkUnitCache::operator()(const cp::ContinuumWorkUnit &wu)
{
   if (itsCacheEnabled) {
       itsCachedWorkUnit = wu;
       // the following modifies itsCachedWorkUnit as appropriate and creates a slipe from the measurement set
       cacheWorkUnit();
       return itsCachedWorkUnit;
   }
   return wu;
}
   
/// @brief take the appropriate slice from the measurement set and setup caching
/// @details This method does the bulk of work. It modifies itsCachedWorkUnit directly.
void WorkUnitCache::cacheWorkUnit()
{
   // MV: the code was taken from ContinuumWorker with minimal changes during refactoring, it could be improved +
   // synchronisation seems to be 1) specific to our use case 2) probably untidy anyway.
  
   const boost::filesystem::path mspath = boost::filesystem::path(itsCachedWorkUnit.get_dataset());
   const string ms = mspath.filename().string();

   std::ostringstream pstr;

   pstr << itsCachePath << "/" << ms << "_chan_" << itsCachedWorkUnit.get_localChannel() + 1 << "_beam_" << itsCachedWorkUnit.get_beam() << ".ms";

   const std::string outms = pstr.str();
   pstr << ".working";

   const std::string outms_flag = pstr.str();

   // MV: this assumes that different groups will access the same dataset roughly at the same time, so only group 0 needs to make the slice
   if (itsComms.inGroup(0)) {

       struct stat buffer;

       // MV: this is ugly and prone to deadlock. Nothing should write the same file unless the high-level logic is stuffed up somehow,
       // perhaps we'd be better with an exception.

       while (stat(outms_flag.c_str(), &buffer) == 0) {
         // flag file exists - someone is writing

         sleep(1);
       }

       if (stat(outms.c_str(), &buffer) == 0) {
           ASKAPLOG_WARN_STR(logger, "Split file already exists");
       } else if (stat(outms.c_str(), &buffer) != 0 && stat(outms_flag.c_str(), &buffer) != 0) {
          // file cannot be read

          // MV: it's a bad way of doing synchronisation through flag files, prone to race condition. 
          // Leave it as it was for now, we shouldn't really rely on this flag.

          // drop trigger
          std::ofstream trigger;
          trigger.open(outms_flag.c_str());
          trigger.close();
          cp::MSSplitter mySplitter(itsBucketSize, itsTileNCorr, itsTileNChan);

          // MV: the configuration of splitter in the original code depended on the parset (which could change). In particular, it could apply
          // filtering based on beams and scans. I don't think we rely on this behaviour anywhere. But there could be some corner cases.
          // On the other hand, the original code never applied filtering based on the beam in the work unit. This is done below. 
          mySplitter.chooseBeams(std::set<unsigned int>({itsCachedWorkUnit.get_beam()}));

          // MV: I am surprised we need channel + 1 here, I thought everything is zero-based. Leave as it was prior to the refactoring
          mySplitter.split(itsCachedWorkUnit.get_dataset(), outms, itsCachedWorkUnit.get_localChannel() + 1, itsCachedWorkUnit.get_localChannel() + 1, 1);
          unlink(outms_flag.c_str());
          itsFilesForCleanup.insert(outms);
       }
   }
   // MV: the original code prior to refactoring has the barrier part below inside the if-statement above (for group 0). I think this is a bug which would lead to an MPI deadlock,
   // changed it to have all ranks to execute the barrier. Note, master also should do this. Perhaps, it would be better to have a communicator which only involves workers here rather than
   // trying to optimise - see also AXA-2889, made similar changes here. In addition, we're probably going to take
   // this code out long term because tmpfs feature is not used any more. 
   ///wait for all groups this rank to get here
   if (itsComms.nGroups() > 1) {
       ASKAPLOG_DEBUG_STR(logger, "Rank " << itsComms.rank() << " at barrier");
       //itsComms.barrier(itsComms.interGroupCommIndex());
       itsComms.barrier(itsComms.theWorkers());
       ASKAPLOG_DEBUG_STR(logger, "Rank " << itsComms.rank() << " passed barrier");
   }
   itsCachedWorkUnit.set_dataset(outms);
   // MV: the original code prior to refactoring didn't set the local channel in the cached work unit. This happened in the logic at the higher level. Setting it here simplifies the logic.
   itsCachedWorkUnit.set_localChannel(0);
}

} // namespace synthesis

} // namespace askap
