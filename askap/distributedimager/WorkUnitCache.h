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

#ifndef ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CACHE_H
#define ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CACHE_H

// own includes
#include <askap/messages/ContinuumWorkUnit.h>
#include <askap/distributedimager/CubeComms.h>

// std includes
#include <set>
#include <string>

// boost includes
#include <boost/noncopyable.hpp>

namespace askap {

namespace synthesis {

/// @brief Class managing caching of a work unit
/// @details This class is created as part of ContinuumWorker refactoring. It provides caching logic where a part of
/// the dataset can be stored on a disk (e.g. a ramdisk) and accessed from there. This class manages corresponding
/// change of the file name, local channel, etc hiding the details from the end user. The code was moved from
/// the original ContinuumWorker pretty much as it was, this class only makes the code more structured. 
/// @ingroup distributedimager
class WorkUnitCache : public boost::noncopyable { 
public:

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
   WorkUnitCache(bool doCache, const std::string& cachePath, cp::CubeComms& comms,
                 casacore::uInt bucketSize = 65536u, casacore::uInt tileNCorr = 4u, casacore::uInt tileNChan = 1u); 

   /// @brief destructor
   ~WorkUnitCache();

   /// @brief revert to the state immediately after construction
   /// @details This method clears all the caches the same way the destructor does and can be used to execute clean up
   /// at a particular time (rather than when the object goes out of scope).
   void reset();

   /// @brief main entry - cache given work units if necessary
   /// @param[in] wu const reference to a work unit to process
   /// @return const reference to either given work unit or the cached one
   const cp::ContinuumWorkUnit& operator()(const cp::ContinuumWorkUnit &wu);

   /// @brief check if caching is enabled
   /// @details This is largely used for cross-check of compatibility with other modes. In general, caching is hidden from
   /// the end user.
   /// @return true if caching is enabled, false otherwise
   inline bool cacheEnabled() const { return itsCacheEnabled; }
   
private:
   /// @brief take the appropriate slice from the measurement set and setup caching
   /// @details This method does the bulk of work. It modifies itsCachedWorkUnit directly.
   void cacheWorkUnit();

   /// @brief cached work unit
   cp::ContinuumWorkUnit itsCachedWorkUnit;

   /// @brief flag whether to do caching
   /// @details Turning this flag to false in the constructor will force operator() to pass
   /// the reference unchanged. Otherwise, it would return the reference to itsCachedWorkUnit
   /// after appropriate slicing of the measurement set
   const bool itsCacheEnabled;

   /// @brief file names for all datasets where cleanup action is required
   std::set<std::string> itsFilesForCleanup;

   /// @brief path to the directory with cached measurement set slices
   const std::string itsCachePath;

   /// @brief reference to communication object (used for synchronisation between other ranks if
   /// caching is done)
   cp::CubeComms& itsComms;

   /// @brief bucket size for the cached MS
   casacore::uInt itsBucketSize;

   /// @brief number of correlations in the tile for the cached MS
   casacore::uInt itsTileNCorr;

   /// @brief number of channels in the tile for the cached MS
   casacore::uInt itsTileNChan;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CACHE_H

