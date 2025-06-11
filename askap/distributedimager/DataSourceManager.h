/// @file
/// @brief Class managing an instance of TableDataSource object
/// @details This class is created as part of ContinuumWorker refactoring. It provides caching logic which
/// allows us to avoid closing and opening the same dataset unnecessarily. There was also some cleanup logic
/// in the code prior to refactoring the necessity of which is not entirely clear to me (MV) as it should've happened
/// in the destructors. It is possible that the cleanup was needed due to an unrelated bug which I fixed as part of 
/// the refactoring (there was an unnecessary copy of the DataSource object). But for now I just move the original
/// clean up code here. If necessary, caching logic could be made more advanced later on (e.g. to allow simultaneous 
/// opening of multiple datasets).
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

#ifndef ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_DATA_SOURCE_MANAGER_H
#define ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_DATA_SOURCE_MANAGER_H

// own includes
#include <askap/dataaccess/TableDataSource.h>

// std includes
#include <set>
#include <string>

// boost includes
#include <boost/noncopyable.hpp>

namespace askap {

namespace synthesis {

/// @brief Class managing an instance of TableDataSource object
/// @details This class is created as part of ContinuumWorker refactoring. It provides caching logic which
/// allows us to avoid closing and opening the same dataset unnecessarily. There was also some cleanup logic
/// in the code prior to refactoring the necessity of which is not entirely clear to me (MV) as it should've happened
/// in the destructors. It is possible that the cleanup was needed due to an unrelated bug which I fixed as part of 
/// the refactoring (there was an unnecessary copy of the DataSource object). But for now I just move the original
/// clean up code here. If necessary, caching logic could be made more advanced later on (e.g. to allow simultaneous 
/// opening of multiple datasets).
/// @ingroup distributedimager
class DataSourceManager : public boost::noncopyable { 
public:

   /// @brief constructor
   /// @param[in] dataColumn the name of the data column in the measurement set to use
   /// @param[in] clearCache if true, reset or destructor will perform cache cleanup for all previously opened files
   /// @param[in] uvwMachineCacheSize the size of the uvw machine cache for newly created data sources
   /// @param[in] uvwMachineCacheTolerance directional tolerance in radians for the uvw machine cache of data sources
   DataSourceManager(const std::string &dataColumn, bool clearCache, size_t uvwMachineCacheSize, double uvwMachineCacheTolerance);

   /// @brief destructor
   ~DataSourceManager();

   /// @brief revert to the state immediately after construction
   /// @details This method destroys the datasource if it has been created and performs cleanup action if necessary.
   /// The same operation happens in the destructor, so it may result in a more clear code if this object is recreated 
   /// rather than reused (as the overheads are small).
   void reset();

   /// @brief get datasource for the given file name
   /// @details It is created on demand. The existing object is returned if the file name is the same. If the requested 
   /// file name is different, the old data source object is destroyed. Note however, that if the cleanup action is
   /// enabled it is not called until explicit reset or destructor call. The size of uvw machine cache and the corresponding tolerance
   /// set in the constructor are passed to each new data source object.
   /// @param[in] name file name of the measurement set to use
   /// @return reference to the data source object
   /// @note Technically, const data source would be sufficient for our use. But there is some technical debt in the code
   /// requiring non-const one throughout. Leave as it was before the refactoring for now. It may be changed in the future.
   accessors::TableDataSource& dataSource(const std::string &name);

   /// @brief force creation of the new data source next time
   /// @details This method disposes the old data source object (if created previously, otherwise - no operation) which 
   /// forces the creation of a new one next time dataSource method is called. This is handy if measurement set itself 
   /// represents a cache and has been replaced since the last call to dataSource method
   void forceNewDataSourceNextTime();

private:
   /// @brief cached data source object
   boost::shared_ptr<accessors::TableDataSource> itsDataSource;

   /// @brief file name for the cached data source object
   std::string itsCachedFileName;

   /// @brief file names for all datasets where cleanup action is required
   std::set<std::string> itsFilesForCleanup;

   /// @brief the name of the data column to use
   const std::string itsDataColumn;

   /// @brief flag whether to clear cache as part of the tear down
   const bool itsClearCache;

   /// @brief number of uvw machines cached for each new datasource
   size_t itsUVWMachineCacheSize;

   /// @brief direction tolerance in radians for the uvw machine cache
   /// @details This is passed to datasource every time a new object is created
   double itsUVWMachineCacheTolerance;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_DATA_SOURCE_MANAGER_H

