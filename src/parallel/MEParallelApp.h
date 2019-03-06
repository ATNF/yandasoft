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

#ifndef ME_PARALLEL_APP_H
#define ME_PARALLEL_APP_H

// ASKAPsoft includes
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <parallel/MEParallel.h>
#include <gridding/IVisGridder.h>


// std includes
#include <string>
#include <vector>

namespace askap {

namespace synthesis {

/// @brief support for parallel applications using the measurement equation classes. 
/// @details This code implements common behavior for imaging, calibration and
/// continuum subtraction. Unlike MEParallel it has some application-specific code in
/// addition to parallelism.
/// @ingroup parallel
class MEParallelApp : public MEParallel
{
public:
   /// @brief constructor 
   /// @details sets communication object and parameter set
   /// @param[in] comms communication object
   /// @param[in] parset parameter set
   MEParallelApp(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset);

protected:
   
   /// @brief obtain data column name
   /// @details to be used in derived classes
   /// @return reference to the string with the data column name
   inline const std::string& dataColumn() const { return itsDataColName;}
   
   /// @brief obtain the array of names for all measurement sets
   /// @details to be used in derived classes
   /// @return reference to the vector of strings with the file names
   inline const std::vector<std::string>& measurementSets() const { return itsMs;}

   /// @brief obtain uvw-machine cache size
   /// @details to be used in derived classes
   /// @return number of uvw-machines in the cache
   inline size_t uvwMachineCacheSize() const { return itsUVWMachineCacheSize; }
   
   /// @brief obtain direction tolerance for uvw-machine cache
   /// @details to be used in derived classes
   /// @return direction tolerance (in radians) for uvw machine cache
   inline double uvwMachineCacheTolerance() const { return itsUVWMachineCacheTolerance; }
   
   /// @brief obtain gridder
   /// @details to be used in derived classes
   /// @return shared pointer to the gridder template
   inline IVisGridder::ShPtr gridder() const { return itsGridder; }
   

protected:
   /// @brief set the list of measurement sets
   /// @details For bandpass calibration we need a different logic - number of measurement sets
   /// is not 1, but less than the number of workers. Using this method allows to bypass
   /// checks that the number of workers is equal to the number of measurement sets specified
   /// in the parset
   /// @param[in] ms vector with measurement sets
   inline void setMeasurementSets(const std::vector<std::string>& ms) { itsMs = ms;}
   
private:   

   /// @brief name of the data column to use.
   std::string itsDataColName;

   /// Names of measurement sets, one per prediffer
   std::vector<std::string> itsMs;	
   
   /// @brief uvw machine cache size
   size_t itsUVWMachineCacheSize;
      
   /// @brief direction tolerance (in radians) for uvw machine cache
   double itsUVWMachineCacheTolerance;

   /// @brief gridder to be used
   IVisGridder::ShPtr itsGridder;		    			  	
}; 

} // namespace synthesis

} // namespace askap


#endif // #ifndef ME_PARALLEL_APP_H

