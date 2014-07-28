/// @file
///
/// Support for parallel statistics accululation to advise on imaging parameters
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

#ifndef SYNTHESIS_ADVISE_PARALLEL_H
#define SYNTHESIS_ADVISE_PARALLEL_H

#include <Common/ParameterSet.h>
#include <parallel/MEParallelApp.h>
#include <measurementequation/VisMetaDataStats.h>
#include <casa/Quanta/MVDirection.h>

#include <boost/shared_ptr.hpp>
#include <string>

namespace askap {

namespace synthesis {

/// @brief parallel helper for the advise utility
/// @details This class does the core operation to run statistics estimators on every measurement set
/// and aggregating the result. Most non-trivial actions happen in the parallel mode.
/// @note It may be a bit untidy to derive this class from MEParallelApp just to reuse a bunch of existing code,
/// but some subtle features like frequency conversion setup may come handy in the future. The goal is that it should
/// work with only the single parameter present in the parset which describes the measurement set(s). 
/// @ingroup parallel
class AdviseParallel : public MEParallelApp 
{
public:
   /// @brief Constructor from ParameterSet
   /// @details The parset is used to construct the internal state. We could
   /// also support construction from a python dictionary (for example).
   /// The command line inputs are needed solely for MPI - currently no
   /// application specific information is passed on the command line.
   /// @param comms communication object 
   /// @param parset ParameterSet for inputs
   AdviseParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset);

   /// @brief make the estimate
   /// @details This method iterates over one or more datasets, accumulates and aggregates statistics. If
   /// tangent point is not defined, two iterations are performed. The first one is to estimate the tangent
   /// point and the second to obtain  
   void estimate();
   
   /// @brief perform the accumulation for the given dataset
   /// @details This method iterates over the given dataset, predicts visibilities according to the
   /// model and subtracts these model visibilities from the original visibilities in the dataset.
   /// This is the core operation of the doSubtraction method, which manages the parallel aspect of it.
   /// All actual calculations are done inside this helper method.
   /// @param[in] ms measurement set name
   void calcOne(const std::string &ms);
      
   /// @brief calculate "normal equations", i.e. statistics for this dataset
   virtual void calcNE();
     
   /// @brief summarise stats into log
   /// @details This method just summarises all stats accumulated in the call to estimate() method
   /// into log. Nothing is done for worker process.
   void summary() const; 

   // stubs for pure virtual methods which we don't use

   /// @brief solve normal equations
   virtual void solveNE() {}

   /// Write the model (runs only in the master)
   /// @param[in] postfix a string to be added to the file name
   virtual void writeModel(const std::string &postfix = std::string()) {}

   /// @brief helper method to get statistics estimator
   /// @return const reference to the statistics estimator
   /// @note An exception is thrown if the estimator is not defined or the method is called from
   /// worker process.
   const VisMetaDataStats& estimator() const;
   
protected:
   
   /// @brief a hopefully temporary method to define missing fields in parset
   /// @details We reuse some code for general synthesis application, but it requires some
   /// parameters (like gridder) to be defined. This method fills the parset with stubbed fields.
   /// Hopefully, it is a temporary approach.
   /// @param parset ParameterSet for inputs
   /// @return new parset 
   static LOFAR::ParameterSet addMissingFields(const LOFAR::ParameterSet& parset);


   /// @brief helper method to broadcast statistics to all workers
   /// @details It seems better conceptually, if all ranks hold the same statistics at the end of
   /// the estimate() call. This allows workers to update their own parsets in the parallel way.
   /// This helper method encapsulates all actions required to broadcast statistics estimator from
   /// the master to all workers.
   void broadcastStatistics();
        
private:
   
   /// @brief optional tangent point
   /// @details Desired tangent point may be given up front. It changes the statistics slightly.
   casa::MVDirection itsTangent;
   
   /// @brief true, if tangent point is defined
   bool itsTangentDefined;
   
   /// @brief w-tolerance for snap-shot imaging
   /// @details Or a negative value if no snap-shot imaging is required.
   double itsWTolerance;
   
   /// @brief statistics estimator
  boost::shared_ptr<VisMetaDataStats> itsEstimator;    
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_ADVISE_PARALLEL_H

