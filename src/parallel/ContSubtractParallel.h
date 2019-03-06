/// @file
///
/// ContSubtractParallel: Support for parallel continuum subtraction using model
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

#ifndef CONT_SUBTRACT_PARALLEL_H
#define CONT_SUBTRACT_PARALLEL_H

// ASKAPsoft includes
#include <Common/ParameterSet.h>
#include <parallel/MEParallelApp.h>

namespace askap {

namespace synthesis {

/// @brief parallel helper for continuum subtraction
/// @details This class does the core operation to subtract continuum by doing visibility
/// prediction from the given model in parallel.
/// @ingroup parallel
class ContSubtractParallel : public MEParallelApp 
{
public:
   /// @brief Constructor from ParameterSet
   /// @details The parset is used to construct the internal state. We could
   /// also support construction from a python dictionary (for example).
   /// The command line inputs are needed solely for MPI - currently no
   /// application specific information is passed on the command line.
   /// @param comms communication object 
   /// @param parset ParameterSet for inputs
   ContSubtractParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset);
 
   /// @brief Initialise continuum subtractor
   /// @details The parameters are taken from the parset file supplied in the constructor.
   /// This method does initialisation which may involve communications in the parallel case
   /// (i.e. distribution of the models between workers). Technically, we could've done this in 
   /// the constructor.
   void init();
   
   /// @brief perform the subtraction
   /// @details This method iterates over one or more datasets, predicts visibilities according to 
   /// the model and subtracts these model visibilities from the original visibilities in the
   /// dataset. The intention is to call this method in a worker.
   void doSubtraction();
 
 protected:
   /// @brief read the models from the parset file
   inline void readModels() const { SynParallel::readModels(itsModel); }
   
   /// @brief initialise measurement equation
   /// @details This method initialises measurement equation
   void initMeasurementEquation();
   
   /// @brief perform the subtraction for the given dataset
   /// @details This method iterates over the given dataset, predicts visibilities according to the
   /// model and subtracts these model visibilities from the original visibilities in the dataset.
   /// This is the core operation of the doSubtraction method, which manages the parallel aspect of it.
   /// All actual calculations are done inside this helper method.
   /// @param[in] ms measurement set name
   void calcOne(const std::string &ms);
      
   // stubs for pure virtual methods which we don't use
   /// @brief calculate normal equations
   inline void calcNE() {}
   
   /// @brief solve normal equations
   inline void solveNE() {}
   
   /// @brief write results
   inline void writeModel(const std::string&) {}
      
 private:
   /// @brief model is read by the master and distributed?
   /// @details Depending on the model file name (containing %w or not), the model
   /// can either be read in the master and distributed across the workers or read
   /// by workers directly. This data member is true, if the model is read by the master
   bool itsModelReadByMaster;    
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef CONT_SUBTRACT_PARALLEL_H


