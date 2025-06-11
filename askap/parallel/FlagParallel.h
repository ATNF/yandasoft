/// @file FlagParallel.h
///
/// FlagParallel: Support for parallel flagging
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
/// @author Mark Wieringa
///

#ifndef FLAG_PARALLEL_H
#define FLAG_PARALLEL_H

// ASKAPsoft includes
#include <Common/ParameterSet.h>
#include <askap/parallel/MEParallelApp.h>
#include <askap/flagging/IFlagger.h>

namespace askap {

namespace synthesis {

/// @brief parallel helper for flagging
/// @details This class does the core operation to flag visibility data
/// @ingroup parallel
class FlagParallel : public MEParallelApp
{
public:
   /// @brief Constructor from ParameterSet
   /// @details The parset is used to construct the internal state.
   /// @param comms communication object
   /// @param parset ParameterSet for inputs
   FlagParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset);

   /// @brief perform the flagging
   /// @details This method iterates over one or more datasets and flag visibilities according to
   /// the parameters given. The intention is to call this method in a worker.
   void doFlag();

 protected:
   /// @brief perform the flagging for the given dataset
   /// @details This method iterates over the given dataset and flags visibilities according to the
   /// parset
   /// @param[in] ms measurement set name
   /// @param[in] distributeByTile set to true to do distribution by tiles if possible
   void flagOne(const std::string &ms, bool distributeByTile = false);

   // stubs for pure virtual methods which we don't use
   /// @brief calculate normal equations
   inline void calcNE() override {}

   /// @brief solve normal equations
   inline void solveNE() override {}

   inline void writeModel(const std::string &postfix = std::string()) override {} 

 private:
     std::vector<std::shared_ptr<IFlagger> > itsFlaggers;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef FLAG_PARALLEL_H
