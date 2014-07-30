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

#ifndef SYNTHESIS_OP_CAL_IMPL_H
#define SYNTHESIS_OP_CAL_IMPL_H

// casa includes
#include <casa/Arrays/Matrix.h>

// ASKAPsoft includes
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <opcal/ScanStats.h>
#include <opcal/GenericCalInfo.h>


namespace askap {

namespace synthesis {

/// @brief Generic operations-specific calibration
/// @details This class is intended for the types of calibration which 
/// cannot follow ASKAP's predict-forward approach, i.e. which require 
/// observations of various fields done in some special way. It is intended
/// for experimentation with calibration as well as some operation-specific
/// tasks like baseline and pointing refinements.  
/// @ingroup opcal
class OpCalImpl {
public:
   /// @brief Constructor from ParameterSet
   /// @details The parset is used to construct the internal state. We could
   /// also support construction from a python dictionary (for example).
   /// @param comms communication object 
   /// @param parset ParameterSet for inputs
   /// @note We don't use parallel aspect at this stage, the code expects a single rank if compiled with MPI.
   OpCalImpl(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset);


   /// @brief main entry point
   void run();
   
protected:
   
   /// @brief gather scan statistics
   /// @details This method iterates over data for all supplied MSs and fills itsScanStats at the end.
   /// Optional parameters describing how to break long observations are taken from the parset.
   void inspectData();
   
   /// @brief perform calibration for every scan
   /// @details This method runs calibration procedure for each scan in itsScanStats, initialises and
   /// fills itsCalData
   void runCalibration(); 
     
   /// @brief obtain configuration parset
   /// @return parset with configuration
   inline const LOFAR::ParameterSet& config() const {return itsConfig;}
   
private:
   /// @brief configuration
   LOFAR::ParameterSet itsConfig;
   
   /// @brief details for every scan
   ScanStats  itsScanStats;
   
   /// @brief calibration data
   /// @details rows are scans (matching itsScanStats.size()), columns are antennas
   casa::Matrix<GenericCalInfo> itsCalData;    
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_OP_CAL_IMPL_H


