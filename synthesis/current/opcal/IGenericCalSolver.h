/// @file
///
/// @brief Interface to generic high-level solver
/// @details This interface is intended to be used together with OpCalImpl class 
/// and corresponding application. The idea is that the OpCalImpl class provides
/// all the infrastructure to bin available data into a number of scans and then
/// execute ccalibrator-like algorithm for each of the scans. Then a high-level
/// solver is called for the results. This approach allows us to reuse the code
/// where most of the complexity is, i.e. doing the actual calibration solution
/// (with the support of non-trivial sky model), and have a pluggable post-processing
/// code doing the actual analysis. Example applications are baseline solutions,
/// pointing corrections and may be even holography. Essentially any non-standard
/// calibration which doesn't conform with the standard ccalibrator or cbpcalibrator
/// applications might benefit from OpCalImpl. Concrete high-level solvers are
/// expected to be derived from this interface. 
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

#ifndef SYNTHESIS_I_GENERIC_CAL_SOLVER_H
#define SYNTHESIS_I_GENERIC_CAL_SOLVER_H

#include <opcal/ScanStats.h>
#include <opcal/GenericCalInfo.h>

#include <casa/Arrays/Matrix.h>

namespace askap {

namespace synthesis {


/// @brief Interface to generic high-level solver
/// @details This interface is intended to be used together with OpCalImpl class 
/// and corresponding application. The idea is that the OpCalImpl class provides
/// all the infrastructure to bin available data into a number of scans and then
/// execute ccalibrator-like algorithm for each of the scans. Then a high-level
/// solver is called for the results. This approach allows us to reuse the code
/// where most of the complexity is, i.e. doing the actual calibration solution
/// (with the support of non-trivial sky model), and have a pluggable post-processing
/// code doing the actual analysis. Example applications are baseline solutions,
/// pointing corrections and may be even holography. Essentially any non-standard
/// calibration which doesn't conform with the standard ccalibrator or cbpcalibrator
/// applications might benefit from OpCalImpl. Concrete high-level solvers are
/// expected to be derived from this interface. 
/// @ingroup opcal
struct IGenericCalSolver {

  /// @brief virtual desctructor to keep the compiler happy
  virtual ~IGenericCalSolver() {}
  
  /// @brief abstract method to be called when the results are ready
  /// @details This method is called (in derived classes) when the results
  /// are ready. Two parameters describe scans and contain calibration data, respectively.
  /// @param[in] scans description of scans, note separate beams are present as separate scans.
  ///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
  ///            not necessarily match the scans known to online system
  /// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
  ///            (column index is antenna ID used in the measurement set).
  virtual void process(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata) = 0;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_I_GENERIC_CAL_SOLVER_H

