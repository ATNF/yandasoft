/// @file
///
/// @brief Solver for pointing corrections
/// @details This is one of the examples of high-level solvers which 
/// refines antenna pointing. The idea behind the algorithm is similar
/// to the pointing corrections done by online ATCA software. 
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

#ifndef SYNTHESIS_POINTING_SOLVER_H
#define SYNTHESIS_POINTING_SOLVER_H

#include <opcal/IGenericCalSolver.h>
#include <Common/ParameterSet.h>

// std includes
#include <string>
#include <utility>

// casa includes
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/casa/Quanta/MVDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/casa/Arrays/Vector.h>

namespace askap {

namespace synthesis {

/// @brief Solver for pointing corrections
/// @details This is one of the examples of high-level solvers which 
/// refines antenna pointing. The idea behind the algorithm is similar
/// to the pointing corrections done by online ATCA software. 
/// @ingroup opcal
class PointingSolver : virtual public IGenericCalSolver {
public:
  /// @brief constructor
  /// @param[in] parset parameters of the solver
  explicit PointingSolver(const LOFAR::ParameterSet& parset);
    
  /// @brief method to be called when the results are ready
  /// @details This method is called (in derived classes) when the results
  /// are ready. Two parameters describe scans and contain calibration data, respectively.
  /// @param[in] scans description of scans, note separate beams are present as separate scans.
  ///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
  ///            not necessarily match the scans known to online system
  /// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
  ///            (column index is antenna ID used in the measurement set).
  virtual void process(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata);

protected:
  /// @brief obtain MRO reference position
  /// @return position measure
  static casa::MPosition mroPosition();

  /// @brief obtain Az and El at MRO for the given direction
  /// @param[in] dir direction of interest in J2000
  /// @param[in] time time in seconds since 0 MJD
  /// @return direction in AzEl
  static casa::MVDirection mroAzEl(const casa::MVDirection &dir, double time);
  
  /// @brief obtain parallactic angle at MRO for the given direction
  /// @param[in] dir direction of interest in J2000
  /// @param[in] time time in seconds since 0 MJD
  /// @return parallactic angle in radians
  static double mroPA(const casa::MVDirection &dir, double time);
  
  
  /// @brief solve for corrections for just one antenna
  /// @param[in] caldata vector of data for all points of the pattern
  /// @return pair with offsets in two coordinates
  std::pair<double, double> solveOne(const casa::Vector<GenericCalInfo> &caldata) const;
  
private:
   
  /// @brief angular tolerance to what is considered a part of the same pattern in radians
  double itsMaxPatternSize;   
     
};

} // namespace synthesis

} // namespace askap


#endif // #ifndef SYNTHESIS_POINTING_SOLVER_H




