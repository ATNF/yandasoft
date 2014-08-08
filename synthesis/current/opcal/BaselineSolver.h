/// @file
///
/// @brief Solver for array baselines
/// @details This is one of the examples of high-level solvers which 
/// refines antenna positions on the ground. Unlike generic algorithms
/// which just throw measurements to a least square solver, this one
/// attempts to unwrap the phase to make it more robust against large
/// errors in antenna positions (or timing). 
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

#ifndef SYNTHESIS_BASELINE_SOLVER_H
#define SYNTHESIS_BASELINE_SOLVER_H

#include <opcal/IGenericCalSolver.h>
#include <Common/ParameterSet.h>


namespace askap {

namespace synthesis {

/// @brief Solver for array baselines
/// @details This is one of the examples of high-level solvers which 
/// refines antenna positions on the ground. Unlike generic algorithms
/// which just throw measurements to a least square solver, this one
/// attempts to unwrap the phase to make it more robust against large
/// errors in antenna positions (or timing). 
/// @ingroup opcal
class BaselineSolver : virtual public IGenericCalSolver {
public:
  /// @brief constructor
  /// @param[in] parset parameters of the solver
  explicit BaselineSolver(const LOFAR::ParameterSet& parset);
    
  /// @brief abstract method to be called when the results are ready
  /// @details This method is called (in derived classes) when the results
  /// are ready. Two parameters describe scans and contain calibration data, respectively.
  /// @param[in] scans description of scans, note separate beams are present as separate scans.
  ///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
  ///            not necessarily match the scans known to online system
  /// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
  ///            (column index is antenna ID used in the measurement set).
  virtual void process(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata);

protected:

  /// @brief solve for dX and dY and populate corrections
  /// @param[in] scans description of scans, note separate beams are present as separate scans.
  ///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
  ///            not necessarily match the scans known to online system
  /// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
  ///            (column index is antenna ID used in the measurement set).
  void solveForXY(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata);

  /// @brief solve for dZ and populate corrections
  /// @param[in] scans description of scans, note separate beams are present as separate scans.
  ///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
  ///            not necessarily match the scans known to online system
  /// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
  ///            (column index is antenna ID used in the measurement set).
  void solveForZ(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata);

protected:
  /// @brief obtain MRO reference position
  /// @return position measure
  static casa::MPosition mroPosition();

private:
   
  /// @brief if true, solve for dX and dY
  bool itsSolveXY;
  
  /// @brief if true, solve for dZ
  bool itsSolveZ;
  
  /// @brief reference antenna
  casa::uInt itsRefAnt;
  
  /// @brief buffer for corrections
  /// @details rows are antennas, columns are X,Y and Z
  casa::Matrix<double> itsCorrections;
  
  /// @brief uncertainties for corrections
  /// @details rows are antennas, columns are X,Y and Z
  casa::Matrix<double> itsErrors;
   
};

} // namespace synthesis

} // namespace askap


#endif // #ifndef SYNTHESIS_BASELINE_SOLVER_H



