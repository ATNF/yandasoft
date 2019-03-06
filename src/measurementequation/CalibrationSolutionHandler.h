/// @file
/// 
/// @brief helper class to manage calibration solution source
/// @details We need a similar functionality in a number of places
/// to update calibration solution accessor if new solution is
/// available. This class encapsulates this functionality.
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

#ifndef CALIBRATION_SOLUTION_HANDLER_H
#define CALIBRATION_SOLUTION_HANDLER_H

// own includes
#include <calibaccess/ICalSolutionConstSource.h>
#include <calibaccess/ICalSolutionConstAccessor.h>
#include <utils/ChangeMonitor.h>

// boost includes
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief helper class to manage calibration solution source
/// @details We need a similar functionality in a number of places
/// to update calibration solution accessor if new solution is
/// available. This class encapsulates this functionality.
/// @ingroup measurementequation
class CalibrationSolutionHandler {
public:
  /// @brief default constructor
  /// @details It constructs the handler class with uninitialised shared pointer to the
  /// solution source
  CalibrationSolutionHandler();

  /// @brief construct with the given solution source
  /// @details
  /// @param[in] css shared pointer to solution source
  explicit CalibrationSolutionHandler(const boost::shared_ptr<accessors::ICalSolutionConstSource> &css);
  
  /// @brief setup handler with the given solution source  
  /// @details
  /// @param[in] css shared pointer to solution source
  void setCalSolutionSource(const boost::shared_ptr<accessors::ICalSolutionConstSource> &css);

  /// @brief helper method to update accessor pointer if necessary
  /// @details This method updates the accessor shared pointer if it is 
  /// uninitialised, or if it has been updated for the given time.
  /// @param[in] time timestamp (seconds since 0 MJD)
  void updateAccessor(const double time) const;
  
  /// @brief helper method to get current solution accessor
  /// @details This method returns a reference to the current solution
  /// accessor or throws an exception if it is uninitialised
  /// (this shouldn't happen if updateAccessor is called first)
  /// @return a const reference to the calibration solution accessor
  const accessors::ICalSolutionConstAccessor& calSolution() const;
  
  /// @brief obtain change monitor
  /// @details This class is handy if one wants to track changes in the
  /// solution accessor without analysing its content in detail (i.e. for
  /// caching of derived information). One just has to compare change monitors
  /// and if they're not equal, obtain a new solution accessor via calSolution
  /// (and a new change monitor)
  inline scimath::ChangeMonitor changeMonitor() const { return itsChangeMonitor;}

private:
  /// @brief solution source to work with
  boost::shared_ptr<accessors::ICalSolutionConstSource> itsCalSolutionSource;
  
  /// @brief shared pointer to the current solution accessor 
  /// @details It is updated every time the time changes.
  mutable boost::shared_ptr<accessors::ICalSolutionConstAccessor> itsCalSolutionAccessor;  
  
  /// @brief solution ID corresponding to the current solution accessor
  mutable long itsCurrentSolutionID;

  /// @brief change monitor (to track changes in the accessor returned by calSolution)
  mutable scimath::ChangeMonitor itsChangeMonitor;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef CALIBRATION_SOLUTION_HANDLER_H

