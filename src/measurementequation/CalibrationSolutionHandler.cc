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

#include <measurementequation/CalibrationSolutionHandler.h>
#include <askap/AskapError.h>

namespace askap {

namespace synthesis {

/// @brief default constructor
/// @details It constructs the handler class with uninitialised shared pointer to the
/// solution source
CalibrationSolutionHandler::CalibrationSolutionHandler() : itsCurrentSolutionID(-1) {}

/// @brief construct with the given solution source
/// @details
/// @param[in] css shared pointer to solution source
CalibrationSolutionHandler::CalibrationSolutionHandler(const boost::shared_ptr<accessors::ICalSolutionConstSource> &css) :
     itsCalSolutionSource(css), itsCurrentSolutionID(-1) 
{
  ASKAPCHECK(itsCalSolutionSource, 
      "An attempt to initialise CalibrationSolutionHandler with a void calibration solution source shared pointer");
}

/// @brief setup handler with the given solution source  
/// @details
/// @param[in] css shared pointer to solution source
void CalibrationSolutionHandler::setCalSolutionSource(const boost::shared_ptr<accessors::ICalSolutionConstSource> &css)
{
  ASKAPCHECK(css, 
      "An attempt to initialise CalibrationSolutionHandler with a void calibration solution source shared pointer");  
  itsCalSolutionSource = css;
  itsCalSolutionAccessor.reset();
  itsCurrentSolutionID = -1;
  itsChangeMonitor.notifyOfChanges();    
}


/// @brief helper method to update accessor pointer if necessary
/// @details This method updates the accessor shared pointer if it is 
/// uninitialised, or if it has been updated for the given time.
/// @param[in] time timestamp (seconds since 0 MJD)
void CalibrationSolutionHandler::updateAccessor(const double time) const
{
  ASKAPDEBUGASSERT(itsCalSolutionSource);
  const long newID = itsCalSolutionSource->solutionID(time);
  if ((newID != itsCurrentSolutionID) || !itsCalSolutionAccessor) {
      itsCalSolutionAccessor = itsCalSolutionSource->roSolution(newID);
      itsCurrentSolutionID = newID;
      itsChangeMonitor.notifyOfChanges();    
  }
}
  
/// @brief helper method to get current solution accessor
/// @details This method returns a reference to the current solution
/// accessor or throws an exception if it is uninitialised
/// (this shouldn't happen if updateAccessor is called first)
/// @return a const reference to the calibration solution accessor
const accessors::ICalSolutionConstAccessor& CalibrationSolutionHandler::calSolution() const
{
  ASKAPASSERT(itsCalSolutionAccessor);
  return *itsCalSolutionAccessor;
}


} // namespace synthesis

} // namespace askap

