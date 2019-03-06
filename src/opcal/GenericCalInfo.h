/// @file
///
/// @brief Generic calibration information for one antenna
/// @details This is essentially a data holding structure which is 
/// intended to be used with OpCalImpl and related solvers. The main
/// OpCalImpl class obtains delays and complex phases for each interval 
/// identified bases on some criteria. These results are aggregated in a
/// single structure per antenna for each interval and then are passed to
/// a high-level "solver" for analysis. The main reason behind these classes
/// is experimentation with calibration which otherwise would require too much
/// scripting and/or hand work if ccalibrator is used directly.
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

#ifndef SYNTHESIS_GENERIC_CAL_INFO_H
#define SYNTHESIS_GENERIC_CAL_INFO_H

// casa includes
#include <casacore/casa/BasicSL/Complex.h>

namespace askap {

namespace synthesis {

/// @brief Generic calibration information for one antenna
/// @details This is essentially a data holding structure which is 
/// intended to be used with OpCalImpl and related solvers. The main
/// OpCalImpl class obtains delays and complex phases for each interval 
/// identified bases on some criteria. These results are aggregated in a
/// single structure per antenna for each interval and then are passed to
/// a high-level "solver" for analysis. The main reason behind these classes
/// is experimentation with calibration which otherwise would require too much
/// scripting and/or hand work if ccalibrator is used directly.
/// @ingroup opcal
struct GenericCalInfo {
  
  /// @brief default constructor - all parameters are undefined
  GenericCalInfo();
  
  /// @brief constructor setting parameters explicitly
  /// @details
  /// @param[in] gain complex gain to set
  /// @param[in] delay delay to set  
  GenericCalInfo(const casa::Complex &gain, double delay);
  
  /// @brief obtain delay
  /// @return delay (in seconds)
  /// @note an exception is thrown if the value is undefined
  double delay() const;
  
  /// @brief obtain complex gain
  /// @return complex gain
  /// @note an exception is thrown if the value is undefined
  const casa::Complex& gain() const;
  
  /// @brief check if delay is defined
  /// @return true, if delay has been defined
  inline bool delayDefined() const { return itsDelayDefined; }
  
  /// @brief check if gain is defined
  /// @return true, if gain has been defined
  inline bool gainDefined() const { return itsGainDefined; }
  
  /// @brief invalidate the values
  void invalidate();
  
  /// @brief set complex gain
  /// @param[in] gain complex gain to set
  void setGain(const casa::Complex &gain);
  
  /// @brief set delay
  /// @param[in] delay delay in seconds
  void setDelay(double delay);  

private:
  /// @brief complex gain
  casa::Complex itsGain;
  
  /// @brief true, if gain is defined
  bool itsGainDefined;
  
  /// @brief delay (in seconds)
  double itsDelay;

  /// @brief true, if delay is defined
  bool itsDelayDefined;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_GENERIC_CAL_INFO_H



