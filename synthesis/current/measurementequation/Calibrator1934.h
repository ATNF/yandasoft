/// @file
///
/// @brief Model of the 1934-638 calibrator observation
/// @details
///     This is a specialised class describting calibration observation of
///     1934-638 to be used as a calibration model. We could've expressed 1934-638
///     via existing classes (like an unpolarised point source) by making the flux model
///     more flexible. However, it seems easier at this stage to just implement a special case
///     than to work on the most flexible solution.  
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

#ifndef ASKAP_SYNTHESIS_CALIBRATOR_1934_H
#define ASKAP_SYNTHESIS_CALIBRATOR_1934_H

#include <measurementequation/IUnpolarizedComponent.h>
#include <measurementequation/IParameterizedComponent.h>

namespace askap {

namespace synthesis {

/// @brief Model of the 1934-638 calibrator observation
/// @details
///     This is a specialised class describting calibration observation of
///     1934-638 to be used as a calibration model. We could've expressed 1934-638
///     via existing classes (like an unpolarised point source) by making the flux model
///     more flexible. However, it seems easier at this stage to just implement a special case
///     than to work on the most flexible solution.  
/// @ingroup measurementequation  
struct Calibrator1934 : virtual public IUnpolarizedComponent,
                        virtual public IParameterizedComponent {
  
  
  /// @brief calculate stokes I visibilities for this component
  /// @details This variant of the method calculates just the visibilities
  /// (without derivatives) for a number of frequencies. This method 
  /// is used in the implementation
  /// of the IComponent interface if stokes I is requested. Otherwise the result
  /// is filled with 0.
  /// @param[in] uvw  baseline spacings (in metres)
  /// @param[in] freq vector of frequencies to do calculations for (in Hz)
  /// @param[out] result an output buffer used to store values
  virtual void calculate(const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    std::vector<double> &result) const;
  
  /// @brief calculate stokes I visibilities and derivatives for this component
  /// @details This variant of the method does simultaneous calculations of
  /// the values and derivatives. The component is assumed to have no free parameters, so
  /// all derivatives are always zero. This method is used in the implementation
  /// of the IComponent interface if stokes I is requested. Otherwise result
  /// is filled with 0.
  /// @param[in] uvw  baseline spacings (in metres)
  /// @param[in] freq vector of frequencies to do calculations for (in Hz)
  /// @param[out] result an output buffer used to store values
  virtual void calculate(const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    std::vector<casa::AutoDiff<double> > &result) const;                    


  // to import full API
  using IUnpolarizedComponent::calculate;
  
  /// @brief get number of parameters
  /// @return a number of parameters  this component depends upon. 
  virtual size_t nParameters() const throw();
  
  /// @brief get the name of the given parameter
  /// @details All parameters are handled in the synthesis code using their
  /// string name, which allows to fix or free any of them easily. This method
  /// allows to obtain this string name using a integer index
  /// @param[in] index an integer index of the parameter (should be less than
  /// nParameters).
  /// @return a const reference to the string name of the parameter 
  virtual const std::string& parameterName(size_t index) const;
  

  /// @brief flux model for 1934-638
  /// @details This method estimates the flux density of 1934-638 for a given frequency using the cm-wavelength 
  /// model of Reynolds et al. (see miriad or query 1934-638 in the ATCA calibrator database for a
  /// reference).   
  /// @param[in] freqInMHz frequency of interest (in MHz)
  /// @return estimated flux density in Jy
  static double fluxDensity(const double freqInMHz);  
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_CALIBRATOR_1934_H
