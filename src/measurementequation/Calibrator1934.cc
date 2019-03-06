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

#include <measurementequation/Calibrator1934.h>

#include <cmath>

namespace askap {

namespace synthesis {

/// @brief calculate stokes I visibilities for this component
/// @details This variant of the method calculates just the visibilities
/// (without derivatives) for a number of frequencies. This method 
/// is used in the implementation
/// of the IComponent interface if stokes I is requested. Otherwise the result
/// is filled with 0.
/// @param[in] freq vector of frequencies to do calculations for (in Hz)
/// @param[out] result an output buffer used to store values
void Calibrator1934::calculate(const casa::RigidVector<casa::Double, 3> &,
        const casa::Vector<casa::Double> &freq, std::vector<double> &result) const
{
   result.resize(2 * freq.nelements());
   std::vector<double>::iterator it = result.begin();
   for (casa::Vector<casa::Double>::const_iterator ci = freq.begin(); ci != freq.end(); ++ci,++it) {
        ASKAPDEBUGASSERT(it != result.end());
        *(it++) = fluxDensity(*ci / 1e6);         
        // calibrator is always assumed to be in the phase centre, so the imaginary part is always 0.
        ASKAPDEBUGASSERT(it != result.end());        
        *it = 0.;
   }
}        

/// @brief calculate stokes I visibilities and derivatives for this component
/// @details This variant of the method does simultaneous calculations of
/// the values and derivatives. The component is assumed to have no free parameters, so
/// all derivatives are always zero. This method is used in the implementation
/// of the IComponent interface if stokes I is requested. Otherwise result
/// is filled with 0.
/// @param[in] uvw  baseline spacings (in metres)
/// @param[in] freq vector of frequencies to do calculations for (in Hz)
/// @param[out] result an output buffer used to store values
void Calibrator1934::calculate(const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    std::vector<casa::AutoDiff<double> > &result) const                    
{
   result.resize(2 * freq.nelements());
   std::vector<casa::AutoDiff<double> >::iterator it = result.begin();
   for (casa::Vector<casa::Double>::const_iterator ci = freq.begin(); ci != freq.end(); ++ci,++it) {
        ASKAPDEBUGASSERT(it != result.end());
        *(it++) = casa::AutoDiff<double>(fluxDensity(*ci / 1e6));         
        // calibrator is always assumed to be in the phase centre, so the imaginary part is always 0.
        ASKAPDEBUGASSERT(it != result.end());        
        *it = 0.;
   }
}

/// @brief flux model for 1934-638
/// @details This method estimates the flux density of 1934-638 for a given frequency using the cm-wavelength 
/// model of Reynolds et al. (see miriad or query 1934-638 in the ATCA calibrator database for a
/// reference).   
/// @param[in] freqInMHz frequency of interest (in MHz)
/// @return estimated flux density in Jy
double Calibrator1934::fluxDensity(const double freqInMHz)
{
  ASKAPCHECK( (freqInMHz > 500) && (freqInMHz < 10000), 
      "The flux model of 1934-638 is only valid from 500 MHz to 10 GHz, you have freq = "<<freqInMHz<<" MHz");
  const double lgF = log(freqInMHz) / log(10.);
  // polynomial fit
  const double lgS = -30.7667 + (26.4908 - (7.0977 - 0.6053334 * lgF) * lgF) * lgF;
  return exp(lgS * log(10.));
}  

/// @brief get number of parameters
/// @return a number of parameters  this component depends upon. 
size_t Calibrator1934::nParameters() const throw()
{
  return 0;
}
  
/// @brief get the name of the given parameter
/// @details All parameters are handled in the synthesis code using their
/// string name, which allows to fix or free any of them easily. This method
/// allows to obtain this string name using a integer index
/// @param[in] index an integer index of the parameter (should be less than
/// nParameters).
// @return a const reference to the string name of the parameter 
const std::string& Calibrator1934::parameterName(size_t index) const
{
  ASKAPTHROW(AskapError, "Calibrator1934 class has no parameters defined, index="<<index);
}


} // namespace synethesis

} // namespace askap

