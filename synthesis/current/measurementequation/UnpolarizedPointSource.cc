/// @file
///
/// @brief A component representing an unpolarized point source with flat spectrum
/// @details
///     This is an implementation of IComponent for the point source model.
///     The point source is assumed unpolarized with a flat spectrum
///     (i.e. spectral index 0).
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
/// 

#include <measurementequation/UnpolarizedPointSource.h>

#include <scimath/Mathematics/RigidVector.h>

namespace askap {

namespace synthesis {


/// @brief actual calculations
/// @details templated method for actual calculations. Made private, because
/// it is used and declared in UnpolarizedPointSource.cc only
/// @param[in] uvw  baseline spacings (in metres)
/// @param[in] freq vector of frequencies to do calculations for
/// @param[in] params RigidVector with parameters
/// @param[out] result an output buffer used to store values
template<typename T>
void UnpolarizedPointSource::calcPoint(
                    const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    const casa::RigidVector<T, 3> &params,
                    std::vector<T> &result)
{
  const T ra=params(1);
  const T dec=params(2);
  const T flux=params(0);
  
  const T n =  casa::sqrt(T(1.0) - (ra*ra+dec*dec));
  const T delay = casa::C::_2pi * (ra * uvw(0) + dec * uvw(1) + 
                                   (n-T(1.0)) * uvw(2))/casa::C::c;
  typename std::vector<T>::iterator it=result.begin();
  for (casa::Vector<casa::Double>::const_iterator ci=freq.begin(); 
       ci!=freq.end();++ci,++it)
      {
        const T phase = delay * (*ci);
        *it = flux * cos(phase)/n;
        *(++it) = flux * sin(phase)/n;
      }
}

/// @brief construct the point source component
/// @details 
/// @param[in] name a name of the component. Will be added to all parameter
///            names (e.g. after direction.ra) 
/// @param[in] flux flux density in Jy
/// @param[in] ra offset in right ascension w.r.t. the current phase 
/// centre (in radians)
/// @param[in] dec offset in declination w.r.t. the current phase
/// centre (in radians)
UnpolarizedPointSource::UnpolarizedPointSource(const std::string &name, 
          double flux, double ra, double dec) : 
          UnpolarizedComponent<3>(casa::RigidVector<double, 3>(flux,ra,dec)) 
{
  parameterNames() = casa::RigidVector<std::string, 3>("flux.i"+name,
            "direction.ra"+name, "direction.dec"+name);
}

              
/// @brief calculate stokes I visibilities for this component
/// @details This variant of the method calculates just the visibilities
/// (without derivatives) for a number of frequencies. This method is 
/// used to in the implementation of the IComponent interface if 
/// stokes I is requested. Otherwise result is filled with 0.
/// @param[in] uvw  baseline spacings (in metres)
/// @param[in] freq vector of frequencies to do calculations for
/// @param[out] result an output buffer used to store values
void UnpolarizedPointSource::calculate(
                    const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    std::vector<double> &result) const
{
  calcPoint(uvw,freq,parameters(),result);
}                    
  

/// @brief calculate stokes I visibilities and derivatives for this component
/// @details This variant of the method does simultaneous calculations of
/// the values and derivatives. This method is used to in the implementation
/// of the IComponent interface if stokes I is requested. Otherwise result
/// is filled with 0.
/// @param[in] uvw  baseline spacings (in metres)
/// @param[in] freq vector of frequencies to do calculations for
/// @param[out] result an output buffer used to store values
void UnpolarizedPointSource::calculate(
                    const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    std::vector<casa::AutoDiff<double> > &result) const
{
  const casa::RigidVector<double, 3> &params = parameters();
  const casa::RigidVector<casa::AutoDiff<double>, 3>  paramsAutoDiff(
                              casa::AutoDiff<double>(params(0),3, 0),
                              casa::AutoDiff<double>(params(1),3, 1),
                              casa::AutoDiff<double>(params(2),3, 2));
  calcPoint(uvw,freq,paramsAutoDiff,result);
}

} // namespace askap

} // namespace synthesis
              
