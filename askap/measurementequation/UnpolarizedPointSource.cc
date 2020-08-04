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

#include <askap/measurementequation/UnpolarizedPointSource.h>

#include <casacore/scimath/Mathematics/RigidVector.h>

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
                    const casacore::RigidVector<casacore::Double, 3> &uvw,
                    const casacore::Vector<casacore::Double> &freq,
                    const casacore::RigidVector<T, 5> &params,
                    std::vector<T> &result)
{
  const T flux0=params(0);
  const T ra=params(1);
  const T dec=params(2);
  const T spectral_index=params(3);
  const T ref_freq=params(4);
  
  const T n =  casacore::sqrt(T(1.0) - (ra*ra+dec*dec));
  const T delay = casacore::C::_2pi * (ra * uvw(0) + dec * uvw(1) + 
                                   (n-T(1.0)) * uvw(2))/casacore::C::c;
  typename std::vector<T>::iterator it=result.begin();
  for (casacore::Vector<casacore::Double>::const_iterator ci=freq.begin(); 
       ci!=freq.end();++ci,++it)
      {
        const casacore::Double f = *ci;
        // cannot use the spectral_index==0 conditional because templated type casacore::AutoDiff doesn't like it
        //const T flux = spectral_index==0 ? flux0 : flux0 * pow(f/ref_freq,spectral_index);
        const T flux = flux0 * pow(f/ref_freq,spectral_index);
        const T phase = delay * (*ci);
        *it = flux * cos(phase)/n;
        *(++it) = flux * sin(phase)/n;
      }
}

/// @brief construct the point source component
/// @details 
/// @param[in] name a name of the component. Will be added to all parameter
///            names (e.g. after direction.ra) 
/// @param[in] flux flux density in Jy at ref_freq
/// @param[in] ra offset in right ascension w.r.t. the current phase 
/// centre (in radians)
/// @param[in] dec offset in declination w.r.t. the current phase
/// centre (in radians)
/// @param[in] spectral_index spectral index in Hz
/// @param[in] ref_freq referece frequency for parameter "flux"
UnpolarizedPointSource::UnpolarizedPointSource(const std::string &name, 
          double flux, double ra, double dec, double spectral_index, double ref_freq) : 
          UnpolarizedComponent<5>(casacore::RigidVector<double, 5>(flux,ra,dec,spectral_index,ref_freq)) 
{
  parameterNames() = casacore::RigidVector<std::string, 5>("flux.i"+name,
            "direction.ra"+name, "direction.dec"+name,
            "flux.spectral_index"+name, "flux.ref_freq"+name);
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
                    const casacore::RigidVector<casacore::Double, 3> &uvw,
                    const casacore::Vector<casacore::Double> &freq,
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
                    const casacore::RigidVector<casacore::Double, 3> &uvw,
                    const casacore::Vector<casacore::Double> &freq,
                    std::vector<casacore::AutoDiff<double> > &result) const
{
  const casacore::RigidVector<double, 5> &params = parameters();
  const casacore::RigidVector<casacore::AutoDiff<double>, 5>  paramsAutoDiff(
                              casacore::AutoDiff<double>(params(0),3, 0),
                              casacore::AutoDiff<double>(params(1),3, 1),
                              casacore::AutoDiff<double>(params(2),3, 2),
                              casacore::AutoDiff<double>(params(1),3, 3),
                              casacore::AutoDiff<double>(params(2),3, 4));
  calcPoint(uvw,freq,paramsAutoDiff,result);
}

} // namespace askap

} // namespace synthesis
              
