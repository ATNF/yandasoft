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
//template<typename T>
void UnpolarizedPointSource::calcPoint(
                    const casacore::RigidVector<casacore::Double, 3> &uvw,
                    const casacore::Vector<casacore::Double> &freq,
                    const casacore::RigidVector<casacore::AutoDiff<double>, 5> &params,
                    std::vector<casacore::AutoDiff<double>> &result)
{
  const casacore::AutoDiff<double> flux0=params(0);
  const casacore::AutoDiff<double> ra=params(1);
  const casacore::AutoDiff<double> dec=params(2);
  const casacore::AutoDiff<double> spectral_index=params(3);
  const casacore::AutoDiff<double> ref_freq=params(4);

  const casacore::AutoDiff<double> n =  casacore::sqrt(casacore::AutoDiff<double>(1.0) - (ra*ra+dec*dec));
  const casacore::AutoDiff<double> delay = casacore::C::_2pi * (ra * uvw(0) + dec * uvw(1) +
                                   (n-casacore::AutoDiff<double>(1.0)) * uvw(2))/casacore::C::c;
  typename std::vector<casacore::AutoDiff<double>>::iterator it=result.begin();
  for (casacore::Vector<casacore::Double>::const_iterator ci=freq.begin();
       ci!=freq.end();++ci,++it)
      {
        const casacore::Double f = *ci;
        casacore::AutoDiff<double> flux = flux0;
        if (spectral_index != casacore::AutoDiff<double>(0)) {
          flux *= pow(f/ref_freq,spectral_index);
        }
        const casacore::AutoDiff<double> phase = delay * (*ci);
        *it = flux * cos(phase);
        *(++it) = flux * sin(phase);
      }
}
void UnpolarizedPointSource::calcPoint(
                    const casacore::RigidVector<casacore::Double, 3> &uvw,
                    const casacore::Vector<casacore::Double> &freq,
                    const casacore::RigidVector<double, 5> &params,
                    std::vector<double> &result)
{
  const double flux0=params(0);
  const double ra=params(1);
  const double dec=params(2);
  const double spectral_index=params(3);
  const double ref_freq=params(4);

  const double n =  casacore::sqrt(double(1.0) - (ra*ra+dec*dec));
  const double delay = casacore::C::_2pi * (ra * uvw(0) + dec * uvw(1) +
                                   (n-double(1.0)) * uvw(2))/casacore::C::c;
  typename std::vector<double>::iterator it=result.begin();
  for (casacore::Vector<casacore::Double>::const_iterator ci=freq.begin();
       ci!=freq.end();++ci,++it)
      {
        const casacore::Double f = *ci;
        double flux = flux0;
        if (spectral_index != 0) {
          flux *= pow(f/ref_freq,spectral_index);
        }
        const double phase = delay * (*ci);
        *it = flux * cos(phase);
        *(++it) = flux * sin(phase);
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
  casacore::RigidVector<casacore::AutoDiff<double>, 5>  paramsAutoDiff;
  for (casacore::uInt i=0; i<5; ++i) {
       paramsAutoDiff(i)=casacore::AutoDiff<double>(params(i), 5, i);
  }
  calcPoint(uvw,freq,paramsAutoDiff,result);
}

} // namespace askap

} // namespace synthesis
