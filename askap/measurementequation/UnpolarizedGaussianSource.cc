/// @file
///
/// @brief A component representing an unpolarized gaussian source with the
//  flat spectrum
/// @details
///     This is an implementation of IComponent for the gaussian source model.
///     The source is assumed unpolarized with a flat spectrum
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


#include <askap/measurementequation/UnpolarizedGaussianSource.h>

#include <casacore/scimath/Mathematics/RigidVector.h>

namespace askap {

namespace synthesis {

/// @brief construct the point source component
/// @details
/// @param[in] name a name of the component. Will be added to all parameter
///            names (e.g. after direction.ra)
/// @param[in] flux flux density in Jy at ref_freq
/// @param[in] ra offset in right ascension w.r.t. the current phase
/// centre (in radians)
/// @param[in] dec offset in declination w.r.t. the current phase
/// centre (in radians)
/// @param[in] maj major axis in radians
/// @param[in] min minor axis in radians
/// @param[in] pa  position angle in radians
/// @param[in] spectral_index spectral index in Hz
/// @param[in] ref_freq referece frequency for parameter "flux"
UnpolarizedGaussianSource::UnpolarizedGaussianSource(const std::string &name,
        double flux, double ra,  double dec,
        double maj, double min, double pa, double spectral_index, double ref_freq)  :
          UnpolarizedComponent<8>(casacore::RigidVector<double, 8>())
{
  casacore::RigidVector<double, 8> &params = parameters();
  params(0)=flux;
  params(1)=ra;
  params(2)=dec;
  params(3)=maj;
  params(4)=min;
  params(5)=pa;
  params(6)=spectral_index;
  params(7)=ref_freq;

  casacore::RigidVector<std::string, 8> &names = parameterNames();
 const char *nameTemplates[] = {"flux.i","direction.ra","direction.dec",
                 "shape.bmaj","shape.bmin","shape.bpa",
                 "flux.spectral_index","flux.ref_freq"};
  for (size_t i=0;i<8;++i) {
       names(i)=std::string(nameTemplates[i])+name;
  }
}

/// @brief calculate stokes I visibilities for this component
/// @details This variant of the method calculates just the visibilities
/// (without derivatives) for a number of frequencies. This method is
/// used to in the implementation of the IComponent interface if
/// stokes I is requested. Otherwise result is filled with 0.
/// @param[in] uvw  baseline spacings (in metres)
/// @param[in] freq vector of frequencies to do calculations for
/// @param[out] result an output buffer used to store values
void UnpolarizedGaussianSource::calculate(
                    const casacore::RigidVector<casacore::Double, 3> &uvw,
                    const casacore::Vector<casacore::Double> &freq,
                    std::vector<double> &result) const
{
  calcGaussian(uvw,freq,parameters(),result);
}

/// @brief calculate stokes I visibilities and derivatives for this component
/// @details This variant of the method does simultaneous calculations of
/// the values and derivatives. This method is used to in the implementation
/// of the IComponent interface if stokes I is requested. Otherwise result
/// is filled with 0.
/// @param[in] uvw  baseline spacings (in metres)
/// @param[in] freq vector of frequencies to do calculations for
/// @param[out] result an output buffer used to store values
void UnpolarizedGaussianSource::calculate(const casacore::RigidVector<casacore::Double, 3> &uvw,
                    const casacore::Vector<casacore::Double> &freq,
                    std::vector<casacore::AutoDiff<double> > &result) const
{
  const casacore::RigidVector<double, 8> &params = parameters();
  casacore::RigidVector<casacore::AutoDiff<double>, 8>  paramsAutoDiff;
  for (casacore::uInt i=0; i<8; ++i) {
       paramsAutoDiff(i)=casacore::AutoDiff<double>(params(i), 8, i);
  }
  calcGaussian(uvw,freq,paramsAutoDiff,result);
}

/// @brief actual calculations
/// @details templated method for actual calculations. Made private, because
/// it is used and declared in UnpolarizedGaussianSource.cc only
/// @param[in] uvw  baseline spacings (in metres)
/// @param[in] freq vector of frequencies to do calculations for
/// @param[in] params RigidVector with parameters
/// @param[out] result an output buffer used to store values
template<typename T>
void UnpolarizedGaussianSource::calcGaussian(
                    const casacore::RigidVector<casacore::Double, 3> &uvw,
                    const casacore::Vector<casacore::Double> &freq,
                    const casacore::RigidVector<T, 8> &params,
                    std::vector<T> &result)
{
  const T ra=params(1);
  const T dec=params(2);
  const T flux0=params(0);
  const T bmaj=params(3);
  const T bmin=params(4);
  const T bpa=params(5);
  const T spectral_index=params(6);
  const T ref_freq=params(7);
  const T n =  casacore::sqrt(T(1.0) - (ra*ra+dec*dec));
  const T delay = casacore::C::_2pi * (ra * uvw(0) + dec * uvw(1) +
                                   (n-T(1.0)) * uvw(2))/casacore::C::c;
  // exp(-a*x^2) transforms to exp(-pi^2*u^2/a)
  // a=4log(2)/FWHM^2 so scaling = pi^2*FWHM/(4log(2))
  const T scale = std::pow(casacore::C::pi,2)/(4*log(2.0));
  const T up=( cos(bpa)*uvw(0) + sin(bpa)*uvw(1))/casacore::C::c;
  const T vp=(-sin(bpa)*uvw(0) + cos(bpa)*uvw(1))/casacore::C::c;
  const T r=(bmaj*bmaj*up*up+bmin*bmin*vp*vp)*scale;

  typename std::vector<T>::iterator it=result.begin();
  for (casacore::Vector<casacore::Double>::const_iterator ci=freq.begin();
       ci!=freq.end();++ci,++it)
      {
        const casacore::Double currentFreq = *ci;
        T flux = flux0;
        if (spectral_index != T(0)) {
          flux *= pow(currentFreq/ref_freq,spectral_index);
        }
        const T phase = delay * currentFreq;
        const T decorr = exp( - r * currentFreq * currentFreq);
        *it = flux * decorr * cos(phase);
        *(++it) = flux *decorr * sin(phase);
      }

}

} // namespace askap

} // namespace synthesis
