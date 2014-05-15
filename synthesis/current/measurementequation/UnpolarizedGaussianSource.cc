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


#include <measurementequation/UnpolarizedGaussianSource.h>

#include <scimath/Mathematics/RigidVector.h>

namespace askap {

namespace synthesis {

/// @brief construct the point source component
/// @details 
/// @param[in] name a name of the component. Will be added to all parameter
///            names (e.g. after direction.ra) 
/// @param[in] flux flux density in Jy
/// @param[in] ra offset in right ascension w.r.t. the current phase 
/// centre (in radians)
/// @param[in] dec offset in declination w.r.t. the current phase
/// centre (in radians)
/// @param[in] maj major axis in radians
/// @param[in] min minor axis in radians
/// @param[in] pa  position angle in radians
UnpolarizedGaussianSource::UnpolarizedGaussianSource(const std::string &name,
        double flux, double ra, 
        double dec, double maj, double min, double pa)  : 
          UnpolarizedComponent<6>(casa::RigidVector<double, 6>()) 
{
  casa::RigidVector<double, 6> &params = parameters();
  params(0)=flux;
  params(1)=ra;
  params(2)=dec;
  params(3)=maj;
  params(4)=min;
  params(5)=pa;
  
  casa::RigidVector<std::string, 6> &names = parameterNames();
  const char *nameTemplates[] = {"flux.i","direction.ra","direction.dec",
                 "shape.bmaj","shape.bmin","shape.bpa"};
  for (size_t i=0;i<6;++i) {
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
                    const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
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
void UnpolarizedGaussianSource::calculate(const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    std::vector<casa::AutoDiff<double> > &result) const
{
  const casa::RigidVector<double, 6> &params = parameters();
  casa::RigidVector<casa::AutoDiff<double>, 6>  paramsAutoDiff;
  for (casa::uInt i=0; i<6; ++i) {
       paramsAutoDiff(i)=casa::AutoDiff<double>(params(i),6, i);
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
                    const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    const casa::RigidVector<T, 6> &params,
                    std::vector<T> &result)
{
  const T ra=params(1);
  const T dec=params(2);
  const T flux=params(0);
  const T bmaj=params(3);
  const T bmin=params(4);
  const T bpa=params(5);
  const T n =  casa::sqrt(T(1.0) - (ra*ra+dec*dec));
  const T delay = casa::C::_2pi * (ra * uvw(0) + dec * uvw(1) + 
                                   (n-T(1.0)) * uvw(2))/casa::C::c;
  // exp(-a*x^2) transforms to exp(-pi^2*u^2/a)
  // a=4log(2)/FWHM^2 so scaling = pi^2*FWHM/(4log(2))
  const T scale = std::pow(casa::C::pi,2)/(4*log(2.0));
  const T up=( cos(bpa)*uvw(0) + sin(bpa)*uvw(1))/casa::C::c;
  const T vp=(-sin(bpa)*uvw(0) + cos(bpa)*uvw(1))/casa::C::c;
  const T r=(bmaj*bmaj*up*up+bmin*bmin*vp*vp)*scale;
  
  typename std::vector<T>::iterator it=result.begin();
  for (casa::Vector<casa::Double>::const_iterator ci=freq.begin(); 
       ci!=freq.end();++ci,++it)
      {
        const casa::Double currentFreq = *ci;
        const T phase = delay * currentFreq;
        const T decorr = exp( - r * currentFreq * currentFreq);
        *it = flux * decorr * cos(phase);
        *(++it) = flux *decorr * sin(phase);
      }

}

} // namespace askap
 
} // namespace synthesis

