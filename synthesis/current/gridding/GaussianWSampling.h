/// @file
///
/// @brief Gaussian w-sampling 
/// @details
/// W-dependent gridders support non-linear sampling in the w-space (through WDependentGridderBase).
/// This class implements IWSampling interface and provides gaussian sampling in the w-space. The
/// class is parameterised with a single parameter being the number of w-planes covering 50% of w-term range.
/// Other parameters of the distribution formula 
///
///     y = sign(x)*A*(1-exp(-x*x/(2*sigma*sigma)))
///
/// are derived from this single parameter under assumption that the whole [-wmax,wmax] interval should be sampled,
/// so the first and the last w-planes should always correspond to -wmax and +wmax, and the middle w-plane should always
/// correspond to zero w-term. The gaussian w-sampling may be helpful if we take into account a typical density of 
/// samples in w-space. Some experimentation is needed to find what values of the free parameter are actually useful.
///
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

#ifndef GAUSSIAN_W_SAMPLING_H
#define GAUSSIAN_W_SAMPLING_H

#include <gridding/IWSampling.h>
#include <cmath>

namespace askap {

namespace synthesis {

/// @brief Gaussian w-sampling 
/// @details
/// W-dependent gridders support non-linear sampling in the w-space (through WDependentGridderBase).
/// This class implements IWSampling interface and provides gaussian sampling in the w-space. The
/// class is parameterised with a single parameter being the number of w-planes covering 50% of w-term range.
/// Other parameters of the distribution formula 
///
///     y = sign(x)*A*(1-exp(-x*x/(2*sigma*sigma)))
///
/// are derived from this single parameter under assumption that the whole [-wmax,wmax] interval should be sampled,
/// so the first and the last w-planes should always correspond to -wmax and +wmax, and the middle w-plane should always
/// correspond to zero w-term. The gaussian w-sampling may be helpful if we take into account a typical density of 
/// samples in w-space. Some experimentation is needed to find what values of the free parameter are actually useful.
/// @ingroup gridding
struct GaussianWSampling : public IWSampling {

  /// @brief initialise the class
  /// @details
  /// @param[in] wplanes50 the fraction of w-planes covering 50% of the w-term range [-wmax,wmax]. The first and the last
  /// w-planes always correspond to -wmax and +wmax, while the mid-plane always corresponds to zero w-term.
  explicit GaussianWSampling(const double wplanes50);
    
  /// @brief plane to w-term conversion (mapping)
  /// @details This is a forward method mapping scaled w-plane to scaled w-term.
  /// @param[in] plane plane number scaled down to interval [-1;1]
  /// @return w-term scaled down to interval [-1;1]
  /// @note The result is unpredictable, if the plane is outside [-1;1] interval
  virtual double map(double plane) const;

  /// @brief w-term to plane conversion (indexing)
  /// @details This is a reverse method indexing dimensionless w-tern to get dimensionless w-plane.
  /// @param[in] wterm scaled down to interval [-1;1]
  /// @return w-term scaled down to interval [-1;1]
  /// @note The result is unpredictable, if the input w-term is outside [-1;1] interval
  virtual double index(double wterm) const;  

protected:

  /// @brief derive distribution paramters
  /// @details This method has been designed to calculate distribution parameters (sigma and amplitude) stored
  /// as members of this class from the input parameter being the number of w-planes containing 50% of the w-term
  /// range [-wmax,wmax] or to be exact [-1,1] as this class works with the normalised w-term. 
  /// @param[in] wplanes50 the fraction of w-planes covering 50% of the w-term range [-wmax,wmax]. The first and the last
  /// w-planes always correspond to -wmax and +wmax, while the mid-plane always corresponds to zero w-term.
  void calculateDistributionParameters(const double wplanes50);    

  /// @brief target function to solve equation numerically
  /// @details This method is called from inside the dichotomy algorithm, it uses 
  /// itsTwoSigmaSquared as an input and does not have any arguments.
  /// @return value of the target function (equal to wplanes50/2 for the match)
  inline double targetFunction() const 
    { return sqrt(-itsTwoSigmaSquared*log(0.5+0.5*exp(-1./itsTwoSigmaSquared))); }
  
private:
  /// @brief two sigma squared (sigma is the parameter of the Gaussian)
  double itsTwoSigmaSquared;
  
  /// @brief amplitude of the Gaussian
  double itsAmplitude;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef GAUSSIAN_W_SAMPLING_H

