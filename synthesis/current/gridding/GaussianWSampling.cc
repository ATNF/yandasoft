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

#include <gridding/GaussianWSampling.h>
#include <askap/AskapError.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.gaussianwsampling");


#include <cmath>

using namespace askap;
using namespace askap::synthesis;

/// @brief initialise the class
/// @details
/// @param[in] wplanes50 the fraction of w-planes covering 50% of the w-term range [-wmax,wmax]. The first and the last
/// w-planes always correspond to -wmax and +wmax, while the mid-plane always corresponds to zero w-term.
GaussianWSampling::GaussianWSampling(const double wplanes50) : itsTwoSigmaSquared(1.), itsAmplitude(0.) 
{
  calculateDistributionParameters(wplanes50);
}

/// @brief plane to w-term conversion (mapping)
/// @details This is a forward method mapping scaled w-plane to scaled w-term.
/// @param[in] plane plane number scaled down to interval [-1;1]
/// @return w-term scaled down to interval [-1;1]
/// @note The result is unpredictable, if the plane is outside [-1;1] interval
double GaussianWSampling::map(double plane) const
{
  ASKAPDEBUGASSERT((plane>=-1.) && (plane<=1.));
  // y = sign(x)*A*(1-exp(-x*x/(2*sigma*sigma)))
  const double absval = itsAmplitude*(1-exp(-plane*plane/itsTwoSigmaSquared));
  if (plane<0) {
      return -absval;
  } else if (plane>0) {
      return absval;
  } 
  return 0.;
}

/// @brief w-term to plane conversion (indexing)
/// @details This is a reverse method indexing dimensionless w-tern to get dimensionless w-plane.
/// @param[in] wterm scaled down to interval [-1;1]
/// @return w-term scaled down to interval [-1;1]
/// @note The result is unpredictable, if the input w-term is outside [-1;1] interval
double GaussianWSampling::index(double wterm) const
{
  ASKAPDEBUGASSERT((wterm>=-1.) && (wterm<=1.));
  ASKAPDEBUGASSERT(itsAmplitude!=0.);
  if ((wterm == -1.) || (wterm == 1.)) {
      return wterm;
  }
  const double expterm = 1.-fabs(wterm)/itsAmplitude;
  ASKAPCHECK(expterm>0., "Unable to invert the formula for w-term sampling, wterm = "<<wterm<<
             " implies infinite plane number; expterm="<<expterm);
  const double planeSquared = -log(expterm)*itsTwoSigmaSquared;
  ASKAPCHECK(planeSquared>=0., "Unable to invert the formula for w-term sampling, wterm = "<<wterm<<
             " implies a negative planeSquared = "<<planeSquared);
  if (wterm<0) {
      return -sqrt(planeSquared);
  } else if (wterm>0) {
      return sqrt(planeSquared);
  }             
  return 0.;
}

/// @brief derive distribution paramters
/// @details This method has been designed to calculate distribution parameters (sigma and amplitude) stored
/// as members of this class from the input parameter being the number of w-planes containing 50% of the w-term
/// range [-wmax,wmax] or to be exact [-1,1] as this class works with the normalised w-term. 
/// @param[in] wplanes50 the fraction of w-planes covering 50% of the w-term range [-wmax,wmax]. The first and the last
/// w-planes always correspond to -wmax and +wmax, while the mid-plane always corresponds to zero w-term.
void GaussianWSampling::calculateDistributionParameters(const double wplanes50)
{
  // we have to solve the equation
  // x=sqrt(-twosigmasquared*log(0.5+0.5*exp(-1/twosigmasquared)))
  // numerically to get twosigmasquared from x=wplanes50 
  // This equation has an interesting asymptotic behavior at both large and small sigmas. 
  // This asymptotic behavior can be used to derive search range for numerical solution of the equation
  // sigma -> infinity, x -> 1/sqrt(2)
  // sigma -> 0, x/sigma -> sqrt(2*log(2))
  // Two estimates for sigma are (x/sqrt(2*log(2)), sqrt(0.05/(1/sqrt(2)-x))) or directly for twosigmasquared are
  // (x*x/log(2), 0.1/(1/sqrt(2)-x)), the first estimate is using small sigma approximation, the second using large
  // sigma approximation.
  
  ASKAPCHECK(wplanes50>0, "Fraction of w-planes containing 50% of -wmax,wmax range of w-terms should be positive, you have normalised wplanes50="
                          << wplanes50);
  ASKAPCHECK(wplanes50*sqrt(2)<1., "Normalised fraction of w-planes containing 50% of -wmax,wmax range of wterms should not exceed 1/sqrt(2)."<< 
         " Otherwise, the solution for gaussian distribution does not exist. You have normalised wplanes50="<<wplanes50);

  double minTwoSigmaSquared = wplanes50*wplanes50 / log(2.) / 2.; // scaled down with factor of 2, just in case
  double maxTwoSigmaSquared = 0.2/(1/sqrt(2) - wplanes50); // scaled up with factor of 2, just in case
  if (maxTwoSigmaSquared < minTwoSigmaSquared) {
      minTwoSigmaSquared = maxTwoSigmaSquared;
      maxTwoSigmaSquared = 5.; // value well outside the ambiguity region
  }
  ASKAPLOG_INFO_STR(logger, "Gaussian w-sampling: searching for twoSigmaSquared in ("<<minTwoSigmaSquared<<","<<
                            maxTwoSigmaSquared<<
                            ") interval for normalised fraction of w-planes containing 50% of -wmax,wmax range of w-terms equal to "<<
                            wplanes50);
  // first verify that interval selection is correct
  itsTwoSigmaSquared = minTwoSigmaSquared;
  ASKAPCHECK(targetFunction() <= wplanes50, "lower end of the search interval corresponds to higher fraction of w-planes. This is a bug!");
  itsTwoSigmaSquared = maxTwoSigmaSquared;
  ASKAPCHECK(targetFunction() >= wplanes50, "lower end of the search interval corresponds to lower fraction of w-planes. This is a bug!");  
  // now do dichotomy                         
  size_t maxNiter = 100;
  const double tolerance = 1e-5;
  for (size_t it = 0; it < maxNiter; ++it) {
       itsTwoSigmaSquared = (minTwoSigmaSquared + maxTwoSigmaSquared) / 2.;
       if (targetFunction() < wplanes50) {
           minTwoSigmaSquared = itsTwoSigmaSquared;
       } else {
           maxTwoSigmaSquared = itsTwoSigmaSquared;
       }
       if (maxTwoSigmaSquared - minTwoSigmaSquared < tolerance) {
           ASKAPLOG_INFO_STR(logger, "Search converged at iteration = "<<it<<" itsTwoSigmaSquared="<<itsTwoSigmaSquared);
           break;
       }
       ASKAPCHECK(it + 1 != maxNiter, "Failed to converge within "<<maxNiter<<
                  " iterations with tolerance="<<tolerance<<", please investigate");
  }    
  itsAmplitude = 1./(1.-exp(-1./itsTwoSigmaSquared));
}

