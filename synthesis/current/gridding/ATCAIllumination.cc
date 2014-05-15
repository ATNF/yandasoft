/// @file 
/// @brief ATCA L-band illumination model
/// @details This class represents a ATCA L-band illumination model.
/// It includes a disk with Jamesian illumination and optionally feed leg shadows.
/// The glish scripts written by Tim Cornwell were used as a guide.
/// Optionally a phase slope can be applied to simulate offset pointing.
///
/// @copyright (c) 2008 CSIRO
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

#include <casa/BasicSL/Constants.h>
#include <casa/Arrays/ArrayMath.h>


#include <gridding/ATCAIllumination.h>
#include <askap/AskapError.h>

#include <profile/AskapProfiler.h>

using namespace askap;
using namespace askap::synthesis;

/// @brief construct the model
/// @param[in] diam disk diameter in metres
/// @param[in] blockage a diameter of the central hole in metres
ATCAIllumination::ATCAIllumination(double diam, double blockage) :
  itsDiameter(diam), itsBlockage(blockage), itsDoTapering(false), itsDoFeedLegs(false),
  itsDoFeedLegWedges(false) 
{
  ASKAPDEBUGASSERT(diam>0);
  ASKAPDEBUGASSERT(blockage>=0);  
  ASKAPDEBUGASSERT(diam > blockage);
}

/// @brief switch on the tapering simulation
/// @details This method assigns defocusing phase and switches on the simulation of tapering.
/// @param[in] maxDefocusingPhase the value of the phase in radians at the dish edge, it will
/// be linearly increased with the radius to simulate defocusing
void ATCAIllumination::simulateTapering(double maxDefocusingPhase)
{
  itsDoTapering = true;
  itsMaxDefocusingPhase = maxDefocusingPhase;
}

/// @brief switch on the feed leg simulation
/// @details This method assigns parameters of the feed leg shadows and allows the
/// simulation of feed legs. Calling this method also makes the pattern asymmetric.
/// @param[in] width width in metres of each feed leg shadow
/// @param[in] rotation angle in radians of the feed lag shadows with respect to u,v axes
/// @param[in] shadowingFactor attenuation of the illumination caused by feed legs, assign
/// zero to get a total blockage.
void ATCAIllumination::simulateFeedLegShadows(double width, double rotation, 
                                              double shadowingFactor)
{
  itsDoFeedLegs = true;
  itsFeedLegsHalfWidth = width/2.;
  itsFeedLegsRotation = rotation;
  itsFeedLegsShadowing = shadowingFactor;
  ASKAPCHECK((shadowingFactor<=1.) && (shadowingFactor>=0.), 
	     "shadowingFactor is supposed to be from [0,1] interval, you have "<<shadowingFactor); 
}  

/// @brief switch on the simulation of feed leg wedges
/// @details This method assigns the parameters of the feed leg wedges and allows
/// their simulation. 
/// @param[in] wedgeShadowingFactor1 additional attenuation inside the wedge for the feed
/// leg which is rotated to u axis by the angle specified in simulateFeedLegShadows.
/// @param[in] wedgeShadowingFactor2 the same as wedgeShadowingFactor1, but for orthogonal
/// feed legs
/// @param[in] wedgeOpeningAngle opening angle of the wedge in radians
/// @param[in] wedgeStartingRadius starting radius in metres of the wedge
/// @note simulateFeedLegShadows should also be called prior to the first use of this object.
void ATCAIllumination::simulateFeedLegWedges(double wedgeShadowingFactor1, double wedgeShadowingFactor2, 
					     double wedgeOpeningAngle, double wedgeStartingRadius)
{
  itsDoFeedLegWedges = true;
  itsFeedLegsWedgeShadowing1 = wedgeShadowingFactor1;
  itsFeedLegsWedgeShadowing2 = wedgeShadowingFactor2;
  itsWedgeOpeningAngle = wedgeOpeningAngle;
  itsWedgeStartingRadius = wedgeStartingRadius;
  ASKAPCHECK((wedgeShadowingFactor1<=1.) && (wedgeShadowingFactor1>=0.), 
	     "wedgeShadowingFactor1 is supposed to be from [0,1] interval, you have "<<wedgeShadowingFactor1); 
  ASKAPCHECK((wedgeShadowingFactor2<=1.) && (wedgeShadowingFactor2>=0.), 
	     "wedgeShadowingFactor2 is supposed to be from [0,1] interval, you have "<<wedgeShadowingFactor2); 
  ASKAPCHECK((wedgeOpeningAngle < M_PI/2.) && (wedgeOpeningAngle > 0), 
	     "Opening angle is supposed to be from the (0, pi/2) interval, you have "<<wedgeOpeningAngle);
}   


/// @brief helper method used inside the calculation of Jamesian illumination
/// @details The pattern is modelled as two stiched functions, this method
/// represents a function applied at small radii.
/// @param[in] fractionalRadius radius given as a fraction of dish radius
/// @return amplitude of illumination
double ATCAIllumination::innerJamesian(double fractionalRadius)
{
  // constants were taken from Tim's glish script
  const double dsub = 1.8/11.0;
  return 0.58+0.27*atan(130*(fractionalRadius-1.1*dsub));
}

/// @brief helper method used inside the calculation of Jamesian illumination
/// @details The pattern is modelled as two stiched functions, this method
/// represents a function applied at large radii.
/// @param[in] fractionalRadius radius given as a fraction of dish radius
/// @return amplitude of illumination
double ATCAIllumination::outerJamesian(double fractionalRadius)
{
  // constants were taken from Tim's glish script
  const double a2 = 0.95;
  const double a1 = 8.;
  return 1.-a2*exp(a1*(casa::square(fractionalRadius)-1.));
}

/// @brief a helper method to return 1D illumination
/// @details There are some build in parameters of the taper. We may need to expose them
/// to the class interface in the future, but currently this method doesn't use the 
/// data members and hence made static.
/// @param[in] fractionalRadius radius given as a fraction of dish radius
/// @return amplitude of illumination
double ATCAIllumination::jamesian(double fractionalRadius)
{
  if (fractionalRadius > 1.) {
    return 0.;
  }
  // constants were taken from Tim's glish script
  const double cd=0.3;
  const double transitionFractionalRadius = 0.4;
  
  const double value = (fractionalRadius < transitionFractionalRadius ? 
			innerJamesian(fractionalRadius) / innerJamesian(transitionFractionalRadius) :
			outerJamesian(fractionalRadius) / outerJamesian(transitionFractionalRadius));
  return value * (1. - cd * (1. - fractionalRadius));
}


/// @brief obtain illumination pattern
/// @details This is the main method which populates the 
/// supplied uv-pattern with the values corresponding to the model
/// represented by this object. It has to be overridden in the 
/// derived classes. An optional phase slope can be applied to
/// simulate offset pointing.
/// @param[in] freq frequency in Hz for which an illumination pattern is required
/// @param[in] pattern a UVPattern object to fill
/// @param[in] l angular offset in the u-direction (in radians)
/// @param[in] m angular offset in the v-direction (in radians)
/// @param[in] pa parallactic angle, or strictly speaking the angle between 
/// uv-coordinate system and the system where the pattern is defined (unused)
void ATCAIllumination::getPattern(double freq, UVPattern &pattern, double l, 
				  double m, double pa) const
{
  ASKAPTRACE("ATCAIllumination::getPattern");

  //std::cout<<"ATCAIllumination::getPattern is called for pa = "<<pa/M_PI*180.<<std::endl;
  const casa::uInt oversample = pattern.overSample();
  const double cellU = pattern.uCellSize()/oversample;
  const double cellV = pattern.vCellSize()/oversample;
  
  ASKAPCHECK(itsDoFeedLegWedges ? itsDoFeedLegs : true, 
	     "you can only switch on simulation of feed leg wedges with the simulation of feed legs");
  ASKAPDEBUGASSERT(!itsDoFeedLegWedges || (itsWedgeStartingRadius < itsDiameter/2.));
  ASKAPDEBUGASSERT(!itsDoFeedLegs || (itsFeedLegsHalfWidth < itsDiameter/2.));  
  
  // scaled l and m to take the calculations out of the loop
  // these quantities are effectively dimensionless 
  const double lScaled = 2.*casa::C::pi*cellU *l;
  const double mScaled = 2.*casa::C::pi*cellV *m;
  
  // zero value of the pattern by default
  pattern.pattern().set(0.);
  
  ASKAPCHECK(std::abs(std::abs(cellU/cellV)-1.)<1e-7, 
	     "Rectangular cells are not supported at the moment");
  
  const double cell = std::abs(cellU*(casa::C::c/freq));
  
  const double dishRadiusInCells = itsDiameter/(2.0*cell);  
  
  // squares of the disk and blockage area radii
  const double rMaxSquared = casa::square(dishRadiusInCells);
  const double rMinSquared = casa::square(itsBlockage/(2.0*cell));     
  
  // feed legs/wedges length parameters
  const double feedLegsHalfWidthInCells = itsFeedLegsHalfWidth / cell;
  const double wedgeStartingRadiusInCells = itsWedgeStartingRadius / cell;
  
  // sizes of the grid to fill with pattern values
  const casa::uInt nU = pattern.uSize();
  const casa::uInt nV = pattern.vSize();
  
  ASKAPCHECK((casa::square(double(nU)) > rMaxSquared) || 
	     (casa::square(double(nV)) > rMaxSquared),
	     "The pattern buffer passed to ATCAIllumination::getPattern is too small for the given model. "
	     "Sizes should be greater than "<<sqrt(rMaxSquared)<<" on each axis, you have "
	     <<nU<<" x "<<nV);
  
  // maximum possible support for this class corresponds to the dish size
  pattern.setMaxSupport(1+2*casa::uInt(dishRadiusInCells)/oversample);
  
  double sum=0.; // normalisation factor
  for (casa::uInt iU=0; iU<nU; ++iU) {
    const double offsetU = double(iU)-double(nU)/2.;
    const double offsetUSquared = casa::square(offsetU);
    for (casa::uInt iV=0; iV<nV; ++iV) {
      const double offsetV = double(iV)-double(nV)/2.;
      const double offsetVSquared = casa::square(offsetV);
      const double radiusSquared = offsetUSquared + offsetVSquared;
      if ( (radiusSquared >= rMinSquared) && (radiusSquared <= rMaxSquared)) {
	// don't need to multiply by wavelength here because we
	// divided the radius (i.e. the illumination pattern is given
	// in a relative coordinates in frequency

	// initial value is just the phase slope related to pointing offset
	casa::DComplex value(1.0);
	if (itsDoTapering) {
	  const double fractionalRadiusSquared = radiusSquared/rMaxSquared;
	  const double extraPhase = itsMaxDefocusingPhase * fractionalRadiusSquared;
	  value *= double(jamesian(sqrt(fractionalRadiusSquared))) *
	    casa::DComplex(cos(extraPhase),sin(extraPhase));
	}
	if (itsDoFeedLegs) {
	  // rotated offsets
	  const double dU = cos(itsFeedLegsRotation+pa)*offsetU + 
	    sin(itsFeedLegsRotation+pa)*offsetV;
	  const double dV = -sin(itsFeedLegsRotation+pa)*offsetU + 
	    cos(itsFeedLegsRotation+pa)*offsetV;				                        
	  // decrease the illumination within the feed leg shadows                  
	  if ((std::abs(dU) < feedLegsHalfWidthInCells) ||
	      (std::abs(dV) < feedLegsHalfWidthInCells)) {
	    value *= double(itsFeedLegsShadowing);
	  }
	  if (itsDoFeedLegWedges) {
	    const double wedgeAngleTan = tan(itsWedgeOpeningAngle);
	    if ((dU > wedgeStartingRadiusInCells) && (std::abs(dV) <
						      wedgeAngleTan*std::abs(dU - wedgeStartingRadiusInCells))) {
	      value *= double(itsFeedLegsShadowing * itsFeedLegsWedgeShadowing1);
	    }
	    
	    if ((dU < -wedgeStartingRadiusInCells) && (std::abs(dV) <
						       wedgeAngleTan*std::abs(dU + wedgeStartingRadiusInCells))) {
	      value *= double(itsFeedLegsShadowing * itsFeedLegsWedgeShadowing1);
	    }
	    
	    if ((dV > wedgeStartingRadiusInCells) && (std::abs(dU) <
						      wedgeAngleTan*std::abs(dV - wedgeStartingRadiusInCells))) {
	      value *= double(itsFeedLegsShadowing * itsFeedLegsWedgeShadowing2);
	    }
	    
	    if ((dV < -wedgeStartingRadiusInCells) && (std::abs(dU) <
						       wedgeAngleTan*std::abs(dV + wedgeStartingRadiusInCells))) {
	      value *= double(itsFeedLegsShadowing * itsFeedLegsWedgeShadowing2);
	    }                                                    
	  }
	}
	
	
	// add up the norm contribution and assign the final value
	//sum += casa::abs(value);
	sum += std::real(std::conj(value)*value);
	// Now apply the phase slope related to pointing offset
	const double phase = lScaled*offsetU + mScaled*offsetV;
	pattern(iU, iV) = value*casa::DComplex(cos(phase), -sin(phase));
      }
    }
  }
  
  // We need to divide by the sum squared so that the convolution is normalized to unity
  ASKAPCHECK(sum > 0., "Integral of the pattern should be non-zero");
  sum = sqrt(sum);
  pattern.pattern() *= casa::DComplex(double(nU)*double(nV)/double(sum),0.);
}

/// @brief check whether the pattern is symmetric
/// @details Some illumination patterns are known a priori to be symmetric.
/// This method returns true if feed legs are not simulated to reflect this
/// @return true if feed legs are not simulated
bool ATCAIllumination::isSymmetric() const
{
  return !itsDoFeedLegs;
}
