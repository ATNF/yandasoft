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

#ifndef ATCA_ILLUMINATION_H
#define ATCA_ILLUMINATION_H

#include <gridding/IBasicIllumination.h>

namespace askap {

namespace synthesis {

/// @brief ATCA L-band illumination model
/// @details This class represents a ATCA L-band illumination model.
/// It includes a disk with Jamesian illumination and optionally feed leg shadows.
/// The glish scripts written by Tim Cornwell were used as a guide.
/// Optionally a phase slope can be applied to simulate offset pointing.
/// @ingroup gridding
struct ATCAIllumination : virtual public IBasicIllumination {

  /// @brief construct the model
  /// @param[in] diam antenna diameter in metres
  /// @param[in] blockage diameter of the central hole in metres
  /// @note by default no taper and no feed leg shadows are simulated. Use
  /// methods like ... to switch these effects off. Default behavior is equivalent
  /// to a disk model.
  ATCAIllumination(double diam, double blockage);
    
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
  virtual void getPattern(double freq, UVPattern &pattern, double l = 0., 
                          double m = 0., double pa = 0.) const;

  /// @brief check whether the pattern is symmetric
  /// @details Some illumination patterns are known a priori to be symmetric.
  /// This method returns true if feed legs are not simulated to reflect this
  /// @return true if feed legs are not simulated
  virtual bool isSymmetric() const;

  // class-specific methods to configure other modes
  // can add methods to switch off these modes if there is a use case
  
  /// @brief switch on the tapering simulation
  /// @details This method assigns defocusing phase and switches on the simulation of tapering.
  /// @param[in] maxDefocusingPhase the value of the phase in radians at the dish edge, it will
  /// be quadraticly increased with the radius to simulate defocusing
  void simulateTapering(double maxDefocusingPhase);
  
  /// @brief switch on the feed leg simulation
  /// @details This method assigns parameters of the feed leg shadows and allows the
  /// simulation of feed legs. Calling this method also makes the pattern asymmetric.
  /// @param[in] width width in metres of each feed leg shadow
  /// @param[in] rotation angle in radians of the feed lag shadows with respect to u,v axes
  /// @param[in] shadowingFactor attenuation of the illumination caused by feed legs, assign
  /// zero to get a total blockage.
  void simulateFeedLegShadows(double width, double rotation, double shadowingFactor);
  
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
  void simulateFeedLegWedges(double wedgeShadowingFactor1, double wedgeShadowingFactor2, 
        double wedgeOpeningAngle, double wedgeStartingRadius);
        
protected:
  /// @brief a helper method to return 1D illumination
  /// @details There are some build in parameters of the taper. We may need to expose them
  /// to the class interface in the future, but currently this method doesn't use the 
  /// data members and hence made static.
  /// @param[in] fractionalRadius radius given as a fraction of dish radius
  /// @return amplitude of illumination
  static double jamesian(double fractionalRadius);

  /// @brief helper method used inside the calculation of Jamesian illumination
  /// @details The pattern is modelled as two stiched functions, this method
  /// represents a function applied at small radii.
  /// @param[in] fractionalRadius radius given as a fraction of dish radius
  /// @return amplitude of illumination
  static double innerJamesian(double fractionalRadius);

  /// @brief helper method used inside the calculation of Jamesian illumination
  /// @details The pattern is modelled as two stiched functions, this method
  /// represents a function applied at large radii.
  /// @param[in] fractionalRadius radius given as a fraction of dish radius
  /// @return amplitude of illumination
  static double outerJamesian(double fractionalRadius);
    
private:
  /// @brief antenna diameter in metres
  double itsDiameter;
  
  /// @brief diameter of the central hole in metres
  double itsBlockage;
  
  /// @brief Apply a jamesian taper, if true
  /// @details If this flag is false, a uniformly illuminated disk is simulated
  bool itsDoTapering;
  
  /// @brief Make feed legs shadows, if true
  /// @details If this flag is false no feed leg shadows are simulated
  bool itsDoFeedLegs;
  
  /// @brief Simulate feed legs wedge, if true
  /// @details If this flag is false, no feed leg wedges are simulated
  bool itsDoFeedLegWedges;
  
  /// @brief phase slope with radius
  /// @details If the taper is applied, the code allows also to apply a quadratic phase
  /// slope with radius to simulate defocusing. This parameter is the value of phase
  /// in radians at the dish edge. It is not used if itsDoTaper is false. Assign zero 
  /// to switch the feature off.
  double itsMaxDefocusingPhase;
  
  /// @brief half width of the feed leg shadows
  /// @details This parameter describes how wide are the shadows from feed legs.
  /// It is not used if itsDoFeedLegs is false.
  double itsFeedLegsHalfWidth;
  
  /// @brief angle in radians of the feed legs with respect to u and v coordinate system
  /// @details Feed leg shadows need not to be aligned with the horizontal and vertical axes.
  /// This parameter controls the rotation of the shadow pattern.
  double itsFeedLegsRotation;
  
  /// @brief shadowing factor for the feed leg shadows (fraction of the attenuation)
  /// @details Due to the diffraction there is no total blockage behind the feed legs. 
  /// This parameter controls the degree of attenuation for the main feed leg shadows.
  /// Reasonable values are less than 1. (as 1. means no feed leg shadows), put 0. for 
  /// the total blockage. 
  double itsFeedLegsShadowing;
  
  /// @brief additional shadowing factor for the first pair of feed leg wedges  
  /// @details Due to multiple reflections in the optical path, the illumination is worse
  /// close to the dish edge. This additional effect can be different for each pair of
  /// feed legs. This parameter is a factor applied on top of itsFeedLegsShadowing for
  /// the pair of feed legs rotated to the angle given by itsFeedLegsRotation.
  double itsFeedLegsWedgeShadowing1;

  /// @brief additional shadowing factor for the first pair of feed leg wedges  
  /// @details Due to multiple reflections in the optical path, the illumination is worse
  /// close to the dish edge. This additional effect can be different for each pair of
  /// feed legs. This parameter is a factor applied on top of itsFeedLegsShadowing for
  /// the pair of feed legs orthogonal to that rotated to the angle given by 
  /// itsFeedLegsRotation
  double itsFeedLegsWedgeShadowing2;
  
  /// @brief opening angle of surface-leg blockage in radians
  /// @details The model has triangular blockage near the edge of each feed leg shadow.
  /// This parameter controls the opening angle of this shadow. It is unused if feed leg
  /// shadows are not simulated.
  double itsWedgeOpeningAngle;
  
  /// @brief radius in metres where the wedge starts
  /// @details The model has triangular blockage near the edge of each feed leg shadow.
  /// This parameter controls the starting point of the wedge.
  double itsWedgeStartingRadius;
};


} // namespace synthesis

} // namespace askap

#endif // #ifndef ATCA_ILLUMINATION_H

