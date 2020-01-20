/// @file 
/// @brief SKA_LOW illumination model
/// @details This class represents a SKA_LOW illumination model, 
/// which is the Fourier transform of the SKA_LOW_PB PrimaryBeam model.
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

#ifndef SKA_LOW_ILLUMINATION_H
#define SKA_LOW_ILLUMINATION_H

#include <askap/gridding/IBasicIllumination.h>
#include <askap/dataaccess/IConstDataAccessor.h>

namespace askap {

namespace synthesis {

/// @brief SKA_LOW illumination model
/// @details This class represents an SKA_LOW illumination model, 
/// which is the Fourier transform of the SKA_LOW_PB PrimaryBeam model.
/// Optionally a phase slope can be applied to simulate offset pointing.
/// @ingroup gridding
struct SKA_LOWIllumination : virtual public IBasicIllumination {

  /// @brief construct the model
  /// @param[in] diam station diameter in metres
  // DAM BLOCKAGE /// @param[in] blockage a diameter of the central hole in metres
  // DAM BLOCKAGE SKA_LOWIllumination(double diam, double blockage);
  SKA_LOWIllumination();
    
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
  virtual void getPattern(double freq, UVPattern &pattern, double l, 
                          double m, double pa) const;

  virtual void getPattern(double freq, UVPattern &pattern,
                          const casacore::MVDirection &imageCentre = {},
                          const casacore::MVDirection &beamCentre = {},
                          const double pa = 0., const bool isPSF = false) const;

  /// @brief check whether the pattern is symmetric
  /// @details Some illumination patterns like this one are trivial and known a priori to
  /// be symmetric. This method always returns true to reflect this
  /// @return always true 
  virtual bool isSymmetric() const;

  /// @brief check whether the output pattern is image-based, rather than an illumination pattern.
  /// @details Some illumination patterns need to be generated in the image domain, and given
  /// the standard usage (FFT to image-domain for combination with other functions) any image
  /// domain function may as well stay in the image domain. So check the state before doing the FFT.
  /// @return false 
  virtual bool isImageBased() const;

  // class-specific methods

  /// @brief  Set pointing parameters
  /// @details Set pointing parameters
  /// @param[in] az pointing azimuth in degrees
  /// @param[in] za pointing zenith angle in degrees
  void setPointing(double az, double za, double diam);

private:
  
  // The SKA_LOW_PB PrimaryBeam model is temporarily using a simple Gaussian
  // beam as a placeholder until a new SKA beam library is available. To remain
  // consistent, the following parameters are liable to change without warning.

  // setDipoleDelays should set steering delays, but since we are
  // temporarily using a Gaussian, just set pointing centre azimuth
  // and zenith angle
  /// @brief pointing azimuth in degrees
  double itsAz0;
  /// @brief pointing zenith angle in degrees
  double itsZa0;
  /// @brief station diameter in metres
  double itsDiameter;

};


} // namespace synthesis

} // namespace askap

#endif // #ifndef SKA_LOW_ILLUMINATION_H

