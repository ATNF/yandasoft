/// @file
/// @brief Interface to a basic illumination pattern
/// @details This class is an abstract base (i.e. an interface) to
/// an hierarchy of classes representing illumination patterns.
/// It provides a method to obtain illumination pattern by populating a
/// pre-defined grid supplied as a UVPattern object.
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

#ifndef I_BASIC_ILLUMINATION_H
#define I_BASIC_ILLUMINATION_H

#include <askap/gridding/UVPattern.h>
#include <askap/dataaccess/IConstDataAccessor.h>

namespace askap {

namespace synthesis {

/// @brief Interface to a basic illumination pattern
/// @details This class is an abstract base (i.e. an interface) to
/// an hierarchy of classes representing illumination patterns.
/// It provides a method to obtain illumination pattern by populating a
/// pre-defined grid supplied as a UVPattern object.
/// @ingroup gridding
struct IBasicIllumination {

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
  /// uv-coordinate system and the system where the pattern is defined
  virtual void getPattern(double freq, UVPattern &pattern, double l,
                          double m, double pa) const = 0;

  /// @brief obtain illumination pattern
  /// @details This is the main method which populates the
  /// supplied uv-pattern with the values corresponding to the model
  /// represented by this object. It has to be overridden in the
  /// derived classes. An optional phase slope can be applied to
  /// simulate offset pointing.
  /// @param[in] freq frequency in Hz for which an illumination pattern is required
  /// @param[in] pattern a UVPattern object to fill
  /// @param[in] imageCentre direction of image
  /// @param[in] beamCentre  direction of beam (pointing)
  /// @param[in] pa parallactic angle, or strictly speaking the angle between
  /// uv-coordinate system and the system where the pattern is defined (unused)
  /// @param[in] isPSF specify if we want the pattern for image or psf (no phase slope)
  /// @param[in] feed  feed number for case where pattern differs between feeds
  virtual void getPattern(double freq, UVPattern &pattern,
                          const casacore::MVDirection &imageCentre = {},
                          const casacore::MVDirection &beamCentre = {},
                          const double pa = 0., const bool isPSF = false,
                          const int feed = 0) const = 0;

  /// @brief check whether the pattern is symmetric
  /// @details Some illumination patterns are trivial and it may be known a priori that
  /// the pattern does not depend on the parallactic angle. This method allows to check
  /// whether such trivial case exists. If true is returned, getPattern ignores pa
  /// parameter.
  /// @return true if the pattern is symmetric, false otherwise
  virtual bool isSymmetric() const = 0;

  /// @brief check whether the output pattern is image-based, rather than an illumination pattern.
  /// @details Some illumination patterns need to be generated in the image domain, and given
  /// the standard usage (FFT to image-domain for combination with other functions) any image
  /// domain function may as well stay in the image domain. So check the state before doing the FFT.
  /// @return false
  virtual bool isImageBased() const = 0;

  /// @brief check whether the output pattern is feed dependent
  /// @details Some illumination patterns vary with feed (number) and no shortcuts can
  /// be taken
  /// @return false
  virtual bool isFeedDependent() const = 0;

  /// @brief empty virtual destructor to keep the compiler happy
  virtual ~IBasicIllumination();
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef I_BASIC_ILLUMINATION_H
