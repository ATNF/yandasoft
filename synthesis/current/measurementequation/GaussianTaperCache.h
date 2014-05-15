/// @file
///
/// @brief mainain an array with cached gaussian taper values
/// @details This is a common code between a number of classes required
/// to apply a gaussian taper (e.g. Gaussian Taper preconditioner)
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

#ifndef GAUSSIAN_TAPER_CACHE_H
#define GAUSSIAN_TAPER_CACHE_H

#include <casa/Arrays/Array.h>
#include <casa/Arrays/IPosition.h>

namespace askap {

namespace synthesis {

/// @brief mainain an array with cached gaussian taper values
/// @details This is a common code between a number of classes required
/// to apply a gaussian taper (e.g. Gaussian Taper preconditioner).
/// @note This class is relatively generic. It can be moved to a higher level,
/// if necessary.
/// @ingroup measurementequation
class GaussianTaperCache {
public:
   /// @brief set up the taper handler
   /// @details This constructor just sets the taper size. The size is full width at
   /// half maximum expressed in pixels.
   /// @param[in] majFWHM full width at half maximum of the major axis in pixels
   /// @param[in] minFWHM full width at half maximum of the minor axis in pixels 
   /// @param[in] pa position angle in radians
   GaussianTaperCache(double majFWHM, double minFWHM, double pa);
   
   /// @brief set up the taper handler with a circularly symmetric taper 
   /// @details This constructor just sets the taper size, same for both axis.
   /// The size is full width at half maximum expressed in pixels 
   /// @param[in] fwhm size in pixels
   explicit GaussianTaperCache(double fwhm);
   
   /// @brief Copy constructor
   /// @param[in] other input object
   /// @note We need a copy constructor because casa arrays use reference semantics
   GaussianTaperCache(const GaussianTaperCache &other);
   
   /// @brief assignment operator
   /// @param[in] other object to assign from
   /// @return reference to itself
   /// @note We need an assignment operator because casa arrays use reference semantics
   GaussianTaperCache& operator=(const GaussianTaperCache &other);   

   /// @brief obtain taper
   /// @details This method returns cached taper for a given shape. The taper
   /// is regenerated if the requested shape does not match the internal cache.
   /// The output is guaranteed to have the requested shape.
   /// @param[in] shape required shape
   /// @return array with the taper (casa arrays use reference semantics)
   casa::Array<casa::Complex> taper(const casa::IPosition &shape) const;
   
   /// @return major axis in pixels
   inline double majorAxis() const {return itsMajorAxis;}

   /// @return minor axis in pixels
   inline double minorAxis() const {return itsMinorAxis;}
   
   /// @return position angle in radians
   inline double posAngle() const {return itsPA;}
   
protected:           
   /// @brief build the cache 
   /// @details This method populates the cache using the values of
   /// data members
   /// @param[in] shape shape of the required array
   void initTaperCache(const casa::IPosition &shape) const;
   
private:
   /// @brief Major axis (sigma, rather than FWHM) in pixels
   double itsMajorAxis;
   /// @brief Minor axis (sigma, rather than FWHM) in pixels
   double itsMinorAxis;
   /// @brief position angle in radians
   double itsPA;
   /// @brief cache of the taper image
   /// @note May be we can make it float?
   mutable casa::Array<casa::Complex> itsTaperCache;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef GAUSSIAN_TAPER_CACHE_H

