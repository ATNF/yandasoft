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

#include <measurementequation/GaussianTaperCache.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.gaussiantapercache");

#include <askap/AskapError.h>

#include <casa/BasicSL/Constants.h>
#include <casa/Arrays/ArrayMath.h>
#include <scimath/Mathematics/SquareMatrix.h>
#include <scimath/Mathematics/RigidVector.h>

#include <utils/PaddingUtils.h>

namespace askap {

namespace synthesis {


/// @brief set up the taper handler
/// @details This constructor just sets the taper size. The size is full width at
/// half maximum expressed in pixels.
/// @param[in] majFWHM full width at half maximum of the major axis in pixels
/// @param[in] minFWHM full width at half maximum of the minor axis in pixels 
/// @param[in] pa position angle in radians
GaussianTaperCache::GaussianTaperCache(double majFWHM, double minFWHM, double pa) :
     itsMajorAxis(majFWHM/sqrt(8.*log(2.))), itsMinorAxis(minFWHM/sqrt(8.*log(2.))),
     itsPA(pa) {}

/// @brief set up the taper handler with a circularly symmetric taper 
/// @details This constructor just sets the taper size, same for both axis.
/// The size is full width at half maximum expressed in pixels 
/// @param[in] fwhm size in pixels
GaussianTaperCache::GaussianTaperCache(double fwhm) : 
     itsMajorAxis(fwhm/sqrt(8.*log(2.))), itsPA(0.) 
{
  itsMinorAxis = itsMajorAxis;
}
           
/// @brief Copy constructor
/// @param[in] other input object
GaussianTaperCache::GaussianTaperCache(const GaussianTaperCache &other) :
   itsMajorAxis(other.itsMajorAxis), itsMinorAxis(other.itsMinorAxis), itsPA(other.itsPA),
   itsTaperCache(other.itsTaperCache.copy()) {}

/// @brief assignment operator
/// @param[in] other object to assign from
/// @return reference to itself
/// @note We need an assignment operator because casa arrays use reference semantics
GaussianTaperCache& GaussianTaperCache::operator=(const GaussianTaperCache &other)
{
  if (this != &other) {
      itsMajorAxis = other.itsMajorAxis;
      itsMinorAxis = other.itsMinorAxis;
      itsPA = other.itsPA;
      itsTaperCache.assign(other.itsTaperCache.copy());
  }
  return *this;
}

/// @brief obtain taper
/// @details This method returns cached taper for a given shape. The taper
/// is regenerated if the requested shape does not match the internal cache.
/// The output is guaranteed to have the requested shape.
/// @param[in] shape required shape
/// @return array with the taper (casa arrays use reference semantics)
casa::Array<casa::Complex> GaussianTaperCache::taper(const casa::IPosition &shape) const
{
  if (!shape.isEqual(itsTaperCache.shape())) {
      initTaperCache(shape);
  }
  return itsTaperCache;
}
           
/// @brief build the cache 
/// @details This method populates the cache using the values of
/// data members
/// @param[in] shape shape of the required array
void GaussianTaperCache::initTaperCache(const casa::IPosition &shape) const
{
  ASKAPDEBUGASSERT(shape.nelements() >= 2);

#ifdef ASKAP_DEBUG
  // if shape is exactly 2, nonDegenerate(2) would throw an exception. Hence, we need
  // a special check to avoid this.
  if (shape.nelements() > 2) {
     ASKAPASSERT(shape.nonDegenerate(2).nelements() == 2);
  }
#endif  

  itsTaperCache.resize(shape);
  const casa::Int nx = shape[0];
  const casa::Int ny = shape[1];
  casa::IPosition index(shape.nelements(),0);
  casa::SquareMatrix<casa::Double, 2> rotation(casa::SquareMatrix<casa::Double, 2>::General);
  // rotation direction is flipped here as we rotate the gaussian, not
  // the coordinate

  rotation(0,0) = rotation(1,1) = sin(itsPA);
  rotation(1,0) = cos(itsPA);
  rotation(0,1) = -rotation(1,0);
  
  // the following formula introduces some error if position angle is not 0
  // may be we need just to sum values?
  //const double normFactor = 2.*M_PI*itsMajorAxis*itsMinorAxis*erf(double(nx)/(2.*sqrt(2.)*itsMajorAxis))*
  //            erf(double(ny)/(2.*sqrt(2.)*itsMinorAxis));
  double sum = 0.;    
  const double maxRadius = double(casa::min(nx,ny)/2);        
  for (index[0] = 0; index[0]<nx; ++index[0]) {
       for (index[1] = 0; index[1]<ny; ++index[1]) {
            casa::RigidVector<casa::Double, 2> offset;
            offset(0) = (double(index[0])-double(nx)/2.);
            offset(1) = (double(index[1])-double(ny)/2.);
            if (sqrt(offset(0)*offset(0)+offset(1)*offset(1)) > maxRadius) {
                // fill on a circular rather than rectangular support
                itsTaperCache(index) = 0.;
                continue;
            } 
            // operator* is commented out in RigidVector due to
            // problems with some compilers. We have to use operator*= instead.
            // according to manual it is equivalent to v=Mv, rather than to v=v*M
            offset *= rotation;
            const double taperingFactor = exp(-casa::square(offset(0)/itsMajorAxis)/2.-
                       casa::square(offset(1)/itsMinorAxis)/2.);
            sum += taperingFactor;
            itsTaperCache(index) = taperingFactor;
       }
  }
  //std::cout<<"normFactor/sum: "<<normFactor/sum<<std::endl;
  //  itsTaperCache /= casa::Complex(sum,0.);
}

} // namespace synthesis

} // namespace askap

