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

#include <askap/measurementequation/GaussianTaperCache.h>

#include <askap/askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.gaussiantapercache");

#include <askap/AskapError.h>

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/scimath/Mathematics/SquareMatrix.h>
#include <casacore/scimath/Mathematics/RigidVector.h>

#include <askap/scimath/utils/PaddingUtils.h>

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
      itsTaperCache.assign(other.itsTaperCache); //.copy());
  }
  return *this;
}

casa::Vector<double> beam2Poly(const casa::Vector<double> & beam) {
    casa::Vector<double> p(3);
    double ct = cos(beam(2));
    double st = sin(beam(2));
    p(0) = casa::square(beam(0)*ct) + casa::square(beam(1)*st);
    p(1) = 2 * (casa::square(beam(1))-casa::square(beam(0)))*st*ct;
    p(2) = casa::square(beam(0)*st) + casa::square(beam(1)*ct);
    return p;
}

casa::Vector<double> poly2Beam(const casa::Vector<double> & poly) {
    casa::Vector<double> bm(3,0.);
    double a = poly(0) + poly(2);
    double b = sqrt(casa::square(poly(0)-poly(2)) + casa::square(poly(1)));
    bm(0) = sqrt((a+b)/2.);
    if (a-b > 0) {
        bm(1) = sqrt((a-b)/2.);
    }
    if (abs(poly(1))+abs(poly(0)-poly(2))>0) {
        bm(2) = 0.5 * atan2(-poly(1),poly(0)-poly(2));
    }
    return bm;
}

/// @brief tune taper parameters based on achieved resolution
/// @param[in] beam - fitted beam fwhm major,minor in image pixels and pos angle in radians
/// @param[in] tolerance - fractional tolerance in fwhm, also tolerance in rad for pa
/// @return true if converged within tolerance
bool GaussianTaperCache::tuneTaper(casa::Vector<double> beam, double tolerance) const
{
    ASKAPASSERT(itsTaper.nelements() == 3 && itsTaperCache.shape().nelements()>=2);
    ASKAPLOG_DEBUG_STR(logger,"Current beam parameters: "<< beam );
    // beam is in image pixels, taper in uvplane pixels
    // The relation between FWHMs in fourier and image plane is
    /// uvFWHM = (Npix / pixFWHM) * (4*log(2)/pi), where Npix is the number of pixels
    /// and pixFWHM is the image-plane FWHM in pixels.
    const double fwhm2sigma = sqrt(8.*log(2.));
    int nx = itsTaperCache.shape()(0);
    int ny = itsTaperCache.shape()(1);
    // convert from beam fwhm in image pixels to taper sigma in uv pixels
    // swap major/minor when going from image to uv plane
    double taper0 = (nx / beam(0)) * (4 * log(2)/casa::C::pi)/fwhm2sigma;
    double taper1 = (ny / beam(1)) * (4 * log(2)/casa::C::pi)/fwhm2sigma;
    double pa = beam(2);
    ASKAPLOG_DEBUG_STR(logger,"Eqv taper parameters   : "<< taper0*fwhm2sigma <<" "<<
    taper1*fwhm2sigma <<" "<< pa);
    ASKAPLOG_DEBUG_STR(logger,"Orig taper parameters  : "<< majorAxis()*fwhm2sigma <<" "<<
    minorAxis()*fwhm2sigma <<" "<< posAngle());
    casa::Vector<double> cbeam(3);
    cbeam(0) = nx / (majorAxis() * fwhm2sigma / (4 * log(2)/casa::C::pi));
    cbeam(1) = ny / (minorAxis() * fwhm2sigma / (4 * log(2)/casa::C::pi));
    cbeam(2) = posAngle();
    ASKAPLOG_DEBUG_STR(logger,"Req beam parameters    : "<< cbeam);

    if (abs(taper0 / majorAxis() - 1) < tolerance &&
        abs(taper1 / minorAxis() - 1) < tolerance &&
        abs((pa - posAngle())*(taper0-taper1)/taper0) < tolerance) return true;
    ASKAPLOG_DEBUG_STR(logger,"Actual taper parameters: "<< itsTaper(0)*fwhm2sigma <<" "<<
    itsTaper(1)*fwhm2sigma <<" "<< itsTaper(2));
    casa::Vector<double> abeam(3);
    abeam(0) = nx / (itsTaper(0) * fwhm2sigma / (4 * log(2)/casa::C::pi));
    abeam(1) = ny / (itsTaper(1) * fwhm2sigma / (4 * log(2)/casa::C::pi));
    abeam(2) = itsTaper(2);
    ASKAPLOG_DEBUG_STR(logger,"Actual beam parameters : "<< abeam);

    // Transform to polynomial coeffs and back for update
    casa::Vector<double> a = beam2Poly(abeam);
    //ASKAPLOG_DEBUG_STR(logger,"Actual beam poly parameters : "<< a);
    casa::Vector<double> b = beam2Poly(beam);
    casa::Vector<double> c = beam2Poly(cbeam);
    a += 0.5 * (c - b);
    abeam = poly2Beam(a);
    // Avoid division by zero - minimum convolving beam size of 0.01 pixels
    abeam(0) = fmax(abeam(0),0.01);
    abeam(1) = fmax(abeam(1),0.01);
    ASKAPLOG_DEBUG_STR(logger,"Updated beam parameters: "<< abeam);// << " poly: "<<a);
    itsTaper(0) = (nx / abeam(0)) * (4 * log(2)/casa::C::pi)/fwhm2sigma;
    itsTaper(1) = (ny / abeam(1)) * (4 * log(2)/casa::C::pi)/fwhm2sigma;
    itsTaper(2) = abeam(2);
    ASKAPLOG_DEBUG_STR(logger,"Updated taper parameters: "<< itsTaper(0)*fwhm2sigma <<" "<<
    itsTaper(1)*fwhm2sigma <<" "<< itsTaper(2));
    initTaperCache(itsTaperCache.shape());
    return false;
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

  if (itsTaper.nelements()!=3) {
      itsTaper.resize(3);
      itsTaper(0) = itsMajorAxis;
      itsTaper(1) = itsMinorAxis;
      itsTaper(2) = itsPA;
  }

  casa::SquareMatrix<casa::Double, 2> rotation(casa::SquareMatrix<casa::Double, 2>::General);
  // rotation direction is flipped here as we rotate the gaussian, not
  // the coordinate

  rotation(0,0) = rotation(1,1) = sin(itsTaper(2));
  rotation(1,0) = cos(itsTaper(2));
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
            const double taperingFactor = exp(-casa::square(offset(0)/itsTaper(0))/2.-
                       casa::square(offset(1)/itsTaper(1))/2.);
            sum += taperingFactor;
            itsTaperCache(index) = taperingFactor;
       }
  }
  //std::cout<<"normFactor/sum: "<<normFactor/sum<<std::endl;
  //  itsTaperCache /= casa::Complex(sum,0.);
}

} // namespace synthesis

} // namespace askap
