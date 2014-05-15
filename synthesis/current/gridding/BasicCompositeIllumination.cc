/// @file 
/// @brief Basic composite illumination pattern
/// @details This class is implements an basic composite illumination pattern corresponding
/// to a given weights and offsets of physical feeds. It can be used for simulation and/or
/// imaging with a synthetic beam. As an implementation of IBasicIllumination interface, 
/// this class provides a method to obtain illumination pattern by populating a pre-defined 
/// grid supplied as a UVPattern object. It looks like handling of illumination patterns
/// inside gridders has to be generalised (i.e. main method should receive a full accessor
/// with all the metadata instead of just the pointing offsets, frequency, etc). Such 
/// transition would definitely require an interface change in this class.
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
#include <scimath/Mathematics/SquareMatrix.h>



#include <gridding/BasicCompositeIllumination.h>

#include <askap/AskapError.h>

namespace askap {

namespace synthesis {

/// @brief construct the pattern using given weights and offsets
/// @param[in] pattern single-feed illumination pattern (assumed the same for all feeds)
/// @param[in] feedOffsets offsets of physical feeds in radians
/// @param[in] weights complex weights for each feed
/// @note The size of two vectors should be the same
BasicCompositeIllumination::BasicCompositeIllumination(const boost::shared_ptr<IBasicIllumination> &pattern,
            const casa::Vector<casa::RigidVector<casa::Double, 2> > &feedOffsets,
            const casa::Vector<casa::Complex> &weights) : itsPattern(pattern),
            itsFeedOffsets(feedOffsets), itsWeights(weights), itsSymmetricFlag(true)
{
   ASKAPDEBUGASSERT(itsPattern);
   ASKAPDEBUGASSERT(itsFeedOffsets.nelements() == itsWeights.nelements());
   if (itsPattern->isSymmetric()) {
      // iterate through all offsets and check whether any of them is non-zero
      for (casa::uInt iFeed = 0; iFeed<feedOffsets.nelements(); ++iFeed) {
           const casa::RigidVector<casa::Double, 2> &offset = feedOffsets[iFeed];
           // use tolerance around 1 nanoarcsecond
           if ((std::abs(offset(0))>5e-15) || (std::abs(offset(1))>5e-15)) {
               itsSymmetricFlag = false;
               break;
           }
      }
   } else {
      // single-feed illumination pattern is asymmetric. Hence, the composite pattern is
      // also asymmetric.
      itsSymmetricFlag = false;
   }
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
/// @param[in] pa parallactic angle (in radians), or strictly speaking the angle between 
/// uv-coordinate system and the system where the pattern is defined
void BasicCompositeIllumination::getPattern(double freq, UVPattern &pattern, double l, 
                           double m, double pa) const
{
   itsPattern->getPattern(freq,pattern,l,m,pa);
   // now apply the phase screen appropriate to the feed configuration/weights
   
   const casa::uInt oversample = pattern.overSample();
   const double cellU = pattern.uCellSize()/oversample;
   const double cellV = pattern.vCellSize()/oversample;
       
   // sizes of the grid to apply the phase screen to
   const casa::uInt nU = pattern.uSize();
   const casa::uInt nV = pattern.vSize();
   // number of feeds
   const casa::uInt nFeeds = itsWeights.nelements();
   
   double sum=0.; // normalisation factor
   
   casa::SquareMatrix<casa::Double, 2> rotation(casa::SquareMatrix<casa::Double, 2>::General);
   if (!itsSymmetricFlag) {
       rotation(0,0) = rotation(1,1) = cos(pa);
       rotation(0,1) = sin(pa);
       rotation(1,0) = -rotation(0,1);
   }
        
   for (casa::uInt iU=0; iU<nU; ++iU) {
	    const double offsetU = double(iU)-double(nU)/2.;
		for (casa::uInt iV=0; iV<nV; ++iV) {
	         const double offsetV = double(iV)-double(nV)/2.;
	         casa::Complex weight(0.,0.);
             for (casa::uInt iFeed = 0; iFeed<nFeeds; ++iFeed) {
                  // operator* is commented out in RigidVector.h due to problems with some
                  // compilers. We have to use operator*= instead (operator*= is
                  // equivalent to v=Mv, rather than v=vM according to inline doc).
                  casa::RigidVector<casa::Double, 2> feedOffset(itsFeedOffsets[iFeed]);
                  if (!itsSymmetricFlag) {
                      feedOffset *= rotation;
                  }
                  // don't need to multiply by wavelength here because the
			      // illumination pattern is given
			      // in a relative coordinates in frequency
                  const double phase = 2.*casa::C::pi*(cellU *feedOffset(0)*offsetU+
                                    cellV *feedOffset(1)*offsetV);
                  weight+=itsWeights[iFeed]*casa::Complex(cos(phase), -sin(phase));
			 }
			 pattern(iU, iV) *= casa::DComplex(weight);
			 sum += std::abs(weight);
		}
	}
	
	
    ASKAPCHECK(sum > 0., "Integral of the synthetic pattern should be non-zero");
    pattern.pattern() *= casa::DComplex(float(nU)*float(nV)/float(sum),0.); 
}

/// @brief check whether the pattern is symmetric
/// @details Some illumination patterns are trivial and it may be known a priori that
/// the pattern does not depend on the parallactic angle. This method allows to check
/// whether such trivial case exists. If true is returned, getPattern ignores pa
/// parameter.
/// @return true if the pattern is symmetric, false otherwise
bool BasicCompositeIllumination::isSymmetric() const
{
  return itsSymmetricFlag;
}

} // namespace synthesis

} // namespace askap
