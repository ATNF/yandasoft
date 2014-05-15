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

#ifndef BASIC_COMPOSITE_ILLUMINATION_H
#define BASIC_COMPOSITE_ILLUMINATION_H

#include <gridding/IBasicIllumination.h>

#include <scimath/Mathematics/RigidVector.h>
#include <casa/Arrays/Vector.h>

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief Basic composite illumination pattern
/// @details This class is implements an basic composite illumination pattern corresponding
/// to a given weights and offsets of physical feeds. It can be used for simulation and/or
/// imaging with a synthetic beam. As an implementation of IBasicIllumination interface, 
/// this class provides a method to obtain illumination pattern by populating a pre-defined 
/// grid supplied as a UVPattern object. 
/// @note It looks like handling of illumination patterns
/// inside gridders has to be generalised (i.e. main method should receive a full accessor
/// with all the metadata instead of just the pointing offsets, frequency, etc). Such 
/// transition would definitely require an interface change in this class.
/// @ingroup gridding
struct BasicCompositeIllumination : public IBasicIllumination {
   /// @brief construct the pattern using given weights and offsets
   /// @param[in] pattern single-feed illumination pattern (assumed the same for all feeds)
   /// @param[in] feedOffsets offsets of physical feeds in radians
   /// @param[in] weights complex weights for each feed
   /// @note The size of two vectors should be the same
   BasicCompositeIllumination(const boost::shared_ptr<IBasicIllumination> &pattern,
            const casa::Vector<casa::RigidVector<casa::Double, 2> > &feedOffsets,
            const casa::Vector<casa::Complex> &weights);   
            
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
   virtual void getPattern(double freq, UVPattern &pattern, double l = 0., 
                           double m = 0., double pa = 0.) const;
 
   /// @brief check whether the pattern is symmetric
   /// @details Some illumination patterns are trivial and it may be known a priori that
   /// the pattern does not depend on the parallactic angle. This method allows to check
   /// whether such trivial case exists. If true is returned, getPattern ignores pa
   /// parameter.
   /// @return true if the pattern is symmetric, false otherwise
   virtual bool isSymmetric() const;
            
private:
   /// @brief single-feed illumination pattern (assumed the same for all feeds)
   boost::shared_ptr<IBasicIllumination> itsPattern;
 
   /// @brief offsets of physical feeds in radians
   casa::Vector<casa::RigidVector<casa::Double, 2> > itsFeedOffsets;
   
   /// @brief complex weights for each physical feed
   casa::Vector<casa::Complex> itsWeights;
   
   /// @brief flag showing that this pattern is symmetric
   /// @details Whether or not it is the case depends on the assigned offsets
   /// (i.e. any non-zero offset means automatically that this pattern is asymmetric,
   /// it is checked in the constructor)
   bool itsSymmetricFlag;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef BASIC_COMPOSITE_ILLUMINATION_H
