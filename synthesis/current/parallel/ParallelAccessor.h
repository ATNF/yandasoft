/// @file 
/// @brief accessor returned by parallel write iterator
/// @details We could've moved this class to accessors as it bears no 
/// parallelism. However, at this stage it is used only here and therefore some methods
/// we don't need here are left without implementation. The main functionality of this class
/// is to provide cached mechanism for rotated uvw and associated delays for the client (worker)
/// side. At the server side this is done by the table iterator.
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

#ifndef ASKAP_SYNTHESIS_PARALLEL_ACCESSOR_H
#define ASKAP_SYNTHESIS_PARALLEL_ACCESSOR_H

#include <dataaccess/DataAccessorStub.h>
#include <dataaccess/UVWRotationHandler.h>

namespace askap {

namespace synthesis {

/// @brief accessor returned by parallel write iterator
/// @details We could've moved this class to accessors as it bears no 
/// parallelism. However, at this stage it is used only here and therefore some methods
/// we don't need here are left without implementation. The main functionality of this class
/// is to provide cached mechanism for rotated uvw and associated delays for the client (worker)
/// side. At the server side this is done by the table iterator.
/// @ingroup parallel
class ParallelAccessor : public accessors::DataAccessorStub {
public:
   /// @brief constructor
   /// @param[in] cacheSize uvw-machine cache size
   /// @param[in] tolerance pointing direction tolerance in radians, exceeding
   /// which leads to initialisation of a new UVW machine and recompute of the rotated uvws/delays
   explicit ParallelAccessor(size_t cacheSize = 1, double tolerance = 1e-6);
   
   // override some stub methods
   
   /// @brief uvw after rotation
   /// @details This method calls UVWMachine to rotate baseline coordinates 
   /// for a new tangent point. Delays corresponding to this correction are
   /// returned by a separate method.
   /// @param[in] tangentPoint tangent point to rotate the coordinates to
   /// @return uvw after rotation to the new coordinate system for each row
   virtual const casa::Vector<casa::RigidVector<casa::Double, 3> >&
                 rotatedUVW(const casa::MDirection &tangentPoint) const;
   
   /// @brief delay associated with uvw rotation
   /// @details This is a companion method to rotatedUVW. It returns delays corresponding
   /// to the baseline coordinate rotation. An additional delay corresponding to the 
   /// translation in the tangent plane can also be applied using the image 
   /// centre parameter. Set it to tangent point to apply no extra translation.
   /// @param[in] tangentPoint tangent point to rotate the coordinates to
   /// @param[in] imageCentre image centre (additional translation is done if imageCentre!=tangentPoint)
   /// @return delays corresponding to the uvw rotation for each row
   virtual const casa::Vector<casa::Double>& uvwRotationDelay(
                 const casa::MDirection &tangentPoint, const casa::MDirection &imageCentre) const;
   
   /// Velocity for each channel
   /// @return a reference to vector containing velocities for each
   ///         spectral channel (vector size is nChannel). Velocities
   ///         are given as Doubles, the frame/units are specified by
   ///         the DataSource object (via IDataConverter).
   virtual const casa::Vector<casa::Double>& velocity() const;
   
private:
   /// @brief handler of uvw rotations
   accessors::UVWRotationHandler itsRotatedUVW;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_PARALLEL_ACCESSOR_H

