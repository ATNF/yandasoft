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

#include <parallel/ParallelAccessor.h>
#include <askap/AskapError.h>

namespace askap {

namespace synthesis {

/// @brief constructor
/// @param[in] cacheSize uvw-machine cache size
/// @param[in] tolerance pointing direction tolerance in radians, exceeding
/// which leads to initialisation of a new UVW machine and recompute of the rotated uvws/delays
ParallelAccessor::ParallelAccessor(size_t cacheSize, double tolerance) : accessors::DataAccessorStub(false),
              itsRotatedUVW(cacheSize, tolerance) {}


/// @brief uvw after rotation
/// @details This method calls UVWMachine to rotate baseline coordinates 
/// for a new tangent point. Delays corresponding to this correction are
/// returned by a separate method.
/// @param[in] tangentPoint tangent point to rotate the coordinates to
/// @return uvw after rotation to the new coordinate system for each row
const casa::Vector<casa::RigidVector<casa::Double, 3> >&
                 ParallelAccessor::rotatedUVW(const casa::MDirection &tangentPoint) const
{
  return itsRotatedUVW.uvw(*this, tangentPoint);
}                 
   
/// @brief delay associated with uvw rotation
/// @details This is a companion method to rotatedUVW. It returns delays corresponding
/// to the baseline coordinate rotation. An additional delay corresponding to the 
/// translation in the tangent plane can also be applied using the image 
/// centre parameter. Set it to tangent point to apply no extra translation.
/// @param[in] tangentPoint tangent point to rotate the coordinates to
/// @param[in] imageCentre image centre (additional translation is done if imageCentre!=tangentPoint)
/// @return delays corresponding to the uvw rotation for each row
const casa::Vector<casa::Double>& ParallelAccessor::uvwRotationDelay(
                 const casa::MDirection &tangentPoint, const casa::MDirection &imageCentre) const
{
  return itsRotatedUVW.delays(*this,tangentPoint,imageCentre);
}
   
/// Velocity for each channel
/// @return a reference to vector containing velocities for each
///         spectral channel (vector size is nChannel). Velocities
///         are given as Doubles, the frame/units are specified by
///         the DataSource object (via IDataConverter).
const casa::Vector<casa::Double>& ParallelAccessor::velocity() const
{
  ASKAPTHROW(AskapError, "ParallelAccessor::velocity() has not been implemented");
}

} // namespace synthesis

} // namespace askap

