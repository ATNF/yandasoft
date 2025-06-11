/// @file
/// @brief Interface class to access UVWeights
/// @details This interface is intended to be used inside gridder to structure access 
/// to uv weights at the stage of application in traditional weighting. In simplest form
/// it can be simply the weight collection. But other options may involve index translation
/// (e.g. if we want to apply the weight obtained for one beam to another or interpolate,
/// etc).
///
/// @copyright (c) 2023 CSIRO
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

// own includes
#include <askap/gridding/IUVWeightAccessor.h>

namespace askap {

namespace synthesis {


/// @brief virtual destructor to keep the compiler happy
IUVWeightAccessor::~IUVWeightAccessor() {}

} // namespace synthesis

} // namespace askap


