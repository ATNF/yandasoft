/// @file
/// @brief Interface class to an object function working on UVWeights
/// @details The derived classes of this interface are the ones which actually calculate
/// robust or other weights from the accumulated grid of weights. All operations are expected
/// to be in situ. The object function is called in the finalise method of the builder class.
/// At this stage we assume that all frequency planes and all indices can be processed independently
/// (and potentially in parallel). Therefore, this method works with the matrix interface rather than
/// the cube (low-level access is needed for performance here).
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

// own include 
#include <askap/gridding/IUVWeightCalculator.h>

namespace askap {

namespace synthesis {


/// @brief virtual destructor to keep the compiler happy
IUVWeightCalculator::~IUVWeightCalculator() {}

} // namespace synthesis

} // namespace askap

