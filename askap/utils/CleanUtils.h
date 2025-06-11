#ifndef SYNTHESIS_CLEANUTILS_H
#define SYNTHESIS_CLEANUTILS_H

/// @file CleanUtils.h
///
/// @brief Utilities to facilitate cleaning multiple images
///
/// @copyright (c) 2024 CSIRO
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
/// @author Mark Wieringa <mark.wieringa@csiro.au>

#include <casacore/casa/Arrays/Matrix.h>
#include <boost/optional.hpp>
#include <askap/scimath/fitting/Params.h>
// Local packages includes

using namespace casa;

namespace askap {
namespace utils {

/// @brief helper method to compute overlap masks for a set of images
/// @details when cleaning a set of images together we want to avoid putting the
/// same source in the model multiple times, this method will create a mask for
/// the biggest image that excludes the regions covered by other fields
/// @param[in] ip Params contains image parameters (starting with "image.")
/// @param[in] taylorMap map with image base names and number of Taylor terms
/// @param[in] extraOversamplingFactor if true, apply the oversampling factor when
/// generating the output mask
/// @return a Matrix with 1 for pixels with no overlap and 0 when there is overlap.
/// If there is only a single image centre present, the Matrix will have shape (0,0)
Matrix<imtype> overlapMask(const scimath::Params& ip, const std::map<std::string,int>& taylorMap,
    boost::optional<float> extraOversamplingFactor);

} // namespace utils
} // namespace askap

#endif
