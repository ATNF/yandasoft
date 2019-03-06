/// @file DeconvolverHelpers.h
///
/// DeconvolverHelpers: Helpers class for deconvolver
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///
#ifndef ASKAP_SYNTHESIS_DECONVOLVERHELPERS_H
#define ASKAP_SYNTHESIS_DECONVOLVERHELPERS_H

#include <Common/ParameterSet.h>
#include <boost/shared_ptr.hpp>

#include <deconvolution/DeconvolverBase.h>

namespace askap {
    namespace synthesis {
        /// @brief Helpers class for deconvolvers
        /// @detail This is a collection of static functions that are of use for
        /// cdeconvolver.
        /// @ingroup deconvolution
        class DeconvolverHelpers {
            public:

                /// @brief Get image as an array
                /// @detail An image will be converted into an array. So for e.g. dirty = gaskap.residual
                /// an array will be created from from the image file gaskap.residual
                /// @param name Name of image in the parset file
                /// @param parset ParameterSet containing description of images
                static casa::Array<casa::Float> getArrayFromImage(const casa::String name,
                                                                  const LOFAR::ParameterSet &parset);

                /// @brief Put an array as an image
                /// @detail An array will be written as an image, cloned from the
                /// file named templateName in the parset file
                /// an array will be created from from the image file gaskap.residual
                /// @param name Name of image in the parset file
                /// @param templateName Name of template image
                /// @param parset ParameterSet containing description of images
                static void putArrayToImage(const casa::Array<casa::Float> imageArray,
                                            const casa::String name, const casa::String templateName,
                                            const LOFAR::ParameterSet &parset);

            private:
                DeconvolverHelpers();
        };

    }
}
#endif
