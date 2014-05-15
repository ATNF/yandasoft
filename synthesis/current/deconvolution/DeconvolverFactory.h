/// @file DeconvolverFactory.h
///
/// DeconvolverFactory: Factory class for deconvolver
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
#ifndef ASKAP_SYNTHESIS_DECONVOLVERFACTORY_H_
#define ASKAP_SYNTHESIS_DECONVOLVERFACTORY_H_

#include <Common/ParameterSet.h>
#include <boost/shared_ptr.hpp>

#include <deconvolution/DeconvolverBase.h>

namespace askap {
    namespace synthesis {
        /// @brief Factory class for deconvolvers
        /// @ingroup deconvolution
        class DeconvolverFactory {
            public:
                /// @brief Factory class for all gridders.
                /// @todo Python version of factory
                DeconvolverFactory();

                /// @brief Make a shared pointer for a deconvolver
                /// @param parset ParameterSet containing description of
                /// deconvolver to be constructed.
                /// Cdeconvolver.dirty = gaskap.dirty
                /// Cdeconvolver.psf = gaskap.psf
                /// Cdeconvolver.mask = gaskap.mask
                /// Cdeconvolver.weight = gaskap.weight
                ///
                /// Cdeconvolver.solver = Clean
                /// Cdeconvolver.Clean.algorithm = Hogbom
                /// Cdeconvolver.Clean.gain = 0.1
                /// Cdeconvolver.Clean.tolerance = 1e-4
                /// Cdeconvolver.Clean.threshold = 0.001
                static DeconvolverBase<casa::Float, casa::Complex>::ShPtr make(const LOFAR::ParameterSet& parset);

            protected:

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
                /// @param templateName Name of template image
                /// @param name Name of image in the parset file
                /// @param parset ParameterSet containing description of images
                static casa::Array<casa::Float> putArrayToImage(const casa::String name,
                                                                const casa::String templateName,
                                                                const LOFAR::ParameterSet &parset);
        };

    }
}
#endif
