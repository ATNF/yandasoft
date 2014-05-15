/// @copyright (c) 1995,1996,1997,1998,1999,2000,2001,2002
/// Associated Universities, Inc. Washington DC, USA.
/// @copyright (c) 2012 CSIRO
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
/// NOTE: This file was imported from casacore as it is scheduled to be removed
/// from the casacore distribution.

#ifndef ASKAP_SYNTHESIS_IMAGECONVOLVER_H
#define ASKAP_SYNTHESIS_IMAGECONVOLVER_H

// ASKAPsoft Includes
#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <images/Images/ImageInterface.h>
#include <lattices/Lattices/Lattice.h>
#include <coordinates/Coordinates/CoordinateSystem.h>

namespace askap {
namespace synthesis {

/// @brief This class does convolution of an image by an Array or Lattice.
/// This class convolves an image by a specified kernel (Array or
/// Lattice).  If the kernel does not have enough dimensions, degenerate
/// ones are added.The class object has no state.  The functions could be static.
/// The convolution is done via FFT.  Thus input pixels which
/// are masked are set to 0 before the convolution.  The mask
/// is transferred to the output image.  No additional scaling
/// of the output image values is done.
template <typename T>
class ImageConvolver {
    public:

        enum ScaleTypes {
            // None; neither autoscaling nor direct scaling
            NONE,

            // Autoscale (normalize kernel to unit sum)
            AUTOSCALE,

            // SCALE (apply given scale factor)
            SCALE,

            // Number
            NTypes
        };

        // Constructor
        ImageConvolver();

        // Destructor
        ~ImageConvolver();

        // Convolve by an Image, Lattice or Array.  If convolving by an image
        // some rudimentary coordinate checks are made and warnings optionally issued
        // (warnOnly) if things are not commensurate.
        // If the output image needs a mask and doesn't have one,
        // it will be given one if possible. The input mask is transferred to
        // the output. The miscInfo, imageInfo, units and logger will be copied
        // from the input to the output unless you indicate not
        // to (copyMiscellaneous).  Any restoring beam is deleted from the
        // output image ImageInfo object.  The input CoordinateSystem
        // is transferred to the output image.  Degenerate  axes are added
        // to the kernel if it does not have enough dimensions.   If "autoScale"
        // is true, the kernel is normalized to have unit sum.  Otherwise,
        // the kernel is scaled (multiplied) by the value "scale".
        void convolve(casa::ImageInterface<T>& imageOut,
                      casa::ImageInterface<T>& imageIn,
                      const casa::ImageInterface<T>& kernel,
                      ScaleTypes scaleType, casa::Double scale,
                      casa::Bool copyMiscellaneous, casa::Bool warnOnly);

        void convolve(casa::ImageInterface<T>& imageOut,
                      casa::ImageInterface<T>& imageIn,
                      const casa::Lattice<T>& kernel,
                      ScaleTypes scaleType,
                      casa::Double scale,
                      casa::Bool copyMiscellaneous);

        void convolve(casa::ImageInterface<T>& imageOut,
                      casa::ImageInterface<T>& imageIn,
                      const casa::Array<T>& kernel,
                      ScaleTypes scaleType,
                      casa::Double scale,
                      casa::Bool copyMiscellaneous);

    private:

        // Make mask for image
        void makeMask(casa::ImageInterface<T>& out) const;

        // Check Coordinates of kernel and image
        void checkCoordinates(const casa::CoordinateSystem& cSysImage,
                              const casa::CoordinateSystem& cSysKernel,
                              casa::Bool warnOnly) const;
};

} // End synthesis namespace
} // End askap namespace

#include <measurementequation/ImageConvolver.tcc>

#endif
