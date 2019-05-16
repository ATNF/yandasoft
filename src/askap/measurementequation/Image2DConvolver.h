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

#ifndef ASKAP_SYNTHESIS_IMAGE2DCONVOLVER_H
#define ASKAP_SYNTHESIS_IMAGE2DCONVOLVER_H

// ASKAPsoft includes
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/Unit.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/images/Images/ImageInfo.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/scimath/Mathematics/VectorKernel.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>

namespace askap {
namespace synthesis {

/// @brief This class does 2D convolution of an image by a functional form.
///
/// Motivation:
/// Convolution is a standard image processing requirement.  The
/// class object has no state.  The functions could be static.
/// The convolution is done via FFT.  Thus input pixels which
/// are masked are set to 0 before the convolution.  The mask
/// is transferred to the output image.  No additional scaling
// of the output image values is done.
template <typename T>
class Image2DConvolver {
    public:

        /// Constructor
        Image2DConvolver();

        /// Destructor
        ~Image2DConvolver();

        /// Convolve. If the output image needs a mask and doesn't have one,
        /// it will be given one if possible.  The miscInfo, imageInfo,
        /// units and logger will be copied from the input to the output
        /// unless you indicate not to (copyMiscellaneous).
        void convolve(casacore::ImageInterface<T>& imageOut,
                      casacore::ImageInterface<T>& imageIn,
                      casacore::VectorKernel::KernelTypes kernelType,
                      const casacore::IPosition& pixelAxes,
                      const casacore::Vector<casacore::Quantum<casacore::Double> >& parameters,
                      casacore::Bool autoScale, casacore::Double scale,
                      casacore::Bool copyMiscellaneous = casacore::True);

    private:

        /// Check kernel parameters
        void checkKernelParameters(casacore::VectorKernel::KernelTypes kernelType,
                                   const casacore::Vector<casacore::Quantum<casacore::Double> >& parameters) const;
        //
        void dealWithRestoringBeam(casacore::String& brightnessUnitOut,
                                   casacore::Vector<casacore::Quantum<casacore::Double> >& beamOut,
                                   casacore::Array<T>& kernelArray,
                                   T kernelVolume,
                                   casacore::VectorKernel::KernelTypes kernelType,
                                   const casacore::Vector<casacore::Quantum<casacore::Double> >& parameters,
                                   const casacore::IPosition& axes,
                                   const casacore::CoordinateSystem& cSys,
                                   const casacore::ImageInfo& imageInfo,
                                   const casacore::Unit& brightnessUnit,
                                   casacore::Bool autoscale, casacore::Double scale) const;
        //
        T fillKernel(casacore::Matrix<T>& kernelMatrix,
                     casacore::VectorKernel::KernelTypes kernelType,
                     const casacore::IPosition& kernelShape,
                     const casacore::IPosition& axes,
                     const casacore::Vector<casacore::Double>& parameters) const;
        //
        void fillGaussian(T& maxVal, T& volume,
                          casacore::Matrix<T>& pixels, T height, T xCentre,
                          T yCentre, T majorAxis, T ratio,
                          T positionAngle) const;
        //
        T makeKernel(casacore::Array<T>& kernel,
                     casacore::VectorKernel::KernelTypes kernelType,
                     const casacore::Vector<casacore::Quantum<casacore::Double> >& parameters,
                     const casacore::IPosition& axes,
                     const casacore::ImageInterface<T>& inImage) const;
        //
        casacore::IPosition shapeOfKernel(casacore::VectorKernel::KernelTypes kernelType,
                                      const casacore::Vector<casacore::Double>& parameters,
                                      casacore::uInt ndim,
                                      const casacore::IPosition& axes) const;
        //
        casacore::uInt sizeOfGaussian(casacore::Double width, casacore::Double nSigma) const;

        //
        // Legacy methods from casacore::ImageUtilities. The functionality of these
        // four methods below is now provided by WCBox/Polygon - but the interface is 
        // sufficiently different to make them incompatible for the moment.
        //
        // Convert a length and position angle in world units (for a non-coupled 
        // coordinate) to pixels. The length is in some 2D plane in the 
        // CoordinateSystem specified  by pixelAxes.

        void worldWidthsToPixel (casacore::LogIO& os, casacore::Vector<casacore::Double>& dParameters,
                                   const casacore::Vector<casacore::Quantum<casacore::Double> >& parameters,
                                   const casacore::CoordinateSystem& cSys,
                                   const casacore::IPosition& pixelAxes,
                                   casacore::Bool doRef=casacore::False) const;



        casacore::Bool pixelWidthsToWorld (casacore::LogIO& os,
                                   casacore::Vector<casacore::Quantum<casacore::Double> >& wParameters,
                                   const casacore::Vector<casacore::Double>& pParameters,
                                   const casacore::CoordinateSystem& cSys,
                                   const casacore::IPosition& pixelAxes,
                                   casacore::Bool doRef=casacore::False) const;

        casacore::Double worldWidthToPixel (casacore::LogIO& os, casacore::Double positionAngle,
                                    const casacore::Quantum<casacore::Double>& length,
                                    const casacore::CoordinateSystem& cSys,
                                    const casacore::IPosition& pixelAxes) const;



        casacore::Quantum<casacore::Double> pixelWidthToWorld (casacore::LogIO& os, casacore::Double positionAngle,
                                             casacore::Double length,
                                             const casacore::CoordinateSystem& cSys,
                                             const casacore::IPosition& pixelAxes) const;

        // Convert 2d sky shape (parameters=major axis, minor axis, position angle)
        // from pixels to world at reference pixel. pixelAxes describes which
        // 2 pixel axes of the coordinate system our 2D shape is in.
        // On input pa is positive for +x -> +y in pixel frame
        // On output pa is positive N->E
        // Returns True if major/minor exchanged themselves on conversion to world.
        casacore::Bool skyPixelWidthsToWorld (casacore::LogIO& os,
                                      casacore::Vector<casacore::Quantum<casacore::Double> >& wParameters,
                                      const casacore::CoordinateSystem& cSys,
                                      const casacore::Vector<casacore::Double>& pParameters,
                                      const casacore::IPosition& pixelAxes, casacore::Bool doRef) const;

};

}
}

#include <measurementequation/Image2DConvolver.tcc>

#endif
