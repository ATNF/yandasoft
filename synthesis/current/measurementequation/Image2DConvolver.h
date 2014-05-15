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
#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/Unit.h>
#include <casa/Arrays/IPosition.h>
#include <images/Images/ImageInfo.h>
#include <images/Images/ImageInterface.h>
#include <scimath/Mathematics/VectorKernel.h>
#include <coordinates/Coordinates/CoordinateSystem.h>

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
        void convolve(casa::ImageInterface<T>& imageOut,
                      casa::ImageInterface<T>& imageIn,
                      casa::VectorKernel::KernelTypes kernelType,
                      const casa::IPosition& pixelAxes,
                      const casa::Vector<casa::Quantum<casa::Double> >& parameters,
                      casa::Bool autoScale, casa::Double scale,
                      casa::Bool copyMiscellaneous = casa::True);

    private:

        /// Check kernel parameters
        void checkKernelParameters(casa::VectorKernel::KernelTypes kernelType,
                                   const casa::Vector<casa::Quantum<casa::Double> >& parameters) const;
        //
        void dealWithRestoringBeam(casa::String& brightnessUnitOut,
                                   casa::Vector<casa::Quantum<casa::Double> >& beamOut,
                                   casa::Array<T>& kernelArray,
                                   T kernelVolume,
                                   casa::VectorKernel::KernelTypes kernelType,
                                   const casa::Vector<casa::Quantum<casa::Double> >& parameters,
                                   const casa::IPosition& axes,
                                   const casa::CoordinateSystem& cSys,
                                   const casa::ImageInfo& imageInfo,
                                   const casa::Unit& brightnessUnit,
                                   casa::Bool autoscale, casa::Double scale) const;
        //
        T fillKernel(casa::Matrix<T>& kernelMatrix,
                     casa::VectorKernel::KernelTypes kernelType,
                     const casa::IPosition& kernelShape,
                     const casa::IPosition& axes,
                     const casa::Vector<casa::Double>& parameters) const;
        //
        void fillGaussian(T& maxVal, T& volume,
                          casa::Matrix<T>& pixels, T height, T xCentre,
                          T yCentre, T majorAxis, T ratio,
                          T positionAngle) const;
        //
        T makeKernel(casa::Array<T>& kernel,
                     casa::VectorKernel::KernelTypes kernelType,
                     const casa::Vector<casa::Quantum<casa::Double> >& parameters,
                     const casa::IPosition& axes,
                     const casa::ImageInterface<T>& inImage) const;
        //
        casa::IPosition shapeOfKernel(casa::VectorKernel::KernelTypes kernelType,
                                      const casa::Vector<casa::Double>& parameters,
                                      casa::uInt ndim,
                                      const casa::IPosition& axes) const;
        //
        casa::uInt sizeOfGaussian(casa::Double width, casa::Double nSigma) const;
};

}
}

#include <measurementequation/Image2DConvolver.tcc>

#endif
