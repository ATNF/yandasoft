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

// Include own header file first
#include <measurementequation/Image2DConvolver.h>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Exceptions/Error.h>
#include <components/ComponentModels/GaussianShape.h>
#include <casacore/coordinates/Coordinates/CoordinateUtil.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/lattices/LatticeMath/Fit2D.h>
#include <casacore/scimath/Functionals/Gaussian2D.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/ImageInfo.h>
#include <casacore/images/Images/ImageUtilities.h>
#include <casacore/scimath/Mathematics/Convolver.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Quanta/Unit.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/iostream.h>
#include <casacore/casa/Logging/LogIO.h>

// Local package includes
#include <measurementequation/ImageConvolver.h>

ASKAP_LOGGER(i2dconvlogger, ".measurementequation.image2dconvolver");

namespace askap {
namespace synthesis {

template <typename T>
Image2DConvolver<T>::Image2DConvolver()
{
}

template <typename T>
Image2DConvolver<T>::~Image2DConvolver()
{
}

template <typename T>
void Image2DConvolver<T>::convolve(casa::ImageInterface<T>& imageOut,
                                   casa::ImageInterface<T>& imageIn,
                                   casa::VectorKernel::KernelTypes kernelType,
                                   const casa::IPosition& pixelAxes,
                                   const casa::Vector<casa::Quantum<casa::Double> >& parameters,
                                   casa::Bool autoScale, casa::Double scale,
                                   casa::Bool copyMiscellaneous)
{
    // Checks
    if (parameters.nelements() != 3) {
        ASKAPTHROW(AskapError, "The world parameters vector must be of length 3");
    }
    if (pixelAxes.nelements() != 2) {
        ASKAPTHROW(AskapError, "You must give two pixel axes to convolve");
    }
    const Int nDim = imageIn.ndim();
    if (pixelAxes(0) < 0 || pixelAxes(0) >= nDim ||
            pixelAxes(1) < 0 || pixelAxes(1) >= nDim) {
        ASKAPTHROW(AskapError, "The pixel axes " << pixelAxes << " are illegal");
    }
    //
    if (nDim < 2) {
        ASKAPTHROW(AskapError, "The image axes must have at least 2 pixel axes");
    }
    //
    const casa::IPosition& inShape = imageIn.shape();
    const casa::IPosition& outShape = imageOut.shape();
    if (!inShape.isEqual(outShape)) {
        ASKAPTHROW(AskapError, "Input and output images must have the same shape");
    }

    // Generate Kernel Array (height unity)
    casa::Array<T> kernel;
    T kernelVolume = makeKernel(kernel, kernelType, parameters, pixelAxes, imageIn);

    // Figure out output image restoring beam (if any), output units and scale
    // factor for convolution kernel array
    casa::Vector<casa::Quantum<casa::Double> > beamOut;
    const CoordinateSystem& cSys = imageIn.coordinates();
    const ImageInfo& imageInfo = imageIn.imageInfo();
    const Unit& brightnessUnit = imageIn.units();
    casa::String brightnessUnitOut;

    dealWithRestoringBeam(brightnessUnitOut, beamOut, kernel, kernelVolume,
                          kernelType, parameters,
                          pixelAxes, cSys, imageInfo, brightnessUnit,
                          autoScale, scale);

    // Convolve.  We have already scaled the convolution kernel (with some
    // trickery cleverer than what ImageConvolver can do) so no more scaling
    casa::Double scale2 = 1.0;
    askap::synthesis::ImageConvolver<T> aic;
    aic.convolve(imageOut, imageIn, kernel, ImageConvolver<T>::NONE,
                 scale2, copyMiscellaneous);

    // Overwrite some bits and pieces in the output image to do with the
    // restoring beam  and image units
    imageOut.setUnits(brightnessUnitOut);
    ImageInfo iiOut = imageOut.imageInfo();
 
    casa::Bool holdsOneSkyAxis;
    casa::Bool hasSky = CoordinateUtil::holdsSky(holdsOneSkyAxis, cSys, pixelAxes.asVector());
    if (hasSky && beamOut.nelements() == 3) {
        iiOut.setRestoringBeam(beamOut);
        imageOut.setImageInfo(iiOut);
    } else {

        // If one of the axes is in the sky plane, we must
        // delete the restoring beam as it is no longer meaningful
        if (holdsOneSkyAxis) {
            ASKAPLOG_WARN_STR(i2dconvlogger, "Because you convolved just one of the sky axes");
            ASKAPLOG_WARN_STR(i2dconvlogger, "The output image does not have a valid spatial restoring beam");
            iiOut.removeRestoringBeam();
            imageOut.setImageInfo(iiOut);
        }
    }
}


// Private functions

template <typename T>
T Image2DConvolver<T>::makeKernel(casa::Array<T>& kernelArray,
                                  casa::VectorKernel::KernelTypes kernelType,
                                  const casa::Vector<casa::Quantum<casa::Double> >& parameters,
                                  const casa::IPosition& pixelAxes,
                                  const casa::ImageInterface<T>& imageIn) const
{
    // Check number of parameters
    checkKernelParameters(kernelType, parameters);

    // Convert kernel widths to pixels from world.  Demands major and minor
    // both in pixels or both in world, else exception
    casa::Vector<casa::Double> dParameters;
    const CoordinateSystem cSys = imageIn.coordinates();

    // Use the reference value for the shape conversion direction
    casa::Vector<casa::Quantum<casa::Double> > wParameters(5);
    for (casa::uInt i = 0; i < 3; i++) {
        wParameters(i + 2) = parameters(i);
    }

    const casa::Vector<casa::Double> refVal = cSys.referenceValue();
    const casa::Vector<String> units = cSys.worldAxisUnits();
    Int wAxis = cSys.pixelAxisToWorldAxis(pixelAxes(0));
    wParameters(0) = casa::Quantum<casa::Double>(refVal(wAxis), units(wAxis));
    wAxis = cSys.pixelAxisToWorldAxis(pixelAxes(1));
    wParameters(1) = casa::Quantum<casa::Double>(refVal(wAxis), units(wAxis));
    LogIO os;
    
    worldWidthsToPixel(os, dParameters, wParameters, cSys, pixelAxes, casa::False);

    // Create n-Dim kernel array shape
    const casa::IPosition kernelShape = shapeOfKernel(kernelType, dParameters, imageIn.ndim(), pixelAxes);

    // Create kernel array. We will fill the n-Dim array (shape non-unity
    // only for pixelAxes) through its 2D casa::Matrix incarnation. Aren't we clever.
    kernelArray.resize(kernelShape);
    casa::Array<T> kernelArray2 = kernelArray.nonDegenerate(pixelAxes);
    casa::Matrix<T> kernelMatrix = static_cast<casa::Matrix<T> >(kernelArray2);

    // Fill kernel casa::Matrix with functional (height unity)
    return fillKernel(kernelMatrix, kernelType, kernelShape, pixelAxes, dParameters);
}


template <typename T>
void Image2DConvolver<T>::dealWithRestoringBeam(
        casa::String& brightnessUnitOut,
        casa::Vector<casa::Quantum<casa::Double> >& beamOut,
        casa::Array<T>& kernelArray,
        T kernelVolume,
        casa::VectorKernel::KernelTypes,
        const casa::Vector<casa::Quantum<casa::Double> >& parameters,
        const casa::IPosition& pixelAxes,
        const casa::CoordinateSystem& cSys,
        const casa::ImageInfo& imageInfo,
        const casa::Unit& brightnessUnitIn,
        casa::Bool autoScale, casa::Double scale) const
{
    // Find out if convolution axes hold the sky.  Scaling from
    // Jy/beam and Jy/pixel only really makes sense if this is True
    casa::Bool holdsOneSkyAxis;
    const casa::Bool hasSky = casa::CoordinateUtil::holdsSky(holdsOneSkyAxis, cSys, pixelAxes.asVector());
    if (hasSky) {
        ASKAPLOG_INFO_STR(i2dconvlogger, "You are convolving the sky");
    } else {
        ASKAPLOG_INFO_STR(i2dconvlogger, "You are not convolving the sky");
    }

    // Generate an array holding the restoring beam if needed
    casa::Vector<casa::Quantum<casa::Double> > beamIn = imageInfo.restoringBeam().toVector();
    beamOut.resize(0);

    // Get brightness units
    const casa::String bUnitIn = upcase(brightnessUnitIn.getName());
    //
    const casa::Vector<casa::Double>& refPix = cSys.referencePixel();
    if (hasSky && bUnitIn == String("JY/PIXEL")) {

        // Easy case.  Peak of convolution kernel must be unity
        // and output units are Jy/beam.  All other cases require
        // numerical convolution of beams
        brightnessUnitOut = String("Jy/beam");
        beamOut.resize(3);

        // Exception already generated if only one of major and minor in pixel units
        if (parameters(0).getFullUnit().getName() == String("pix")) {
            casa::Vector<casa::Double> pixelParameters(5);
            pixelParameters(0) = refPix(pixelAxes(0));
            pixelParameters(1) = refPix(pixelAxes(1));
            pixelParameters(2) = parameters(0).getValue();
            pixelParameters(3) = parameters(1).getValue();
            pixelParameters(4) = parameters(2).getValue(Unit("rad"));
            casa::Vector<casa::Quantum<casa::Double> > worldParameters;
            //
            casa::LogIO os;
            pixelWidthsToWorld(os, worldParameters, pixelParameters,
                                               cSys, pixelAxes, casa::False);
            //
            beamOut(0) = worldParameters(0);
            beamOut(1) = worldParameters(1);
        } else {
            beamOut(0) = parameters(0);
            beamOut(1) = parameters(1);
        }
        beamOut(2) = parameters(2);                // Input p.a. is positive N->E
        //
        if (!autoScale) {
            const T t = static_cast<T>(scale);
            kernelArray *= t;
            ASKAPLOG_WARN_STR(i2dconvlogger, "Autoscaling is recommended for Jy/pixel convolution");
        }
    } else {

        // Is there an input restoring beam and are we convolving the sky to which it
        // pertains ?  If not, all we can do is use user scaling or normalize the convolution
        // kernel to unit volume.  There is no point to convolving the input beam either as it pertains
        // only to the sky
        if (hasSky && bUnitIn == String("JY/BEAM") && beamIn.nelements() == 3) {

            // Convert restoring beam parameters to pixels.  Output pa is pos +x -> +y in pixel frame.
            casa::Vector<casa::Quantum<casa::Double> > wParameters(5);
            const casa::Vector<casa::Double> refVal = cSys.referenceValue();
            const casa::Vector<String> units = cSys.worldAxisUnits();
            Int wAxis = cSys.pixelAxisToWorldAxis(pixelAxes(0));
            wParameters(0) = casa::Quantum<casa::Double>(refVal(wAxis), units(wAxis));
            wAxis = cSys.pixelAxisToWorldAxis(pixelAxes(1));
            wParameters(1) = casa::Quantum<casa::Double>(refVal(wAxis), units(wAxis));
            for (casa::uInt i = 0; i < 3; i++) {
                wParameters(i + 2) = beamIn(i);
            }
            casa::Vector<casa::Double> dParameters;
            casa::LogIO os;
            worldWidthsToPixel(os, dParameters, wParameters, cSys, pixelAxes, False);

            // Create 2-D beam array shape
            casa::IPosition dummyAxes(2, 0, 1);
            casa::IPosition beamShape = shapeOfKernel(casa::VectorKernel::GAUSSIAN,
                                        dParameters, 2, dummyAxes);

            // Create beam casa::Matrix and fill with height unity
            casa::Matrix<T> beamMatrixIn(beamShape(0), beamShape(1));
            fillKernel(beamMatrixIn, casa::VectorKernel::GAUSSIAN, beamShape,
                       dummyAxes, dParameters);

            // Get 2-D version of convolution kenrel
            casa::Array<T> kernelArray2 = kernelArray.nonDegenerate(pixelAxes);
            casa::Matrix<T> kernelMatrix = static_cast<casa::Matrix<T> >(kernelArray2);

            // Convolve input restoring beam array by convolution kernel array
            casa::Matrix<T> beamMatrixOut;
            Convolver<T> conv(beamMatrixIn, kernelMatrix.shape());   // matrixIn     = input restoring beam
            conv.linearConv(beamMatrixOut, kernelMatrix);

            // Scale kernel
            const T maxValOut = max(beamMatrixOut);
            if (autoScale) {
                brightnessUnitOut = String("Jy/beam");
                kernelArray /= maxValOut;
            } else {
                T t = static_cast<T>(scale);
                kernelArray *= t;
            }

            // Fit output beam matrix with a Gaussian, for better or worse
            // Fit2D is not templated.  So all our templating is useless
            // other than for Float until I template Fit2D
            casa::Fit2D fitter(os);
            const casa::uInt n = beamMatrixOut.shape()(0);
            //
            casa::Vector<casa::Double> bParameters = fitter.estimate(Fit2D::GAUSSIAN, beamMatrixOut);
            casa::Vector<casa::Bool> bParameterMask(bParameters.nelements(), True);
            bParameters(1) = (n - 1) / 2;      // x centre
            bParameters(2) = bParameters(1);    // y centre
            /*
               bParameterMask(1) = False;         // dont allow centre to float in fit
               bParameterMask(2) = False;
               */

            // Set range so we don't include too many pixels in fit which will make it very slow
            fitter.addModel(casa::Fit2D::GAUSSIAN, bParameters, bParameterMask);
            casa::Array<casa::Float> sigma;
            fitter.setIncludeRange(maxValOut / 10.0, maxValOut + 0.1);
            Fit2D::ErrorTypes error = fitter.fit(beamMatrixOut, sigma);
            if (error == Fit2D::NOCONVERGE ||
                    error == Fit2D::FAILED ||
                    error == Fit2D::NOGOOD) {
                ASKAPTHROW(AskapError, "Failed to fit the output beam");
            }
            casa::Vector<casa::Double> bSolution = fitter.availableSolution();

            // Convert to world units. Ho hum.
            casa::Vector<casa::Double> pixelParameters(5);
            pixelParameters(0) = refPix(pixelAxes(0));
            pixelParameters(1) = refPix(pixelAxes(1));
            pixelParameters(2) = bSolution(3);
            pixelParameters(3) = bSolution(4);
            pixelParameters(4) = bSolution(5);
            //
            pixelWidthsToWorld(os, beamOut, pixelParameters, cSys, pixelAxes, False);
        } else {
            if (autoScale) {
                // Conserve flux is the best we can do
                kernelArray /= kernelVolume;
            } else {
                const T t = static_cast<T>(scale);
                kernelArray *= t;
            }
        }
    }

    // Put beam position angle into range +/- 180 in case it has eluded us so far
    if (beamOut.nelements() == 3) {
        casa::MVAngle pa(beamOut(2).getValue(Unit("rad")));
        pa();
        beamOut(2) = casa::Quantum<casa::Double>(pa.degree(), Unit("deg"));
    }
}


template <typename T>
void Image2DConvolver<T>::checkKernelParameters(
        casa::VectorKernel::KernelTypes kernelType,
        const casa::Vector<casa::Quantum<casa::Double> >& parameters) const
{
    if (kernelType == casa::VectorKernel::BOXCAR) {
        ASKAPTHROW(AskapError, "Boxcar kernel not yet implemented");
        //
        if (parameters.nelements() != 3) {
            ASKAPTHROW(AskapError, "Boxcar kernels require 3 parameters");
        }
    } else if (kernelType == casa::VectorKernel::GAUSSIAN) {
        if (parameters.nelements() != 3) {
            ASKAPTHROW(AskapError, "Gaussian kernels require 3 parameters");
        }
    } else {
        ASKAPTHROW(AskapError, "The kernel type "
            << casa::VectorKernel::fromKernelType(kernelType) << " is not allowed");
    }
}

// Work out how big the array holding the kernel should be.
// Simplest algorithm possible. Shape is presently square.
template <typename T>
casa::IPosition Image2DConvolver<T>::shapeOfKernel(
        casa::VectorKernel::KernelTypes kernelType,
        const casa::Vector<casa::Double>& parameters,
        casa::uInt ndim,
        const casa::IPosition& axes) const
{
    // Find 2D shape
    casa::uInt n;
    if (kernelType == casa::VectorKernel::GAUSSIAN) {
        const casa::uInt n1 = sizeOfGaussian(parameters(0), 5.0);
        const casa::uInt n2 = sizeOfGaussian(parameters(1), 5.0);
        n = max(n1, n2);
        if (n % 2 == 0) n++;                                 // Make shape odd so centres well
    } else if (kernelType == casa::VectorKernel::BOXCAR) {
        n = 2 * Int(max(parameters(0), parameters(1)) + 0.5);
        if (n % 2 == 0) n++;                                 // Make shape odd so centres well
    } else {
        ASKAPTHROW(AskapError, "Unrecognized kernel type");  // Earlier checking should prevent this
    }

    // Now find the shape for the image and slot the 2D shape in
    // in the correct axis locations
    casa::IPosition shape(ndim, 1);
    shape(axes(0)) = n;
    shape(axes(1)) = n;

    return shape;
}

template <typename T>
casa::uInt Image2DConvolver<T>::sizeOfGaussian(casa::Double width, casa::Double nSigma) const
{
    // +/- 5sigma is a volume error of less than 6e-5%
    casa::Double sigma = width / sqrt(casa::Double(8.0) * C::ln2);
    return (Int(nSigma*sigma + 0.5) + 1) * 2;
}

template <typename T>
T Image2DConvolver<T>::fillKernel(casa::Matrix<T>& kernelMatrix,
                                  casa::VectorKernel::KernelTypes kernelType,
                                  const casa::IPosition& kernelShape,
                                  const casa::IPosition& axes,
                                  const casa::Vector<casa::Double>& parameters) const
{
    // Centre functional in array (shape is odd)
    // Need to think about these T castes for Complex images
    const T xCentre = static_cast<T>((kernelShape(axes(0)) - 1) / 2.0);
    const T yCentre = static_cast<T>((kernelShape(axes(1)) - 1) / 2.0);
    const T height = static_cast<T>(1.0);

    // Create functional.  We only have gaussian2d functionals
    // at this point.  Later the filling code can be moved out
    // of the if statement
    T maxValKernel, volumeKernel;
    const T pa = static_cast<T>(parameters(2));
    const T ratio = static_cast<T>(parameters(1) / parameters(0));
    const T major = static_cast<T>(parameters(0));
    if (kernelType == casa::VectorKernel::GAUSSIAN) {
        fillGaussian(maxValKernel, volumeKernel, kernelMatrix, height,
                     xCentre, yCentre, major, ratio, pa);
    } else if (kernelType == casa::VectorKernel::BOXCAR) {
        // TODO: Unsure why fillBoxcar is commented out. Rather than an uninitialised
        // volumeKernel, throw an exception until this is fixed.
        ASKAPTHROW(AskapError, "VectorKernel::BOXCAR not supported");
        /*
           fillBoxcar (maxValKernel, volumeKernel, kernelMatrix, height,
           xCentre, yCentre, major, ratio, pa);
           */
    } else {
        // Earlier checking should prevent this
        ASKAPTHROW(AskapError, "Unrecognized kernel type");
    }

    return volumeKernel;
}

template <typename T>
void Image2DConvolver<T>::fillGaussian(T& maxVal, T& volume,
                                       casa::Matrix<T>& pixels, T height, T xCentre,
                                       T yCentre, T majorAxis, T ratio,
                                       T positionAngle) const
//
// pa positive in +x ->+y pixel coordinate frame
//
{
    const casa::uInt n1 = pixels.shape()(0);
    const casa::uInt n2 = pixels.shape()(1);
    AlwaysAssert(n1 == n2, AipsError);
    positionAngle += C::pi_2;        // +y -> -x
    const Gaussian2D<T> g2d(height, xCentre, yCentre, majorAxis,
                            ratio, positionAngle);
    maxVal = -1.0e30;
    volume = 0.0;
    casa::Vector<T> pos(2);
    for (casa::uInt j = 0; j < n1; j++) {
        pos(1) = static_cast<T>(j);
        for (casa::uInt i = 0; i < n1; i++) {
            pos(0) = static_cast<T>(i);
            T val = g2d(pos);
            pixels(i, j) = val;
            //
            maxVal = max(val, maxVal);
            volume += val;
        }
    }
}
template < typename T >
void Image2DConvolver<T>::worldWidthsToPixel (casa::LogIO& os,
        casa::Vector<casa::Double>& dParameters,
        const casa::Vector<casa::Quantum<casa::Double> >& wParameters,
        const casa::CoordinateSystem& cSys,
        const casa::IPosition& pixelAxes,
        casa::Bool doRef) const
//
// world parameters: x, y, major, minor, pa
// pixel parameters: major, minor, pa (rad)
//
{
    if (pixelAxes.nelements()!=2) {
        os << "You must give two pixel axes" << LogIO::EXCEPTION;
    }
    if (wParameters.nelements()!=5) {
        os << "The world parameters vector must be of length 5" << LogIO::EXCEPTION;
    }
    //
    dParameters.resize(3);
    Int c0, c1, axisInCoordinate0, axisInCoordinate1;
    cSys.findPixelAxis(c0, axisInCoordinate0, pixelAxes(0));
    cSys.findPixelAxis(c1, axisInCoordinate1, pixelAxes(1));

    // Find units

    String majorUnit = wParameters(2).getFullUnit().getName();
    String minorUnit = wParameters(3).getFullUnit().getName();

    // This saves me trying to handle mixed pixel/world units which is a pain for coupled coordinates

    if ( (majorUnit==String("pix") && minorUnit!=String("pix"))  ||
            (majorUnit!=String("pix") && minorUnit==String("pix")) ) {
        os << "If pixel units are used, both major and minor axes must have pixel units" << LogIO::EXCEPTION;
    }

    // Some checks

    Coordinate::Type type0 = cSys.type(c0);
    Coordinate::Type type1 = cSys.type(c1);
    if (type0 != type1) {
        if (majorUnit!=String("pix") || minorUnit!=String("pix")) {
            os << "The coordinate types for the convolution axes are different" << endl;
            os << "Therefore the units of the major and minor axes of " << endl;
            os << "the convolution kernel widths must both be pixels" << LogIO::EXCEPTION;
        }
    }
    if (type0==Coordinate::DIRECTION && type1==Coordinate::DIRECTION &&  c0!=c1) {
        os << "The given axes do not come from the same Direction coordinate" << endl;
        os << "This situation requires further code development" << LogIO::EXCEPTION;
    }
    if (type0==Coordinate::STOKES || type1==Coordinate::STOKES) {
        os << "Cannot convolve Stokes axes" << LogIO::EXCEPTION;
    }

    // Deal with pixel units separately.    Both are in pixels if either is in pixels.

    if (majorUnit==String("pix")) {
        dParameters(0) = max(wParameters(2).getValue(), wParameters(3).getValue());
        dParameters(1) = min(wParameters(2).getValue(), wParameters(3).getValue());
        //
        if (type0==Coordinate::DIRECTION && type1==Coordinate::DIRECTION) {
            const DirectionCoordinate& dCoord = cSys.directionCoordinate (c0);

            // Use GaussianShape to get the position angle right. Use the specified
            // direction or the reference direction

            MDirection world;
            if (doRef) {
                dCoord.toWorld(world, dCoord.referencePixel());
            } else {
                world = MDirection(wParameters(0), wParameters(1), dCoord.directionType());
            }
            //
            Quantum<casa::Double> tmpMaj(1.0, Unit("arcsec"));
            GaussianShape gaussShape(world, tmpMaj, dParameters(1)/dParameters(0),
                    wParameters(4));                              // pa is N->E
            Vector<casa::Double> pars = gaussShape.toPixel (dCoord);
            dParameters(2) = pars(4);                                              // pa: +x -> +y
        } else {

            // Some 'mixed' plane; the pa is already +x -> +y

            dParameters(2) = wParameters(4).getValue(Unit("rad"));                  // pa
        }
        return;
    }

    // Continue on if non-pixel units

    if (type0==Coordinate::DIRECTION && type1==Coordinate::DIRECTION) {

        // Check units are angular

        Unit rad(String("rad"));
        if (!wParameters(2).check(rad.getValue())) {
            os << "The units of the major axis must be angular" << LogIO::EXCEPTION;
        }
        if (!wParameters(3).check(rad.getValue())) {
            os << "The units of the minor axis must be angular" << LogIO::EXCEPTION;
        }

        // Make a Gaussian shape to convert to pixels at specified location

        const DirectionCoordinate& dCoord = cSys.directionCoordinate (c0);
        //
        MDirection world;
        if (doRef) {
            dCoord.toWorld(world, dCoord.referencePixel());
        } else {
            world = MDirection(wParameters(0), wParameters(1), dCoord.directionType());
        }
        GaussianShape gaussShape(world, wParameters(2), wParameters(3), wParameters(4));
        Vector<casa::Double> pars = gaussShape.toPixel (dCoord);
        dParameters(0) = pars(2);
        dParameters(1) = pars(3);
        dParameters(2) = pars(4);      // radians; +x -> +y
    } else {

        // The only other coordinates currently available are non-coupled
        // ones and linear except for Tabular, which can be non-regular.
        // Urk.

        // Find major and minor axes in pixels

        dParameters(0) = worldWidthToPixel (os, dParameters(2), wParameters(2),
                cSys, pixelAxes);
        dParameters(1) = worldWidthToPixel (os, dParameters(2), wParameters(3),
                cSys, pixelAxes);
        dParameters(2) = wParameters(4).getValue(Unit("rad"));                // radians; +x -> +y
    }

    // Make sure major > minor

    Double tmp = dParameters(0);
    dParameters(0) = max(tmp, dParameters(1));
    dParameters(1) = min(tmp, dParameters(1));
}  
template < typename T >
casa::Bool Image2DConvolver<T>::pixelWidthsToWorld (casa::LogIO& os,
        casa::Vector<casa::Quantum<casa::Double> >& wParameters,
        const casa::Vector<casa::Double>& pParameters,
        const casa::CoordinateSystem& cSys,
        const casa::IPosition& pixelAxes,
        casa::Bool doRef) const
//
// pixel parameters: x, y, major, minor, pa (rad)
// world parameters: major, minor, pa
//
{
    if (pixelAxes.nelements()!=2) {
        os << "You must give two pixel axes" << LogIO::EXCEPTION;
    }
    if (pParameters.nelements()!=5) {
        os << "The parameters vector must be of length 5" << LogIO::EXCEPTION;
    }
    //
    Int c0, axis0, c1, axis1;
    cSys.findPixelAxis(c0, axis0, pixelAxes(0));
    cSys.findPixelAxis(c1, axis1, pixelAxes(1));
    Bool flipped = False;
    if (cSys.type(c1)==Coordinate::DIRECTION  && cSys.type(c0)==Coordinate::DIRECTION) {
        if (c0==c1) {
            flipped = skyPixelWidthsToWorld(os, wParameters, cSys, pParameters, pixelAxes, doRef);
        } else {
            os << "Cannot yet handle axes from different DirectionCoordinates" << LogIO::EXCEPTION;
        }
    } else {
        wParameters.resize(3);

        // Major/minor

        Quantum<casa::Double> q0 = pixelWidthToWorld (os, pParameters(4), pParameters(2),
                cSys, pixelAxes);
        Quantum<casa::Double> q1 = pixelWidthToWorld (os, pParameters(4), pParameters(3),
                cSys, pixelAxes);
        //
        if (q0.getValue() < q1.getValue(q0.getFullUnit())) {
            flipped = True;
            wParameters(0) = q1;
            wParameters(1) = q0;
        } else {
            wParameters(0) = q0;
            wParameters(1) = q1;
        }

        // Position angle; radians; +x -> +y

        wParameters(2).setValue(pParameters(4));
        wParameters(2).setUnit(Unit("rad"));
    }
    return flipped;
}
template < typename T>
casa::Double Image2DConvolver<T>::worldWidthToPixel (casa::LogIO& os, casa::Double positionAngle,
        const casa::Quantum<casa::Double>& length,
        const casa::CoordinateSystem& cSys,
        const casa::IPosition& pixelAxes) const
{
    //
    Int worldAxis0 = cSys.pixelAxisToWorldAxis(pixelAxes(0));
    Int worldAxis1 = cSys.pixelAxisToWorldAxis(pixelAxes(1));

    // Units of the axes must be consistent for now.
    // I will be able to relax this criterion when I get the time

    casa::Vector<casa::String> units = cSys.worldAxisUnits();
    casa::Unit unit0(units(worldAxis0));
    casa::Unit unit1(units(worldAxis1));
    if (unit0 != unit1) {
        os << "Units of the two axes must be conformant" << LogIO::EXCEPTION;
    }
    casa::Unit unit(unit0);

    // Check units are ok

    if (!length.check(unit.getValue())) {
        ostringstream oss;
        oss << "The units of the world length (" << length.getFullUnit().getName()
            << ") are not consistent with those of Coordinate System ("
            << unit.getName() << ")";
        casa::String s(oss);
        os << s << LogIO::EXCEPTION;
    }
    //
    casa::Double w0 = cos(positionAngle) * length.getValue(unit);
    casa::Double w1 = sin(positionAngle) * length.getValue(unit);

    // Find pixel coordinate of tip of axis  relative to reference pixel

    casa::Vector<casa::Double> world = cSys.referenceValue().copy();
    world(worldAxis0) += w0;
    world(worldAxis1) += w1;
    //
    casa::Vector<casa::Double> pixel;
    if (!cSys.toPixel (pixel, world)) {
        os << cSys.errorMessage() << LogIO::EXCEPTION;
    }
    // 
    casa::Double lengthInPixels = hypot(pixel(pixelAxes(0)), pixel(pixelAxes(1)));
    return lengthInPixels;
}
template < typename T>

casa::Quantum<casa::Double> Image2DConvolver<T>::pixelWidthToWorld (casa::LogIO& os,
        casa::Double positionAngle,
        casa::Double length,
        const casa::CoordinateSystem& cSys2,
        const casa::IPosition& pixelAxes) const
{
    casa::CoordinateSystem cSys(cSys2);
    Int worldAxis0 = cSys.pixelAxisToWorldAxis(pixelAxes(0));
    Int worldAxis1 = cSys.pixelAxisToWorldAxis(pixelAxes(1));

    // Units of the axes must be consistent for now.
    // I will be able to relax this criterion when I get the time

    casa::Vector<casa::String> units = cSys.worldAxisUnits().copy();
    casa::Unit unit0(units(worldAxis0));
    casa::Unit unit1(units(worldAxis1));
    if (unit0 != unit1) {
        os << "Units of the axes must be conformant" << LogIO::EXCEPTION;
    }

    // Set units to be the same for both axes

    units(worldAxis1) = units(worldAxis0);
    if (!cSys.setWorldAxisUnits(units)) {
        os << cSys.errorMessage() << LogIO::EXCEPTION;
    }
    //
    casa::Double p0 = cos(positionAngle) * length;
    casa::Double p1 = sin(positionAngle) * length;

    // Find world coordinate of tip of length relative to reference pixel

    casa::Vector<casa::Double> pixel= cSys.referencePixel().copy();
    pixel(pixelAxes(0)) += p0;
    pixel(pixelAxes(1)) += p1;
    //
    casa::Vector<casa::Double> world;
    if (!cSys.toWorld(world, pixel)) {
        os << cSys.errorMessage() << LogIO::EXCEPTION;
    }
    // 
    casa::Double lengthInWorld = hypot(world(worldAxis0), world(worldAxis1));
    casa::Quantum<casa::Double> q(lengthInWorld, Unit(units(worldAxis0)));
    //
    return q;
}
template <typename T>

casa::Bool Image2DConvolver<T>::skyPixelWidthsToWorld (casa::LogIO& os,
        casa::Vector<casa::Quantum<casa::Double> >& wParameters,
        const casa::CoordinateSystem& cSys,
        const casa::Vector<Double>& pParameters,
        const casa::IPosition& pixelAxes, casa::Bool doRef) const
//
// pixel parameters: x, y, major, minor, pa (rad)
// world parameters: major, minor, pa
//
{

    // What coordinates are these axes ?

    Int c0, c1, axisInCoordinate0, axisInCoordinate1;
    cSys.findPixelAxis(c0, axisInCoordinate0, pixelAxes(0));
    cSys.findPixelAxis(c1, axisInCoordinate1, pixelAxes(1));

    // See what sort of coordinates we have. Make sure it is called
    // only for the Sky.  More development needed otherwise.

    casa::Coordinate::Type type0 = cSys.type(c0);
    casa::Coordinate::Type type1 = cSys.type(c1);
    if (type0!=Coordinate::DIRECTION || type1!=Coordinate::DIRECTION) {
        os << "Can only be called for axes holding the sky" << LogIO::EXCEPTION;
    }
    if (c0!=c1) {
        os << "The given axes do not come from the same Direction coordinate" << endl;
        os << "This situation requires further code development" << LogIO::EXCEPTION;
    }

    // Is the 'x' (first axis) the Longitude or Latitude ?

    casa::Vector<casa::Int> dirPixelAxes = cSys.pixelAxes(c0);
    casa::Bool xIsLong = dirPixelAxes(0)==pixelAxes(0) && dirPixelAxes(1)==pixelAxes(1);
    casa::uInt whereIsX = 0;
    casa::uInt whereIsY = 1;
    if (!xIsLong) {
        whereIsX = 1;
        whereIsY = 0;
    }

    // Encode a pretend GaussianShape from these values as a means
    // of converting to world.

    const casa::DirectionCoordinate& dCoord = cSys.directionCoordinate(c0);
    GaussianShape gaussShape;
    casa::Vector<casa::Double> cParameters(pParameters.copy());
    //
    if (doRef) {
        cParameters(0) = dCoord.referencePixel()(whereIsX);     // x centre
        cParameters(1) = dCoord.referencePixel()(whereIsY);     // y centre
    } else {
        if (xIsLong) {
            cParameters(0) = pParameters(0);
            cParameters(1) = pParameters(1);
        } else {
            cParameters(0) = pParameters(1);
            cParameters(1) = pParameters(0);
        }
    }
    // 
    casa::Bool flipped = gaussShape.fromPixel (cParameters, dCoord);
    wParameters.resize(3);
    wParameters(0) = gaussShape.majorAxis();
    wParameters(1) = gaussShape.minorAxis();
    wParameters(2) = gaussShape.positionAngle();
    //
    return flipped;
} 



} // end namespace synthesis
} // end namespace askap
