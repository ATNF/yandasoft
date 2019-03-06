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
#include <measurementequation/ImageConvolver.h>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/images/Regions/RegionHandler.h>
#include <casacore/images/Regions/ImageRegion.h>
#include <casacore/images/Images/ImageUtilities.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/LatticeMath/LatticeConvolver.h>
#include <casacore/lattices/Lattices/LatticeUtilities.h>
#include <casacore/lattices/LEL/LatticeExpr.h>
#include <casacore/lattices/LEL/LatticeExprNode.h>
#include <casacore/coordinates/Coordinates/CoordinateUtil.h>
#include <casacore/casa/iostream.h>

ASKAP_LOGGER(imageconvlogger, ".measurementequation.imageconvolver");

namespace askap {
namespace synthesis {

template <typename T>
ImageConvolver<T>::ImageConvolver()
{
}

template <typename T>
ImageConvolver<T>::~ImageConvolver()
{
}

using namespace casa;

template <typename T>
void ImageConvolver<T>::convolve(casa::ImageInterface<T>& imageOut,
                                 casa::ImageInterface<T>& imageIn,
                                 const casa::ImageInterface<T>& kernel,
                                 ScaleTypes scaleType, casa::Double scale,
                                 casa::Bool copyMiscellaneous, casa::Bool warnOnly)
{
    // Check Coordinates
    checkCoordinates(imageIn.coordinates(), kernel.coordinates(), warnOnly);

    // Convolve
    convolve(imageOut, imageIn, kernel, scaleType, scale, copyMiscellaneous);
}

template <typename T>
void ImageConvolver<T>::convolve(casa::ImageInterface<T>& imageOut,
                                 casa::ImageInterface<T>& imageIn,
                                 const casa::Array<T>& kernel,
                                 ScaleTypes scaleType, casa::Double scale,
                                 casa::Bool copyMiscellaneous)
{
    casa::ArrayLattice<T> kernelLattice(kernel);
    convolve(imageOut, imageIn, kernelLattice, scaleType, scale, copyMiscellaneous);
}


template <typename T>
void ImageConvolver<T>::convolve(casa::ImageInterface<T>& imageOut,
                                 casa::ImageInterface<T>& imageIn,
                                 const Lattice<T>& kernel,
                                 ScaleTypes scaleType, casa::Double scale,
                                 casa::Bool copyMiscellaneous)
{
    // Check
    const casa::IPosition& inShape = imageIn.shape();
    const casa::IPosition& outShape = imageOut.shape();
    if (!inShape.isEqual(outShape)) {
        ASKAPTHROW(AskapError, "Input and output images must have same shape");
    }

    if (kernel.ndim() > imageIn.ndim()) {
        ASKAPTHROW(AskapError, "Kernel lattice has more axes than the image!");
    }

    // Add degenerate axes if needed
    casa::Lattice<T>* pNewKernel = 0;
    casa::LatticeUtilities::addDegenerateAxes(pNewKernel, kernel, inShape.nelements());

    // Normalize kernel.
    casa::LatticeExprNode node;
    if (scaleType == AUTOSCALE) {
        node = casa::LatticeExprNode((*pNewKernel) / sum(*pNewKernel));
    } else if (scaleType == SCALE) {
        T t = static_cast<T>(scale);
        node = casa::LatticeExprNode(t * (*pNewKernel));
    } else if (scaleType == NONE) {
        node = casa::LatticeExprNode(*pNewKernel);
    }
    casa::LatticeExpr<T> kernelExpr(node);

    // Create convolver
    casa::LatticeConvolver<T> lc(kernelExpr, imageIn.shape(),  casa::ConvEnums::LINEAR);

    if (imageIn.isMasked()) {
        // Generate output mask if needed
        makeMask(imageOut);

        // Copy input mask to output.  Copy input pixels
        // and set to zero where masked
        casa::LogIO os;
        casa::LatticeUtilities::copyDataAndMask(os, imageOut, imageIn, True);

        // Convolve in situ
        lc.convolve(imageOut);
    } else {
        // Convolve to output
        lc.convolve(imageOut, imageIn);
    }

    // Overwrite output CoordinateSystem
    imageOut.setCoordinateInfo(imageIn.coordinates());

    // Copy miscellaneous things across as required
    if (copyMiscellaneous) casa::ImageUtilities::copyMiscellaneous(imageOut, imageIn);

    // Delete the restoring beam (should really check that the beam is in the
    // plane of convolution)
    casa::ImageInfo ii = imageOut.imageInfo();
    ii.removeRestoringBeam();
    imageOut.setImageInfo(ii);

    // Clean up
    delete pNewKernel;
}

template <typename T>
void ImageConvolver<T>::makeMask(casa::ImageInterface<T>& out)  const
{
    if (out.canDefineRegion()) {
        // Generate mask name
        casa::String maskName = out.makeUniqueRegionName(String("mask"), 0);

        // Make the mask if it does not exist
        if (!out.hasRegion(maskName, casa::RegionHandler::Masks)) {
            out.makeMask(maskName, true, true, false, true);
            ASKAPLOG_DEBUG_STR(imageconvlogger, "Created mask `" << maskName << "'");
        }
    } else {
        ASKAPLOG_WARN_STR(imageconvlogger, "Cannot make requested mask for this type of image");
    }
}

template <typename T>
void ImageConvolver<T>::checkCoordinates(const casa::CoordinateSystem& cSysImage,
        const casa::CoordinateSystem& cSysKernel,
        casa::Bool warnOnly) const
{
    const casa::uInt nPixelAxesK = cSysKernel.nPixelAxes();
    const casa::uInt nPixelAxesI = cSysImage.nPixelAxes();
    if (nPixelAxesK > nPixelAxesI) {
        ASKAPTHROW(AskapError, "Kernel has more pixel axes than the image");
    }

    const casa::uInt nWorldAxesK = cSysKernel.nWorldAxes();
    const casa::uInt nWorldAxesI = cSysImage.nWorldAxes();
    if (nWorldAxesK > nWorldAxesI) {
        ASKAPTHROW(AskapError, "Kernel has more world axes than the image");
    }

    const casa::Vector<casa::Double>& incrI = cSysImage.increment();
    const casa::Vector<casa::Double>& incrK = cSysKernel.increment();
    const casa::Vector<String>& unitI = cSysImage.worldAxisUnits();
    const casa::Vector<String>& unitK = cSysKernel.worldAxisUnits();

    for (casa::uInt i = 0; i < nWorldAxesK; i++) {
        // Compare Coordinate types and reference
        if (casa::CoordinateUtil::findWorldAxis(cSysImage, i) !=
                casa::CoordinateUtil::findWorldAxis(cSysKernel, i)) {
            if (warnOnly) {
                ASKAPLOG_WARN_STR(imageconvlogger, "Coordinate types are not the same for axis " << i + 1);
            } else {
                ASKAPTHROW(AskapError, "Coordinate types are not the same for axis " << i + 1);
            }
        }

        // Compare units
        const casa::Unit u1(unitI[i]);
        const casa::Unit u2(unitK[i]);
        if (u1 != u2) {
            if (warnOnly) {
                ASKAPLOG_WARN_STR(imageconvlogger, "Axis units are not consistent for axis " << i + 1);
            } else {
                ASKAPTHROW(AskapError, "Axis units are not consistent for axis " << i + 1);
            }
        }

        // Compare increments ; this is not a very correct test as there may be
        // values in the LinearTransform inside the coordinate.  Should really
        // convert some values...  See how we go.
        casa::Quantum<casa::Double> q2(incrK[i], u2);
        casa::Double val2 = q2.getValue(u1);
        if (!near(incrI[i], val2, 1.0e-6)) {
            if (warnOnly) {
                ASKAPLOG_WARN_STR(imageconvlogger, "Axis increments are not consistent for axis " << i + 1);
            } else {
                ASKAPTHROW(AskapError,"Axis increments are not consistent for axis " << i + 1);
            }
        }
    }
}

} // End namespace synthesis
} // End namespace askap
