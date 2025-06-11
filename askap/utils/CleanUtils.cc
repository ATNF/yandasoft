/// @file CleanUtils.cc
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

// Package level header file
#include "askap/askap_synthesis.h"

#include <casacore/lattices/LRegions/LCPolygon.h>

// Local packages includes
#include <askap/utils/CleanUtils.h>
#include <askap/askap/AskapLogging.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/imagemath/linmos/LinmosAccumulator.h>
#include <askap/measurementequation/ImageParamsHelper.h>

ASKAP_LOGGER(logger, "utils.cleanutils");


using namespace casacore;
using namespace askap;
using namespace synthesis;

/// @brief helper method to compute overlap masks for a set of images
/// @details when cleaning a set of images together we want to avoid putting the
/// same source in the model multiple times, this method will create a mask for
/// the largest image that excludes the regions covered by other fields
/// @param[in] ip Params contains image parameters (starting with "image.")
/// @param[in] taylorMap map with image base names and number of Taylor terms
/// @param[in] extraOversamplingFactor if true, apply the oversampling factor when
/// generating the output mask
/// @return a Matrix with 1 for pixels with no overlap and 0 when there is overlap.
/// If there is only a single image centre present, the Matrix will have shape (0,0)
Matrix<imtype> askap::utils::overlapMask(const scimath::Params& ip, const std::map<std::string,int>& taylorMap,
    boost::optional<float> extraOversamplingFactor)
{
    // make list of unique image names (.taylor0 only), their centres and their sizes
    std::vector<std::string> names;
    std::vector<DirectionCoordinate> DCs;
    std::vector<IPosition> shapes;
    for (const auto&  tmIt : taylorMap) {
        ImageParamsHelper iph(tmIt.first);
        if(tmIt.second > 1) {
            iph.makeTaylorTerm(0);
        }
        const std::string name = iph.paramName();
        const DirectionCoordinate dirCoord = ip.axes(name).directionAxis();
        const IPosition shape = ip.shape(name).getFirst(2);
        // Check we haven't seen this one before
        bool seen = false;
        for (int i=0; i < DCs.size(); i++) {
            if (allEQ(dirCoord.referenceValue(), DCs[i].referenceValue()))
            {
                // seen before, store this one if it is bigger

                if (shape.product() > shapes[i].product()) {
                    shapes[i] = shape;
                    names[i] = name;
                    DCs[i] = dirCoord;
                }
                seen = true;
            }
        }
        if (!seen) {

            names.push_back(name);
            DCs.push_back(dirCoord);
            shapes.push_back(shape);
        }
    }

    // Nothing to do if we have no, or a single field
    if (names.size() < 2) {
        return Matrix<imtype>();
    }

    // find biggest image
    size_t maxSize = 0;
    int mainImage = 0;
    for (int i=0; i < shapes.size(); i++) {
        // Apply oversampling
        if (extraOversamplingFactor) {
            shapes[i] =
                scimath::PaddingUtils::paddedShape(shapes[i],*extraOversamplingFactor);
            // fix coord systems
            Vector<double> refPix({shapes[i][0]/2.,shapes[i][1]/2.});
            DCs[i].setReferencePixel(refPix);
            DCs[i].setIncrement(DCs[i].increment()/double(*extraOversamplingFactor));
        }
        const size_t size = shapes[i].product();
        if ( size > maxSize) {
            maxSize = size;
            mainImage = i;
        }
    }
    ASKAPASSERT(maxSize > 0);
    // Create default mask
    Matrix<imtype> mask(shapes[mainImage](0),shapes[mainImage](1),static_cast<imtype>(1));

    // Work out overlap for each image and set pixels to zero
    const DirectionCoordinate& refDC = DCs[mainImage];
    bool anyOverlap = false;
    for (int i=0; i<names.size(); i++) {
        ASKAPLOG_DEBUG_STR(logger,"Field "<<i<<" : "<<names[i]);
        // apply mask for all images smaller than the biggest one
        // (if we have multiple images of the biggest size we are probably faceting)
        if (shapes[i].product() < shapes[mainImage].product()) {
            Vector<float> x,y;
            const DirectionCoordinate& inDC = DCs[i];
            const Vector<IPosition> edgePoints = imagemath::LinmosAccumulator<float>::
                convertImageEdgePointsToRef(inDC,shapes[i],refDC, false, true, &x, &y);
            bool overlap = false;
            for (const IPosition & point : edgePoints) {
                if (point >= 0  && point < shapes[mainImage]) {
                    overlap = true;
                    anyOverlap = true;
                    break;
                }
            }
            if (overlap) {
                ASKAPLOG_INFO_STR(logger,names[i]<<" overlaps "<<names[mainImage]<<", setting a mask on the main image");
                LCPolygon poly(x,y,shapes[mainImage]);
                mask(poly.boundingBox())(poly.maskArray()) = 0;
            }
        }
    }
    // if there is no overlap, we don't need the mask
    if (!anyOverlap) {
        mask.resize(0,0);
    }
    return mask;

}
