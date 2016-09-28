/// @file linmosUtils.cc
///
/// @brief combine a number of images as a linear mosaic
/// @details This is a utility to merge images into a mosaic. Images can be set
/// explicitly or found automatically based on input tags. 
///
/// @copyright (c) 2012,2014,2015 CSIRO
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
/// @author Daniel Mitchell <daniel.mitchell@csiro.au>

// Package level header file
#include "askap_synthesis.h"

// System includes
#include <sstream>
#include <typeinfo>
#include <iostream>

// other 3rd party
#include <Common/ParameterSet.h>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/images/Images/ImageRegrid.h>

// Local packages includes
#include <measurementequation/SynthesisParamsHelper.h>
#include <linmos/LinmosAccumulator.h>

ASKAP_LOGGER(logger, ".linmos");

using namespace casa;
using namespace askap;
using namespace askap::synthesis;
/// @brief helper method to load beam offsets from the parset file
/// @details shares the same format as csimulator feed definition. This is needed to support ASKAP BETA,
///    which initially uses the same image centre for all beams, leaving beam offsets unspecified.
///    Therefore, this information has to be supplied by other means. Copied from testlinmos.
/// @param[in] const LOFAR::ParameterSet &parset : parset containing spacing and offset parameters
/// @param[in] const Vector<std::string> beamNames : which offsets to get from the parset
/// @param[in] MVDirection centre : the pointing centre, which all offsets are relative to
/// @return Vector<MVDirection> : a MVDirection for each name in beamNames
Vector<MVDirection> loadBeamOffsets(const LOFAR::ParameterSet &parset,
                                    const Vector<std::string> beamNames,
                                    MVDirection centre) {

    Vector<MVDirection> centres (beamNames.size(), centre);

    ASKAPLOG_INFO_STR(linmoslogger, " -> looking for the feed spacing");
    Quantity qspacing = askap::asQuantity(parset.getString("feeds.spacing"));
    double spacing = qspacing.getValue("rad");
    ASKAPLOG_INFO_STR(linmoslogger, "    beam spacing set to " << qspacing);       

    ASKAPLOG_INFO_STR(linmoslogger, " -> looking for a feed offset for each image");
    for (uint beam = 0; beam < beamNames.size(); ++beam) {
         const string parName = "feeds." + beamNames[beam];
         const Vector<double> xy(parset.getDoubleVector(parName));
         //ASKAPCHECK(xy.size() == 2, "Expect two elements for each offset");
         // the shift appears to be positive in HA, so multiply by -1. Simulator.cc states:
         // "x direction is flipped to convert az-el type frame to ra-dec"
         centres[beam].shift(-xy[0]*spacing, xy[1]*spacing, casa::True);
         ASKAPLOG_INFO_STR(linmoslogger, " -> " << parName << " centre: " << centres[beam] );
    }
    return centres;
}

/// @brief helper method to get beam centres from parset and/or image metadata
/// @details separate from loadParset to allow metadata to be read from input images
/// @param[in] const LOFAR::ParameterSet &parset: linmos parset
/// @param[in] const accessors::IImageAccess &iacc: image accessor
/// @param[in] const string outImgName: current mosaic name
/// @return bool true=success, false=fail
Vector<MVDirection> loadBeamCentres(const LOFAR::ParameterSet &parset,
                                    const accessors::IImageAccess &iacc,
                                    const vector<string> &inImgNames) {

    // if setting weights using beam models, check the input for extra information

    ASKAPLOG_INFO_STR(linmoslogger, "Looking for parset options associated with primary-beam models");

    MVDirection centre;
    bool centreDefined = false;

    // set the centre of the "feeds" offset parameters (e.g. the boresight of the PAF)
    if (parset.isDefined("feeds.centre")) {
        ASKAPLOG_INFO_STR(linmoslogger, "Found centre of the feeds to use in beam models:");
        const vector<string> feedsCentre(parset.getStringVector("feeds.centre"));
        ASKAPCHECK(feedsCentre.size()==2, " -> the feeds.centre vector should have 2 elements");
        centre = convertDir(feedsCentre[0], feedsCentre[1]);
        ASKAPLOG_INFO_STR(linmoslogger, " -> "<<feedsCentre<<", = "<<centre);
        centreDefined = true;
    }
    else if (parset.isDefined("feeds.centreref")) {
        uint centreref = parset.getInt("feeds.centreref");
        if (centreref<inImgNames.size()) {
            ASKAPLOG_INFO_STR(linmoslogger, "Using the reference pixel of input image "<<centreref<<
                " as the centre of the feeds to use in beam models");
            const CoordinateSystem coordSys = iacc.coordSys(inImgNames[centreref]);
            const int DCpos = coordSys.findCoordinate(Coordinate::DIRECTION,-1);
            const DirectionCoordinate DC = coordSys.directionCoordinate(DCpos);
            DC.toWorld(centre,DC.referencePixel());
            ASKAPLOG_INFO_STR(linmoslogger, " -> "<<centre);
            centreDefined = true;
        }
        else {
            ASKAPLOG_WARN_STR(linmoslogger, "Found unsuitable centreref parameter: "<<centreref);
        }
    }

    // centres for each beam
    if (centreDefined) {

        if (parset.isDefined("feeds.offsetsfile")) {

            ASKAPLOG_INFO_STR(linmoslogger,  "Loading beam offsets from " << parset.getString("feeds.offsetsfile"));
            LOFAR::ParameterSet feed_parset(parset.getString("feeds.offsetsfile"));

            vector<string> beamNames;
            ASKAPLOG_INFO_STR(linmoslogger, " -> looking for feed names");
            if (parset.isDefined("feeds.names")) {
                beamNames = parset.getStringVector("feeds.names", true);
                ASKAPLOG_INFO_STR(linmoslogger,  "    using names given in the main parset");
            } else if (feed_parset.isDefined("feeds.names")) {
                beamNames = feed_parset.getStringVector("feeds.names", true);
                ASKAPLOG_INFO_STR(linmoslogger,  "    using names given in the feed-offset parset");
            }
            ASKAPCHECK(beamNames.size() > 0, "No beams specified");
            ASKAPCHECK(beamNames.size() == inImgNames.size(),
               "Number of beams does not match number of input files");

            return loadBeamOffsets(feed_parset, beamNames, centre);

            if (parset.isDefined("feeds.spacing")) {
                ASKAPLOG_WARN_STR(linmoslogger, "Feed info specified in parset but ignored. Using offset file");
            }

        } else {

            ASKAPCHECK(parset.isDefined("names"), "No names specified in parset");
            return loadBeamOffsets(parset, parset.getStringVector("names", true), centre);

        }

    } else {
        ASKAPLOG_WARN_STR(linmoslogger, "Centre of the feeds not found. Setting beam centres to input ref. pixels");
    }

    return Vector<MVDirection>();

}

