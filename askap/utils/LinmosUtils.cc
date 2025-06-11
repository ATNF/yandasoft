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
#include "askap/askap_synthesis.h"

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
#include <casacore/casa/IO.h>
#include <casacore/coordinates/Coordinates/Coordinate.h>

// Local packages includes
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <linmos/LinmosAccumulator.h>

ASKAP_LOGGER(logger, ".linmos");

using namespace casacore;
using namespace askap::synthesis;

namespace askap {

/// @brief helper method to load beam offsets from the parset file
/// @details shares the same format as csimulator feed definition. This is needed to support ASKAP BETA,
///    which initially uses the same image centre for all beams, leaving beam offsets unspecified.
///    Therefore, this information has to be supplied by other means. Copied from testlinmos.
/// @param[in] parset parset containing spacing and offset parameters
/// @param[in] beamNames which offsets to get from the parset
/// @param[in] centre the pointing centre, which all offsets are relative to
/// @return a MVDirection for each name in beamNames
Vector<MVDirection> loadBeamOffsets(const LOFAR::ParameterSet &parset,
                                    const std::vector<std::string> &beamNames,
                                    MVDirection centre) {

    Vector<MVDirection> centres (beamNames.size(), centre);

    ASKAPLOG_INFO_STR(linmoslogger, " -> looking for the feed spacing");
    Quantity qspacing = askap::asQuantity(parset.getString("feeds.spacing"));
    const double spacing = qspacing.getValue("rad");
    ASKAPLOG_INFO_STR(linmoslogger, "    beam spacing set to " << qspacing);

    ASKAPLOG_INFO_STR(linmoslogger, " -> looking for a feed offset for each image");
    for (uint beam = 0; beam < beamNames.size(); ++beam) {
         const string parName = "feeds." + beamNames[beam];
         const std::vector<double> xy(parset.getDoubleVector(parName));
         ASKAPCHECK(xy.size() >= 2, "Expect two elements for each offset");
         // the shift appears to be positive in HA, so multiply by -1. Simulator.cc states:
         // "x direction is flipped to convert az-el type frame to ra-dec"
         //centres[beam].shift(-xy[0]*spacing, xy[1]*spacing, True);
         // Note: with ASKAP-36 the x shift as taken from the footprint file doesn't need this flip
         centres[beam].shift(xy[0]*spacing, xy[1]*spacing, True);
         ASKAPLOG_INFO_STR(linmoslogger, " -> " << parName << " centre: " << centres[beam] );
    }
    return centres;
}

/// @brief helper method to get beam centres from parset and/or image metadata
/// @details separate from loadParset to allow metadata to be read from input images
/// @param[in] parset linmos parset
/// @param[in] iacc image accessor
/// @param[in] outImgName current mosaic name
/// @return bool true=success, false=fail
Vector<MVDirection> loadBeamCentres(const LOFAR::ParameterSet &parset,
                                    const accessors::IImageAccess<Float> &iacc,
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

            //if (parset.isDefined("feeds.spacing")) {
            //    ASKAPLOG_WARN_STR(linmoslogger, "Feed info specified in parset but ignored. Using offset file");
            //}

        } else {

            ASKAPCHECK(parset.isDefined("names"), "No names specified in parset");
            return loadBeamOffsets(parset, parset.getStringVector("names", true), centre);

        }

    } else {
        ASKAPLOG_WARN_STR(linmoslogger, "Centre of the feeds not found. Setting beam centres to input ref. pixels");
        // no other information, so set the centre of the beam to be the reference pixel
        Vector<MVDirection> centres(inImgNames.size());
        for (int i = 0; i < centres.size(); i++) {
            const CoordinateSystem coordSys = iacc.coordSys(inImgNames[i]);
            const int dcPos = coordSys.findCoordinate(Coordinate::DIRECTION,-1);
            const DirectionCoordinate inDC = coordSys.directionCoordinate(dcPos);
            inDC.toWorld(centres[i],inDC.referencePixel());
        }
        return centres;
    }

    return Vector<MVDirection>();

}

/// @brief copy selected keywords from the reference image to the output
/// @param[in] string outName : output image name
/// @param[in] string inName : input image name
/// @param[in] vector<string> keywords : list of keyword names to copy
void copyKeywords(const string & outName, const string& inName, const vector<string> & keywords) {

    accessors::IImageAccess<Float>& iacc = SynthesisParamsHelper::imageHandler();

    for (const string& key : keywords) {
        if (key.size() > 0) {
            pair<string,string> valueAndComment = iacc.getMetadataKeyword(inName, key);
            if (valueAndComment.first.size() > 0) {
                iacc.setMetadataKeyword(outName, key, valueAndComment.first,valueAndComment.second);
            }
        }
    }
}


/// @brief save a table with the image containing the beamcentres and other information
/// @param[in] string outImgName : output image name, where the table will be saved
/// @param vector<string> inImgNames : input images
/// @param Vector<MVDirection> beamCentres : list of beam centres, must be same number as input images
void saveMosaicTable(const string & outImgName,const vector<string> & inImgNames,
                     const Vector<MVDirection> & beamCentres)
{
    ASKAPCHECK(inImgNames.size()==beamCentres.size(),"inImgNames and beamCentres must be the same length");
    int n = inImgNames.size();
    ASKAPLOG_INFO_STR(linmoslogger,"Saving the mosaic table");
    accessors::IImageAccess<Float>& iacc = SynthesisParamsHelper::imageHandler();
    Record record;
    // We want columns for beam number, inImgName, beamCentre (RA/Dec/J2000?)
    Record subRecord;
    Vector<Int> colBeamNumber(n);
    Vector<String> colImgName(n);
    Vector<Double> colRA(n);
    Vector<Double> colDEC(n);
    Vector<Double> colBMaj(n,0.);
    Vector<Double> colBMin(n,0.);
    Vector<Double> colBPA(n,0.);
    for (int i=0; i<n; i++) {
        colBeamNumber(i) = imagemath::LinmosAccumulator<Float>::getBeamFromImageName(inImgNames[i]);
        colImgName(i) = inImgNames[i];
        colRA(i) = beamCentres[i].getLong("deg").getValue();
        colDEC(i) = beamCentres[i].getLat("deg").getValue();
        Vector<Quantity> beam = iacc.beamInfo(inImgNames[i]);
        if (beam.nelements()==3) {
          colBMaj(i) = beam[0].getValue("arcsec");
          colBMin(i) = beam[1].getValue("arcsec");
          colBPA(i) = beam[2].getValue("deg");
        }
    }

    subRecord.define("BEAM",colBeamNumber);
    subRecord.define("NAME",colImgName);
    subRecord.define("RA",colRA);
    subRecord.define("DEC",colDEC);
    // Add the reference psf parameters
    subRecord.define("BMAJ",colBMaj);
    subRecord.define("BMIN",colBMin);
    subRecord.define("BPA",colBPA);
    Vector<String> units(7);
    units(0) = "";
    units(1) = "";
    units(2) = "deg"; //? FITS unit?
    units(3) = "deg";
    units(4) = "arcsec";
    units(5) = "arcsec";
    units(6) = "deg";
    subRecord.define("Units",units);
    record.defineRecord("MOSAIC",subRecord);
    // write it out
    iacc.setInfo(outImgName,record);
}

Vector<Float> readWeightsTable(const string& inImgName) {
    Vector<Float> wtVec;
    accessors::IImageAccess<Float>& iacc = SynthesisParamsHelper::imageHandler();
    const std::pair<std::string,std::string> imWtKw = iacc.getMetadataKeyword(inImgName,"IMWEIGHT");
    if (imWtKw.first!="") {
      try {
        float wt = std::stod(imWtKw.first);
        wtVec.resize(1);
        wtVec(0) = wt;
      } catch (const std::invalid_argument&) {
        ASKAPLOG_WARN_STR(logger, "Invalid float value for header keyword IMWEIGHT : "<<imWtKw.first);
      } catch (const std::out_of_range&) {
        ASKAPLOG_WARN_STR(logger, "Out of range double value for header keyword IMWEIGHT : "<<imWtKw.first);
      }
    } else {
      Record info;
      iacc.getInfo(inImgName,"WEIGHTS",info);
      if (!info.empty()) {
        Record wtsRec = info.asRecord("WEIGHTS");
        // get subrecord with arrays
        wtsRec = wtsRec.asRecord("WEIGHTS");
        ASKAPCHECK(wtsRec.isDefined("WEIGHT"),"WEIGHTS table does not have WEIGHT column");
        wtVec.assign(wtsRec.asArrayFloat("WEIGHT"));
      }
    }
    return wtVec;
}


} // namespace askap
