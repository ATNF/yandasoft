/// @file linmos.cc
///
/// @brief combine a number of images as a linear mosaic
/// @details This is a standalone utility to merge images into
///     a mosaic. Some code/functionality can later be moved into cimager,
///     but for now it is handy to have it separate. 
///
/// @copyright (c) 2012 CSIRO
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
/// @author Daniel Mitchell <daniel.mitchell@csiro.au> (2014)

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
#include <casa/Arrays/Array.h>
#include <images/Images/ImageRegrid.h>

#include <images/Images/PagedImage.h>

// ASKAPsoft includes
#include <askap/Application.h>
#include <askap/AskapError.h>
#include <askap/AskapLogging.h>
#include <askap/AskapUtil.h>
#include <askap/StatReporter.h>
#include <fitting/Params.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <utils/MultiDimArrayPlaneIter.h>
#include <scimath/Mathematics/Interpolate2D.h>

#include <fitting/ParamsCasaTable.h>

ASKAP_LOGGER(logger, ".linmos");

using namespace casa;
using namespace askap;
using namespace askap::synthesis;

// See loadParset for these options.
enum weight_types {FROM_WEIGHT_IMAGES=0, FROM_BP_MODEL};
// FROM_WEIGHT_IMAGES   Obtain pixel weights from weight images (parset "weights" entries)
// FROM_BP_MODEL        Generate pixel weights using a Gaussian primary-beam model
enum weight_states {CORRECTED=0, INHERENT, WEIGHTED};
// CORRECTED            Direction-dependent beams/weights have been divided out of input images
// INHERENT             Input images retain the natural primary-beam weighting of the visibilities
// WEIGHTED             Input images have full primary-beam-squared weighting

/// @brief the top-level class supporting linmos
class LinmosAccumulator {

  public:

    LinmosAccumulator();

    /// @brief check parset parameters
    /// @details check parset parameters for consistency and set any dependent variables
    ///     weighttype: FromWeightImages or FromPrimaryBeamModel. No default.
    ///     weightstate: Corrected, Inherent or Weighted. Default: Corrected.
    /// @param[in] const LOFAR::ParameterSet &parset: linmos parset
    /// @return bool true=success, false=fail
    bool loadParset(const LOFAR::ParameterSet &parset);

    /// @brief get beam centres from parset and/or image metadata
    /// @details separate from loadParset to allow metadata to be read from input images
    /// @param[in] const LOFAR::ParameterSet &parset: linmos parset
    /// @param[in] const accessors::IImageAccess &iacc: image accessor
    /// @param[in] const string outImgName: current mosaic name
    /// @return bool true=success, false=fail
    bool loadBeamCentres(const LOFAR::ParameterSet &parset,
                         const accessors::IImageAccess &iacc,
                         const string outImgName);

    /// @brief set up a single mosaic
    /// @param[in] const vector<string> &inImgNames : vector of images to mosaic
    /// @param[in] const vector<string> &inWgtNames : vector of weight images, if required
    /// @param[in] const string &outImgName : output mosaic image name
    /// @param[in] const string &outWgtName : output mosaic weight image name
    void setSingleMosaic(const vector<string> &inImgNames, const vector<string> &inWgtNames,
                         const string &outImgName, const string &outWgtName);


    /// @brief set up a single mosaic for each taylor term
    /// @param[in] const vector<string> &inImgNames : vector of images to mosaic
    /// @param[in] const vector<string> &inWgtNames : vector of weight images, if required
    /// @param[in] const string &outImgName : output mosaic image name
    /// @param[in] const string &outWgtName : output mosaic weight image name
    void findAndSetTaylorTerms(const vector<string> &inImgNames, const vector<string> &inWgtNames,
                               const string &outImgName, const string &outWgtName);

    /// @brief search the current directory for suitable mosaics
    /// @details based on a vector of image tags, look for sets of images with names
    ///     that contain all tags but are otherwise equal and contain an allowed prefix.
    /// @param[in] const vector<string> &imageTags : vector of image tags
    void findAndSetMosaics(const vector<string> &imageTags);

    /// @brief test whether the output buffers are empty and need initialising
    /// @return bool
    bool outputBufferSetupRequired(void);

    /// @brief set the input coordinate system and shape
    /// @param[in] const string& inImgName: name of the input image
    /// @param[in] const accessors::IImageAccess& iac
    /// @param[in] const int n: input order of the input file
    void setInputParameters(const string& inImgName, const accessors::IImageAccess& iacc, const int n);

    /// @brief set the output coordinate system and shape, based on the overlap of input images
    /// @details This method is based on the SynthesisParamsHelper::add and
    ///     SynthesisParamsHelper::facetSlicer. It has been reimplemented here
    ///     so that images can be read into memory separately.
    /// @param[in] const vector<string>& inImgNames: names of the input images (those given for parset key 'names')
    /// @param[in] const accessors::IImageAccess& iac
    void setOutputParameters(const vector<string>& inImgNames, const accessors::IImageAccess& iacc);

    /// @brief set up any 2D temporary output image buffers required for regridding
    void initialiseOutputBuffers(void);

    /// @brief set up any 2D temporary input image buffers required for regridding
    void initialiseInputBuffers(void);

    /// @brief set up regridder
    void initialiseRegridder(void);

    /// @brief load the temporary image buffers with the current plane of the current input image
    /// @param[in] const scimath::MultiDimArrayPlaneIter& planeIter: current plane id
    /// @param[in] Array<float>& inPix: image buffer
    /// @param[in] Array<float>& inWgtPix: weight image buffer
    void loadInputBuffers(const scimath::MultiDimArrayPlaneIter& planeIter,
                          Array<float>& inPix, Array<float>& inWgtPix, Array<float>& inSenPix);

    /// @brief call the regridder for the buffered plane
    void regrid(void);

    /// @brief add the current plane to the accumulation arrays
    /// @details This method adds from the regridded buffers
    /// @param[out] Array<float>& outPix: accumulated weighted image pixels
    /// @param[out] Array<float>& outWgtPix: accumulated weight pixels
    /// @param[in] const IPosition& curpos: indices of the current plane
    void accumulatePlane(Array<float>& outPix, Array<float>& outWgtPix,
                         Array<float>& outSenPix, const IPosition& curpos);

    /// @brief add the current plane to the accumulation arrays
    /// @details This method adds directly from the input arrays
    /// @param[out] Array<float>& outPix: accumulated weighted image pixels
    /// @param[out] Array<float>& outWgtPix: accumulated weight pixels
    /// @param[in] const Array<float>& inPix: input image pixels
    /// @param[in] const Array<float>& inWgtPix: input weight pixels
    /// @param[in] const IPosition& curpos: indices of the current plane
    void accumulatePlane(Array<float>& outPix, Array<float>& outWgtPix, Array<float>& outSenPix,
                         const Array<float>& inPix, const Array<float>& inWgtPix,
                         const Array<float>& inSenPix, const IPosition& curpos);

    /// @brief divide the weighted pixels by the weights for the current plane
    /// @param[in,out] Array<float>& outPix: accumulated deweighted image pixels
    /// @param[in] const Array<float>& outWgtPix: accumulated weight pixels
    /// @param[in] const IPosition& curpos: indices of the current plane
    void deweightPlane(Array<float>& outPix, const Array<float>& outWgtPix,
                       Array<float>& outSenPix, const IPosition& curpos);

    /// @brief check to see if the input and output coordinate grids are equal
    /// @return bool: true if they are equal
    bool coordinatesAreEqual(void);

    // return metadata for the current input image
    IPosition inShape(void) {return itsInShape;}
    CoordinateSystem inCoordSys(void) {return itsInCoordSys;}

    // return metadata for the output image
    IPosition outShape(void) {return itsOutShape;}
    CoordinateSystem outCoordSys(void) {return itsOutCoordSys;}

    int weightType(void) {return itsWeightType;}
    int weightState(void) {return itsWeightState;}
    int numTaylorTerms(void) {return itsNumTaylorTerms;}
    bool doSensitivity(void) {return itsDoSensitivity;}
    void doSensitivity(bool value) {itsDoSensitivity = value;}
    string taylorTag(void) {return itsTaylorTag;}

    map<string,string> outWgtNames(void) {return itsOutWgtNames;}
    map<string,string> outSenNames(void) {return itsOutSenNames;}
    map<string,vector<string> > inImgNameVecs(void) {return itsInImgNameVecs;}
    map<string,vector<string> > inWgtNameVecs(void) {return itsInWgtNameVecs;}
    map<string,vector<string> > inSenNameVecs(void) {return itsInSenNameVecs;}
    map<string,bool> outWgtDuplicates(void) {return itsOutWgtDuplicates;}
    map<string,bool> genSensitivityImage(void) {return itsGenSensitivityImage;}

  private:

    /// @brief convert the current input shape and coordinate system to the reference (output) system
    /// @param[in] const DirectionCoordinate& refDC: reference direction coordinate
    /// @return IPosition vector containing BLC and TRC of the current input image, relative to another coord. system
    Vector<IPosition> convertImageCornersToRef(const DirectionCoordinate& refDC);

    /// @brief check to see if the input coordinate system is consistent enough with the reference system to merge
    /// @param[in] const CoordinateSystem& refCoordSys: reference coordinate system
    /// @return bool: true if they are consistent
    bool coordinatesAreConsistent(const CoordinateSystem& refCoordSys);

    // regridding options
    ImageRegrid<float> itsRegridder;
    IPosition itsAxes;
    String itsMethod;
    Int itsDecimate;
    Bool itsReplicate;
    Bool itsForce;
    Interpolate2D::Method itsEmethod;
    // regridding buffers
    TempImage<float> itsInBuffer, itsInWgtBuffer, itsInSenBuffer, itsInSnrBuffer;
    TempImage<float> itsOutBuffer, itsOutWgtBuffer, itsOutSnrBuffer;
    // metadata objects
    IPosition itsInShape;
    CoordinateSystem itsInCoordSys;
    IPosition itsOutShape;
    CoordinateSystem itsOutCoordSys;
    // options
    int itsWeightType;
    int itsWeightState;
    int itsNumTaylorTerms;
    bool itsDoSensitivity;

    float itsCutoff;

    // 
    Vector<MVDirection> itsCentres;
    MVDirection itsInCentre;

    // Set some objects to support multiple mosaics.
    const string itsMosaicTag;
    const string itsTaylorTag;

    map<string,string> itsOutWgtNames;
    map<string,string> itsOutSenNames;
    map<string,vector<string> > itsInImgNameVecs;
    map<string,vector<string> > itsInWgtNameVecs;
    map<string,vector<string> > itsInSenNameVecs;
    map<string,bool> itsOutWgtDuplicates;
    map<string,bool> itsGenSensitivityImage;

};

LinmosAccumulator::LinmosAccumulator() : itsMethod("linear"), itsDecimate(3), itsReplicate(false), itsForce(false),
                                         itsWeightType(-1), itsWeightState(-1), itsNumTaylorTerms(-1),
                                         itsCutoff(0.01), itsMosaicTag("linmos"), itsTaylorTag("taylor.0") {}


// functions used by the linmos accumulator class

/// @brief function to convert RA and Dec strings to a MVDirection
/// @param[in] const std::string &ra
/// @param[in] const std::string &dec
/// @return MVDirection
MVDirection convertDir(const std::string &ra, const std::string &dec) {
  Quantity tmpra,tmpdec;
  Quantity::read(tmpra, ra);
  Quantity::read(tmpdec,dec);
  return MVDirection(tmpra,tmpdec);  
}

/// @brief helper method to load beam offsets from the parset file
/// @details shares the same format as csimulator feed definition. This is needed to support ASKAP BETA,
///    which initially uses the same image centre for all beams, leaving beam offsets unspecified.
///    Therefore, this information has to be supplied by other means. Copied from testlinmos.
/// @param[in] const LOFAR::ParameterSet &parset : parset containing spacing and offset parameters
/// @param[in] const Vector<std::string> beamNames : which offsets to get from the parset
/// @param[in] MVDirection centre : the pointing centre, which all offsets are relative to
/// @return Vector<MVDirection> : a MVDirection for each name in beamNames
Vector<MVDirection> loadBeamOffsets(const LOFAR::ParameterSet &parset, const Vector<std::string> beamNames,
                                          MVDirection centre) {
    
    Vector<MVDirection> centres (beamNames.size(), centre);

    ASKAPLOG_INFO_STR(logger, " -> looking for the feed spacing");
    Quantity qspacing = asQuantity(parset.getString("feeds.spacing"));
    double spacing = qspacing.getValue("rad");
    ASKAPLOG_INFO_STR(logger, "    beam spacing set to " << qspacing);       

    ASKAPLOG_INFO_STR(logger, " -> looking for a feed offset for each image");
    for (uint beam = 0; beam < beamNames.size(); ++beam) {
         const string parName = "feeds." + beamNames[beam];
         const Vector<double> xy(parset.getDoubleVector(parName));
         ASKAPCHECK(xy.size() == 2, "Expect two elements for each offset");
         // the shift appears to be positive in HA, so multiply by -1. Simulator.cc states:
         // "x direction is flipped to convert az-el type frame to ra-dec"
         centres[beam].shift(-xy[0]*spacing, xy[1]*spacing, casa::True);
         ASKAPLOG_INFO_STR(logger, " -> " << parName << " centre: " << centres[beam] );
    }
    return centres;
}


// functions in the linmos accumulator class

bool LinmosAccumulator::loadParset(const LOFAR::ParameterSet &parset) {

    const vector<string> inImgNames = parset.getStringVector("names", true);
    const vector<string> inWgtNames = parset.getStringVector("weights", vector<string>(), true);
    const string weightTypeName = parset.getString("weighttype");
    const string weightStateName = parset.getString("weightstate", "Corrected");

    const bool findMosaics = parset.getBool("findmosaics", false);

    // Check the input images
    ASKAPCHECK(inImgNames.size()>0, "Number of input images should be greater than 0");

    // Check weighting options. One of the following must be set:
    //  - weightTypeName==FromWeightImages: get weights from input weight images
    //    * the number of weight images and their shapes must match the input images
    //  - weightTypeName==FromPrimaryBeamModel: set weights using a Gaussian beam model 
    //    * the direction coordinate centre will be used as beam centre, unless ...
    //    * an output weight image will be written, so an output file name is required

    if (boost::iequals(weightTypeName, "FromWeightImages")) {
        itsWeightType = FROM_WEIGHT_IMAGES;
        ASKAPLOG_INFO_STR(logger, "Weights are coming from weight images");
    } else if (boost::iequals(weightTypeName, "FromPrimaryBeamModel")) {
        itsWeightType = FROM_BP_MODEL;
        ASKAPLOG_INFO_STR(logger, "Weights to be set using a Gaussian primary-beam models");
    } else {
        ASKAPLOG_ERROR_STR(logger, "Unknown weighttype " << weightTypeName);
        return false;
    }

    if (findMosaics) {

        ASKAPLOG_INFO_STR(logger, "Image names to be automatically generated. Searching...");
        // check for useless parameters
        if (parset.isDefined("outname") || parset.isDefined("outweight")) {
            ASKAPLOG_WARN_STR(logger, "  - output file names are specified in parset but ignored.");
        }
        if (parset.isDefined("nterms")) {
            ASKAPLOG_WARN_STR(logger, "  - nterms is specified in parset but ignored.");
        }

        findAndSetMosaics(inImgNames);

        ASKAPCHECK(itsInImgNameVecs.size() > 0, "No suitable mosaics found.");
        ASKAPLOG_INFO_STR(logger, itsInImgNameVecs.size() << " suitable mosaics found.");

    } else {

        string outImgName = parset.getString("outname");
        string outWgtName = parset.getString("outweight");

        // If reading weights from images, check the input for those
        if (itsWeightType == FROM_WEIGHT_IMAGES) {
            ASKAPCHECK(inImgNames.size()==inWgtNames.size(), "# weight images should equal # images");
        }
 
        // Check for taylor terms
 
        if (parset.isDefined("nterms")) {
 
            itsNumTaylorTerms = parset.getInt32("nterms");
            findAndSetTaylorTerms(inImgNames, inWgtNames, outImgName, outWgtName);
 
        } else {

            setSingleMosaic(inImgNames, inWgtNames, outImgName, outWgtName);

        }

    }

    if (itsWeightType == FROM_WEIGHT_IMAGES) {
 
        // if reading weights from images, check for inputs associated with other kinds of weighting
        if (parset.isDefined("feeds.centre") ||
            parset.isDefined("feeds.centreref") ||
            parset.isDefined("feeds.offsetsfile") ||
            parset.isDefined("feeds.names") ||
            parset.isDefined("feeds.spacing") ) {
            ASKAPLOG_WARN_STR(logger, "Beam information specified in parset but ignored. Using weight images");
        }
 
    } else if (itsWeightType == FROM_BP_MODEL) {

        // check for inputs associated with other kinds of weighting
        if (inWgtNames.size()>0) {
            ASKAPLOG_WARN_STR(logger, "Weight images specified in parset but ignored. Using a primary-beam model");
        }

    }

    // Check the initial weighting state of the input images

    if (boost::iequals(weightStateName, "Corrected")) {
        ASKAPLOG_INFO_STR(logger, "Input image state: Direction-dependent beams/weights have been divided out");
        itsWeightState = CORRECTED;
    } else if (boost::iequals(weightStateName, "Inherent")) {
        ASKAPLOG_INFO_STR(logger, "Input image state: natural primary-beam weighting of the visibilities is retained");
        itsWeightState = INHERENT;
    } else if (boost::iequals(weightStateName, "Weighted")) {
        ASKAPLOG_INFO_STR(logger, "Input image state: full primary-beam-squared weighting");
        itsWeightState = WEIGHTED;
    } else {
        ASKAPLOG_ERROR_STR(logger, "Unknown weightstyle " << weightStateName);
        return false;
    }

    if (parset.isDefined("regrid.method")) itsMethod = parset.getString("regrid.method");
    if (parset.isDefined("regrid.decimate")) itsDecimate = parset.getInt("regrid.decimate");
    if (parset.isDefined("regrid.replicate")) itsReplicate = parset.getBool("regrid.replicate");
    if (parset.isDefined("regrid.force")) itsForce = parset.getBool("regrid.force");

    if (parset.isDefined("psfref")) {
        ASKAPCHECK(parset.getUint("psfref")<inImgNames.size(), "PSF reference-image number is too large");
    }

    return true;

}

bool LinmosAccumulator::loadBeamCentres(const LOFAR::ParameterSet &parset,
                                        const accessors::IImageAccess &iacc,
                                        const string outImgName) {

    const vector<string> inImgNames = itsInImgNameVecs[outImgName];

    if (itsWeightType == FROM_BP_MODEL) {

        // if setting weights using beam models, check the input for extra information

        ASKAPLOG_INFO_STR(logger, "Looking for parset options associated with primary-beam models");

        MVDirection centre;
        bool centreDefined = false;

        // set the centre of the "feeds" offset parameters (e.g. the boresight of the PAF)
        if (parset.isDefined("feeds.centre")) {
            ASKAPLOG_INFO_STR(logger, "Found centre of the feeds to use in beam models:");
            const vector<string> feedsCentre(parset.getStringVector("feeds.centre"));
            ASKAPCHECK(feedsCentre.size()==2, " -> the feeds.centre vector should have 2 elements");
            centre = convertDir(feedsCentre[0], feedsCentre[1]);
            ASKAPLOG_INFO_STR(logger, " -> "<<feedsCentre<<", = "<<centre);
            centreDefined = true;
        }
        else if (parset.isDefined("feeds.centreref")) {
            uint centreref = parset.getInt("feeds.centreref");
            if ((centreref>=0) && (centreref<inImgNames.size())) {
                ASKAPLOG_INFO_STR(logger, "Using the reference pixel of input image "<<centreref<<
                    " as the centre of the feeds to use in beam models");
                const CoordinateSystem coordSys = iacc.coordSys(inImgNames[centreref]);
                const int DCpos = coordSys.findCoordinate(Coordinate::DIRECTION,-1);
                const DirectionCoordinate DC = coordSys.directionCoordinate(DCpos);
                DC.toWorld(centre,DC.referencePixel());
                ASKAPLOG_INFO_STR(logger, " -> "<<centre);
                centreDefined = true;
            }
            else {
                ASKAPLOG_WARN_STR(logger, "Found unsuitable centreref parameter: "<<centreref);
            }
        }

        // centres for each beam
        if (centreDefined) {

            if (parset.isDefined("feeds.offsetsfile")) {

                ASKAPLOG_INFO_STR(logger,  "Loading beam offsets from " << parset.getString("feeds.offsetsfile"));
                LOFAR::ParameterSet feed_parset(parset.getString("feeds.offsetsfile"));

                vector<string> beamNames;
                ASKAPLOG_INFO_STR(logger, " -> looking for feed names");
                if (parset.isDefined("feeds.names")) {
                    beamNames = parset.getStringVector("feeds.names", true);
                    ASKAPLOG_INFO_STR(logger,  "    using names given in the main parset");
                } else if (feed_parset.isDefined("feeds.names")) {
                    beamNames = feed_parset.getStringVector("feeds.names", true);
                    ASKAPLOG_INFO_STR(logger,  "    using names given in the feed-offset parset");
                }
                ASKAPCHECK(beamNames.size() > 0, "No beams specified");
                ASKAPCHECK(beamNames.size() == inImgNames.size(),
                   "Number of beams does not match number of input files");

                itsCentres = loadBeamOffsets(feed_parset, beamNames, centre);

                if (parset.isDefined("feeds.spacing")) {
                    ASKAPLOG_WARN_STR(logger, "Feed info specified in parset but ignored. Using offset file");
                }

            } else {

                ASKAPCHECK(parset.isDefined("names"), "No names specified in parset");
                itsCentres = loadBeamOffsets(parset, parset.getStringVector("names", true), centre);

            }

        } else {
            ASKAPLOG_WARN_STR(logger, "Centre of the feeds not found. Setting beam centres to input ref. pixels");
        }

    }

    return true;

}

void LinmosAccumulator::setSingleMosaic(const vector<string> &inImgNames, const vector<string> &inWgtNames,
                                        const string &outImgName, const string &outWgtName) {

    // set some variables for the sensitivity image searches
    string image_tag = "image", restored_tag = ".restored", tmpName;
    itsDoSensitivity = true; // set false if any sensitivity images are missing or if not an image* mosaic

    // Check the input images
    for (size_t img = 0; img < inImgNames.size(); ++img) {

        // make sure the output image will not be overwritten
        ASKAPCHECK(inImgNames[img]!=outImgName, "Output image, "<<outImgName<<", is present among the inputs");

        // if this is an "image*" file, see if there is an appropriate sensitivity image 
        if (itsDoSensitivity) {
            tmpName = inImgNames[img];
            size_t image_pos = tmpName.find(image_tag);
            // if the file starts with image_tag, look for a sensitivity image
            if (image_pos == 0) {
                tmpName.replace(image_pos, image_tag.length(), "sensitivity");
                // remove any ".restored" sub-string from the file name
                size_t restored_pos = tmpName.find(restored_tag);
                if (restored_pos != string::npos) {
                    tmpName.replace(restored_pos, restored_tag.length(), "");
                }
                if (boost::filesystem::exists(tmpName)) {
                    itsInSenNameVecs[outImgName].push_back(tmpName);
                } else {
                    ASKAPLOG_WARN_STR(logger, "Cannot find file "<<tmpName<<" . Ignoring sensitivities.");
                    itsDoSensitivity = false;
                }
            } else {
                ASKAPLOG_WARN_STR(logger, "Input not an image* file. Ignoring sensitivities.");
                itsDoSensitivity = false;
            }
        }

    }

    // set a single key for the various file-name maps
    itsOutWgtNames[outImgName] = outWgtName;
    itsInImgNameVecs[outImgName] = inImgNames;
    if (itsWeightType == FROM_WEIGHT_IMAGES) {
        itsInWgtNameVecs[outImgName] = inWgtNames;
    }
    if (itsDoSensitivity) {
        itsGenSensitivityImage[outImgName] = true;
        // set an output sensitivity file name
        tmpName = outImgName;
        tmpName.replace(0, image_tag.length(), "sensitivity");
        // remove any ".restored" sub-string from the weights file name
        size_t restored_pos = tmpName.find(restored_tag);
        if (restored_pos != string::npos) {
            tmpName.replace(restored_pos, restored_tag.length(), "");
        }
        itsOutSenNames[outImgName] = tmpName;
    } else {
        itsGenSensitivityImage[outImgName] = false;
        // if some but not all sensitivity images were found, remove this key from itsInSenNameVecs
        if (itsInSenNameVecs.find(outImgName)!=itsInSenNameVecs.end()) {
            itsInSenNameVecs.erase(outImgName);
        }
    }

} // LinmosAccumulator::setSingleMosaic()

void LinmosAccumulator::findAndSetTaylorTerms(const vector<string> &inImgNames, const vector<string> &inWgtNames,
                                              const string &outImgNameOrig, const string &outWgtNameOrig) {

    ASKAPLOG_INFO_STR(logger, "Looking for "<<itsNumTaylorTerms<<" taylor terms");
    ASKAPCHECK(itsNumTaylorTerms>=0, "Number of taylor terms should be greater than or equal to 0");

    size_t pos0, pos1;
    pos0 = outImgNameOrig.find(itsTaylorTag);
    ASKAPCHECK(pos0!=string::npos, "Cannot find "<<itsTaylorTag<<" in output file "<<outImgNameOrig);
    pos1 = outImgNameOrig.find(itsTaylorTag, pos0+1); // make sure there aren't multiple entries.
    ASKAPCHECK(pos1==string::npos, "There are multiple  "<<itsTaylorTag<<" strings in output file "<<outImgNameOrig);

    // set some variables for the sensitivity image searches
    string image_tag = "image", restored_tag = ".restored", tmpName;
    itsDoSensitivity = true; // set false if any sensitivity images are missing or if not an image* mosaic

    for (int n = 0; n < itsNumTaylorTerms; ++n) {

        string outImgName = outImgNameOrig;
        string outWgtName = outWgtNameOrig;
        const string taylorN = "taylor." + boost::lexical_cast<string>(n);

        // set a new key for the various output file-name maps
        outImgName.replace(outImgName.find(itsTaylorTag), itsTaylorTag.length(), taylorN);
        outWgtName.replace(outWgtName.find(itsTaylorTag), itsTaylorTag.length(), taylorN);
        itsOutWgtNames[outImgName] = outWgtName;

        for (uint img = 0; img < inImgNames.size(); ++img) {

            // do some tests
            string inImgName = inImgNames[img]; // short cut
            pos0 = inImgName.find(itsTaylorTag);
            ASKAPCHECK(pos0!=string::npos, "Cannot find "<<itsTaylorTag<<" in input file "<<inImgName);
            pos1 = inImgName.find(itsTaylorTag, pos0+1); // make sure there aren't multiple entries.
            ASKAPCHECK(pos1==string::npos, "There are multiple "<<itsTaylorTag<<" strings in input file "<<inImgName);

            // set a new key for the input file-name-vector map
            inImgName.replace(pos0, itsTaylorTag.length(), taylorN);
            itsInImgNameVecs[outImgName].push_back(inImgName);

            // Check the input image
            ASKAPCHECK(inImgName!=outImgName, "Output image, "<<outImgName<<", is present among the inputs");

            if (itsWeightType == FROM_WEIGHT_IMAGES) {
                // do some tests
                string inWgtName = inWgtNames[img]; // short cut
                pos0 = inWgtName.find(itsTaylorTag);
                ASKAPCHECK(pos0!=string::npos, "Cannot find "<<itsTaylorTag<< " in input weight file "<<inWgtName);
                pos1 = inWgtName.find(itsTaylorTag, pos0+1); // make sure there aren't multiple entries.
                ASKAPCHECK(pos1==string::npos, "There are multiple " << itsTaylorTag <<
                                               " strings in input file "<<inWgtName);

                // set a new key for the input weights file-name-vector map
                inWgtName.replace(pos0, itsTaylorTag.length(), taylorN);
                itsInWgtNameVecs[outImgName].push_back(inWgtName);

                // Check the input weights image
                ASKAPCHECK(inWgtName!=outWgtName, "Output wgt image, "<<outWgtName<<", is among the inputs");
            }

            // if this is an "image*" file, see if there is an appropriate sensitivity image 
            if (itsDoSensitivity) {
                tmpName = inImgName;
                size_t image_pos = tmpName.find(image_tag);
                // if the file starts with image_tag, look for a sensitivity image
                if (image_pos == 0) {
                    tmpName.replace(image_pos, image_tag.length(), "sensitivity");
                    // remove any ".restored" sub-string from the file name
                    size_t restored_pos = tmpName.find(restored_tag);
                    if (restored_pos != string::npos) {
                        tmpName.replace(restored_pos, restored_tag.length(), "");
                    }
                    if (boost::filesystem::exists(tmpName)) {
                        itsInSenNameVecs[outImgName].push_back(tmpName);
                    } else {
                        ASKAPLOG_WARN_STR(logger, "Cannot find file "<<tmpName<<" . Ignoring sensitivities.");
                        itsDoSensitivity = false;
                    }
                } else {
                    ASKAPLOG_WARN_STR(logger, "Input not an image* file. Ignoring sensitivities.");
                    itsDoSensitivity = false;
                }
            }

        } // img loop (input image)

        // check whether any sensitivity images were found
        if (itsDoSensitivity) {
            itsGenSensitivityImage[outImgName] = true;
            // set an output sensitivity file name
            tmpName = outImgName;
            tmpName.replace(0, image_tag.length(), "sensitivity");
            // remove any ".restored" sub-string from the weights file name
            size_t restored_pos = tmpName.find(restored_tag);
            if (restored_pos != string::npos) {
                tmpName.replace(restored_pos, restored_tag.length(), "");
            }
            itsOutSenNames[outImgName] = tmpName;
        } else {
            itsGenSensitivityImage[outImgName] = false;
            // if some but not all sensitivity images were found, remove this key from itsInSenNameVecs
            if (itsInSenNameVecs.find(outImgName)!=itsInSenNameVecs.end()) {
                itsInSenNameVecs.erase(outImgName);
            }
        }

    } // n loop (taylor term)

} // void LinmosAccumulator::findAndSetTaylorTerms()

void LinmosAccumulator::findAndSetMosaics(const vector<string> &imageTags) {

    vector<string> prefixes;
    prefixes.push_back("image");
    prefixes.push_back("residual");
    //prefixes.push_back("weights"); // these need to be handled separately
    //prefixes.push_back("sensitivity"); // these need to be handled separately
    //prefixes.push_back("mask");

    // if this directory name changes from "./", the erase call below may also need to change
    boost::filesystem::path p (".");

    typedef vector<boost::filesystem::path> path_vec;
    path_vec v;

    copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(v));

    // find mosaics by looking for images that contain one of the tags. Then see which of those contain all tags.
    const string searchTag = imageTags[0];

    for (path_vec::const_iterator it (v.begin()); it != v.end(); ++it) {

        // set name of the current file name and remove "./"
        string name = it->string();
        name.erase(0,2);

        // make sure this is a directory
        // a sym link to a directory will pass this test
        if (!boost::filesystem::is_directory(*it)) {
            //ASKAPLOG_INFO_STR(logger, name << " is not a directory. Ignoring.");
            continue;
        }

        // see if the name contains the desired tag (i.e., contains the first tag in "names")
        size_t pos = name.find(searchTag);
        if (pos == string::npos) {
            //ASKAPLOG_INFO_STR(logger, name << " is not a match. Ignoring.");
            continue;
        }

        // set some variables for problem sub-strings
        string restored_tag = ".restored";
        size_t restored_pos;

        // see if the name contains a desired prefix, and if so, check the other input names and weights
        int full_set = 0, full_wgt_set = 0;
        string mosaicName = name, nextName = name, tmpName;
        for (vector<string>::const_iterator pre (prefixes.begin()); pre != prefixes.end(); ++pre) {
            if (name.find(*pre) == 0) {

                // both of these must remain set to 1 for this mosaic to be established
                full_set = 1;
                full_wgt_set = 1;

                // set the output mosaic name
                mosaicName.replace(pos, searchTag.length(), itsMosaicTag);

                // file seems good, but check that it is present in all input images
                for (uint img = 0; img < imageTags.size(); ++img) {

                    // name is initially set for image 0
                    nextName = name;
                    // replace the image 0 tag with the current image's tag
                    if (img > 0) {
                        nextName.replace(pos, searchTag.length(), imageTags[img]);
                        // check that the file exists
                        if (!boost::filesystem::exists(nextName)) {
                            full_set = -1;
                            break;
                        }
                    }
                    // add the image to this mosaics inputs
                    itsInImgNameVecs[mosaicName].push_back(nextName);

                    // see if there is an appropriate sensitivity image
                    tmpName = nextName;
                    tmpName.replace(0, (*pre).length(), "sensitivity");
                    // remove any ".restored" sub-string from the weights file name
                    restored_pos = tmpName.find(restored_tag);
                    if (restored_pos != string::npos) {
                        tmpName.replace(restored_pos, restored_tag.length(), "");
                    }
                    if (boost::filesystem::exists(tmpName)) {
                        itsInSenNameVecs[mosaicName].push_back(tmpName);
                    }

                    // look for weights image if required (weights are not needed when combining sensitivity images)
                    if (itsWeightType == FROM_WEIGHT_IMAGES) {
                        // replace the prefix with "weights"
                        nextName.replace(0, (*pre).length(), "weights");
                        // remove any ".restored" sub-string from the weights file name
                        restored_pos = nextName.find(restored_tag);
                        if (restored_pos != string::npos) {
                            nextName.replace(restored_pos, restored_tag.length(), "");
                        }
                        // check that the file exists
                        if (!boost::filesystem::exists(nextName)) {
                            full_wgt_set = -1;
                            break;
                        }
                        // add the file to this mosaics inputs
                        itsInWgtNameVecs[mosaicName].push_back(nextName);
                    }

                }

                // set the output weights image name
                // replace the mosaic prefix with "weights"
                nextName = mosaicName;
                nextName.replace(0, (*pre).length(), "weights");
                // remove any ".restored" sub-string from the weights file name
                restored_pos = nextName.find(restored_tag);
                if (restored_pos != string::npos) {
                    nextName.replace(restored_pos, restored_tag.length(), "");
                }
                itsOutWgtNames[mosaicName] = nextName;

                itsGenSensitivityImage[mosaicName] = false;
                if (itsInSenNameVecs.find(mosaicName)!=itsInSenNameVecs.end()) { // if key is found
                    if (itsInImgNameVecs[mosaicName].size()==itsInSenNameVecs[mosaicName].size()) {
                        itsGenSensitivityImage[mosaicName] = true;
                        // set an output sensitivity file name
                        tmpName = mosaicName;
                        tmpName.replace(0, (*pre).length(), "sensitivity");
                        // remove any ".restored" sub-string from the weights file name
                        restored_pos = tmpName.find(restored_tag);
                        if (restored_pos != string::npos) {
                            tmpName.replace(restored_pos, restored_tag.length(), "");
                        }
                        itsOutSenNames[mosaicName] = tmpName;
                    } else {
                        itsInSenNameVecs.erase(mosaicName);
                    }
                }

                break; // found the prefix, so leave the loop

            }
        }

        if (full_set==0) {
            // this file did not have a relevant prefix, so just move on
            continue;
        }

        if ((full_set == -1) || ((itsWeightType == FROM_WEIGHT_IMAGES) && (full_wgt_set == -1))) {
            // this file did have a relevant prefix, but failed
            if (full_set == -1) {
                ASKAPLOG_INFO_STR(logger, mosaicName << " does not have a full set of input files. Ignoring.");
            }
            if ((itsWeightType == FROM_WEIGHT_IMAGES) && (full_wgt_set == -1)) {
                ASKAPLOG_INFO_STR(logger, mosaicName << " does not have a full set of weights files. Ignoring.");
            }

            // if any of these were started for the current failed key, clean up and move on
            if (itsOutWgtNames.find(mosaicName)!=itsOutWgtNames.end()) itsOutWgtNames.erase(mosaicName);
            if (itsOutSenNames.find(mosaicName)!=itsOutSenNames.end()) itsOutSenNames.erase(mosaicName);
            if (itsInImgNameVecs.find(mosaicName)!=itsInImgNameVecs.end()) itsInImgNameVecs.erase(mosaicName);
            if (itsInWgtNameVecs.find(mosaicName)!=itsInWgtNameVecs.end()) itsInWgtNameVecs.erase(mosaicName);
            if (itsInSenNameVecs.find(mosaicName)!=itsInSenNameVecs.end()) itsInSenNameVecs.erase(mosaicName);

            continue;
        }

        // double check the size of the various maps and vectors. These should have been caught already
        ASKAPCHECK(itsInImgNameVecs.size()==itsOutWgtNames.size(), mosaicName << "Inconsistent name maps.");
        if (itsWeightType == FROM_WEIGHT_IMAGES) {
            ASKAPCHECK(itsInImgNameVecs.size()==itsInWgtNameVecs.size(), mosaicName <<
                       "Something has gone wrong with automatic mosaic search. Inconsistent name maps.");
            ASKAPCHECK(itsInImgNameVecs[mosaicName].size()==itsInWgtNameVecs[mosaicName].size(), mosaicName <<
                       "Something has gone wrong with automatic mosaic search. Inconsistent name vectors.");
        }

        ASKAPLOG_INFO_STR(logger, mosaicName << " seems complete. Mosaicking.");

        // it is possible that there may be duplicate itsOutWgtNames/itsOutSenNames (e.g. for image.* and residual.*)
        // check that the input is the same for these duplicates, and then only write once
        // if this is common, we should be avoiding more that just duplicate output
        string mosaicOrig;
        itsOutWgtDuplicates[mosaicName] = false;
        for(map<string,string>::iterator ii=itsOutWgtNames.begin(); ii!=itsOutWgtNames.find(mosaicName); ++ii) {
            if (itsOutWgtNames[mosaicName].compare((*ii).second) == 0) {
                itsOutWgtDuplicates[mosaicName] = true;
                mosaicOrig = (*ii).first;
                break;
            }
        }

        // if this is a duplicate, just remove it. Can't with weights because we need them unaveraged
        if (itsOutSenNames.find(mosaicName)!=itsOutSenNames.end()) {
            for(map<string,string>::iterator ii=itsOutSenNames.begin(); ii!=itsOutSenNames.find(mosaicName); ++ii) {
                if (itsOutSenNames[mosaicName].compare((*ii).second) == 0) {
                    ASKAPLOG_INFO_STR(logger, "  - sensitivity image done in an earlier mosaic. Will not redo here.");
                    itsGenSensitivityImage[mosaicName] = false;
                    itsOutSenNames.erase(mosaicName);
                    itsInSenNameVecs.erase(mosaicName);
                    break;
                }
            }
        }
        if (itsOutSenNames.find(mosaicName)!=itsOutSenNames.end()) {
            ASKAPLOG_INFO_STR(logger, "  - sensitivity images found. Generating mosaic sens. image.");
        }

    } // it loop (over potential images in this directory)

} // void LinmosAccumulator::findAndSetMosaics()

bool LinmosAccumulator::outputBufferSetupRequired(void) {
    return ( itsOutBuffer.shape().nelements() == 0 );
}

void LinmosAccumulator::setInputParameters(const string& inImgName, const accessors::IImageAccess& iacc, const int n) {
    // set the input coordinate system and shape
    itsInCoordSys = iacc.coordSys(inImgName);
    itsInShape = iacc.shape(inImgName);

    if (itsWeightType == FROM_BP_MODEL) {
        // set the centre of the beam
        if ( itsCentres.size() > 0 ) {
            itsInCentre = itsCentres[n];
        } else {
            // no other information, so set the centre of the beam to be the reference pixel
            const int dcPos = itsInCoordSys.findCoordinate(Coordinate::DIRECTION,-1);
            const DirectionCoordinate inDC = itsInCoordSys.directionCoordinate(dcPos);
            inDC.toWorld(itsInCentre,inDC.referencePixel());
        }
    }
}

void LinmosAccumulator::setOutputParameters(const vector<string>& inImgNames, const accessors::IImageAccess& iacc) {

    ASKAPLOG_INFO_STR(logger, "Determining output image properties based on the overlap of input images");
    ASKAPCHECK(inImgNames.size()>0, "Number of input images should be greater that 0");

    const IPosition refShape = iacc.shape(inImgNames[0]);
    ASKAPDEBUGASSERT(refShape.nelements() >= 2);
    const CoordinateSystem refCS = iacc.coordSys(inImgNames[0]);
    const int dcPos = refCS.findCoordinate(Coordinate::DIRECTION,-1);
    const DirectionCoordinate refDC = refCS.directionCoordinate(dcPos);
    IPosition refBLC(refShape.nelements(),0);
    IPosition refTRC(refShape);       
    for (uInt dim=0; dim<refShape.nelements(); ++dim) {
        --refTRC(dim); // these are added back later. Is this just to deal with degenerate axes?
    }
    ASKAPDEBUGASSERT(refBLC.nelements() >= 2);
    ASKAPDEBUGASSERT(refTRC.nelements() >= 2);
 
    IPosition tempBLC = refBLC;
    IPosition tempTRC = refTRC;

    // Loop over input images, converting their image bounds to the ref system 
    // and expanding the new overlapping image bounds where appropriate.
    for (uint img = 1; img < inImgNames.size(); ++img ) {

        // short cuts
        const string inImgName = inImgNames[img];

        // 
        itsInShape = iacc.shape(inImgName);
        itsInCoordSys = iacc.coordSys(inImgName);

        // test to see if the loaded coordinate system is close enough to the reference system for merging
        ASKAPCHECK(coordinatesAreConsistent(refCS), "Input images have inconsistent coordinate systems");
        // could also test whether they are equal and set a regrid tag to false if all of them are

        Vector<IPosition> corners = convertImageCornersToRef(refDC);

        const IPosition newBLC = corners[0];
        const IPosition newTRC = corners[1];
        ASKAPDEBUGASSERT(newBLC.nelements() >= 2);
        ASKAPDEBUGASSERT(newTRC.nelements() >= 2);
        for (casa::uInt dim=0; dim<2; ++dim) {
            if (newBLC(dim) < tempBLC(dim)) {
                tempBLC(dim) = newBLC(dim);
            }
            if (newTRC(dim) > tempTRC(dim)) {
                tempTRC(dim) = newTRC(dim);
            }
        }
 
    }

    itsOutShape = refShape;
    itsOutShape(0) = tempTRC(0) - tempBLC(0) + 1;
    itsOutShape(1) = tempTRC(1) - tempBLC(1) + 1;
    ASKAPDEBUGASSERT(itsOutShape(0) > 0);
    ASKAPDEBUGASSERT(itsOutShape(1) > 0);       
    Vector<Double> refPix = refDC.referencePixel();
    refPix[0] -= Double(tempBLC(0) - refBLC(0));
    refPix[1] -= Double(tempBLC(1) - refBLC(1));
    DirectionCoordinate newDC(refDC);
    newDC.setReferencePixel(refPix);

    // set up a coord system for the merged images
    itsOutCoordSys = refCS;
    itsOutCoordSys.replaceCoordinate(newDC, dcPos);

}

void LinmosAccumulator::initialiseOutputBuffers(void) {
    // set up temporary images needed for regridding (which is done on a plane-by-plane basis so ignore other dims)

    // set up the coord. sys.
    int dcPos = itsOutCoordSys.findCoordinate(Coordinate::DIRECTION,-1);
    ASKAPCHECK(dcPos>=0, "Cannot find the directionCoordinate");
    const DirectionCoordinate dcTmp = itsOutCoordSys.directionCoordinate(dcPos);
    CoordinateSystem cSysTmp;
    cSysTmp.addCoordinate(dcTmp);

    // set up the shape
    Vector<Int> shapePos = itsOutCoordSys.pixelAxes(dcPos);
    // check that the length is equal to 2 and the both elements are >= 0
    ASKAPCHECK(shapePos.nelements()>=2, "Cannot find the directionCoordinate");
    ASKAPCHECK((shapePos[0]==0 && shapePos[1]==1) || (shapePos[1]==0 && shapePos[0]==1),
               "Linmos currently requires the direction coordinates to come before any others");

    IPosition shape = IPosition(2,itsOutShape(shapePos[0]),itsOutShape(shapePos[1]));

    // apparently the +100 forces it to use the memory
    double maxMemoryInMB = double(shape.product()*sizeof(float))/1024./1024.+100;
    itsOutBuffer = TempImage<float>(shape, cSysTmp, maxMemoryInMB);
    ASKAPCHECK(itsOutBuffer.shape().nelements()>0, "Output buffer does not appear to be set");
    if (itsWeightType == FROM_WEIGHT_IMAGES) {
        itsOutWgtBuffer = TempImage<float>(shape, cSysTmp, maxMemoryInMB);
        ASKAPCHECK(itsOutWgtBuffer.shape().nelements()>0, "Output weights buffer does not appear to be set");
    }
    if (itsDoSensitivity) {
        itsOutSnrBuffer = TempImage<float>(shape, cSysTmp, maxMemoryInMB);
        ASKAPCHECK(itsOutSnrBuffer.shape().nelements()>0, "Output sensitivity buffer does not appear to be set");
    }
}

void LinmosAccumulator::initialiseInputBuffers() {
    // set up temporary images needed for regridding (which is done on a plane-by-plane basis so ignore other dims)

    // set up a coord. sys. the planes
    int dcPos = itsInCoordSys.findCoordinate(Coordinate::DIRECTION,-1);
    ASKAPCHECK(dcPos>=0, "Cannot find the directionCoordinate");
    const DirectionCoordinate dc = itsInCoordSys.directionCoordinate(dcPos);
    CoordinateSystem cSys;
    cSys.addCoordinate(dc);

    // set up the shape
    Vector<Int> shapePos = itsInCoordSys.pixelAxes(dcPos);
    // check that the length is equal to 2 and the both elements are >= 0

    IPosition shape = IPosition(2,itsInShape(shapePos[0]),itsInShape(shapePos[1]));

    double maxMemoryInMB = double(shape.product()*sizeof(float))/1024./1024.+100;
    itsInBuffer = TempImage<float>(shape,cSys,maxMemoryInMB);
    ASKAPCHECK(itsInBuffer.shape().nelements()>0, "Input buffer does not appear to be set");
    if (itsWeightType == FROM_WEIGHT_IMAGES) {
        itsInWgtBuffer = TempImage<float>(shape,cSys,maxMemoryInMB);       
        ASKAPCHECK(itsInWgtBuffer.shape().nelements()>0, "Input weights buffer does not appear to be set");
    }
    if (itsDoSensitivity) {
        itsInSenBuffer = TempImage<float>(shape,cSys,maxMemoryInMB);       
        itsInSnrBuffer = TempImage<float>(shape,cSys,maxMemoryInMB);       
        ASKAPCHECK(itsInSnrBuffer.shape().nelements()>0, "Input sensitivity buffer does not appear to be set");
    }

}

void LinmosAccumulator::initialiseRegridder() {
    ASKAPLOG_INFO_STR(logger, "Initialising regridder for " << itsMethod << " interpolation");
    itsAxes = IPosition::makeAxisPath(itsOutBuffer.shape().nelements());
    itsEmethod = Interpolate2D::stringToMethod(itsMethod);
}

void LinmosAccumulator::loadInputBuffers(const scimath::MultiDimArrayPlaneIter& planeIter,
                                         Array<float>& inPix, Array<float>& inWgtPix, Array<float>& inSenPix) {
    itsInBuffer.put(planeIter.getPlane(inPix));
    if (itsWeightType == FROM_WEIGHT_IMAGES) {
        itsInWgtBuffer.put(planeIter.getPlane(inWgtPix));
    }
    if (itsDoSensitivity) {
        // invert sensitivities before regridding to avoid artefacts at sharp edges in the sensitivity image
        itsInSenBuffer.put(planeIter.getPlane(inSenPix));
        float sensitivity;
        float minVal, maxVal;
        IPosition minPos, maxPos;
        minMax(minVal,maxVal,minPos,maxPos,inSenPix);
        float senCutoff = itsCutoff * maxVal; // inSenPix is prop. to image sigma/gain
        IPosition pos(2);
        for (int x=0; x<inSenPix.shape()[0];++x) {
            for (int y=0; y<inSenPix.shape()[1];++y) {
                pos[0] = x;
                pos[1] = y;
                sensitivity = itsInSenBuffer.getAt(pos);
                if (sensitivity>senCutoff) {
                    itsInSnrBuffer.putAt(1.0 / (sensitivity * sensitivity), pos);
                } else {
                    itsInSnrBuffer.putAt(0.0, pos);
                }
            }
        }
    }
}

void LinmosAccumulator::regrid() {
    ASKAPLOG_INFO_STR(logger, " - regridding with dec="<<itsDecimate<<" rep="<<itsReplicate<<" force="<<itsForce);
    ASKAPCHECK(itsOutBuffer.shape().nelements()>0, "Output buffer does not appear to be set");
    itsRegridder.regrid(itsOutBuffer, itsEmethod, itsAxes, itsInBuffer, itsReplicate, itsDecimate, false, itsForce);
    if (itsWeightType == FROM_WEIGHT_IMAGES) {
        itsRegridder.regrid(itsOutWgtBuffer, itsEmethod, itsAxes, itsInWgtBuffer, itsReplicate, itsDecimate, false,
                            itsForce);
    }
    if (itsDoSensitivity) {
        itsRegridder.regrid(itsOutSnrBuffer, itsEmethod, itsAxes, itsInSnrBuffer, itsReplicate, itsDecimate, false,
                            itsForce);
    }
cout << "finished regrid. ";
}

void LinmosAccumulator::accumulatePlane(Array<float>& outPix, Array<float>& outWgtPix,
                                        Array<float>& outSenPix, const IPosition& curpos) {

    // copy the pixel iterator containing all dimensions
    IPosition fullpos(curpos);
    // set a pixel iterator that does not have the higher dimensions
    IPosition pos(2);

    // set the weights, either to those read in or using the primary-beam model
    TempImage<float> wgtBuffer;
    if (itsWeightType == FROM_WEIGHT_IMAGES) {
        wgtBuffer = itsOutWgtBuffer;
    } else {

        Vector<double> pixel(2,0.);
        MVDirection world0, world1;
        float offsetBeam, pb;

        // get coordinates of the spectral axis and the current frequency
        const int scPos = itsInCoordSys.findCoordinate(Coordinate::SPECTRAL,-1);
        const SpectralCoordinate inSC = itsInCoordSys.spectralCoordinate(scPos);
        int chPos = itsInCoordSys.pixelAxes(scPos)[0];
        const float freq = inSC.referenceValue()[0] + (curpos[chPos] - inSC.referencePixel()[0]) * inSC.increment()[0];

        // set FWHM for the current beam
        // Removing the factor of 1.22 gives a good match to the simultation weight images
        //const float fwhm = 1.22*3e8/freq/12;
        const float fwhm = 3e8/freq/12;

        // get coordinates of the direction axes
        const int dcPos = itsInCoordSys.findCoordinate(Coordinate::DIRECTION,-1);
        const DirectionCoordinate inDC = itsInCoordSys.directionCoordinate(dcPos);
        const DirectionCoordinate outDC = itsOutCoordSys.directionCoordinate(dcPos);

        // set the centre of the input beam (needs to be more flexible -- and correct...)
        inDC.toWorld(world0,inDC.referencePixel());

        // apparently the +100 forces it to use the memory
        double maxMemoryInMB = double(itsOutBuffer.shape().product()*sizeof(float))/1024./1024.+100;
        wgtBuffer = TempImage<float>(itsOutBuffer.shape(), itsOutBuffer.coordinates(), maxMemoryInMB);

        // step through the pixels, setting the weights (power primary beam squared)
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                pos[0] = x;
                pos[1] = y;

                // get the current pixel location and distance from beam centre
                pixel[0] = double(x);
                pixel[1] = double(y);
                outDC.toWorld(world1,pixel);
                offsetBeam = world0.separation(world1);

                // set the weight
                pb = exp(-offsetBeam*offsetBeam*4.*log(2.)/fwhm/fwhm);
                wgtBuffer.putAt(pb * pb, pos);

            }
        }

    }

    float minVal, maxVal;
    IPosition minPos, maxPos;
    minMax(minVal,maxVal,minPos,maxPos,wgtBuffer);
    float wgtCutoff = itsCutoff * itsCutoff * maxVal; // wgtBuffer is prop. to image (gain/sigma)^2

    // Accumulate the pixels of this slice.
    // Could restrict it (and the regrid) to a smaller region of interest.
    if (itsWeightState == CORRECTED) {
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                pos[0] = x;
                pos[1] = y;
                if (wgtBuffer.getAt(pos)>=wgtCutoff) {
                    outPix(fullpos)    = outPix(fullpos)    + itsOutBuffer.getAt(pos) * wgtBuffer.getAt(pos);
                    outWgtPix(fullpos) = outWgtPix(fullpos) + wgtBuffer.getAt(pos);
                }
            }
        }
    } else if (itsWeightState == INHERENT) {
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                pos[0] = x;
                pos[1] = y;
                if (wgtBuffer.getAt(pos)>=wgtCutoff) {
                    outPix(fullpos)    = outPix(fullpos)    + itsOutBuffer.getAt(pos) * sqrt(wgtBuffer.getAt(pos));
                    outWgtPix(fullpos) = outWgtPix(fullpos) + wgtBuffer.getAt(pos);
                }
            }
        }
    } else if (itsWeightState == WEIGHTED) {
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                pos[0] = x;
                pos[1] = y;
                if (wgtBuffer.getAt(pos)>=wgtCutoff) {
                    outPix(fullpos)    = outPix(fullpos)    + itsOutBuffer.getAt(pos);
                    outWgtPix(fullpos) = outWgtPix(fullpos) + wgtBuffer.getAt(pos);
                }
            }
        }
    }
    // Accumulate sensitivity for this slice.
    if (itsDoSensitivity) {
        float invVariance;
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                pos[0] = x;
                pos[1] = y;
                invVariance = itsOutSnrBuffer.getAt(pos);
                if (wgtBuffer.getAt(pos)>=wgtCutoff) {
                    outSenPix(fullpos) = outSenPix(fullpos) + invVariance;
                }
            }
        }
    }

}

void LinmosAccumulator::accumulatePlane(Array<float>& outPix, Array<float>& outWgtPix, Array<float>& outSenPix,
                                        const Array<float>& inPix, const Array<float>& inWgtPix,
                                        const Array<float>& inSenPix, const IPosition& curpos) {

    ASKAPASSERT(inPix.shape() == outPix.shape());

    // copy the pixel iterator containing all dimensions
    IPosition fullpos(curpos);
    // set up an indexing vector for the weights. If weight images are used, these are as in the image.
    IPosition wgtpos(curpos);

    Array<float> wgtPix;

    // set the weights, either to those read in or using the primary-beam model
    if (itsWeightType == FROM_WEIGHT_IMAGES) {
        wgtPix.reference(inWgtPix);
    } else {

        Vector<double> pixel(2,0.);
        MVDirection world;
        float offsetBeam, pb;

        // get coordinates of the spectral axis and the current frequency
        const int scPos = itsInCoordSys.findCoordinate(Coordinate::SPECTRAL,-1);
        const SpectralCoordinate inSC = itsInCoordSys.spectralCoordinate(scPos);
        int chPos = itsInCoordSys.pixelAxes(scPos)[0];
        const float freq = inSC.referenceValue()[0] + (curpos[chPos] - inSC.referencePixel()[0]) * inSC.increment()[0];

        // set FWHM for the current beam
        // Removing the factor of 1.22 gives a good match to the simultation weight images
        //const float fwhm = 1.22*3e8/freq/12;
        const float fwhm = 3e8/freq/12;

        // get coordinates of the direction axes
        const int dcPos = itsInCoordSys.findCoordinate(Coordinate::DIRECTION,-1);
        const DirectionCoordinate outDC = itsOutCoordSys.directionCoordinate(dcPos);

        // set the higher-order dimension to zero, as weights are on a 2D plane
        for (uInt dim=0; dim<curpos.nelements(); ++dim) {
            wgtpos[dim] = 0;
        }
        // set the array
        wgtPix = Array<float>(itsInShape);

        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                wgtpos[0] = x;
                wgtpos[1] = y;

                // get the current pixel location and distance from beam centre
                pixel[0] = double(x);
                pixel[1] = double(y);
                outDC.toWorld(world,pixel);
                offsetBeam = itsInCentre.separation(world);

                // set the weight
                pb = exp(-offsetBeam*offsetBeam*4.*log(2.)/fwhm/fwhm);
                wgtPix(wgtpos) = pb * pb;

            }
        }
    }

    float minVal, maxVal;
    IPosition minPos, maxPos;
    minMax(minVal,maxVal,minPos,maxPos,wgtPix);
    float wgtCutoff = itsCutoff * itsCutoff * maxVal; // wgtPix is prop. to image (gain/sigma)^2

    if (itsWeightState == CORRECTED) {
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                wgtpos[0] = x;
                wgtpos[1] = y;
                if (wgtPix(wgtpos)>=wgtCutoff) {
                    outPix(fullpos)    = outPix(fullpos)    + inPix(fullpos) * wgtPix(wgtpos);
                    outWgtPix(fullpos) = outWgtPix(fullpos) + wgtPix(wgtpos);
                }
            }
        }
    } else if (itsWeightState == INHERENT) {
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                wgtpos[0] = x;
                wgtpos[1] = y;
                if (wgtPix(wgtpos)>=wgtCutoff) {
                    outPix(fullpos)    = outPix(fullpos)    + inPix(fullpos) * sqrt(wgtPix(wgtpos));
                    outWgtPix(fullpos) = outWgtPix(fullpos) + wgtPix(wgtpos);
                }
            }
        }
    } else if (itsWeightState == WEIGHTED) {
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                wgtpos[0] = x;
                wgtpos[1] = y;
                if (wgtPix(wgtpos)>=wgtCutoff) {
                    outPix(fullpos)    = outPix(fullpos)    + inPix(fullpos);
                    outWgtPix(fullpos) = outWgtPix(fullpos) + wgtPix(wgtpos);
                }
            }
        }
    }
    // Accumulate sensitivity for this slice.
    if (itsDoSensitivity) {
        double sensitivity;
        minMax(minVal,maxVal,minPos,maxPos,inSenPix);
        float senCutoff = itsCutoff * maxVal; // inSenPix is prop. to image sigma/gain
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                sensitivity = inSenPix(fullpos);
                if (wgtPix(wgtpos)>=wgtCutoff && sensitivity>senCutoff) {
                    outSenPix(fullpos) = outSenPix(fullpos) + 1.0 / (sensitivity * sensitivity);
                }
            }
        }
    }

}

void LinmosAccumulator::deweightPlane(Array<float>& outPix, const Array<float>& outWgtPix,
                                      Array<float>& outSenPix, const IPosition& curpos) {

    float minVal, maxVal;
    IPosition minPos, maxPos;
    minMax(minVal,maxVal,minPos,maxPos,outWgtPix);
    float wgtCutoff = itsCutoff * itsCutoff * maxVal; // outWgtPix is prop. to image (gain/sigma)^2

    // copy the pixel iterator containing all dimensions
    IPosition fullpos(curpos);

    for (int x=0; x<outPix.shape()[0];++x) {
        for (int y=0; y<outPix.shape()[1];++y) {
            fullpos[0] = x;
            fullpos[1] = y;
            if (outWgtPix(fullpos)>=wgtCutoff) {
                outPix(fullpos) = outPix(fullpos) / outWgtPix(fullpos);
            } else {
                outPix(fullpos) = 0.0;
            }
        }
    }

    if (itsDoSensitivity) {
        minMax(minVal,maxVal,minPos,maxPos,outSenPix);
        float varCutoff = itsCutoff * maxVal; // outSenPix is prop. to image (gain/sigma)^2 but squaring is too much
        for (int x=0; x<outPix.shape()[0];++x) {
            for (int y=0; y<outPix.shape()[1];++y) {
                fullpos[0] = x;
                fullpos[1] = y;
                if (outWgtPix(fullpos)>=wgtCutoff && outSenPix(fullpos)>varCutoff) {
                    outSenPix(fullpos) = sqrt(1.0 / outSenPix(fullpos));
                } else {
                    outSenPix(fullpos) = 0.0;
                }
            }
        }
    }

}

Vector<IPosition> LinmosAccumulator::convertImageCornersToRef(const DirectionCoordinate& refDC) {
    // based on SynthesisParamsHelper::facetSlicer, but don't want to load every input image into a scimath::Param

    ASKAPDEBUGASSERT(itsInShape.nelements() >= 2);
    // add more checks

    const int coordPos = itsInCoordSys.findCoordinate(Coordinate::DIRECTION,-1);
    const DirectionCoordinate inDC = itsInCoordSys.directionCoordinate(coordPos);

    IPosition blc(itsInShape.nelements(),0);
    IPosition trc(itsInShape);
    for (uInt dim=0; dim<itsInShape.nelements(); ++dim) {
         --trc(dim); // these are added back later. Is this just to deal with degenerate axes?
    }
    // currently blc,trc describe the whole input image; convert coordinates
    Vector<Double> pix(2);
    
    // first process BLC
    pix[0] = Double(blc[0]);
    pix[1] = Double(blc[1]);
    MDirection tempDir;
    Bool success = inDC.toWorld(tempDir, pix);
    ASKAPCHECK(success, "Pixel to world coordinate conversion failed for input BLC: "<<inDC.errorMessage());
    success = refDC.toPixel(pix,tempDir);
    ASKAPCHECK(success, "World to pixel coordinate conversion failed for output BLC: "<<refDC.errorMessage());
    blc[0] = casa::Int(round(pix[0]));
    blc[1] = casa::Int(round(pix[1]));
 
    // now process TRC
    pix[0] = Double(trc[0]);
    pix[1] = Double(trc[1]);
    success = inDC.toWorld(tempDir, pix);
    ASKAPCHECK(success, "Pixel to world coordinate conversion failed for input TRC: "<<inDC.errorMessage());
    success = refDC.toPixel(pix,tempDir);
    ASKAPCHECK(success, "World to pixel coordinate conversion failed for output TRC: "<<refDC.errorMessage());
    trc[0] = casa::Int(round(pix[0]));
    trc[1] = casa::Int(round(pix[1]));

    Vector<IPosition> corners(2);
    corners[0] = blc;
    corners[1] = trc;

    return corners;

}

bool LinmosAccumulator::coordinatesAreConsistent(const CoordinateSystem& refCoordSys) {
    // Check to see if it makes sense to combine images with these coordinate systems.
    // Could get more tricky, but right now make sure any extra dimensions, such as frequency
    // and polarisation, are equal in the two systems.
    if ( itsInCoordSys.nCoordinates() != refCoordSys.nCoordinates() ) {
        //ASKAPLOG_INFO_STR(logger, "Coordinates are not consistent: shape mismatch");
        return false;
    }
    if (!allEQ(itsInCoordSys.worldAxisNames(), refCoordSys.worldAxisNames())) {
        //ASKAPLOG_INFO_STR(logger, "Coordinates are not consistent: axis name mismatch");
        return false;
    }
    if (!allEQ(itsInCoordSys.worldAxisUnits(), refCoordSys.worldAxisUnits())) {
        //ASKAPLOG_INFO_STR(logger, "Coordinates are not consistent: axis unit mismatch");
        return false;
    }
    return true;
}

bool LinmosAccumulator::coordinatesAreEqual(void) {
    // Check to see if regridding is required. If they are equal there is no need.

    // Test the these things are set up...

    // Does something better already exist?
    double thresh = 1.0e-12;

    // Check that the input dimensionality is the same as that of the output
    if ( !coordinatesAreConsistent(itsOutCoordSys) ) return false;

    // Also check that the size and centre of each dimension is the same.
    if ( itsInShape != itsOutShape ) {
        //ASKAPLOG_INFO_STR(logger, "Input and output coordinates are not equal: shape mismatch");
        return false;
    }
    // test that the grid properties of each dimension are equal
    for (casa::uInt dim=0; dim<itsInCoordSys.nCoordinates(); ++dim) {

        if ( (itsInCoordSys.referencePixel()[dim] != itsOutCoordSys.referencePixel()[dim]) ||
             (fabs(itsInCoordSys.increment()[dim] - itsOutCoordSys.increment()[dim]) > thresh) ||
             (fabs(itsInCoordSys.referenceValue()[dim] - itsOutCoordSys.referenceValue()[dim]) > thresh) ) {
            //ASKAPLOG_INFO_STR(logger, "Coordinates are not equal: coord system mismatch for dim " << dim);
            return false;
        }
    }
    return true;
}

/// @brief do the merge
/// @param[in] parset subset with parameters
static void merge(const LOFAR::ParameterSet &parset) {

    // initialise an image accumulator
    LinmosAccumulator accumulator;

    // load the parset
    if ( !accumulator.loadParset(parset) ) return;

    // initialise an image accessor
    accessors::IImageAccess& iacc = SynthesisParamsHelper::imageHandler();

    // loop over the mosaics, reading each in an adding to the output pixel arrays
    vector<string> inImgNames, inWgtNames, inSenNames;
    string outImgName, outWgtName, outSenName;
    map<string,string> outWgtNames = accumulator.outWgtNames();
    for(map<string,string>::iterator ii=outWgtNames.begin(); ii!=outWgtNames.end(); ++ii) {

        // get output files for this mosaic
        outImgName = (*ii).first;
        outWgtName = accumulator.outWgtNames()[outImgName];
        ASKAPLOG_INFO_STR(logger, "++++++++++++++++++++++++++++++++++++++++++");
        ASKAPLOG_INFO_STR(logger, "Preparing mosaic " << outImgName);
        if (!accumulator.outWgtDuplicates()[outImgName]) {
            ASKAPLOG_INFO_STR(logger, " - also weights image " << outWgtName);
        }
        accumulator.doSensitivity(false);
        if (accumulator.genSensitivityImage()[outImgName]) {
            outSenName = accumulator.outSenNames()[outImgName];
            accumulator.doSensitivity(true);
            ASKAPLOG_INFO_STR(logger, " - also sensitivity image " << outSenName);
        }

        // get input files for this mosaic
        inImgNames = accumulator.inImgNameVecs()[outImgName];
        ASKAPLOG_INFO_STR(logger, " - input images: "<<inImgNames);
        if (accumulator.weightType() == FROM_WEIGHT_IMAGES) {
            inWgtNames = accumulator.inWgtNameVecs()[outImgName];
            ASKAPLOG_INFO_STR(logger, " - input weights images: " << inWgtNames);
        }
        else if (accumulator.weightType() == FROM_BP_MODEL) {
            accumulator.loadBeamCentres(parset,iacc,outImgName);
        }
        if (accumulator.doSensitivity()) {
            inSenNames = accumulator.inSenNameVecs()[outImgName];
            ASKAPLOG_INFO_STR(logger, " - input sensitivity images: " << inSenNames);
        }

        // set the output coordinate system and shape, based on the overlap of input images
        accumulator.setOutputParameters(inImgNames, iacc);

        // set up the output pixel arrays
        Array<float> outPix(accumulator.outShape(),0.);
        Array<float> outWgtPix(accumulator.outShape(),0.);
        Array<float> outSenPix;
        if (accumulator.doSensitivity()) {
            outSenPix = Array<float>(accumulator.outShape(),0.);
        }

        // set up an indexing vector for the arrays
        IPosition curpos(outPix.shape());
        ASKAPASSERT(curpos.nelements()>=2);
        for (uInt dim=0; dim<curpos.nelements(); ++dim) {
            curpos[dim] = 0;
        }

        // loop over the input images, reading each in an adding to the output pixel arrays
        for (uInt img = 0; img < inImgNames.size(); ++img ) {

            // short cuts
            string inImgName = inImgNames[img];
            string inWgtName, inSenName;

            ASKAPLOG_INFO_STR(logger, "Processing input image " << inImgName);
            if (accumulator.weightType() == FROM_WEIGHT_IMAGES) {
                inWgtName = inWgtNames[img];
                ASKAPLOG_INFO_STR(logger, " - and input weight image " << inWgtName);
            }
            if (accumulator.doSensitivity()) {
                inSenName = inSenNames[img];
                ASKAPLOG_INFO_STR(logger, " - and input sensitivity image " << inSenName);
            }

            // set the input coordinate system and shape
            accumulator.setInputParameters(inImgName, iacc, img);

            Array<float> inPix = iacc.read(inImgName);
            Array<float> inWgtPix;
            Array<float> inSenPix;
            if (accumulator.weightType() == FROM_WEIGHT_IMAGES) {
                inWgtPix = iacc.read(inWgtName);
                ASKAPASSERT(inPix.shape() == inWgtPix.shape());
            }
            if (accumulator.doSensitivity()) {
                inSenPix = iacc.read(inSenName);
                ASKAPASSERT(inPix.shape() == inSenPix.shape());
            }

            // set up an iterator for all directionCoordinate planes in the input image
            scimath::MultiDimArrayPlaneIter planeIter(accumulator.inShape());

            // test whether to simply add weighted pixels, or whether a regrid is required
            bool regridRequired = !accumulator.coordinatesAreEqual();

            // if regridding is required, set up buffer some images
            if ( regridRequired ) {

                ASKAPLOG_INFO_STR(logger, " - regridding -- input pixel grid is different from the output");

                // currently all output planes have full-size, so only initialise once
                // would be faster if this was reduced to the size of the current input image
                if ( accumulator.outputBufferSetupRequired() ) {
                    ASKAPLOG_INFO_STR(logger, " - initialising output buffers and the regridder");
                    // set up temp images required for regridding
                    //accumulator.initialiseOutputBuffers();
                    // set up regridder
                    accumulator.initialiseRegridder();
                }

                // set up temp images required for regridding
                // need to do this here if some do and some do not have sensitivity images
                accumulator.initialiseOutputBuffers();

                // set up temp images required for regridding
                // are those of the previous iteration correctly freed?
                accumulator.initialiseInputBuffers();

            } else {
                ASKAPLOG_INFO_STR(logger, " - not regridding -- input pixel grid is the same as the output");
            }

            // iterator over planes (e.g. freq & polarisation), regridding and accumulating weights and weighted images
            for (; planeIter.hasMore(); planeIter.next()) {

                // set the indices of any higher-order dimensions for this slice
                curpos = planeIter.position();

                ASKAPLOG_INFO_STR(logger, " - slice " << curpos);

                if ( regridRequired ) {

                    // load input buffer for the current plane
                    accumulator.loadInputBuffers(planeIter, inPix, inWgtPix, inSenPix);
                    // call regrid for any buffered images
                    accumulator.regrid();
                    // update the accululation arrays for this plane
                    accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, curpos);

                } else {

                    // Update the accululation arrays for this plane.
                    accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, inPix, inWgtPix, inSenPix, curpos);

                }

            }

        } // img loop (over input images)

        // deweight the image pixels
        // use another iterator to loop over planes
        ASKAPLOG_INFO_STR(logger, "Deweighting accumulated images");
        scimath::MultiDimArrayPlaneIter deweightIter(accumulator.outShape());
        for (; deweightIter.hasMore(); deweightIter.next()) {
            curpos = deweightIter.position();
            accumulator.deweightPlane(outPix, outWgtPix, outSenPix, curpos);
        }

        // set one of the input images as a reference for metadata (the first by default)
        uint psfref = 0;
        if (parset.isDefined("psfref")) psfref = parset.getUint("psfref");
        ASKAPLOG_INFO_STR(logger, "Getting PSF beam info for the output image from input number " << psfref);
        // get pixel units from the selected reference image
        Table tmpTable(inImgNames[psfref]);
        string units = tmpTable.keywordSet().asString("units");
        // get psf beam information from the selected reference image
        Vector<Quantum<double> > psf = iacc.beamInfo(inImgNames[psfref]);
        if (psf.nelements()<3) 
            ASKAPLOG_WARN_STR(logger, inImgNames[psfref] << ": beamInfo needs at least 3 elements. Not writing PSF");

        // write accumulated images and weight images
        ASKAPLOG_INFO_STR(logger, "Writing accumulated image to " << outImgName);
        iacc.create(outImgName, accumulator.outShape(), accumulator.outCoordSys());
        iacc.write(outImgName,outPix);
        iacc.setUnits(outImgName,units);
        if (psf.nelements()>=3) 
            iacc.setBeamInfo(outImgName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));

        if (accumulator.outWgtDuplicates()[outImgName]) {
            ASKAPLOG_INFO_STR(logger, "Accumulated weight image " << outWgtName << " already written");
        } else {
            ASKAPLOG_INFO_STR(logger, "Writing accumulated weight image to " << outWgtName);
            iacc.create(outWgtName, accumulator.outShape(), accumulator.outCoordSys());
            iacc.write(outWgtName,outWgtPix);
            iacc.setUnits(outWgtName,units);
            if (psf.nelements()>=3) 
                iacc.setBeamInfo(outWgtName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
        }

        if (accumulator.doSensitivity()) {
            ASKAPLOG_INFO_STR(logger, "Writing accumulated sensitivity image to " << outSenName);
            iacc.create(outSenName, accumulator.outShape(), accumulator.outCoordSys());
            iacc.write(outSenName,outSenPix);
            iacc.setUnits(outSenName,units);
            if (psf.nelements()>=3) 
                iacc.setBeamInfo(outSenName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
        }

    } // ii loop (separate mosaics for different image types)

};

class LinmosApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
	    StatReporter stats;
            LOFAR::ParameterSet subset(config().makeSubset("linmos."));
            SynthesisParamsHelper::setUpImageHandler(subset);
            merge(subset);
	    stats.logSummary();
            return 0;
        }
};

int main(int argc, char *argv[])
{
    LinmosApp app;
    return app.main(argc, argv);
}
