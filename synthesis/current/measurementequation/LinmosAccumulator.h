/// @file LinmosAccumulator.h
///
/// @brief combine a number of images as a linear mosaic
/// @details This is a standalone utility to merge images into
///     a mosaic. Some code/functionality can later be moved into cimager,
///     but for now it is handy to have it separate. 
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

#ifndef ASKAP_LINMOS_H
#define ASKAP_LINMOS_H

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

using namespace casa;

namespace askap {

    namespace synthesis {

        /// @brief Base class supporting linear mosaics (linmos)

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
                void setSingleMosaic(const vector<string> &inImgNames,
                                     const vector<string> &inWgtNames,
                                     const string &outImgName, const string &outWgtName);


                /// @brief set up a single mosaic for each taylor term
                /// @param[in] const vector<string> &inImgNames : vector of images to mosaic
                /// @param[in] const vector<string> &inWgtNames : vector of weight images, if required
                /// @param[in] const string &outImgName : output mosaic image name
                /// @param[in] const string &outWgtName : output mosaic weight image name
                void findAndSetTaylorTerms(const vector<string> &inImgNames,
                                           const vector<string> &inWgtNames,
                                           const string &outImgName,
                                           const string &outWgtName);

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
                void setInputParameters(const string& inImgName,
                                        const accessors::IImageAccess& iacc,
                                        const int n);

                /// @brief set the output coordinate system and shape, based on the overlap of input images
                /// @details This method is based on the SynthesisParamsHelper::add and
                ///     SynthesisParamsHelper::facetSlicer. It has been reimplemented here
                ///     so that images can be read into memory separately.
                /// @param[in] const vector<string>& inImgNames: names of the input images
                ///     (those given for parset key 'names')
                /// @param[in] const accessors::IImageAccess& iac
                void setOutputParameters(const vector<string>& inImgNames,
                                         const accessors::IImageAccess& iacc);

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
                                      Array<float>& inPix,
                                      Array<float>& inWgtPix,
                                      Array<float>& inSenPix);

                /// @brief call the regridder for the buffered plane
                void regrid(void);

                /// @brief add the current plane to the accumulation arrays
                /// @details This method adds from the regridded buffers
                /// @param[out] Array<float>& outPix: accumulated weighted image pixels
                /// @param[out] Array<float>& outWgtPix: accumulated weight pixels
                /// @param[in] const IPosition& curpos: indices of the current plane
                void accumulatePlane(Array<float>& outPix,
                                     Array<float>& outWgtPix,
                                     Array<float>& outSenPix,
                                     const IPosition& curpos);

                /// @brief add the current plane to the accumulation arrays
                /// @details This method adds directly from the input arrays
                /// @param[out] Array<float>& outPix: accumulated weighted image pixels
                /// @param[out] Array<float>& outWgtPix: accumulated weight pixels
                /// @param[in] const Array<float>& inPix: input image pixels
                /// @param[in] const Array<float>& inWgtPix: input weight pixels
                /// @param[in] const IPosition& curpos: indices of the current plane
                void accumulatePlane(Array<float>& outPix,
                                     Array<float>& outWgtPix,
                                     Array<float>& outSenPix,
                                     const Array<float>& inPix,
                                     const Array<float>& inWgtPix,
                                     const Array<float>& inSenPix,
                                     const IPosition& curpos);

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

    } // namespace synthesis

} // namespace askap

#include <measurementequation/LinmosAccumulator.tcc>

#endif

