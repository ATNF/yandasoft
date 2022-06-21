#ifndef SYNTHESIS_LINMOSUTILS_H
#define SYNTHESIS_LINMOSUTILS_H

/// @file linmosUtils.h
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
// other 3rd party

#include <Common/ParameterSet.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/coordinates/Coordinates/Coordinate.h>
#include <casacore/casa/Quanta/MVDirection.h>
// Local packages includes
#include <askap/imageaccess/IImageAccess.h>

using namespace casa;

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
                                    MVDirection centre);

/// @brief helper method to get beam centres from parset and/or image metadata
/// @details separate from loadParset to allow metadata to be read from input images
/// @param[in] parset linmos parset
/// @param[in] iacc image accessor
/// @param[in] outImgName current mosaic name
/// @return bool true=success, false=fail
Vector<MVDirection> loadBeamCentres(const LOFAR::ParameterSet &parset,
                                    const accessors::IImageAccess<casacore::Float> &iacc,
                                    const vector<string> &inImgNames);

/// @brief copy selected keywords from the reference image to the output
/// @param[in] string outName : output image name
/// @param[in] string inName : input image name
/// @param[in] vector<string> keywords : list of keyword names to copy
void copyKeywords(const string & outName, const string& inName, const vector<string> & keywords);

/// @brief save a table with the image containing the beamcentres and other information
/// @param[in] string outImgName : output image name, where the table will be saved
/// @param vector<string> inImgNames : input images
/// @param Vector<MVDirection> beamCentres : list of beam centres, must be same number as input images
void saveMosaicTable(const string & outImgName,const vector<string> & inImgNames,
                     const casacore::Vector<casacore::MVDirection> & beamCentres);

/// @brief read weights table or keyword from image header
/// @param[in] string inImgName : image name, where the table is read from
/// @return Vector<Float> : vector with weights, size is 0 if none found,
///    for continuum case (header keyword) size is 1, for spectral case size is nchan
casacore::Vector<casacore::Float> readWeightsTable(const string& inImgName);

} // namespace askap

#endif
