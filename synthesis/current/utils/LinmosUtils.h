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
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/images/Images/ImageRegrid.h>

// Local packages includes
#include <measurementequation/SynthesisParamsHelper.h>
#include <linmos/LinmosAccumulator.h>


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
                                    MVDirection centre);
/// @brief helper method to get beam centres from parset and/or image metadata
/// @details separate from loadParset to allow metadata to be read from input images
/// @param[in] const LOFAR::ParameterSet &parset: linmos parset
/// @param[in] const accessors::IImageAccess &iacc: image accessor
/// @param[in] const string outImgName: current mosaic name
/// @return bool true=success, false=fail
Vector<MVDirection> loadBeamCentres(const LOFAR::ParameterSet &parset,
                                    const accessors::IImageAccess &iacc,
                                    const vector<string> &inImgNames);

#endif


