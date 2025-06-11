/// @file ASKAPIllumination.cc
/// @brief ASKAP illumination model
/// @details This class represents a ASKAP illumination model,
/// represented in the image domain via the ASKAP_PB PrimaryBeam model.
///
/// @copyright (c) 2021 CSIRO
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

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <askap/askap/AskapError.h>
#include <askap/askap/AskapLogging.h>

#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/gridding/ASKAPIllumination.h>
#include <askap/imagemath/primarybeam/PrimaryBeam.h>
#include <askap/imagemath/primarybeam/PrimaryBeamFactory.h>
#include <askap/imagemath/primarybeam/ASKAP_PB.h>

ASKAP_LOGGER(logger, ".gridding.askapillumination");

#include <askap/profile/AskapProfiler.h>

// for debugging - to export intermediate images
#include <askap/scimath/utils/ImageUtils.h>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::imagemath;

/// @brief construct the model
ASKAPIllumination::ASKAPIllumination(const LOFAR::ParameterSet &parset):
itsParset(parset)
{
    ASKAPCHECK(itsParset.getString("primarybeam","ASKAP_PB") == "ASKAP_PB",
               "ASKAPIllumination needs a primarybeam.ASKAP_PB specification");
}

/// @brief obtain illumination pattern
/// @details This is the main method which populates the
/// supplied uv-pattern with the values corresponding to the model
/// represented by this object. It has to be overridden in the
/// derived classes. An optional phase slope can be applied to
/// simulate offset pointing.
/// @param[in] freq frequency in Hz for which an illumination pattern is required
/// @param[in] pattern a UVPattern object to fill
/// @param[in] l angular offset in the u-direction (in radians)
/// @param[in] m angular offset in the v-direction (in radians)
/// @param[in] pa parallactic angle, or strictly speaking the angle between
/// uv-coordinate system and the system where the pattern is defined (unused)
void ASKAPIllumination::getPattern(double freq, UVPattern &pattern, double l,
                          double m, double pa) const
{
    // zero value of the pattern by default
    pattern.pattern().set(0.);
}

/// @param[in] freq frequency in Hz for which an illumination pattern is required
/// @param[in] pattern a UVPattern object to fill
/// @param[in] imageCentre ra & dec of the image centre
/// @param[in] beamCentre ra & dec of the beam pointing centre
/// @param[in] pa polarisation position angle (in radians)
/// @param[in] isPSF bool indicting if this gridder is for a PSF
/// @param[in] feed  feed number for case where pattern differs between feeds
void ASKAPIllumination::getPattern(double freq, UVPattern &pattern,
                          const casacore::MVDirection &imageCentre,
                          const casacore::MVDirection &beamCentre,
                          const double pa, const bool isPSF, const int feed) const
{

    ASKAPTRACE("ASKAPIllumination::getPattern");

    ASKAPCHECK(std::abs(std::abs(pattern.uCellSize()/pattern.vCellSize())-1.)<1e-7,
               "Rectangular cells are not supported, you have "<<pattern.uCellSize()<<","<<pattern.vCellSize());

    // set image and beam centres
    const double ra0 = imageCentre.getLong();
    const double dec0 = imageCentre.getLat();
    double raB, decB;
    if (isPSF) {
   	    ASKAPLOG_INFO_STR(logger, "PB: using centred beam pointing (for PSF):");
        raB = ra0;
        decB = dec0;
    } else {
   	    ASKAPLOG_INFO_STR(logger, "PB: using FEED table for beam pointing:");
        raB = beamCentre.getLong();
        decB = beamCentre.getLat();
    }
    ASKAPLOG_INFO_STR(logger, "PB: image - RA = " << raB*12./casacore::C::pi << " hours"<<
        " Dec = " << decB/casacore::C::degree << " degrees");
    ASKAPLOG_INFO_STR(logger, "PB: beam  - RA = " << raB*12./casacore::C::pi << " hours"<<
        " Dec = " << decB/casacore::C::degree << " degrees");

    // initialise the primary beam
    LOFAR::ParameterSet PBparset=itsParset.makeSubset("");
    PBparset.replace("primarybeam","ASKAP_PB");
    PBparset.replace("primarybeam.ASKAP_PB.beam_number",std::to_string(feed));
    // Do we need to add the pa to the specified alpha?
    const double alpha = itsParset.getDouble("primarybeam.ASKAP_PB.alpha",0.0);
    PBparset.replace("primarybeam.ASKAP_PB.alpha",std::to_string(alpha+pa));
    PBparset.replace("primarybeam.ASKAP_PB.pointing.ra",std::to_string(raB));
    PBparset.replace("primarybeam.ASKAP_PB.pointing.dec",std::to_string(decB));

    PrimaryBeam::ShPtr PB = PrimaryBeamFactory::make(PBparset);

    // sizes of the grid to fill with pattern values
    const casacore::uInt nU = pattern.uSize();
    const casacore::uInt nV = pattern.vSize();

    // Find the actual cellsizes in x and y (radians)
    // corresponding to the limited support
    const double ccellx = 1.0 / (double(nU / pattern.overSample()) * pattern.uCellSize());
    const double ccelly = 1.0 / (double(nV / pattern.overSample()) * pattern.vCellSize());

    // zero pad pixels outside 1/(n * cellsize) to avoid aliasing? This isn't right...
    // const double cellmax = double(nU)/2 / (double(nU) * pattern.uCellSize());

    // zero value of the pattern by default
    pattern.pattern().set(0.);

    double ra, dec;
    //double sum=0.; // normalisation factor
    casacore::Matrix<casacore::Complex> Jones(2,2,0.0);
    for (casacore::uInt iy = 0; iy < nV; ++iy) {
        const double m = (double(iy) - double(nV) / 2) * ccelly;
        for (casacore::uInt ix = 0; ix < nU; ++ix) {
            const double l = (double(ix) - double(nU) / 2) * ccellx;
            // should there be a check of whether x, y and xy_dist are valid?
            const double n = sqrt(1. - l*l - m*m);
            ra  = ra0 + atan2( l, n*cos(dec0) - m*sin(dec0) );
            dec = asin( m*cos(dec0) + n*sin(dec0) );
            Jones = PB->getJones(ra,dec,freq);
            // not sure about polarisation. Need to check useage elsewhere (and in AWProjectVisGridder)
            pattern(ix, iy) = Jones(0,0);
            //sum += abs(pattern(ix, iy));

        }
    }

    //ASKAPLOG_INFO_STR(logger,"Integral of ASKAP illumination pattern ="<<sum);
    //ASKAPCHECK(sum > 0., "Integral of the aperture should be non-zero");
    //pattern.pattern() *= imtypeComplex(1.0/float(sum),0.);

}

/// @brief check whether the pattern is symmetric
/// @details Some illumination patterns are known a priori to be symmetric.
/// @return false
bool ASKAPIllumination::isSymmetric() const
{
  return false;
}

/// @brief check whether the output pattern is image-based, rather than an illumination pattern.
/// @details Some illumination patterns need to be generated in the image domain, and given
/// the standard usage (FFT to image-domain for combination with other functions) any image
/// domain function may as well stay in the image domain. So check the state before doing the FFT.
/// @return true
bool ASKAPIllumination::isImageBased() const
{
  return true;
}

/// @brief check whether the output pattern is feed dependent
/// @details Some illumination patterns vary with feed (number) and no shortcuts can
/// be taken
/// @return false
bool ASKAPIllumination::isFeedDependent() const
{
    return true;
}
