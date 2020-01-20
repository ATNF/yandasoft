/// @file 
/// @brief SKA_LOW illumination model
/// @details This class represents a SKA_LOW illumination model, 
/// which is the Fourier transform of the SKA_LOW_PB PrimaryBeam model.
/// Optionally a phase slope can be applied to simulate offset pointing.
///
/// @copyright (c) 2008 CSIRO
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

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <askap/AskapError.h>
#include <askap/AskapLogging.h>

#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/gridding/SKA_LOWIllumination.h>
#include <askap/imagemath/primarybeam/PrimaryBeam.h>
#include <askap/imagemath/primarybeam/PrimaryBeamFactory.h>
#include <askap/imagemath/primarybeam/SKA_LOW_PB.h>
#include <askap/scimath/fft/FFTWrapper.h>

ASKAP_LOGGER(logger, ".gridding.ska-lowllumination");

#include <profile/AskapProfiler.h>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::imagemath;

/// @brief construct the model
SKA_LOWIllumination::SKA_LOWIllumination()
{
}

/// @brief  Set pointing parameters
/// @details Set pointing parameters
/// @param[in] az pointing azimuth in degrees
/// @param[in] za pointing zenith angle in degrees
void SKA_LOWIllumination::setPointing(double az, double za, double diam)
{
    ASKAPCHECK((za>=0) && (za<=casacore::C::pi_2), "pointing zenith angle must be in range [0,pi/2] rad, not "<<za);
    ASKAPCHECK(diam>0, "aperture diameter much be greater than zero, not "<<diam);

   	ASKAPLOG_INFO_STR(logger, "Aperture diam = "<<diam<<" m, pointing az,za = "<<az<<", "<<za<<" rad");

    itsAz0 = az;
    itsZa0 = za;
    itsDiameter = diam;
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
void SKA_LOWIllumination::getPattern(double freq, UVPattern &pattern, double l, 
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
void SKA_LOWIllumination::getPattern(double freq, UVPattern &pattern,
                          const casacore::MVDirection &imageCentre,
                          const casacore::MVDirection &beamCentre,
                          const double pa, const bool isPSF) const
{ 

    ASKAPTRACE("SKA_LOWIllumination::getPattern");
    const casacore::uInt oversample = pattern.overSample();
    const double cellU = pattern.uCellSize()/oversample;
    const double cellV = pattern.vCellSize()/oversample;

    // calculate beam offset
    double l = 0.0;
    double m = 0.0;
    if (!isPSF) {
        l = sin(beamCentre.getLong() - imageCentre.getLong()) * cos(beamCentre.getLat());
        m = sin(beamCentre.getLat()) * cos(imageCentre.getLat()) -
            cos(beamCentre.getLat()) * sin(imageCentre.getLat()) *
            cos(beamCentre.getLong() - imageCentre.getLong());
    }

    // scaled l and m to take the calculations out of the loop
    // these quantities are effectively dimensionless 
    const double lScaled = 2.*casacore::C::pi*cellU * l;
    const double mScaled = 2.*casacore::C::pi*cellV * m;
    
    // zero value of the pattern by default
    pattern.pattern().set(0.);
    
    ASKAPCHECK(std::abs(std::abs(cellU/cellV)-1.)<1e-7, 
               "Rectangular cells are not supported at the moment, you have ("<<cellU<<" , "<<cellV<<")");
    
    const double cell = std::abs(cellU*(casacore::C::c/freq));
    
    const double apertureRadiusInCells = itsDiameter/(2.0*cell);  
    
    // square of the station radius
    const double rMaxSquared = casacore::square(apertureRadiusInCells);
	
    // sizes of the grid to fill with pattern values
    const casacore::uInt nU = pattern.uSize();
    const casacore::uInt nV = pattern.vSize();
	
    ASKAPCHECK((casacore::square(double(nU)) > rMaxSquared) &&
               (casacore::square(double(nV)) > rMaxSquared),
               "The pattern buffer passed to SKA_LOWIllumination::getPattern is too small for the given model. "
               "Sizes should be greater than "<<sqrt(rMaxSquared)<<" on each axis, you have "
                <<nU<<" x "<<nV);
	
    // maximum possible support for this class corresponds to the dish size
    pattern.setMaxSupport(1+2*casacore::uInt(apertureRadiusInCells)/oversample);

    // initialise the primary beam
    LOFAR::ParameterSet PBparset;
    PBparset.add("primarybeam","SKA_LOW_PB");
    PBparset.add("primarybeam.SKA_LOW_PB.pointing.az",std::to_string(itsAz0));
    PBparset.add("primarybeam.SKA_LOW_PB.pointing.za",std::to_string(itsZa0));
    PrimaryBeam::ShPtr PB = PrimaryBeamFactory::make(PBparset);

    //std::cout << "DAM PBparset:" << std::endl << PBparset << std::endl;
    std::cout << "DAM image size = " << pattern.uSize() << " / " << pattern.overSample() << std::endl;
    std::cout << "DAM apertureRadiusInCells = " << apertureRadiusInCells << std::endl;
    std::cout << "DAM 1+2*apertureRadiusInCells = " << 1+2*casacore::uInt(apertureRadiusInCells) << std::endl;
    std::cout << "DAM rMaxSquared = " << rMaxSquared << std::endl;
    std::cout << "DAM rMaxSquared x OS = " << rMaxSquared*pattern.overSample() << std::endl;
    std::cout << "DAM freq = " << freq << " Hz" << std::endl;
    std::cout << "DAM lambda = " << casacore::C::c/freq << " metres" << std::endl;
    std::cout << "DAM itsDiameter = " << itsDiameter << " metres" << std::endl;
    std::cout << "DAM itsAz0 = " << itsAz0/casacore::C::degree << " degrees" << std::endl;
    std::cout << "DAM itsZa0 = " << itsZa0/casacore::C::degree << " degrees" << std::endl;

    casacore::Matrix<casacore::Complex> Jones(2,2,0.0);
    double az = 90*casacore::C::degree, za;
    za = 10*casacore::C::degree + 0.0 * casacore::C::c/freq / itsDiameter;
    Jones = PB->getJonesAtOffset(az,za,freq);
    std::cout << "DAM gain at offset of " << (za - 10*casacore::C::degree)/casacore::C::degree << " degrees";
    std::cout << " = " << Jones(0,0) << std::endl;
    za = 10*casacore::C::degree + 0.5 * casacore::C::c/freq / itsDiameter;
    Jones = PB->getJonesAtOffset(az,za,freq);
    std::cout << "DAM gain at offset of " << (za - 10*casacore::C::degree)/casacore::C::degree << " degrees";
    std::cout << " = " << Jones(0,0) << std::endl;
    za = 10*casacore::C::degree + 1.0 * casacore::C::c/freq / itsDiameter;
    Jones = PB->getJonesAtOffset(az,za,freq);
    std::cout << "DAM gain at offset of " << (za - 10*casacore::C::degree)/casacore::C::degree << " degrees";
    std::cout << " = " << Jones(0,0) << std::endl;

    // Find the actual cellsizes in x and y (radians)
    // corresponding to the limited support
    const double ccellx = 1.0 / (double(nU / oversample) * pattern.uCellSize());
    const double ccelly = 1.0 / (double(nV / oversample) * pattern.vCellSize());

    double sum=0.; // normalisation factor
    for (casacore::uInt iy = 0; iy < nV; ++iy) {
        const double y = (double(iy) - double(nV) / 2) * ccelly;
        for (casacore::uInt ix = 0; ix < nU; ++ix) {
            const double x = (double(ix) - double(nU) / 2) * ccellx;
            // should there be a check of whether x, y and xy_dist are valid?
            const double az = casacore::C::pi/2.0 - atan2(y,x);
            const double za = asin(sqrt(x*x+y*y));
            Jones = PB->getJonesAtOffset(az,za,freq);
            // note sure about polarisation. Need to check useage elsewhere (and in AWProjectVisGridder)
            pattern(ix, iy) = Jones(0,0);
            sum += abs(pattern(ix, iy));
        }
    }

    //scimath::fft2d(pattern.pattern(), true);

    std::cout << "DAM pixel sum = " << sum << std::endl;

// DAM normalising to give peak of 1. Probably should be the integral...

    //ASKAPCHECK(sum > 0., "Integral of the aperture should be non-zero");
    //pattern.pattern() *= casacore::DComplex(1.0/float(sum),0.);

for (casacore::uInt iU=0; iU<nU; ++iU) {
    const double offsetU = double(iU)-double(nU)/2.;
    if (abs(offsetU) < 10) std::cout << "DAM " << iU << ", " << offsetU << ": " << pattern(iU, nU/2) << std::endl;
}

}

/// @brief check whether the pattern is symmetric
/// @details Some illumination patterns are known a priori to be symmetric.
/// Need to check this one
/// @return false until more is known
bool SKA_LOWIllumination::isSymmetric() const
{
  return false;
}

/// @brief check whether the output pattern is image-based, rather than an illumination pattern.
/// @details Some illumination patterns need to be generated in the image domain, and given
/// the standard usage (FFT to image-domain for combination with other functions) any image
/// domain function may as well stay in the image domain. So check the state before doing the FFT.
/// @return true 
bool SKA_LOWIllumination::isImageBased() const
{
  return true;
}

