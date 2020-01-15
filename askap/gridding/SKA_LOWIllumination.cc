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


#include <askap/gridding/SKA_LOWIllumination.h>
#include <askap/AskapError.h>
// temporary
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.ska-lowllumination");
//

#include <profile/AskapProfiler.h>

using namespace askap;
using namespace askap::synthesis;

/// @brief construct the model
/// @param[in] diam station diameter in metres
// DAM BLOCKAGE /// @param[in] blockage a diameter of the central hole in metres
// DAM BLOCKAGE SKA_LOWIllumination::SKA_LOWIllumination(double diam, double blockage) :
// DAM BLOCKAGE    itsDiameter(diam), itsBlockage(blockage) 
SKA_LOWIllumination::SKA_LOWIllumination()
{
  // DAM BLOCKAGE ASKAPDEBUGASSERT(diam>0);
  // DAM BLOCKAGE ASKAPDEBUGASSERT(blockage>=0);  
  // DAM BLOCKAGE ASKAPDEBUGASSERT(diam > blockage);
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
                          double m, double) const
{
    ASKAPTRACE("SKA_LOWIllumination::getPattern");
    const casacore::uInt oversample = pattern.overSample();
    const double cellU = pattern.uCellSize()/oversample;
    const double cellV = pattern.vCellSize()/oversample;
    
    // scaled l and m to take the calculations out of the loop
    // these quantities are effectively dimensionless 
    const double lScaled = 2.*casacore::C::pi*cellU *l;
    const double mScaled = 2.*casacore::C::pi*cellV *m;
    
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


    double sum=0.; // normalisation factor
    #ifdef _OPENMP_WORKING_WORKING
    #pragma omp parallel default(shared)
    {
        #pragma omp for reduction(+:sum)
    #endif
        for (casacore::uInt iU=0; iU<nU; ++iU) {
             const double offsetU = double(iU)-double(nU)/2.;
             const double offsetUSquared = casacore::square(offsetU);
             for (casacore::uInt iV=0; iV<nV; ++iV) {
                  const double offsetV = double(iV)-double(nV)/2.;
                  const double offsetVSquared = casacore::square(offsetV);
                  const double radiusSquared = offsetUSquared + offsetVSquared;
                  if (radiusSquared <= rMaxSquared) {
                       // don't need to multiply by wavelength here because we
                       // divided the radius (i.e. the illumination pattern is given
                       // in a relative coordinates in frequency
                       const double phase = lScaled*offsetU + mScaled*offsetV;
                       pattern(iU, iV) = casacore::DComplex(cos(phase), -sin(phase));
                       sum += 1.;
                  }
             }
	}
    #ifdef _OPENMP_WORKING_WORKING
    }
    #endif

    ASKAPCHECK(sum > 0., "Integral of the aperture should be non-zero");
    pattern.pattern() *= casacore::DComplex(float(nU)*float(nV)/float(sum),0.);
}

/// @brief check whether the pattern is symmetric
/// @details Some illumination patterns are known a priori to be symmetric.
/// Need to check this one
/// @return false until more is known
bool SKA_LOWIllumination::isSymmetric() const
{
  return false;
}

