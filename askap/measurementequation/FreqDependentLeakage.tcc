/// @file
///
/// @brief Calibration effect: frequency-dependent antenna gains without cross-pol
/// @details This effect represents bandpass calibration parameters. Note, an
/// exception is thrown if the number of spectral channels in the accessor doesn't match
/// the dimension of the bandpass parameter
///
/// @copyright (c) 2007 CSIRO
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

#ifndef FREQ_DEPENDENT_LEAKAGE_TCC
#define FREQ_DEPENDENT_LEAKAGE_TCC

// own includes
#include <askap/scimath/utils/PolConverter.h>
#include <askap/calibaccess/CalParamNameHelper.h>

// casa includes
#include <casacore/measures/Measures/Stokes.h>

// std includes
#include <string>
#include <utility>

namespace askap {

namespace synthesis {


/// @brief main method returning Mueller matrix and derivatives
/// @details This method has to be overloaded (in the template sense) for
/// all classes representing various calibration effects. CalibrationME
/// template will call it when necessary. It returns
/// @param[in] chunk accessor to work with
/// @param[in] row row of the chunk to work with
/// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
/// this effect
inline scimath::ComplexDiffMatrix FreqDependentLeakage::get(const accessors::IConstDataAccessor &chunk,
                                      casacore::uInt row) const
{
   const casacore::uInt chanOffset = static_cast<casacore::uInt>(parameters()->has("chan_offset") ? parameters()->scalarValue("chan_offset") : 0);
   const casacore::uInt nPol = chunk.nPol();
   ASKAPDEBUGASSERT(nPol != 0);
   const casacore::uInt nChan = chunk.nChannel();
   const casacore::Vector<casacore::Stokes::StokesTypes> stokes = chunk.stokes();
   ASKAPDEBUGASSERT(stokes.nelements() == nPol);
   ASKAPDEBUGASSERT(!scimath::PolConverter::isStokes(stokes));

   const casacore::uInt ant1 = chunk.antenna1()[row];
   const casacore::uInt ant2 = chunk.antenna2()[row];

   const casacore::uInt beam1 = chunk.feed1()[row];
   const casacore::uInt beam2 = chunk.feed2()[row];

   // flag showing that the polarisation products are present
   // in the canonic form (e.g. XX,XY,YX,YY for linears)
   bool canonicPolOrder = (nPol == 4);
   ASKAPCHECK(canonicPolOrder, "Need 4 polarisations for leakage calibration");
   for (casacore::uInt pol=0; pol<nPol; ++pol) {

        if (pol<nPol) {
            const casacore::uInt polIndex = scimath::PolConverter::getIndex(stokes[pol]);
            ASKAPDEBUGASSERT(polIndex<4);
            // polIndex is index in the polarisation frame, i.e.
            // XX is 0, XY is 1, YX is 2 and YY is 3
            // we need an index into matrix
            if (polIndex != pol) {
                canonicPolOrder = false;
            }
        } else {
              canonicPolOrder = false;
        }
    }
    ASKAPCHECK(canonicPolOrder, "Only canonic order of polarisation products (e.g. XX,XY,YX,YY) is currently supported");

   // main diagonal is always 1.
   scimath::ComplexDiffMatrix calFactor(nPol, nPol * nChan, 1.);

   for (casacore::uInt pol=0; pol<nPol; ++pol) {

        // leakage terms for antenna 1, polarisation X if XX or XY, or Y if YX or YY
        const std::string d1yx = accessors::CalParamNameHelper::paramName(ant1, beam1, casacore::Stokes::YX, true);
        const std::string d1xy = accessors::CalParamNameHelper::paramName(ant1, beam1, casacore::Stokes::XY, true);

        // gains for antenna 2, polarisation X if XX or YX, or Y if XY or YY
        const std::string d2yx = accessors::CalParamNameHelper::paramName(ant2, beam2, casacore::Stokes::YX, true);
        const std::string d2xy = accessors::CalParamNameHelper::paramName(ant2, beam2, casacore::Stokes::XY, true);

        for (casacore::uInt chan = 0; chan < nChan; ++chan) {
             // we need to think of how to deal with distributed problem on the cluster (i.e. adding
             // some base to the channel number and propagating it through the framework)
             casacore::uInt cOff = chan * nPol;
             calFactor(3, 1 + cOff) = -1.*getParameter(accessors::CalParamNameHelper::addChannelInfo(d1yx,chan+chanOffset));
             calFactor(1, 3 + cOff) =     getParameter(accessors::CalParamNameHelper::addChannelInfo(d1xy,chan+chanOffset));
             calFactor(3, 2 + cOff) = -1.*conj(getParameter(accessors::CalParamNameHelper::addChannelInfo(d2yx,chan+chanOffset)));
             calFactor(2, 3 + cOff) =     conj(getParameter(accessors::CalParamNameHelper::addChannelInfo(d2xy,chan+chanOffset)));
             calFactor(0, 1 + cOff) = calFactor(2, 3 + cOff);
             calFactor(1, 0 + cOff) = calFactor(3, 2 + cOff);
             calFactor(0, 2 + cOff) = calFactor(1, 3 + cOff);
             calFactor(2, 0 + cOff) = calFactor(3, 1 + cOff);

             // cross-diagonal
             calFactor(pol, (3 - pol) + cOff) = (pol % 3 == 0 ? 1. : -1.) * getParameter(
                accessors::CalParamNameHelper::addChannelInfo(pol < 2 ? d1xy : d1yx, chan+chanOffset)) *
                conj(getParameter(accessors::CalParamNameHelper::addChannelInfo(pol % 2 == 0 ? d2xy : d2yx, chan+chanOffset)));
        }
   }
   return calFactor;
}

} // namespace synthesis

} // namespace askap


#endif // #ifndef FREQ_DEPENDENT_LEAKAGE_TCC
