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

#ifndef NO_X_POL_FREQ_DEPENDENT_GAIN_TCC
#define NO_X_POL_FREQ_DEPENDENT_GAIN_TCC

// own includes
#include <utils/PolConverter.h>
#include <calibaccess/CalParamNameHelper.h>

// casa includes
#include <measures/Measures/Stokes.h>

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
inline scimath::ComplexDiffMatrix NoXPolFreqDependentGain::get(const accessors::IConstDataAccessor &chunk, 
                                      casa::uInt row) const
{
   const casa::uInt chanOffset = static_cast<casa::uInt>(parameters()->has("chan_offset") ? parameters()->scalarValue("chan_offset") : 0);
   const casa::uInt nPol = chunk.nPol();
   ASKAPDEBUGASSERT(nPol != 0);   
   const casa::uInt nChan = chunk.nChannel();
   const casa::Vector<casa::Stokes::StokesTypes> stokes = chunk.stokes();   
   ASKAPDEBUGASSERT(stokes.nelements() == nPol);
   ASKAPDEBUGASSERT(!scimath::PolConverter::isStokes(stokes));
   
   const casa::uInt ant1 = chunk.antenna1()[row];
   const casa::uInt ant2 = chunk.antenna2()[row];
   
   const casa::uInt beam1 = chunk.feed1()[row];
   const casa::uInt beam2 = chunk.feed2()[row];
   
   
   scimath::ComplexDiffMatrix calFactor(nPol, nPol * nChan, 0.);

   for (casa::uInt pol=0; pol<nPol; ++pol) {
        
        const casa::uInt polIndex = scimath::PolConverter::getIndex(stokes[pol]);
        // polIndex is index in the polarisation frame, i.e.
        // XX is 0, XY is 1, YX is 2 and YY is 3
        // we need an index into matrix 

        // gains for antenna 1, polarisation X if XX or XY, or Y if YX or YY
        const std::string g1name = accessors::CalParamNameHelper::paramName(ant1, beam1, 
                      polIndex / 2 == 0 ? casa::Stokes::XX : casa::Stokes::YY, true);
            
        // gains for antenna 2, polarisation X if XX or YX, or Y if XY or YY
        const std::string g2name = accessors::CalParamNameHelper::paramName(ant2, beam2, 
                      polIndex % 2 == 0 ? casa::Stokes::XX : casa::Stokes::YY, true);
                                                      
        for (casa::uInt chan = 0; chan < nChan; ++chan) {
             // we need to think of how to deal with distributed problem on the cluster (i.e. adding
             // some base to the channel number and propagating it through the framework)    
             calFactor(pol, pol + chan * nPol) = getParameter(accessors::CalParamNameHelper::addChannelInfo(g1name,chan+chanOffset)) *
                        conj(getParameter(accessors::CalParamNameHelper::addChannelInfo(g2name,chan+chanOffset)));
        }            
   }
   return calFactor;
}

} // namespace synthesis

} // namespace askap


#endif // #ifndef NO_X_POL_FREQ_DEPENDENT_GAIN_TCC
