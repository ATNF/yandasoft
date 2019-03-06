/// @file 
/// @brief Mapping between frequency channels and image planes
/// @details This class provides mapping between image (grid) planes and frequency
/// channels. One image plane can correspond to a number of accessor planes 
/// (multi-frequency synthesis). This class is used inside TableVisGridder.
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

#include <gridding/FrequencyMapper.h>
#include <askap/AskapError.h>

using namespace askap;
using namespace synthesis;

/// @brief default constructor
/// @details The class is left in an uninitialised state
FrequencyMapper::FrequencyMapper() : itsImageNChan(-1) {}
   
/// @brief constructor doing initialisation
/// @details
/// @param[in] axes axes object containing spectral axis of the image cube
/// @param[in] nchan number of frequency channels in the image cube
/// @note an exception is thrown if axes object doesn't contain the spectral axis
FrequencyMapper::FrequencyMapper(const scimath::Axes &axes, int nchan) : itsImageNChan(nchan) 
{
  setupImage(axes, nchan);
}
   
/// @brief setup image
/// @details 
/// @param[in] axes axes object containing spectral axis of the image cube
/// @param[in] nchan number of frequency channels in the image cube
/// @note an exception is thrown if axes object doesn't contain the spectral axis
void FrequencyMapper::setupImage(const scimath::Axes &axes, int nchan)
{
   ASKAPCHECK(axes.has("FREQUENCY"), "FREQUENCY axis is missing in axes object passed to FrequencyMapper:setupImage");
   ASKAPASSERT(nchan>0);
   itsStartFreq = axes.start("FREQUENCY");
   itsEndFreq = axes.end("FREQUENCY");
   itsImageNChan = nchan;
}
   
/// @brief setup mapping 
/// @details 
/// This method sets up actual mapping between image and accessor channels. Only 
/// vector returned by accessor's frequency method is required.    
/// @param[in] freqs vector with frequencies
/// @note The current assumption is that no regridding is required. Therefore, it is
/// expected that no fractional channel offset can occur.
void FrequencyMapper::setupMapping(const casa::Vector<casa::Double> &freqs)
{
   if (itsImageNChan == -2) {
       // special case of the single plane mapper, do nothing
       return;
   }
   ASKAPCHECK(itsImageNChan>0, "An attempt to call setupMapping for uninitialised FrequencyMapper");
   itsMap.resize(freqs.nelements());
   if (itsImageNChan == 1) {
       // special case (tbd - check that assumed definition is what was intended)
       // ignore frequency axis and map everything to a single plane (same behavior as we used to have)
       for (vector<int>::iterator it = itsMap.begin(); it!=itsMap.end(); ++it) {
            *it = 0;
       }            
   } else {
      const double increment = (itsEndFreq - itsStartFreq)/double(itsImageNChan-1);
      ASKAPCHECK(fabs(increment)>0, "Frequency axis in the image has the same start and end frequency "<<itsStartFreq);
      for (casa::uInt chan=0; chan<freqs.nelements(); ++chan) {
           double imgChanDouble = (freqs[chan]-itsStartFreq)/increment;
           imgChanDouble += (imgChanDouble < 0.) ? -0.5 : 0.5;
           const int imgChan = int(imgChanDouble);
           if (imgChan<0 || imgChan>=itsImageNChan) {
               itsMap[chan] = -1;
           } else {
               itsMap[chan] = imgChan;
           }
           /*
           std::cout<<"chan="<<chan<<" freqs[chan]="<<freqs[chan]<<" imgChan="<<imgChan<<
                      " imgChanDouble="<<imgChanDouble<<" mapped to "<<itsMap[chan]<<std::endl;
           */
      }
   }
}

/// @brief test whether the given channel is mapped 
/// @details The measurement does not necessarily contribute to the cube which is being imaged.
/// This method allows to check whether some mapping exists. Operator() throws the exception if
/// it is called for a channel without a mapping.
/// @param[in] chan accessor channel
/// @return true, if the given channel has a mapping
bool FrequencyMapper::isMapped(casa::uInt chan)  const
{
  if (itsImageNChan == -2) {
      return true;
  }
  ASKAPDEBUGASSERT(chan<itsMap.size());
  return itsMap[chan]>=0;
}

/// @brief setup an image where everything is gridded into single plane
/// @details Current unit tests are written without frequency axis. This method was
/// added instead of patching all unit tests (and it may be quite useful for debugging as well).
/// If this method is called, all subsequent calls to operator() would return 0.
/// Calling setupImage would revert operations back to normal.
void FrequencyMapper::setupSinglePlaneGridding()
{
  itsImageNChan = -2;
}

   
/// @brief map accessor channel to image channel
/// @details
/// @param[in] chan accessor channel
/// @note the output is guaranteed to be from [0,itsImageNChan-1] interval.
casa::uInt FrequencyMapper::operator()(casa::uInt chan) const
{
  if (itsImageNChan == -2) {
      return 0;
  }
  ASKAPDEBUGASSERT(chan<itsMap.size());
  ASKAPDEBUGASSERT(itsImageNChan>0);
  ASKAPCHECK(itsMap[chan]>=0, "An attempt to call FrequencyMapper::operator() for unmapped channel "<<chan);
  return casa::uInt(itsMap[chan]);
}

