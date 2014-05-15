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

#ifndef FREQUENCY_MAPPER_H
#define FREQUENCY_MAPPER_H

#include <fitting/Axes.h>
#include <casa/Arrays/Vector.h>

#include <vector>

namespace askap {

namespace synthesis {

/// @brief Mapping between frequency channels and image planes
/// @details This class provides mapping between image (grid) planes and frequency
/// channels. One image plane can correspond to a number of accessor planes 
/// (multi-frequency synthesis). This class is used inside TableVisGridder.
/// @ingroup gridding
struct FrequencyMapper {
   /// @brief default constructor
   /// @details The class is left in an uninitialised state
   FrequencyMapper();
   
   /// @brief constructor doing initialisation
   /// @details
   /// @param[in] axes axes object containing spectral axis of the image cube
   /// @param[in] nchan number of frequency channels in the image cube
   /// @note an exception is thrown if axes object doesn't contain the spectral axis
   FrequencyMapper(const scimath::Axes &axes, int nchan);
   
   /// @brief setup image
   /// @details 
   /// @param[in] axes axes object containing spectral axis of the image cube
   /// @param[in] nchan number of frequency channels in the image cube
   /// @note an exception is thrown if axes object doesn't contain the spectral axis
   void setupImage(const scimath::Axes &axes, int nchan);
   
   /// @brief setup an image where everything is gridded into single plane
   /// @details Current unit tests are written without frequency axis. This method was
   /// added instead of patching all unit tests (and it may be quite useful for debugging as well).
   /// If this method is called, all subsequent calls to operator() would return 0.
   /// Calling setupImage would revert operations back to normal.
   void setupSinglePlaneGridding();
   
   /// @brief setup mapping 
   /// @details 
   /// This method sets up actual mapping between image and accessor channels. Only 
   /// vector returned by accessor's frequency method is required.    
   /// @param[in] freqs vector with frequencies
   /// @note The current assumption is that no regridding is required. Therefore, it is
   /// expected that no fractional channel offset can occur.
   void setupMapping(const casa::Vector<casa::Double> &freqs);
   
   /// @brief test whether the given channel is mapped 
   /// @details The measurement does not necessarily contribute to the cube which is being imaged.
   /// This method allows to check whether some mapping exists. Operator() throws the exception if
   /// it is called for a channel without a mapping.
   /// @param[in] chan accessor channel
   /// @return true, if the given channel has a mapping
   bool isMapped(casa::uInt chan)  const;
   
   /// @brief map accessor channel to image channel
   /// @details
   /// @param[in] chan accessor channel
   /// @note the output is guaranteed to be from [0,itsImageNChan-1] interval.
   casa::uInt operator()(casa::uInt chan) const;
   
private:
   /// @brief start frequency of the image cube
   double itsStartFreq;
   /// @brief end frequency of the image cube
   double itsEndFreq;
   /// @brief number of channels in the image cube
   /// @details A negative value means that the class has not been initialised (-1) or is in a 
   /// single plane image mode (-2).
   int itsImageNChan;
   /// @brief map of accessor channels to image channels
   /// @details the value is negative if no mapping exists for a particular channel
   std::vector<int> itsMap;   
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef FREQUENCY_MAPPER_H

