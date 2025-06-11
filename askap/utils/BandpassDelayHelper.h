/// @file
///
/// Helper class to manage bandpass importing/exporting as well as 
/// delay estimates and some operations to aggregate several solutions.
///
///
/// @copyright (c) 2025 CSIRO
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
///

#ifndef ASKAP_UTILITIES_BANDPASS_DELAY_HELPER_H
#define ASKAP_UTILITIES_BANDPASS_DELAY_HELPER_H

// boost
#include <boost/noncopyable.hpp>

// casa
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/BasicSL/Complex.h>

// own
#include <askap/calibaccess/ICalSolutionConstAccessor.h>
#include <askap/calibaccess/ICalSolutionAccessor.h>
//#include <askap/scimath/utils/DelayEstimator.h>

namespace askap {

namespace utils {

// forward declaration for unit-testing
class BandpassDelayHelperTest;

/// @brief class for delay manipulation in bandpass solution
/// @details This class implements operations on bandpass tables including solving for delays. At this stage it is a proof of concept implementation,
/// we may need to refactor it into several classes later on (read/write operations, actual app, buffer management)
/// @ingroup utils
class BandpassDelayHelper : public boost::noncopyable {
public:
   /// @brief constructor
   /// @param[in] nAnt number of antennas
   /// @param[in] nBeam number of beams
   /// @param[in] nChan number of spectral channels
   /// @note We need these parameters to setup the buffer appropriately. In the current calibration accessor interfaces
   /// there is no way to get this information (because it was intended that the setup won't change). For table-based accessor
   /// which we use, the intention is to set them up from the same parset values as for the output table. In principle, this tool can
   /// be extended in the future to expand the data (e.g. casting to more channels) or to select a subset of the data. But currently it is out of
   /// scope.
   BandpassDelayHelper(casacore::uInt nAnt, casacore::uInt nBeam, casacore::uInt nChan);
    
   /// @brief calculate delays
   /// @details This method estimates delays for all antennas, beams and polarisations for the currently
   /// buffered bandpass solution. Note, it must be read or obtained some other way. For now, the simplest apporach
   /// of phase-slope fitting is used. It is possible to extend it into two stage approach (FFT-spectral averaging-phase slope)
   /// similar to the one used by the delay solver later on. 
   void calcDelays(); 
    
   /// @brief add delays from another object
   /// @details This method adds delays stored in another object of this type (which maps to another calibration table). It is a bit messy to do it
   /// this way but handy to encapsulate all operations with validity flags (this logic can be extracted into a separate class which we can return instead of
   /// having such data as data members of this class).
   /// @param[in] other other instance of the class with calculated delays (which are added to the delays stored in this class)
   void addDelays(const BandpassDelayHelper &other);

   /// @brief flip the sign of stored delays
   /// @details This method essentially multiplies all stored delays by -1. It can be used together with addDelays to implement the subtraction or,
   /// together with calcDelays and applyDelays, to compensate the delays present in the bandpass solution
   void negateDelays();

   /// @brief zero all delays preserving validity flags
   /// @details This method just sets delays to zero, but doesn't touch validity flags. This is handy if one wants to summarise delays 
   /// in the log using the summary() method or add delays from another instance.
   void zeroDelays();

   /// @brief apply currently stored delays to the currently stored bandpass
   /// @details This method modifies currently stored bandpass by applying the phase slope implied by the appropriate delay term. So, combined with
   /// calcDelays and negateDelays calls, this can be used to remove the best fit delay from the loaded bandpass.
   void applyDelays();
    
   /// @brief load bandpass via calibration accessor
   /// @details This method loads bandpass parameters from the calibration solution accessor for the number of
   /// antennas, beams and channels (two polarisations are always expected) this class is initialised with. We
   /// also need spectral resolution (which is not stored with the calibration table as the design assumed a fixed configuration)
   /// to be able to compute (or apply) delays later on in physical units. Working with physical units of delay (rather than per channel) allows
   /// us to combine results obtained with different resolution (e.g. one bandpass table is in full resolution and another is averaged to 1 MHz)
   /// @param[in] acc calibration solution accessor to read bandpass from
   /// @param[in] resolution spectral resolution in Hz
   void loadBandpass(const accessors::ICalSolutionConstAccessor &acc, double resolution);

   /// @brief store bandpass into a calibration accessor
   /// @details This method does the opposite of loadBandpass and stores the current bandpass parameters using the supplied calibration 
   /// solution accessor. Only valid parameters are stored (which ensures that they remain flagged according to the general calibration accessor logic),
   /// but an exception is expected if, say, the number of antennas or beams in the current bandpass with valid data is larger than the given accessor 
   /// can handle. Note, this cannot occur if the output accessor and this class are initialised with the same antennas, beams, etc numbers or if 
   /// the bandpass of the same array has been loaded via loadBandpass first before any manipulation.
   /// @param[in] acc calibration accessor to store the bandpass into
   void storeBandpass(accessors::ICalSolutionAccessor &acc) const;

   /// @brief initialise the class with ideal bandpass
   /// @details This method initialises the class with the amplitude of 1. and zero phase. It doesn't affect the delays if they're calculated earlier.
   /// Thus, this method can be used to export a fake bandpass containing just the delay slope or the inverted one providing the correction term.
   /// This option exists largely for experiments. All antennas, beams, channels and polarisations this class has been setup to handle are marked as valid.
   /// However, the validity flags of delays will be respected if delays are applied to this ideal bandpass.
   /// @param[in] resolution spectral resolution in Hz to be used in subsequent calcDelays and applyDelays calls (as for the bandpass loaded via the accessor)
   void setIdealBandpass(double resolution);


   /// @brief summarise delays in the log
   /// @details This method adds a summary of the current delays in the log. If only one beam has unflagged data, the corresponding delay is printed. 
   /// Otherwise, the range is shown. In the future we can add methods which export delays, e.g. into an ascii file. The output of this method only makes sense if
   /// delays are calculated (e.g. with calcDelays).
   void summary() const;

private:
    /// @brief spectral resolution in Hz of the current bandpass
    /// @details See loadBandpass, this resolution is used by all subsequent calls to calcDelays and applyDelays.
    double itsResolution;
 
    /// @brief bandpass for all antennas, beams, channels and polarisations
    /// @note shape of this array is [nChan x nPol x nBeam x nAnt]. It is initialised in the constructor and used throughout.
    casacore::Array<casacore::Complex> itsBandpass;

    /// @brief bandpass validity flags (true means good)
    /// @note The same shape as itsBandpass
    casacore::Array<bool> itsBandpassValid;

    /// @brief delays for all antennas, beams and polarisations
    /// @note the shape is [nPol x nBeam x nAnt]
    casacore::Cube<float> itsDelay;

    /// @brief validity flags (true means good)
    /// @note The same  shape as itsDelays
    casacore::Cube<bool> itsDelayValid;

    // for unit-testing
    friend BandpassDelayHelperTest;
};

} // namespace utils

} // namespace askap

#endif // #ifndef ASKAP_UTILITIES_BANDPASS_DELAY_HELPER_H
