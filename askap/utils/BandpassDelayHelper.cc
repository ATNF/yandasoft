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

// casa
#include <casacore/casa/Arrays/ArrayLogical.h>

// own
#include <askap/utils/BandpassDelayHelper.h>
#include <askap/scimath/utils/DelayEstimator.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>

// std
#include <set>
#include <sstream>
#include <iomanip>

ASKAP_LOGGER(logger, ".BandpassDelayHelper");

namespace askap {

namespace utils {

/// @brief constructor
/// @param[in] nAnt number of antennas
/// @param[in] nBeam number of beams
/// @param[in] nChan number of spectral channels
/// @note We need these parameters to setup the buffer appropriately. In the current calibration accessor interfaces
/// there is no way to get this information (because it was intended that the setup won't change). For table-based accessor
/// which we use, the intention is to set them up from the same parset values as for the output table. In principle, this tool can
/// be extended in the future to expand the data (e.g. casting to more channels) or to select a subset of the data. But currently it is out of
/// scope.
BandpassDelayHelper::BandpassDelayHelper(casacore::uInt nAnt, casacore::uInt nBeam, casacore::uInt nChan) :
     itsBandpass(casacore::IPosition(4, nChan, 2, nBeam, nAnt)), itsBandpassValid(casacore::IPosition(4, nChan, 2, nBeam, nAnt), false),
     itsDelay(2, nBeam, nAnt), itsDelayValid(2, nBeam, nAnt, false) {}
    
/// @brief calculate delays
/// @details This method estimates delays for all antennas, beams and polarisations for the currently
/// buffered bandpass solution. Note, it must be read or obtained some other way. For now, the simplest apporach
/// of phase-slope fitting is used. It is possible to extend it into two stage approach (FFT-spectral averaging-phase slope)
/// similar to the one used by the delay solver later on. 
void BandpassDelayHelper::calcDelays() 
{
   scimath::DelayEstimator de(itsResolution);
   const casacore::IPosition shape = itsBandpass.shape();
   ASKAPDEBUGASSERT(shape.nelements() == 4);
   // note, implicitly cast to unsigned types as they're accepted by the interface
   const casacore::uInt nChan = shape[0];
   const casacore::uInt nBeam = shape[2];
   const casacore::uInt nAnt = shape[3];
   ASKAPDEBUGASSERT(shape[1] == 2);
   for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
        for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
             for (casacore::uInt pol = 0; pol < 2u; ++pol) {
                  // taking the slice of bandpass and the corresponding array of validity flags for the given antenna and beam
                  // due to the selected order of axes, the resulting 2D array should be contiguous
                  const casacore::IPosition start(4, 0, pol, beam, ant);
                  const casacore::IPosition length(4, nChan, 1, 1, 1);
                  const casacore::Vector<casacore::Complex> bpVec = itsBandpass(casacore::Slicer(start, length)).nonDegenerate(1u);
                  const casacore::Vector<bool> bpValidVec = itsBandpassValid(casacore::Slicer(start, length)).nonDegenerate(1u);
                  ASKAPDEBUGASSERT(bpVec.nelements() == bpValidVec.nelements());
                  ASKAPDEBUGASSERT(bpVec.nelements() == nChan);
                  // take a copy into the working buffer to be able to modify flagged channels (there is no good way to deal with flags here,
                  // although some improvements could be done like interpolation across flagged channels or better delay fitting algorithms)
                  casacore::Vector<casacore::Complex> buf(bpVec.copy());
                  casacore::uInt numValidChan = 0u;
                  for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                       if (bpValidVec[chan]) {
                           buf[chan] = bpVec[chan];
                           ++numValidChan;
                       } else {
                           // some more clever interpolation can be done here, but leave it up to bandpass solver / preprocessor
                           // essentially the same approach delay solver is using (and we can add two stage solution and averaging later on)
                           buf[chan] = chan == 0u ? 0.f : buf[chan - 1];
                       }
                  }
                  if (numValidChan > 1u) {
                      itsDelayValid(pol, beam, ant) = true;
                      // can use more advanced methods (e.g. two-stage approach) here
                      itsDelay(pol, beam, ant) = de.getDelay(buf);
                  } else {
                      itsDelayValid(pol, beam, ant) = false;
                  }
             }
        }
   }
}
    
/// @brief add delays from another object
/// @details This method adds delays stored in another object of this type (which maps to another calibration table). It is a bit messy to do it
/// this way but handy to encapsulate all operations with validity flags (this logic can be extracted into a separate class which we can return instead of
/// having such data as data members of this class).
/// @param[in] other other instance of the class with calculated delays (which are added to the delays stored in this class)
void BandpassDelayHelper::addDelays(const BandpassDelayHelper &other)
{
   const casacore::IPosition shape = itsDelay.shape();
   const casacore::IPosition otherShape = other.itsDelay.shape();
   ASKAPDEBUGASSERT(shape.nelements() == 3u);
   ASKAPDEBUGASSERT(otherShape.nelements() == 3u);
   ASKAPDEBUGASSERT(shape[0] == otherShape[0]);
   ASKAPDEBUGASSERT(shape[0] == 2);
   if (otherShape > shape) {
       ASKAPLOG_INFO_STR(logger, "adding delays computed for a bandpass table with larger number of antennas or beams - some data are going to be ignored");
   }
   if (otherShape < shape) {
       ASKAPLOG_INFO_STR(logger, "adding delays computed for a bandpass table with less antennas or beams - some good data may be flagged as invalid");
   }
   for (int ant = 0; ant < shape[2]; ++ant) {
        for (int beam = 0; beam < shape[1]; ++beam) {
             const bool otherIsOutOfBounds = (ant >= otherShape[2]) || (beam >= otherShape[1]);
             for (int pol = 0; pol < shape[0]; ++pol) {
                  itsDelayValid(pol, beam, ant) &= !otherIsOutOfBounds && other.itsDelayValid(pol, beam, ant);
                  if (itsDelayValid(pol, beam, ant)) {
                      itsDelay(pol, beam, ant) += other.itsDelay(pol, beam, ant);
                  }
             }
        }
   }
}

/// @brief flip the sign of stored delays
/// @details This method essentially multiplies all stored delays by -1. It can be used together with addDelays to implement the subtraction or,
/// together with calcDelays and applyDelays, to compensate the delays present in the bandpass solution
void BandpassDelayHelper::negateDelays()
{
   // just flip the sign for everything, validity flags remain the same
   itsDelay *= -1.f;
}

/// @brief zero all delays preserving validity flags
/// @details This method just sets delays to zero, but doesn't touch validity flags. This is handy if one wants to summarise delays 
/// in the log using the summary() method or add delays from another instance.
void BandpassDelayHelper::zeroDelays()
{
   itsDelay.set(0.f);
}

/// @brief apply currently stored delays to the currently stored bandpass
/// @details This method modifies currently stored bandpass by applying the phase slope implied by the appropriate delay term. So, combined with
/// calcDelays and negateDelays calls, this can be used to remove the best fit delay from the loaded bandpass.
void BandpassDelayHelper::applyDelays()
{
   const casacore::IPosition shape = itsDelay.shape();
   ASKAPDEBUGASSERT(itsBandpass.shape().nelements() == shape.nelements() + 1);
   // keep the type signed because IPosition class has signed coordinates
   const int nChan = itsBandpass.shape()[0];
   ASKAPDEBUGASSERT(shape == itsBandpass.shape().getLast(shape.nelements()));

   for (int ant = 0; ant < shape[2]; ++ant) {
        for (int beam = 0; beam < shape[1]; ++beam) {
             for (int pol = 0; pol < shape[0]; ++pol) {
                  const bool delayValid = itsDelayValid(pol, beam, ant);
                  // taking the slice to get bandpass for the given pol, beam, ant and benefit from reference semantics
                  // (note, due to the order of axes the resulting vector should be contiguous in memory and hence, it should be efficient)
                  const casacore::IPosition start(4, 0, pol, beam, ant);
                  const casacore::IPosition length(4, nChan, 1, 1, 1);
                  casacore::Vector<casacore::Complex> bpVec = itsBandpass(casacore::Slicer(start, length)).nonDegenerate(1u);
                  casacore::Vector<bool> bpValidVec = itsBandpassValid(casacore::Slicer(start, length)).nonDegenerate(1u);
                  if (delayValid) {
                      // apply the slope, keep the original bandpass validity flags
                      // note - check the sign
                      const float phaseStepPerChannel = 2.*casacore::C::pi * itsResolution * itsDelay(pol, beam, ant);
                      for (int chan = 0; chan < nChan; ++chan) {
                           bpVec[chan] *= casacore::polar(1.f, phaseStepPerChannel * chan);
                      }
                  } else {
                      // invalidate all bandpass points as the corresponding delay term is invalid
                      bpValidVec.set(false);
                  }
             }
        }
   }  
}
    
/// @brief load bandpass via calibration accessor
/// @details This method loads bandpass parameters from the calibration solution accessor for the number of
/// antennas, beams and channels (two polarisations are always expected) this class is initialised with. We
/// also need spectral resolution (which is not stored with the calibration table as the design assumed a fixed configuration)
/// to be able to compute (or apply) delays later on in physical units. Working with physical units of delay (rather than per channel) allows
/// us to combine results obtained with different resolution (e.g. one bandpass table is in full resolution and another is averaged to 1 MHz)
/// @param[in] acc calibration solution accessor to read bandpass from
/// @param[in] resolution spectral resolution in Hz
void BandpassDelayHelper::loadBandpass(const accessors::ICalSolutionConstAccessor &acc, double resolution)
{
   itsResolution = resolution;
   const casacore::IPosition shape = itsBandpass.shape();
   ASKAPDEBUGASSERT(shape.nelements() == 4);
   // note, implicitly cast to unsigned types as they're accepted by the interface
   const casacore::uInt nChan = shape[0];
   const casacore::uInt nBeam = shape[2];
   const casacore::uInt nAnt = shape[3];
   ASKAPDEBUGASSERT(shape[1] == 2);
   for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
        for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
             // taking the slice of bandpass and the corresponding array of validity flags for the given antenna and beam
             // due to the selected order of axes, the resulting 2D array should be contiguous
             const casacore::IPosition start(4, 0, 0, beam, ant);
             const casacore::IPosition length(4, nChan, 2, 1, 1);
             casacore::Matrix<casacore::Complex> bpMtr = itsBandpass(casacore::Slicer(start, length)).nonDegenerate(2u);
             casacore::Matrix<bool> bpValidMtr = itsBandpassValid(casacore::Slicer(start, length)).nonDegenerate(2u);
             ASKAPDEBUGASSERT(bpMtr.shape() == bpValidMtr.shape());
             ASKAPDEBUGASSERT(bpMtr.nrow() == nChan);
             ASKAPDEBUGASSERT(bpMtr.ncolumn() == 2u);

             const askap::scimath::JonesIndex index(ant, beam);
             for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                  // MV: I don't think it is a great design choice by throwing an exception if the requested element
                  // is out of range of the stored table, especially because this behaviour is not really documented.
                  // Moreover, it looks like this is a side-effect of reusing the memory buffer class for the calibration accessor
                  // (where doing checks with ASKAPCHECK seems appropriate). Perhaps, we need to add a logic to the table-based
                  // calibration accessor to return invalid gain, etc if antenna or beam indices are wrong. Having said that, I am not
                  // sure it is hugely important from the practical side of things as we don't normally work with different number
                  // of channels, beams or antennas to cause the problem.
                  try {
                     const askap::accessors::JonesJTerm jTerm = acc.bandpass(index, chan); 
                     bpMtr(chan,0) = jTerm.g1();
                     bpValidMtr(chan,0) = jTerm.g1IsValid();
                     bpMtr(chan,1) = jTerm.g2();
                     bpValidMtr(chan,1) = jTerm.g2IsValid();
                  }
                  catch (const CheckError&) {
                     bpValidMtr.row(chan).set(false);
                  }
             }
        }
   }
}

/// @brief store bandpass into a calibration accessor
/// @details This method does the opposite of loadBandpass and stores the current bandpass parameters using the supplied calibration 
/// solution accessor. Only valid parameters are stored (which ensures that they remain flagged according to the general calibration accessor logic),
/// but an exception is expected if, say, the number of antennas or beams in the current bandpass with valid data is larger than the given accessor 
/// can handle. Note, this cannot occur if the output accessor and this class are initialised with the same antennas, beams, etc numbers or if 
/// the bandpass of the same array has been loaded via loadBandpass first before any manipulation.
/// @param[in] acc calibration accessor to store the bandpass into
void BandpassDelayHelper::storeBandpass(accessors::ICalSolutionAccessor &acc) const
{
   const casacore::IPosition shape = itsBandpass.shape();
   ASKAPDEBUGASSERT(shape.nelements() == 4);
   // note, implicitly cast to unsigned types as they're accepted by the interface
   const casacore::uInt nChan = shape[0];
   const casacore::uInt nBeam = shape[2];
   const casacore::uInt nAnt = shape[3];
   ASKAPDEBUGASSERT(shape[1] == 2);
   for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
        for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
             // taking the slice of bandpass and the corresponding array of validity flags for the given antenna and beam
             // due to the selected order of axes, the resulting 2D array should be contiguous
             const casacore::IPosition start(4, 0, 0, beam, ant);
             const casacore::IPosition length(4, nChan, 2, 1, 1);
             const casacore::Matrix<casacore::Complex> bpMtr = itsBandpass(casacore::Slicer(start, length)).nonDegenerate(2u);
             const casacore::Matrix<bool> bpValidMtr = itsBandpassValid(casacore::Slicer(start, length)).nonDegenerate(2u);
             ASKAPDEBUGASSERT(bpMtr.shape() == bpValidMtr.shape());
             ASKAPDEBUGASSERT(bpMtr.nrow() == nChan);
             ASKAPDEBUGASSERT(bpMtr.ncolumn() == 2u);

             const askap::scimath::JonesIndex index(ant, beam);
             for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                  const bool g1IsValid = bpValidMtr(chan,0);
                  const bool g2IsValid = bpValidMtr(chan,1);
                  // ignore points without valid data in at least one of the polarisations
                  if (g1IsValid || g2IsValid) {
                     const askap::accessors::JonesJTerm jTerm(bpMtr(chan,0), g1IsValid, bpMtr(chan,1), g2IsValid);
                     acc.setBandpass(index, jTerm, chan); 
                  }
             }
        }
   }
}

/// @brief initialise the class with ideal bandpass
/// @details This method initialises the class with the amplitude of 1. and zero phase. It doesn't affect the delays if they're calculated earlier.
/// Thus, this method can be used to export a fake bandpass containing just the delay slope or the inverted one providing the correction term.
/// This option exists largely for experiments. All antennas, beams, channels and polarisations this class has been setup to handle are marked as valid.
/// However, the validity flags of delays will be respected if delays are applied to this ideal bandpass.
/// @param[in] resolution spectral resolution in Hz to be used in subsequent calcDelays and applyDelays calls (as for the bandpass loaded via the accessor)
void BandpassDelayHelper::setIdealBandpass(double resolution)
{
   itsResolution = resolution;
   itsBandpass.set(casacore::Complex(1.f, 0.f));
   itsBandpassValid.set(true);
}

/// @brief summarise delays in the log
/// @details This method adds a summary of the current delays in the log. If only one beam has unflagged data, the corresponding delay is printed. 
/// Otherwise, the range is shown. In the future we can add methods which export delays, e.g. into an ascii file. The output of this method only makes sense if
/// delays are calculated (e.g. with calcDelays).
void BandpassDelayHelper::summary() const
{
   const casacore::uInt nBeam = itsDelay.ncolumn();
   const casacore::uInt nAnt = itsDelay.nplane();
   ASKAPDEBUGASSERT(itsDelay.nrow() == 2u);
   ASKAPDEBUGASSERT(itsDelay.shape() == itsDelayValid.shape());
   ASKAPDEBUGASSERT(itsBandpass.shape().nelements() == 4u);
   const casacore::uInt nChan = itsBandpass.shape()[0];
   // assess channel flag statistics
   casacore::uInt chanWithData = 0u;
   for (casacore::uInt chan = 0u; chan < nChan; ++chan) {
        // not a very effective access pattern, but ok for now
        const casacore::IPosition start(4, chan, 0, 0, 0);
        const casacore::IPosition length(4, 1, 2, nBeam, nAnt);
        if (casacore::anyTrue(itsBandpassValid(casacore::Slicer(start, length)))) {
            ++chanWithData;
        } 
   }

   ASKAPLOG_INFO_STR(logger, "Bandpass helper has been configured for "<<nAnt<<" antennas, "<<nBeam<<" beams and "<<nChan<<
                             " channels ("<<chanWithData<<" valid ones)");
   ASKAPLOG_INFO_STR(logger, "Assumed spectral resolution: "<<std::setprecision(5)<<itsResolution / 1e3<<" kHz");
   // first assess flagging per beam to know if we have the results for single beam only (or no valid data at all)
   std::set<casacore::uInt> beamsWithData;
   for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
        // this is not the best access pattern (it's a strided slice in memory), but we do it only once
        const casacore::Matrix<bool> thisBeamDelayValid = itsDelayValid.xzPlane(beam);
        if (casacore::anyTrue(thisBeamDelayValid)) {
            beamsWithData.insert(beam);
        }
   }
   if (beamsWithData.size() == 0) {
       ASKAPLOG_WARN_STR(logger, "All delays appear to be flagged as invalid");
       return;
   }
   const bool oneBeamCase = (beamsWithData.size() == 1);
   // for one beam case the following variable will contain the only beam width data, otherwise the one with smallest index
   const casacore::uInt validBeam = *beamsWithData.begin();
   if (oneBeamCase) {
       ASKAPLOG_INFO_STR(logger, "Delays appear to be valid for one beam only (0-based index "<<validBeam<<")");
   }
   ASKAPDEBUGASSERT(validBeam < nBeam);
   // now print the delays into the log
   for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
        const casacore::Matrix<bool> delayValidMtr = itsDelayValid.xyPlane(ant);
        const casacore::Matrix<float> delayMtr = itsDelay.xyPlane(ant);
        for (casacore::uInt pol = 0; pol < 2u; ++pol) {
             std::ostringstream ss;
             if (pol == 0) {
                 ss<<"antenna (1-based) "<<std::setw(2)<<(ant + 1)<<" XX: ";
             } else {
                 ss<<"                     YY: ";
             }
             const casacore::Vector<bool> delayValidVec = delayValidMtr.row(pol);
             const casacore::Vector<float> delayVec = delayMtr.row(pol);
             if (oneBeamCase) {
                 ASKAPDEBUGASSERT(validBeam < delayValidVec.nelements());
                 if (delayValidVec[validBeam]) {
                    ss<<std::setw(9)<<std::setprecision(2)<<delayVec[validBeam]*1e9<<" ns";
                 } else {
                    ss<<"     flagged";
                 }
             } else {
                 // for multiple beams with valid data only print smallest and largest delay 
                 if (casacore::anyTrue(delayValidVec)) {
                     // use the fact that validBeam is the beam with valid data which has the smallest index
                     // (std::set is by default sorted in the ascending order as per C++ standard)
                     ASKAPDEBUGASSERT(validBeam < delayValidVec.nelements());
                     float minDelay = delayVec[validBeam];
                     float maxDelay = delayVec[validBeam];
                     casacore::uInt minDelayBeam = validBeam;
                     casacore::uInt maxDelayBeam = validBeam;
                     for (casacore::uInt beam = validBeam + 1; beam < nBeam; ++beam) {
                          ASKAPDEBUGASSERT(beam < delayValidVec.nelements());
                          if (delayValidVec[beam]) {
                              const float thisDelay = delayVec[beam];
                              if (minDelay > thisDelay) {
                                  minDelay = thisDelay;
                                  minDelayBeam = beam;
                              }
                              if (maxDelay < thisDelay) {
                                  maxDelay = thisDelay;
                                  maxDelayBeam = beam;
                              }
                          }
                     }
                     ss<<"from "<<std::setw(9)<<std::setprecision(2)<<minDelay*1e9<<" ns (beam "<<std::setw(2)<<minDelayBeam<<") to "<<
                                  std::setw(9)<<std::setprecision(2)<<maxDelay*1e9<<" ns (beam "<<std::setw(2)<<maxDelayBeam<<")";
                 } else {
                   // all beams are flagged for this antenna
                   ss<<std::setw(53)<<std::right<<"all beams flagged";
                 }
             } // one beam vs many beams 
             ASKAPLOG_INFO_STR(logger, "    "<<ss.str());
        } // loop by polarisations
   } // loop by antennas
}

} // namespace utils

} // namespace askap
