/// @file
///
/// @brief An estimator of statistics for metadata associated with visibilities
/// @details Some configuration parameters depend on the metadata, for example
/// cell size depends on the largest baseline. The ASKAP approach is to set
/// all parameters like this a priori to avoid an additional iteration over data.
/// For BETA we could afford iteration over the dataset and, therefore, an "advise"
/// utility could be written. This class handles basic statistics to assist with this.
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
/// @author Max Voronkov <maxim.voronkov@csiro.au
///

#include <askap/measurementequation/VisMetaDataStats.h>
#include <askap/askap/AskapError.h>
#include <Blob/BlobArray.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

ASKAP_LOGGER(logger, ".measurementequation.vismetadatastats");
#include <iomanip>

namespace askap {

namespace synthesis {

/// @brief constructor, initialise class
/// @param[in] useFlagged boolean indicating we take flagged data into account in the stats
/// @param[in] wpercentile, percentile value to calculate for W distribution
VisMetaDataStats::VisMetaDataStats(bool useFlagged, double wpercentile) : itsTangentSet(false), itsUseFlagged(useFlagged),itsAccessorAdapter(-1.),
     itsNVis(0ul), itsMaxU(0.), itsMaxV(0.), itsMaxW(0.), itsMaxResidualW(0.), itsMinFreq(0.), itsMaxFreq(0.),
     itsMaxAntennaIndex(0u), itsMaxBeamIndex(0u), itsRefDirValid(false), itsFieldBLC(0.,0.), itsFieldTRC(0.,0.), itsWPercentileCalculator(wpercentile) {}

/// @brief constructor with explicitly given tangent point
/// @details We need to know tangent point to estimate the w-term correctly
/// (tangent point is required for uvw-rotation). Unless the tangent point
/// is chosen in advance, a two-pass iteration over the data is required.
/// The first iteration is used to find out the centre of the field which
/// can be used as a tangent point during imaging. The second pass determines
/// actual stats on the w-term. In the second pass, this class is initialised
/// with either this version of the constructor or the version specific for
/// the snap-shot imaging.
/// @param[in] tangent tangent point to be used with snap-shot imaging (for uvw-rotation)
/// @param[in] useFlagged boolean indicating we take flagged data into account in the stats
/// @param[in] wpercentile, percentile value to calculate for W distribution
VisMetaDataStats::VisMetaDataStats(const casacore::MVDirection &tangent, bool useFlagged, double wpercentile) : itsTangent(tangent), itsTangentSet(true), itsUseFlagged(useFlagged),
     itsAccessorAdapter(-1.), itsNVis(0ul), itsMaxU(0.), itsMaxV(0.), itsMaxW(0.), itsMaxResidualW(0.), itsMinFreq(0.), itsMaxFreq(0.),
     itsMaxAntennaIndex(0u), itsMaxBeamIndex(0u), itsReferenceDir(tangent), itsRefDirValid(true), itsFieldBLC(0.,0.), itsFieldTRC(0.,0.), itsWPercentileCalculator(wpercentile) {}


/// @brief constructor specific to snap-shot imaging
/// @details For the snap-shot imaging we need to do two passes unless the desired tangent point
/// can be specified up front. The first pass can be used to find out the centre of the field
/// which can be used as a tangent point during imaging. The second pass, where the class is setup
/// with this version of the constructor, can determine the largest residual w-term for the
/// given tangent point and w-tolerance.
/// @note For a coplanar array the largest residual w-term will always be less than the w-tolerance
/// which is a threshold for the fitting of a new plane. For non-coplanar array it is not always the
/// case. This is why a complex two-pass estimation procedure is required.
/// @param[in] tangent tangent point to be used with snap-shot imaging (for uvw-rotation)
/// @param[in] wtolerance threshold triggering fitting of a new plane for snap-shot imaging (wavelengths)
/// @param[in] useFlagged boolean indicating we take flagged data into account in the stats
/// @param[in] wpercentile, percentile value to calculate for W distribution
VisMetaDataStats::VisMetaDataStats(const casacore::MVDirection &tangent, double wtolerance, bool useFlagged, double wpercentile) : itsTangent(tangent),
     itsTangentSet(true), itsUseFlagged(useFlagged), itsAccessorAdapter(wtolerance,false),
     itsNVis(0ul), itsMaxU(0.), itsMaxV(0.), itsMaxW(0.), itsMaxResidualW(0.), itsMinFreq(0.), itsMaxFreq(0.),
     itsMaxAntennaIndex(0u), itsMaxBeamIndex(0u), itsReferenceDir(tangent), itsRefDirValid(true), itsFieldBLC(0.,0.), itsFieldTRC(0.,0.),
     itsWPercentileCalculator(wpercentile)
     {}


/// @brief copy constructor
/// @details An explicit copy constructor is required because accessor adapter is a non-copyable non-trivial type.
/// This causes problems in the parallel mode only when different mpi ranks may need to clone the statistics
/// estimator as part of the reduction process.
/// @param[in] other const reference to the object to copy from
VisMetaDataStats::VisMetaDataStats(const VisMetaDataStats &other) : itsTangent(other.itsTangent), itsTangentSet(other.itsTangentSet),
   itsUseFlagged(other.itsUseFlagged), itsAccessorAdapter(other.itsAccessorAdapter.tolerance(), false), itsNVis(other.itsNVis), itsMaxU(other.itsMaxU),
   itsMaxV(other.itsMaxV), itsMaxW(other.itsMaxW), itsMaxResidualW(other.itsMaxResidualW), itsMinFreq(other.itsMinFreq),
   itsMaxFreq(other.itsMaxFreq), itsMaxAntennaIndex(other.itsMaxAntennaIndex),
   itsMaxBeamIndex(other.itsMaxBeamIndex), itsReferenceDir(other.itsReferenceDir), itsRefDirValid(other.itsRefDirValid),
   itsFieldBLC(other.itsFieldBLC), itsFieldTRC(other.itsFieldTRC), itsWPercentileCalculator(other.itsWPercentileCalculator)
{
  ASKAPCHECK(!other.itsAccessorAdapter.isAssociated(), "An attempt to copy VisMetaDataStats which has the accessor adapter in the attached state");
}

/// @brief reset all accumulated statistics
/// @details After this method, the object will be reset to a pristine state preserving only
/// parameters passed in the constructor, i.e. tangent and wtolerance (if they're defined).
/// All accumulated statistics are reset.
void VisMetaDataStats::reset()
{
   itsNVis = 0ul;
   itsMaxU = 0.;
   itsMaxV = 0.;
   itsMaxW = 0.;
   itsMaxResidualW = 0.;
   itsMaxAntennaIndex = 0u;
   itsMaxBeamIndex = 0u;
   itsFieldBLC = std::pair<double,double>(0.,0.);
   itsFieldTRC = std::pair<double,double>(0.,0.);
   if (itsTangentSet) {
       itsReferenceDir = itsTangent;
       itsRefDirValid  = true;
   } else {
       itsRefDirValid = false;
   }
   itsWPercentileCalculator.init();
}

/// @brief aggregate statistics with that accumulated by another instance
/// @details This class will be run in parallel if the measurement set is distributed.
/// This method is intended to combine statistics as part of reduction.
/// @param[in] other an instance of the estimator to take data from
void VisMetaDataStats::merge(const VisMetaDataStats &other)
{
   ASKAPCHECK(itsTangentSet == other.itsTangentSet, "Different tangent settings detected during VisMetaDataStats merge");
   if (itsTangentSet) {
       ASKAPCHECK(itsTangent.separation(other.itsTangent)<1e-6, "Different tangent directions are used by merged VisMetaDataStats instances");
   }
   ASKAPCHECK(fabs(itsAccessorAdapter.tolerance() - other.itsAccessorAdapter.tolerance()) < 1e-6,
           "Different configuration of w-tolerance detected during VisMetaDataStats merge");
   if (other.nVis() != 0ul) {
       ASKAPDEBUGASSERT(other.itsRefDirValid);
       if (nVis() == 0ul) {
           // just copy all data fields from other
           itsNVis = other.itsNVis;
           itsMaxU = other.itsMaxU;
           itsMaxV = other.itsMaxV;
           itsMaxW = other.itsMaxW;
           itsMaxResidualW = other.itsMaxResidualW;
           itsMinFreq = other.itsMinFreq;
           itsMaxFreq = other.itsMaxFreq;
           itsMaxAntennaIndex = other.itsMaxAntennaIndex;
           itsMaxBeamIndex = other.itsMaxBeamIndex;
           itsRefDirValid = true;
           itsReferenceDir = other.itsReferenceDir;
           itsFieldBLC = other.itsFieldBLC;
           itsFieldTRC = other.itsFieldTRC;
       } else {
           // need to merge properly
           itsNVis += other.itsNVis;
           if (other.itsMaxU > itsMaxU) {
               itsMaxU = other.itsMaxU;
           }
           if (other.itsMaxV > itsMaxV) {
               itsMaxV = other.itsMaxV;
           }
           if (other.itsMaxW > itsMaxW) {
               itsMaxW = other.itsMaxW;
           }
           if (other.itsMaxResidualW > itsMaxResidualW) {
               itsMaxResidualW = other.itsMaxResidualW;
           }
           if (other.itsMaxFreq > itsMaxFreq) {
               itsMaxFreq = other.itsMaxFreq;
           }
           if (other.itsMinFreq < itsMinFreq) {
               itsMinFreq = other.itsMinFreq;
           }
           if (other.itsMaxAntennaIndex > itsMaxAntennaIndex) {
               itsMaxAntennaIndex = other.itsMaxAntennaIndex;
           }
           if (other.itsMaxBeamIndex > itsMaxBeamIndex) {
               itsMaxBeamIndex = other.itsMaxBeamIndex;
           }
           // adjust direction stats taking into account that the reference direction
           // may be different in these two classes
           ASKAPDEBUGASSERT(itsRefDirValid);
           const casacore::MVDirection otherBLCDir = other.getOffsetDir(other.itsFieldBLC);
           const casacore::MVDirection otherTRCDir = other.getOffsetDir(other.itsFieldTRC);
           std::pair<double,double> otherBLC = getOffsets(otherBLCDir);
           std::pair<double,double> otherTRC = getOffsets(otherTRCDir);

           // a hack related to ASKAPSDP-1741. Due to accuracy of MVDirection::shift, referencing
           // of offsets from another object to the same reference may result in inconsistent blc/trc
           if ((otherTRC.first < otherBLC.first) && ((otherBLC.first - otherTRC.first) < 1e-15)) {
               const double temp = otherTRC.first;
               otherTRC.first = otherBLC.first;
               otherBLC.first = temp;
           }
           if ((otherTRC.second < otherBLC.second) && ((otherBLC.second - otherTRC.second) < 1e-15)) {
               const double temp = otherTRC.second;
               otherTRC.second = otherBLC.second;
               otherBLC.second = temp;
           }
           //

           /*
           ASKAPLOG_DEBUG_STR(logger, "otherBLCDir="<<printDirection(otherBLCDir)<<" -> "<<std::setprecision(15)<<otherBLC.first<<" "<<otherBLC.second);
           ASKAPLOG_DEBUG_STR(logger, "otherTRCDir="<<printDirection(otherTRCDir)<<" -> "<<std::setprecision(15)<<otherTRC.first<<" "<<otherTRC.second);
           ASKAPLOG_DEBUG_STR(logger, "trc-blc="<<std::setprecision(15)<<otherTRC.first-otherBLC.first<<" "<<otherTRC.second - otherBLC.second);
           ASKAPLOG_DEBUG_STR(logger, "other(trc-blc)="<<std::setprecision(15)<<other.itsFieldTRC.first-other.itsFieldBLC.first<<" "<<other.itsFieldTRC.second - other.itsFieldBLC.second);
           ASKAPLOG_DEBUG_STR(logger, "refDir="<<printDirection(itsReferenceDir)<<" vs other.refDir="<<printDirection(other.itsReferenceDir));
           */

           // cross checks just in case
           ASKAPDEBUGASSERT(otherBLC.first <= otherTRC.first);
           ASKAPDEBUGASSERT(otherBLC.second <= otherTRC.second);
           //
           if (itsFieldBLC.first > otherBLC.first) {
               itsFieldBLC.first = otherBLC.first;
           }
           if (itsFieldTRC.first < otherTRC.first) {
               itsFieldTRC.first = otherTRC.first;
           }
           if (itsFieldBLC.second > otherBLC.second) {
               itsFieldBLC.second = otherBLC.second;
           }
           if (itsFieldTRC.second < otherTRC.second) {
               itsFieldTRC.second = otherTRC.second;
           }
       }
       itsWPercentileCalculator.merge(other.itsWPercentileCalculator);
   }
}

/// @brief helper method to apply an offset to the current reference direction
/// @details
/// @param[in] offsets pair of offsets to apply
/// @return direction measure
casacore::MVDirection VisMetaDataStats::getOffsetDir(const std::pair<double,double> &offsets) const
{
  ASKAPCHECK(itsRefDirValid, "getOffsetDir() called before any visibility has been processed, nvis="<<nVis());
  casacore::MVDirection result(itsReferenceDir);
  result.shift(offsets.first,offsets.second,casacore::True);
  return result;
  // Note, the accuracy of the shift method doesn't seem to be good enough for some reason.
  // (found while debugging ASKAPSDP-1741)
}

/// @brief helper method to compute offsets of the given direction w.r.t. the reference direction
/// @details
/// @param[in] dir direction measure
/// @return pair with offsets w.r.t. the reference direction
std::pair<double,double> VisMetaDataStats::getOffsets(const casacore::MVDirection &dir) const
{
  ASKAPCHECK(itsRefDirValid, "getOffsets() called before any visibility has been processed, nvis="<<nVis());
  const double offset1 = asin(sin(dir.getLong() - itsReferenceDir.getLong()) * cos(dir.getLat()));
  const double offset2 = asin(sin(dir.getLat()) * cos(itsReferenceDir.getLat()) - cos(dir.getLat()) * sin(itsReferenceDir.getLat())
                                           * cos(dir.getLong() - itsReferenceDir.getLong()));
  return std::pair<double,double>(offset1,offset2);
}


/// @brief process one accessor of data updating statistics
/// @details
/// @param[in] acc read-only accessor with data
void VisMetaDataStats::process(const accessors::IConstDataAccessor &acc)
{
  if (acc.nRow() == 0) {
      return; // no data - nothing to do
  }
  // for now ignore flagging. Technically, some metadata may be ignored if all corresponding data are flagged, but
  // it seems to be too much of the complication now.

  // well the time to use the flags may have come...
  // Note that this is a lot slower, but nearly all the extra time is spent in the acc.flag() call, not here
  const casacore::Cube<casacore::Bool> flag(itsUseFlagged ? casacore::Cube<casacore::Bool>(0,0,0) : acc.flag());
  casacore::Vector<casacore::Bool> rowFlag(acc.nRow(),casacore::False);

  // fill row flags - accessor doesn't provide MS row flag, so create them
  int unflaggedRows = 0;
  if (!itsUseFlagged) {
        for (casacore::uInt row = 0; row < acc.nRow(); row++) {
            if (casacore::allEQ(flag(row,casacore::Slice(),casacore::Slice()),casacore::True)) {
                rowFlag(row) = casacore::True;
            } else {
                unflaggedRows++;
            }
        }
        if (casacore::allEQ(rowFlag,casacore::True)) {
            // no unflagged data
            return;
        }
  }

  // fill channel flags
  casacore::uInt nChan = acc.nChannel();
  casacore::Vector<casacore::Bool> chanFlag(nChan, casacore::False);
  int unflaggedChannels = 0;
  if (!itsUseFlagged) {
      for (casacore::uInt chan = 0; chan < nChan; chan++) {
          if (casacore::allEQ(flag(casacore::Slice(),chan,casacore::Slice()),casacore::True)) {
              chanFlag(chan) = casacore::True;
          } else {
              unflaggedChannels++;
          }
      }
  }

  if (!itsRefDirValid) {
      // estimate reference direction as the first encountered dish pointing
      int i = 0;
      if (!itsUseFlagged) {
          for (casacore::uInt row = 0; row < acc.nRow(); row++) {
              if (!rowFlag(row)) {
                  i = row;
                  break;
              }
          }
      }
      itsReferenceDir = acc.dishPointing1()[i];
      itsRefDirValid = true;
  }

  const casacore::Vector<casacore::MVDirection> &pointingDir = acc.pointingDir1();
  //ASKAPLOG_DEBUG_STR(logger, "referenceDir = "<<printDirection(itsReferenceDir));
  bool first = true;
  for (casacore::uInt row=0; row < acc.nRow(); ++row) {
      if (!itsUseFlagged && rowFlag(row)) {
          continue;
      }
       const std::pair<double,double> offsets = getOffsets(pointingDir[row]);
       if ( (itsNVis == 0ul) && first ) {
            itsFieldBLC = itsFieldTRC = offsets;
            first = false;
       } else {
            if (itsFieldBLC.first > offsets.first) {
                itsFieldBLC.first = offsets.first;
            }
            if (itsFieldTRC.first < offsets.first) {
                itsFieldTRC.first = offsets.first;
            }
            if (itsFieldBLC.second > offsets.second) {
                itsFieldBLC.second = offsets.second;
            }
            if (itsFieldTRC.second < offsets.second) {
                itsFieldTRC.second = offsets.second;
            }
       }
  }
  //ASKAPLOG_DEBUG_STR(logger, "after iteration over "<<acc.nRow()<<" rows blc: "<<std::setprecision(15)<<itsFieldBLC.first<<" "<<itsFieldBLC.second<<" trc: "<<itsFieldTRC.first<<" "<<itsFieldTRC.second);
  double currentMaxFreq = 0;
  double currentMinFreq = -1;
  for (casacore::uInt chan = 0; chan < nChan; chan++) {
      if (!itsUseFlagged && chanFlag(chan)) {
          continue;
      }
      const double freq = acc.frequency()(chan);
      if (freq > currentMaxFreq) {
          currentMaxFreq = freq;
      }
      if (currentMinFreq < 0 || freq < currentMinFreq ) {
          currentMinFreq = freq;
      }
  }

  casacore::uInt currentMaxAntennaIndex = 0;
  casacore::uInt currentMaxBeamIndex = 0;
  for (casacore::uInt row = 0; row < acc.nRow(); row++) {
      if (!itsUseFlagged && rowFlag(row)) {
          continue;
      }
      const casacore::uInt ant = casacore::max(acc.antenna1()(row),acc.antenna2()(row));
      if (ant > currentMaxAntennaIndex) {
          currentMaxAntennaIndex = ant;
      }
      const casacore::uInt beam = casacore::max(acc.feed1()(row),acc.feed2()(row));
      if (beam > currentMaxBeamIndex) {
          currentMaxBeamIndex = beam;
      }
  }

  if (itsNVis == 0ul) {
      itsMinFreq = currentMinFreq;
      itsMaxFreq = currentMaxFreq;
      itsMaxAntennaIndex = currentMaxAntennaIndex;
      itsMaxBeamIndex = currentMaxBeamIndex;
  } else {
      if (itsMinFreq > currentMinFreq) {
          itsMinFreq = currentMinFreq;
      }
      if (itsMaxFreq < currentMaxFreq) {
          itsMaxFreq = currentMaxFreq;
      }
      if (itsMaxAntennaIndex < currentMaxAntennaIndex) {
          itsMaxAntennaIndex = currentMaxAntennaIndex;
      }
      if (itsMaxBeamIndex < currentMaxBeamIndex) {
          itsMaxBeamIndex = currentMaxBeamIndex;
      }
  }

  const double reciprocalToShortestWavelength = currentMaxFreq / casacore::C::c;

  if (itsAccessorAdapter.tolerance() >=0.) {
      ASKAPCHECK(itsTangentSet, "wtolerance has to be set together with the tangent point!")
  }
  unsigned long nvis = 0;
  if (itsTangentSet) {
      const casacore::Vector<casacore::RigidVector<casacore::Double, 3> > &origUVW = acc.rotatedUVW(itsTangent);

      if (itsAccessorAdapter.tolerance() >= 0.) {
          itsAccessorAdapter.associate(acc);
          ASKAPDEBUGASSERT(acc.nRow() == itsAccessorAdapter.nRow());
      }

      const casacore::Vector<casacore::Double> & freq = acc.frequency();
      if (itsUseFlagged) {
          for (casacore::uInt row=0; row < acc.nRow(); ++row) {
              const double currentU = casacore::abs(origUVW[row](0)) * reciprocalToShortestWavelength;
              const double currentV = casacore::abs(origUVW[row](1)) * reciprocalToShortestWavelength;
              const double currentW = casacore::abs(origUVW[row](2)) * reciprocalToShortestWavelength;

              if ((itsNVis == 0ul) && (row == 0)) {
                  itsMaxU = currentU;
                  itsMaxV = currentV;
                  itsMaxW = currentW;
              } else {
                  if (itsMaxU < currentU) {
                      itsMaxU = currentU;
                  }
                  if (itsMaxV < currentV) {
                      itsMaxV = currentV;
                  }
                  if (itsMaxW < currentW) {
                      itsMaxW = currentW;
                  }
              }
              // Have to do this channel by channel to get distribution of W
              for (casacore::uInt chan=0; chan<nChan; chan++) {
                  const double reciprocalWavelength = freq(chan) / casacore::C::c;
                  itsWPercentileCalculator.add(casacore::abs(origUVW[row](2)) * reciprocalWavelength);
              }
          }
          if (itsAccessorAdapter.tolerance() >= 0.) {
              const casacore::Vector<casacore::RigidVector<casacore::Double, 3> > &uvw = itsAccessorAdapter.rotatedUVW(itsTangent);
              for (casacore::uInt row=0; row < itsAccessorAdapter.nRow(); ++row) {
                   const double currentResidualW = casacore::abs(uvw[row](2)) * reciprocalToShortestWavelength;
                   if ((itsNVis == 0ul) && (row == 0)) {
                       itsMaxResidualW = currentResidualW;
                   } else {
                       if (itsMaxResidualW < currentResidualW) {
                           itsMaxResidualW = currentResidualW;
                       }
                   }
              }
              itsAccessorAdapter.detach();
          }
      } else {
          for (casacore::uInt row=0; row < acc.nRow(); ++row) {
              if (rowFlag(row)) {
                  continue;
              }
              for (casacore::uInt chan=0; chan<nChan; chan++) {
                  if (chanFlag(chan)) {
                      continue;
                  }
                  nvis++;
                  const double reciprocalWavelength = freq(chan) / casacore::C::c;
                  const double currentU = casacore::abs(origUVW[row](0)) * reciprocalWavelength;
                  const double currentV = casacore::abs(origUVW[row](1)) * reciprocalWavelength;
                  const double currentW = casacore::abs(origUVW[row](2)) * reciprocalWavelength;

                  if ((itsNVis == 0ul) && first) {
                      itsMaxU = currentU;
                      itsMaxV = currentV;
                      itsMaxW = currentW;
                      first = false;
                  } else {
                      if (itsMaxU < currentU) {
                          itsMaxU = currentU;
                      }
                      if (itsMaxV < currentV) {
                          itsMaxV = currentV;
                      }
                      if (itsMaxW < currentW) {
                          itsMaxW = currentW;
                      }
                  }
                  itsWPercentileCalculator.add(currentW);
              }
          }
          if (itsAccessorAdapter.tolerance() >= 0.) {
              const casacore::Vector<casacore::RigidVector<casacore::Double, 3> > &uvw = itsAccessorAdapter.rotatedUVW(itsTangent);
              bool first = true;
              for (casacore::uInt row=0; row < itsAccessorAdapter.nRow(); ++row) {
                  if (rowFlag(row)) {
                      continue;
                  }
                  for (casacore::uInt chan=0; chan<nChan; chan++) {
                      if (chanFlag(chan)) {
                          continue;
                      }
                      const double reciprocalWavelength = freq(chan) / casacore::C::c;
                      const double currentResidualW = casacore::abs(uvw[row](2)) * reciprocalWavelength;
                       if ((itsNVis == 0ul) && first) {
                           itsMaxResidualW = currentResidualW;
                           first = false;
                       } else {
                           if (itsMaxResidualW < currentResidualW) {
                               itsMaxResidualW = currentResidualW;
                           }
                       }
                   }
              }
              itsAccessorAdapter.detach();
          }
      }
  } else {
      // this is the first pass, do the best effort job as exact tangent point is unknown
      const casacore::Vector<casacore::RigidVector<casacore::Double, 3> > &uvw = acc.uvw();
      bool first = true;
      for (casacore::uInt row=0; row < acc.nRow(); ++row) {
           if (!itsUseFlagged && rowFlag(row)) {
               continue;
           }
           nvis += unflaggedChannels;
           const double currentU = casacore::abs(uvw[row](0)) * reciprocalToShortestWavelength;
           const double currentV = casacore::abs(uvw[row](1)) * reciprocalToShortestWavelength;
           const double currentW = casacore::abs(uvw[row](2)) * reciprocalToShortestWavelength;
           if ((itsNVis == 0ul) && first) {
               itsMaxU = currentU;
               itsMaxV = currentV;
               itsMaxW = currentW;
               first = false;
           } else {
               if (itsMaxU < currentU) {
                   itsMaxU = currentU;
               }
               if (itsMaxV < currentV) {
                   itsMaxV = currentV;
               }
               if (itsMaxW < currentW) {
                   itsMaxW = currentW;
               }
           }
      }
  }

  if (!itsUseFlagged) {
      itsNVis += nvis;
  } else {
      itsNVis += acc.nRow() * acc.nChannel();
  }
}

/// @brief largest residual w-term (for snap-shotting)
/// @return largest value of residual w in wavelengths
double VisMetaDataStats::maxResidualW() const
{
  ASKAPCHECK(itsAccessorAdapter.tolerance()>=0., "maxResidualW() called for an object not configured for snap-shot imaging");
  return itsMaxResidualW;
}

/// @brief most central direction of the observed field
/// @return direction of the centre in the frame used by the accessor
casacore::MVDirection VisMetaDataStats::centre() const {
  ASKAPCHECK(itsRefDirValid, "centre() called before any visibility has been processed, nvis="<<nVis());
  const std::pair<double,double>  cnt((itsFieldTRC.first + itsFieldBLC.first) / 2, (itsFieldTRC.second + itsFieldBLC.second) / 2);
  return getOffsetDir(cnt);
}

/// @brief largest separation of individual pointing from the centre
/// @return largest offsets from the centre() in radians (measure of the field size)
std::pair<double,double> VisMetaDataStats::maxOffsets() const
{
  ASKAPCHECK(itsRefDirValid, "maxOffset() called before any visibility has been processed, nvis="<<nVis());
  const std::pair<double,double>  result((itsFieldTRC.first - itsFieldBLC.first) / 2, (itsFieldTRC.second - itsFieldBLC.second) / 2);
  ASKAPDEBUGASSERT(result.first >= 0.);
  ASKAPDEBUGASSERT(result.second >= 0.);
  return result;
}

// helper operators, we can move them to Base if they're found useful somewhere else
LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream &os, const casacore::MVDirection &dir) {
  os<<dir.get();
  return os;
}

LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream &is, casacore::MVDirection &dir) {
  casacore::Vector<casacore::Double> angles;
  is>>angles;
  ASKAPCHECK(angles.nelements() == 2, "Expect two-element array with angles for a direction measure");
  dir.setAngle(angles[0],angles[1]);
  return is;
}

// helper operators, we can move them to Base if they're found useful somewhere else
LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream &os, const PercentileCalculator &wpc) {
  wpc.writeToBlob(os);
  return os;
}

LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream &is, PercentileCalculator &wpc) {
  wpc.readFromBlob(is);
  return is;
}

// increment the number when format changes
#define VIS_META_DATA_STATS_STREAM_VERSION 2

/// @brief write the object to a blob stream
/// @param[in] os the output stream
void VisMetaDataStats::writeToBlob(LOFAR::BlobOStream& os) const
{
  ASKAPCHECK(!itsAccessorAdapter.isAssociated(), "An attempt to serialise VisMetaDataStats with accessor adapter in the attached state");
  ASKAPCHECK(itsFieldBLC.first <= itsFieldTRC.first, "Inconsistent blc and trc prior to send, difference.first="<<std::setprecision(15)<<itsFieldTRC.first - itsFieldBLC.first);
  ASKAPCHECK(itsFieldBLC.second <= itsFieldTRC.second, "Inconsistent blc and trc prior to send, difference.second="<<std::setprecision(15)<<itsFieldTRC.second - itsFieldBLC.second);
  os.putStart("VisMetaDataStats",VIS_META_DATA_STATS_STREAM_VERSION);
  os<<itsTangent<<itsTangentSet<<itsAccessorAdapter.tolerance()<<(LOFAR::TYPES::uint64)itsNVis<<itsMaxU<<itsMaxV<<itsMaxW<<itsMaxResidualW<<
  itsWPercentileCalculator<<itsMinFreq<<itsMaxFreq<<itsMaxAntennaIndex<<itsMaxBeamIndex<<itsReferenceDir<<itsRefDirValid<<
    itsFieldBLC.first<<itsFieldBLC.second<<itsFieldTRC.first<<itsFieldTRC.second;
  os.putEnd();
}

/// @brief read the object from a blob stream
/// @param[in] is the input stream
void VisMetaDataStats::readFromBlob(LOFAR::BlobIStream& is)
{
  ASKAPCHECK(!itsAccessorAdapter.isAssociated(), "An attempt to de-serialise VisMetaDataStats with accessor adapter in the attached state");
  const int version = is.getStart("VisMetaDataStats");
  ASKAPCHECK(version == VIS_META_DATA_STATS_STREAM_VERSION,
       "Attempting to read from a blob stream an object of the wrong version, expected "<<VIS_META_DATA_STATS_STREAM_VERSION<<
       "got "<<version);
  double wtolerance = -1;
  LOFAR::TYPES::uint64 nVisBuf = 0;
  is >> itsTangent >> itsTangentSet >> wtolerance >> nVisBuf;
  itsAccessorAdapter = accessors::BestWPlaneDataAccessor(wtolerance, false);
  itsNVis = (unsigned long)nVisBuf;
  is >> itsMaxU >> itsMaxV >> itsMaxW >> itsMaxResidualW >> itsWPercentileCalculator >> itsMinFreq >> itsMaxFreq >> itsMaxAntennaIndex >> itsMaxBeamIndex >>
        itsReferenceDir >> itsRefDirValid >> itsFieldBLC.first >> itsFieldBLC.second >> itsFieldTRC.first >> itsFieldTRC.second;
  is.getEnd();
  // consistency check
  ASKAPCHECK(itsFieldBLC.first <= itsFieldTRC.first, "Inconsistent blc and trc received, difference.first="<<std::setprecision(15)<<itsFieldTRC.first - itsFieldBLC.first);
  ASKAPCHECK(itsFieldBLC.second <= itsFieldTRC.second, "Inconsistent blc and trc received, difference.second="<<std::setprecision(15)<<itsFieldTRC.second - itsFieldBLC.second);
}

/// @brief estimate of the field size
/// @details This method uses maxOffsets, centre and itsTangent to estimate the field size applying
/// current knowledge on the guard band around the edge pointing (hard coded for ASKAP).
/// @param[in] forceCentreAtTangent if true, the tangent point is assumed to be in the image centre
/// @return square field size in degrees
double VisMetaDataStats::squareFieldSize(bool forceCentreAtTangent) const
{

  std::pair<double,double> offsets = maxOffsets();

  if (itsTangentSet) {
      // do extra checks to ensure maxOffsets give offsets w.r.t. the tangent point
      ASKAPCHECK(itsRefDirValid, "Reference direction is not valid! There likely to be a logic error.");
      ASKAPCHECK(itsTangent.separation(itsReferenceDir)<1e-6, "Tangent point looks sufficiently different from the reference direction! There likely to be a logic error.");
      if (forceCentreAtTangent) {
          offsets.first = casacore::max(itsFieldBLC.first, itsFieldTRC.first);
          offsets.second = casacore::max(itsFieldBLC.second, itsFieldTRC.second);
      }
  }

  // primary beam fwhm for a 12m antenna
  const double longestWavelength = casacore::C::c / minFreq(); // in metres
  const double pbFWHM = 1.2 * longestWavelength / 12; // in radians
  // the guard band (both sides together) is 1.7*FWHM (roughly to the first null)
  const double sizeInRad =  2. * casacore::max(offsets.first, offsets.second) + 1.7 * pbFWHM;
  return sizeInRad / casacore::C::pi * 180.;
}

/// @brief estimate cell size
/// @details This method uses maxU and maxV to estimate the largest (square) image cell size in arcsec
/// @return square cell size in arcsec
double VisMetaDataStats::squareCellSize() const
{
  const double largestSpacing = casacore::max(maxU(), maxV());
  // Nyquist sampling corresponds to 1/2, 1/6 is the minumum used in practice to achieve a reasonable image quality
  const double cellSizeInRad = 1./largestSpacing/6.;
  return cellSizeInRad / casacore::C::pi * 6.48e5;
}


} // namespace synthesis

} // namespace askap
