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

#include <measurementequation/VisMetaDataStats.h>
#include <askap/AskapError.h>
#include <Blob/BlobArray.h>

namespace askap {

namespace synthesis {

/// @brief constructor, initialise class 
VisMetaDataStats::VisMetaDataStats() : itsTangentSet(false), itsAccessorAdapter(-1.), itsNVis(0ul), itsMaxU(0.), 
     itsMaxV(0.), itsMaxW(0.), itsMaxResidualW(0.), itsMinFreq(0.), itsMaxFreq(0.),
     itsMaxAntennaIndex(0u), itsMaxBeamIndex(0u), itsRefDirValid(false), itsFieldBLC(0.,0.), itsFieldTRC(0.,0.) {}

/// @brief constructor with explicitly given tangent point
/// @details We need to know tangent point to estimate the w-term correctly
/// (tangent point is required for uvw-rotation). Unless the tangent point
/// is chosen in advance, a two-pass iteration over the data is required. 
/// The first iteration is used to find out the centre of the field which
/// can be used as a tangent point during imaging. The second pass determines
/// actual stats on the w-term. In the second pass, this class is initialised 
/// with either this version of the constructor or the version specific for 
/// the snap-shot imaging.
/// @param[in] wtolerance threshold triggering fitting of a new plane for snap-shot imaging (wavelengths)      
VisMetaDataStats::VisMetaDataStats(const casa::MVDirection &tangent) : itsTangent(tangent), itsTangentSet(true), itsAccessorAdapter(-1.),
     itsNVis(0ul), itsMaxU(0.), itsMaxV(0.), itsMaxW(0.), itsMaxResidualW(0.), itsMinFreq(0.), itsMaxFreq(0.), 
     itsMaxAntennaIndex(0u), itsMaxBeamIndex(0u), itsReferenceDir(tangent), itsRefDirValid(true), itsFieldBLC(0.,0.), itsFieldTRC(0.,0.) {} 

   
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
VisMetaDataStats::VisMetaDataStats(const casa::MVDirection &tangent, double wtolerance) : itsTangent(tangent), itsTangentSet(true), 
     itsAccessorAdapter(wtolerance,false), 
     itsNVis(0ul), itsMaxU(0.), itsMaxV(0.), itsMaxW(0.), itsMaxResidualW(0.), itsMinFreq(0.), itsMaxFreq(0.),
     itsMaxAntennaIndex(0u), itsMaxBeamIndex(0u), itsReferenceDir(tangent), itsRefDirValid(true), itsFieldBLC(0.,0.), itsFieldTRC(0.,0.)
     {}


/// @brief copy constructor
/// @details An explicit copy constructor is required because accessor adapter is a non-copyable non-trivial type.
/// This causes problems in the parallel mode only when different mpi ranks may need to clone the statistics
/// estimator as part of the reduction process. 
/// @param[in] other const reference to the object to copy from
VisMetaDataStats::VisMetaDataStats(const VisMetaDataStats &other) : itsTangent(other.itsTangent), itsTangentSet(other.itsTangentSet),
   itsAccessorAdapter(other.itsAccessorAdapter.tolerance(), false), itsNVis(other.itsNVis), itsMaxU(other.itsMaxU),
   itsMaxV(other.itsMaxV), itsMaxW(other.itsMaxW), itsMaxResidualW(other.itsMaxResidualW), itsMinFreq(other.itsMinFreq), 
   itsMaxFreq(other.itsMaxFreq), itsMaxAntennaIndex(other.itsMaxAntennaIndex),
   itsMaxBeamIndex(other.itsMaxBeamIndex), itsReferenceDir(other.itsReferenceDir), itsRefDirValid(other.itsRefDirValid), 
   itsFieldBLC(other.itsFieldBLC), itsFieldTRC(other.itsFieldTRC)
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
           const casa::MVDirection otherBLCDir = other.getOffsetDir(other.itsFieldBLC);
           const casa::MVDirection otherTRCDir = other.getOffsetDir(other.itsFieldTRC);
           const std::pair<double,double> otherBLC = getOffsets(otherBLCDir);
           const std::pair<double,double> otherTRC = getOffsets(otherTRCDir);
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
   }
}

/// @brief helper method to apply an offset to the current reference direction
/// @details
/// @param[in] offsets pair of offsets to apply
/// @return direction measure
casa::MVDirection VisMetaDataStats::getOffsetDir(const std::pair<double,double> &offsets) const
{
  ASKAPCHECK(itsRefDirValid, "getOffsetDir() called before any visibility has been processed, nvis="<<nVis());
  casa::MVDirection result(itsReferenceDir);
  result.shift(offsets.first,offsets.second,casa::True);
  return result;
}
   
/// @brief helper method to compute offsets of the given direction w.r.t. the reference direction
/// @details
/// @param[in] dir direction measure
/// @return pair with offsets w.r.t. the reference direction
std::pair<double,double> VisMetaDataStats::getOffsets(const casa::MVDirection &dir) const
{
  ASKAPCHECK(itsRefDirValid, "getOffsets() called before any visibility has been processed, nvis="<<nVis());
  const double offset1 = sin(dir.getLong() - itsReferenceDir.getLong()) * cos(dir.getLat());
  const double offset2 = sin(dir.getLat()) * cos(itsReferenceDir.getLat()) - cos(dir.getLat()) * sin(itsReferenceDir.getLat())
                                           * cos(dir.getLong() - itsReferenceDir.getLong());
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
  
  if (!itsRefDirValid) {
      // estimate reference direction as the first encountered dish pointing
      itsReferenceDir = acc.dishPointing1()[0];
      itsRefDirValid = true;
  }
  
  const casa::Vector<casa::MVDirection> &pointingDir = acc.pointingDir1();
  for (casa::uInt row=0; row < acc.nRow(); ++row) {
       const std::pair<double,double> offsets = getOffsets(pointingDir[row]);
       if ( (itsNVis == 0ul) && (row == 0) ) {
            itsFieldBLC = itsFieldTRC = offsets;
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
  
  const double currentMaxFreq = casa::max(acc.frequency());
  const double currentMinFreq = casa::min(acc.frequency());
  const casa::uInt currentMaxAntennaIndex = casa::max(casa::max(acc.antenna1()), casa::max(acc.antenna2()));
  const casa::uInt currentMaxBeamIndex = casa::max(casa::max(acc.feed1()), casa::max(acc.feed2()));
  
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

  const double reciprocalToShortestWavelength = currentMaxFreq / casa::C::c;
  
  if (itsAccessorAdapter.tolerance() >=0.) {
      ASKAPCHECK(itsTangentSet, "wtolerance has to be set together with the tangent point!")
  } 
  
  if (itsTangentSet) {
      const casa::Vector<casa::RigidVector<casa::Double, 3> > &origUVW = acc.rotatedUVW(itsTangent);
      
      if (itsAccessorAdapter.tolerance() >= 0.) {
          itsAccessorAdapter.associate(acc);
          ASKAPDEBUGASSERT(acc.nRow() == itsAccessorAdapter.nRow());
      }
                     
      for (casa::uInt row=0; row < acc.nRow(); ++row) {
           const double currentU = casa::abs(origUVW[row](0)) * reciprocalToShortestWavelength;
           const double currentV = casa::abs(origUVW[row](1)) * reciprocalToShortestWavelength;
           const double currentW = casa::abs(origUVW[row](2)) * reciprocalToShortestWavelength;
           
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
      } 
      if (itsAccessorAdapter.tolerance() >= 0.) {
          const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw = itsAccessorAdapter.rotatedUVW(itsTangent);
          for (casa::uInt row=0; row < itsAccessorAdapter.nRow(); ++row) {
               const double currentResidualW = casa::abs(uvw[row](2)) * reciprocalToShortestWavelength;
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
      // this is the first pass, do the best effort job as exact tangent point is unknown
      const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw = acc.uvw();
      for (casa::uInt row=0; row < acc.nRow(); ++row) {
           const double currentU = casa::abs(uvw[row](0)) * reciprocalToShortestWavelength;
           const double currentV = casa::abs(uvw[row](1)) * reciprocalToShortestWavelength;
           const double currentW = casa::abs(uvw[row](2)) * reciprocalToShortestWavelength;
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
      }
  }

  itsNVis += acc.nRow() * acc.nChannel();
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
casa::MVDirection VisMetaDataStats::centre() const {
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
LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream &os, const casa::MVDirection &dir) {
  os<<dir.get();
  return os;
} 

LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream &is, casa::MVDirection &dir) {
  casa::Vector<casa::Double> angles;
  is>>angles;
  ASKAPCHECK(angles.nelements() == 2, "Expect two-element array with angles for a direction measure");
  dir.setAngle(angles[0],angles[1]);
  return is;
}

// increment the number when format changes
#define VIS_META_DATA_STATS_STREAM_VERSION 1

/// @brief write the object to a blob stream
/// @param[in] os the output stream
void VisMetaDataStats::writeToBlob(LOFAR::BlobOStream& os) const
{
  ASKAPCHECK(!itsAccessorAdapter.isAssociated(), "An attempt to serialise VisMetaDataStats with accessor adapter in the attached state");
  os.putStart("VisMetaDataStats",VIS_META_DATA_STATS_STREAM_VERSION);
  os<<itsTangent<<itsTangentSet<<itsAccessorAdapter.tolerance()<<(LOFAR::TYPES::uint64)itsNVis<<itsMaxU<<itsMaxV<<itsMaxW<<itsMaxResidualW<<
    itsMinFreq<<itsMaxFreq<<itsMaxAntennaIndex<<itsMaxBeamIndex<<itsReferenceDir<<itsRefDirValid<<
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
  is >> itsMaxU >> itsMaxV >> itsMaxW >> itsMaxResidualW >> itsMinFreq >> itsMaxFreq >> itsMaxAntennaIndex >> itsMaxBeamIndex >>
        itsReferenceDir >> itsRefDirValid >> itsFieldBLC.first >> itsFieldBLC.second >> itsFieldTRC.first >> itsFieldTRC.second;     
  is.getEnd();       
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
          offsets.first = casa::max(itsFieldBLC.first, itsFieldTRC.first);
          offsets.second = casa::max(itsFieldBLC.second, itsFieldTRC.second);
      }
  }
    
  // primary beam fwhm for a 12m antenna
  const double longestWavelength = casa::C::c / minFreq(); // in metres
  const double pbFWHM = 1.2 * longestWavelength / 12; // in radians
  // the guard band (both sides together) is 1.7*FWHM (roughly to the first null)  
  const double sizeInRad =  2. * casa::max(offsets.first, offsets.second) + 1.7 * pbFWHM;
  return sizeInRad / casa::C::pi * 180.;
}

/// @brief estimate cell size
/// @details This method uses maxU and maxV to estimate the largest (square) image cell size in arcsec
/// @return square cell size in arcsec
double VisMetaDataStats::squareCellSize() const
{
  const double largestSpacing = casa::max(maxU(), maxV());
  // Nyquist sampling corresponds to 1/2, 1/6 is the minumum used in practice to achieve a reasonable image quality 
  const double cellSizeInRad = 1./largestSpacing/6.;
  return cellSizeInRad / casa::C::pi * 6.48e5; 
}


} // namespace synthesis

} // namespace askap

