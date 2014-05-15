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

#ifndef SYNTHESIS_VIS_METADATA_STATS_H
#define SYNTHESIS_VIS_METADATA_STATS_H

#include <dataaccess/IConstDataAccessor.h>
#include <dataaccess/BestWPlaneDataAccessor.h>
#include <measures/Measures/MDirection.h>
#include <fitting/ISerializable.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>

#include <utility>

namespace askap {

namespace synthesis {

/// @brief An estimator of statistics for metadata associated with visibilities
/// @details Some configuration parameters depend on the metadata, for example
/// cell size depends on the largest baseline. The ASKAP approach is to set
/// all parameters like this a priori to avoid an additional iteration over data.
/// For BETA we could afford iteration over the dataset and, therefore, an "advise"
/// utility could be written. This class handles basic statistics to assist with this.
/// @note This class is in the measurement equation directory for now, although it doesn't 
/// quite fit with the rest of the stuff up there
/// @ingroup measurementequation
class VisMetaDataStats : public ISerializable {
public:
   /// @brief constructor, initialise class 
   VisMetaDataStats();
   
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
   explicit VisMetaDataStats(const casa::MVDirection &tangent); 

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
   VisMetaDataStats(const casa::MVDirection &tangent, double wtolerance);

   /// @brief copy constructor
   /// @details An explicit copy constructor is required because accessor adapter is a non-copyable non-trivial type.
   /// This causes problems in the parallel mode only when different mpi ranks may need to clone the statistics
   /// estimator as part of the reduction process. 
   /// @param[in] other const reference to the object to copy from
   VisMetaDataStats(const VisMetaDataStats &other);
   
   
   /// @brief reset all accumulated statistics
   /// @details After this method, the object will be reset to a pristine state preserving only
   /// parameters passed in the constructor, i.e. tangent and wtolerance (if they're defined).
   /// All accumulated statistics are reset.
   void reset();
   
   /// @brief aggregate statistics with that accumulated by another instance
   /// @details This class will be run in parallel if the measurement set is distributed. 
   /// This method is intended to combine statistics as part of the reduction.
   /// @param[in] other an instance of the estimator to take data from
   void merge(const VisMetaDataStats &other);
   
   /// @brief process one accessor of data updating statistics
   /// @details 
   /// @param[in] acc read-only accessor with data
   void process(const accessors::IConstDataAccessor &acc);
   
   // access to the data
   
   /// @brief total number of visibility points processed
   /// @details This method counts all visibility points. One spectral channel is one
   /// visibility point (but polarisations are not counted separately).
   /// @return number of visibility points processed so far
   inline unsigned long nVis() const { return itsNVis;}
         
   /// @brief longest baseline spacing in wavelengths
   /// @return largest absolute value of u in wavelengths
   inline double maxU() const { return itsMaxU; }
   
   /// @brief longest baseline spacing in wavelengths
   /// @return largest absolute value of v in wavelengths
   inline double maxV() const { return itsMaxV; }
      
   /// @brief largest w-term without snap-shotting 
   /// @return largest absolute value of w in wavelengths
   /// @note If the class has been initialised with the default constructor and 
   /// no tangent point is set, then this field returns the largest w-term before
   /// the uvw-rotation.
   inline double maxW() const { return itsMaxW; }
   
   /// @brief largest residual w-term (for snap-shotting)
   /// @return largest value of residual w in wavelengths
   double maxResidualW() const;
   
   /// @brief smallest frequency (same units as in the accessor, usually Hz)
   /// @return frequency
   inline double minFreq() const { return itsMinFreq; }
   
   /// @brief largest frequency (same units as in the accessor, usually Hz)
   /// @return frequency
   inline double maxFreq() const { return itsMaxFreq; }
          
   /// @brief number of antennas
   /// @return largest encountered antenna index + 1
   inline casa::uInt nAntennas() const { return itsMaxAntennaIndex + 1; } 
   
   /// @brief number of beams
   /// @return largest encountered beam index + 1
   inline casa::uInt nBeams() const { return itsMaxBeamIndex + 1; }
   
   /// @brief most central direction of the observed field
   /// @return direction of the centre in the frame used by the accessor
   casa::MVDirection centre() const;
   
   /// @brief largest separation of individual pointing from the centre
   /// @return largest offsets from the centre() in radians (measure of the field size)
   std::pair<double,double> maxOffsets() const;

   // derived statistics
   
   /// @brief estimate of the field size
   /// @details This method uses maxOffsets, centre and itsTangent to estimate the field size applying
   /// current knowledge on the guard band around the edge pointing (hard coded for ASKAP).
   /// @param[in] forceCentreAtTangent if true, the tangent point is assumed to be in the image centre
   /// @return square field size in degrees
   double squareFieldSize(bool forceCentreAtTangent = false) const;
   
   /// @brief estimate cell size
   /// @details This method uses maxU and maxV to estimate the largest (square) image cell size in arcsec
   /// @return square cell size in arcsec
   double squareCellSize() const;
   
   // serialization 
   
   /// @brief write the object to a blob stream
   /// @param[in] os the output stream
   virtual void writeToBlob(LOFAR::BlobOStream& os) const;

   /// @brief read the object from a blob stream
   /// @param[in] is the input stream
   virtual void readFromBlob(LOFAR::BlobIStream& is); 
   
protected:

   /// @brief helper method to apply an offset to the current reference direction
   /// @details
   /// @param[in] offsets pair of offsets to apply
   /// @return direction measure
   casa::MVDirection getOffsetDir(const std::pair<double,double> &offsets) const;
   
   /// @brief helper method to compute offsets of the given direction w.r.t. the reference direction
   /// @details
   /// @param[in] dir direction measure
   /// @return pair with offsets w.r.t. the reference direction
   std::pair<double,double> getOffsets(const casa::MVDirection &dir) const;
   
private:
   /// @brief tangent point for imaging
   /// @details Is not initialised, if itsTangentSet is false
   casa::MVDirection itsTangent;
   
   /// @brief flag that itsTangent is initialised
   bool itsTangentSet;      
   
   /// @brief adapter dealing with plane fitting
   /// @note This adapter is only used when w-tolerance and tangent points are set. Otherwise,
   /// we set tolerance to a negative value used as a flag to work with the original accessor.
   accessors::BestWPlaneDataAccessor itsAccessorAdapter;
   
   /// @brief number of visibilities processed
   unsigned long itsNVis;
   
   /// @brief largest absolute value of u
   double itsMaxU;
   
   /// @brief largest absolute value of v
   double itsMaxV;
   
   /// @brief largest absolute value of w
   /// @note It can be both before and after uvw-rotation depending on the value of itsTangentSet
   double itsMaxW;
   
   /// @brief largest value of residual w for snap-shot imaging
   double itsMaxResidualW;
   
   /// @brief smallest frequency (same units as in the accessor)
   double itsMinFreq;
   
   /// @brief largest frequency (same units as in the accessor)
   double itsMaxFreq;
   
   /// @brief largest antenna index
   casa::uInt itsMaxAntennaIndex;
   
   /// @brief largest beam index
   casa::uInt itsMaxBeamIndex;
   
   // data members related to direction and field size estimates
   
   /// @brief reference direction serving as the origin for the offsets
   /// @details This is initialised to either the tangent point or to the first encountered 
   /// phase centre. For all other phase centres, the offsets w.r.t. this reference
   /// position are calculated and the largest (in both directions) are found. These are
   /// used to estimate the centre and the field size. 
   casa::MVDirection itsReferenceDir;
   
   /// @brief true, if the reference direction has been initialised
   bool itsRefDirValid; 
   
   /// @brief offsets of the "bottom left" corner w.r.t. reference direction (in radians)
   /// @details True angle extreme offsets are calculated during the iteration. These offsets
   /// are used to estimate the centre and the field size.
   std::pair<double,double> itsFieldBLC;
   
   /// @brief offsets of the "top right" corner w.r.t. reference direction (in radians)
   /// @details True angle extreme offsets are calculated during the iteration. These offsets
   /// are used to estimate the centre and the field size.
   std::pair<double,double> itsFieldTRC;
   
   
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_VIS_METADATA_STATS_H

