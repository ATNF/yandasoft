/// @file
/// @brief A structural class for solvers doing cleaning
/// @details There are specific parameters for solvers doing cleaning via LatticeCleaner.
/// It seems that at this stage most of these specialized parameters are handled by
/// the scimath::Solveable class (is it a fat interface?), however a fractional threshold
/// required by just multiscale solver and MSMF solver is defined here.
///
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

#ifndef IMAGE_CLEANING_SOLVER_H
#define IMAGE_CLEANING_SOLVER_H

#include <askap/measurementequation/ImageSolver.h>

namespace askap {

namespace synthesis {

/// @brief A structural class for solvers doing cleaning
/// @details There are specific parameters for solvers doing cleaning via LatticeCleaner.
/// It seems that at this stage most of these specialized parameters are handled by
/// the scimath::Solveable class (is it a fat interface?), however a fractional threshold
/// required by just multiscale solver and MSMF solver is defined here.
/// @ingroup measurementequation
class ImageCleaningSolver : public ImageSolver {
public:
   /// @brief default constructor
   ImageCleaningSolver();

   /// @brief access to a fractional threshold
   /// @return current fractional threshold
   double fractionalThreshold() const;

   /// @brief set a new fractional threshold
   /// @param[in] fThreshold new fractional threshold
   /// @note Assign 0. to switch this option off.
   void setFractionalThreshold(double fThreshold);

   /// @brief access to a deep threshold
   /// @return current deep threshold
   double deepThreshold() const;

   /// @brief set a new deep threshold
   /// @param[in] dThreshold new deep threshold
   /// @note Assign 0. to switch this option off.
   void setDeepThreshold(double dThreshold);

   /// @brief access to a noise threshold
   /// @return current noise threshold
   double noiseThreshold() const
   { return itsNoiseThreshold; }

   /// @brief set a new noise threshold
   /// @param[in] sThreshold new noise threshold
   /// @note Assign 0. to switch this option off.
   void setNoiseThreshold(double sThreshold)
   { itsNoiseThreshold = sThreshold; }

   /// @brief access to a deep noise threshold
   /// @return current deep noise threshold
   double deepNoiseThreshold() const
   { return itsDeepNoiseThreshold; }

   /// @brief set a new deep noise threshold
   /// @param[in] dsThreshold new deep noise threshold
   /// @note Assign 0. to switch this option off.
   void setDeepNoiseThreshold(double dsThreshold)
   { itsDeepNoiseThreshold = dsThreshold; }

   /// @brief access to a masking threshold
   /// @return current masking threshold
   double maskingThreshold() const;

   /// @brief set a new masking threshold
   /// @param[in] mThreshold new masking threshold
   /// @note Assign -1. or any negative number to revert to a default behavior of the
   /// S/N based cleaning. The masking threshold value, which used to be hardcoded in
   /// the casacore when signal-based cleaning was the only available option, equals to 0.9.
   void setMaskingThreshold(double mThreshold);

   /// @brief how to apply thresholds for multiple images
   /// @details When cleaning multiple images, especially a main field and some offset fields
   /// you can avoid overcleaning the offset fields by treating the clean thresholds as if it
   /// is one combined image
   /// @return current setting - true = use first image, false = independent
   bool useFirstImageForThresholds() const {return itsUseFirstImageForThresholds;}

   /// @brief set how to apply thresholds for multiple images
   /// @details When cleaning multiple images, especially a main field and some offset fields
   /// you can avoid overcleaning the offset fields by treating the clean thresholds as if it
   /// is one combined image
   /// @param[in] useFirst : true = use first image, false = independent
   void setUseFirstImageForThresholds(bool useFirst)
   { itsUseFirstImageForThresholds = useFirst;}

   /// @brief use a noise map for a spatially varying noise cutoff?
   /// @details When there are sources with dynamic range issues in the image
   /// a noise map can keep clean from picking up sidelobes around strong sources while
   /// still cleaning deeply in other areas of the image. If non-zero, use the specified
   /// box size for the noise calculation
   casacore::uInt noiseBoxSize() const { return itsNoiseBoxSize;}

   /// @brief set the box size for the noise map for a spatially varying noise cutoff
   /// @details When there are sources with dynamic range issues in the image
   /// a noise map can keep clean from picking up sidelobes around strong sources while
   /// still cleaning deeply in other areas of the image. 
   /// @param[in] noiseBoxSize : 0 = no noise map, >0 : specifies box size
   void setNoiseBoxSize(casacore::uInt noiseBoxSize) { itsNoiseBoxSize = noiseBoxSize;}

   /// @brief set padding factor for this solver
   /// @details Because cleaning usually utilises FFT for performance (to calculate convolution with PSF),
   /// some padding is necessary to avoid wrap around. This parameter controlls the amount of padding.
   /// 1.0 means no padding, the value should be greater than or equal to 1.0.
   /// @param[in] padding padding factor
   void setPaddingFactor(float padding = 1.0);

   /// @brief obtain padding factor
   /// @return current padding factor
   inline float paddingFactor() const {return itsPaddingFactor;}

    /// @brief configure basic parameters of the solver
    /// @details This method encapsulates extraction of basic solver parameters from the parset.
    /// @param[in] parset parset's subset (should have solver.Clean or solver.Dirty removed)
    virtual void configure(const LOFAR::ParameterSet &parset);

    void setInverseCouplingMatrix(casa::Matrix<casa::Double> &InverseMatrix);
    casa::Matrix<casa::Double> getInverseCouplingMatrix();


protected:
   /// @brief helper method to pad an image
   /// @details This method encapsulates all padding logic. In addition double to float conversion happens
   /// here.
   /// @param[in] image input image (to be padded, with double precision at the moment)
   /// @return padded image converted to floats
   casacore::Array<float> padImage(const casacore::Array<imtype> &image) const;

   /// @brief helper method to pad an image
   /// @details This method encapsulates all padding logic. In addition double to float conversion happens
   /// here.
   /// @param[in] image input padded image (with single precision at the moment)
   /// @return image of original (unpadded) shape converted to double precision
   casacore::Array<imtype> unpadImage(const casacore::Array<float> &image) const;

   /// @brief helper method to clip the edges of padded image
   /// @details This method fills the edges of a padded image with 0 (original subimage is left intact, so
   /// unpadImage would return the same result before and after this method). This operation is required
   /// after non-linear transformation of an image in the other domain (i.e. some types of preconditioning).
   /// @param[in] img input padded image to be clipped
   void clipImage(casacore::Array<float> &img) const;

   /// @brief helper method to pad diagonal
   /// @details The difference from padImage is that we don't need double to float conversion for diagonal and
   /// the output array is flattened into a 1D vector
   /// @param[in] diag diagonal array
   /// @return flattened padded vector
   casacore::Vector<imtype> padDiagonal(const casacore::Array<imtype> &diag) const;


private:
  /// @brief Fractional cleaning threshold
  /// @details it means it is defined with respect to the peak flux (i.e. 0.2 from the peak flux).
  /// Assign zero (default) if you don't want any fractional threshold applied.
  double itsFractionalThreshold;

  /// @brief Deep cleaning threshold
  /// @details Threshold for a second level of deep cleaning with model as mask
  /// Assign zero (default) if you don't want any deep threshold applied.
  double itsDeepThreshold;

  /// @brief noise cleaning threshold
  /// @details Threshold defined as multiple of residual rms
  /// Assign zero (default) if you don't want any noise threshold applied.
  double itsNoiseThreshold;

  /// @brief Deep noise cleaning threshold
  /// @details noise threshold for a second level of deep cleaning with model as mask
  /// Assign zero (default) if you don't want any deep threshold applied.
  double itsDeepNoiseThreshold;

  /// @brief specify the noise box size
  /// @details specify the box size to use when calculating the 2d noise map, 
  /// default (=0) is no noise map, just a single value for the image
  casacore::uInt itsNoiseBoxSize;

   /// @brief threshold for cleaning mask
   /// @details This value is passed in the setMask call to the LatticeCleaner. Negative
   /// value (default) means that no thresholding is used to the mask and the mask array
   /// is treated as a weight providing a capability to do S/N-based cleaning. This is the
   /// default approach. If a positive value is assigned, it is treated as a threshold.
   /// If the mask value is lower than this threshold, the corresponding pixel is not cleaned.
   /// This is a classical signal-based cleaning.
   double itsMaskingThreshold;

   /// @brief how to apply thresholds for multiple images
   bool itsUseFirstImageForThresholds;

   /// @brief padding factor for solver
   /// @details Because cleaning usually utilises FFT for performance (to calculate convolution with PSF),
   /// some padding is necessary to avoid wrap around. This parameter controlls the amount of padding.
   /// 1.0 means no padding, the value should be greater than or equal to 1.0.
   float itsPaddingFactor;

   //static casa::Matrix<casa::Double> itsInverseCouplingMatrixCache;
};

} // namespace synthesis

} // namespace askap


#endif // #ifndef IMAGE_CLEANING_SOLVER_H
