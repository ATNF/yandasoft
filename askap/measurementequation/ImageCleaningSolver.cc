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

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>

ASKAP_LOGGER(logger, ".measurementequation.imagecleaningsolver");

#include <askap/measurementequation/ImageCleaningSolver.h>
#include <askap/askap/AskapError.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/scimath/utils/MultiDimArrayPlaneIter.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/LatticeMath/LatticeFFT.h>


namespace askap {

namespace synthesis {

/// @brief default constructor
ImageCleaningSolver::ImageCleaningSolver() :
   itsFractionalThreshold(0.), itsDeepThreshold(0.), itsMaskingThreshold(-1.), itsPaddingFactor(1.) {}

/// @brief access to a fractional threshold
/// @return current fractional threshold
double ImageCleaningSolver::fractionalThreshold() const
{
  return itsFractionalThreshold;
}

/// @brief set a new fractional threshold
/// @param[in] fThreshold new fractional threshold
/// @note Assign 0. to switch this option off.
void ImageCleaningSolver::setFractionalThreshold(double fThreshold)
{
  itsFractionalThreshold = fThreshold;
}

/// @brief access to the deep threshold
/// @return current deep threshold
double ImageCleaningSolver::deepThreshold() const
{
  return itsDeepThreshold;
}

/// @brief set a new deep threshold
/// @param[in] dThreshold new deep threshold
/// @note Assign 0. to switch this option off.
void ImageCleaningSolver::setDeepThreshold(double dThreshold)
{
  itsDeepThreshold = dThreshold;
}

/// @brief access to a masking threshold
/// @return current masking threshold
double ImageCleaningSolver::maskingThreshold() const
{
  return itsMaskingThreshold;
}

/// @brief set a new masking threshold
/// @param[in] mThreshold new masking threshold
/// @note Assign -1. or any negative number to revert to a default behavior of the
/// S/N based cleaning. The masking threshold value, which used to be hardcoded in
/// the casacore when signal-based cleaning was the only available option, equals to 0.9.
void ImageCleaningSolver::setMaskingThreshold(double mThreshold)
{
  itsMaskingThreshold = mThreshold;
}

/// @brief set padding factor for this solver
/// @details Because cleaning usually utilises FFT for performance (to calculate convolution with PSF),
/// some padding is necessary to avoid wrap around. This parameter controlls the amount of padding.
/// 1.0 means no padding, the value should be greater than or equal to 1.0.
/// @param[in] padding padding factor
void ImageCleaningSolver::setPaddingFactor(float padding)
{
  ASKAPCHECK(padding>=1.0, "Padding in the solver is supposed to be greater than or equal to 1.0, you have "<<padding);
  itsPaddingFactor = padding;
}

/// @brief helper method to pad an image
/// @details This method encapsulates all padding logic. In addition double to float conversion happens
/// here.
/// @param[in] image input image (to be padded, with double precision at the moment)
/// @return padded image converted to floats
casacore::Array<float> ImageCleaningSolver::padImage(const casacore::Array<imtype> &image) const
{
  casacore::Array<float> result(scimath::PaddingUtils::paddedShape(image.shape(),paddingFactor()),0.);
  casacore::Array<float> subImage = scimath::PaddingUtils::extract(result,paddingFactor());
  #ifdef ASKAP_FLOAT_IMAGE_PARAMS
  subImage = image;
  #else
  casacore::convertArray<float, double>(subImage, image);
  #endif
  return result;
}

/// @brief helper method to clip the edges of padded image
/// @details This method fills the edges of a padded image with 0 (original subimage is left intact, so
/// unpadImage would return the same result before and after this method). This operation is required
/// after non-linear transformation of an image in the other domain (i.e. some types of preconditioning).
/// @param[in] img input padded image to be clipped
void ImageCleaningSolver::clipImage(casacore::Array<float> &img) const
{
  const casacore::IPosition origShape = scimath::PaddingUtils::unpadShape(img.shape(),paddingFactor());
  if (origShape != img.shape()) {
      scimath::PaddingUtils::clip(img,origShape);
  }
}


/// @brief helper method to pad diagonal
/// @details The difference from padImage is that we don't need double to float conversion for diagonal and
/// the output array is flattened into a 1D vector
/// @param[in] diag diagonal array
/// @return flattened padded vector
casacore::Vector<imtype> ImageCleaningSolver::padDiagonal(const casacore::Array<imtype> &diag) const
{
  if (scimath::PaddingUtils::paddedShape(diag.shape(),paddingFactor()) == diag.shape()) {
      return casacore::Vector<imtype>(diag.reform(casacore::IPosition(1,diag.nelements())));
  }
  casacore::Array<imtype> result(scimath::PaddingUtils::paddedShape(diag.shape(),paddingFactor()),0.);
  scimath::PaddingUtils::extract(result,paddingFactor()) = diag;
  return casacore::Vector<imtype>(result.reform(casacore::IPosition(1,result.nelements())));
}


/// @brief helper method to pad an image
/// @details This method encapsulates all padding logic. In addition double to float conversion happens
/// here.
/// @param[in] image input padded image (with single precision at the moment)
/// @return image of original (unpadded) shape converted to double precision
casacore::Array<imtype> ImageCleaningSolver::unpadImage(const casacore::Array<float> &image) const
{
  casacore::Array<float> wrapper(image);
  const casacore::Array<float> subImage = scimath::PaddingUtils::extract(wrapper,paddingFactor());
  casacore::Array<imtype> result(subImage.shape());
  #ifdef ASKAP_FLOAT_IMAGE_PARAMS
  result = subImage;
  #else
  casacore::convertArray<imtype,float>(result,subImage);
  #endif
  return result;
}


/// @brief zero-pad in the Fourier domain to increase resolution before cleaning
/// @param[in] osfactor extra oversampling factor
/// @param[in] image input image to be oversampled
/// @return oversampled image
/// @todo move osfactor to itsOsFactor to enforce consistency between oversample() & downsample()?
/// @todo use scimath::PaddingUtils::fftPad? Works with imtype rather than float so template there or here?
void ImageCleaningSolver::oversample(casacore::Array<float> &image, const float osfactor, const bool norm) const
{

    ASKAPCHECK(osfactor >= 1.0,
        "Oversampling factor in the solver is supposed to be greater than or equal to 1.0, you have "<<osfactor);
    if (osfactor == 1.) return;

    casacore::Array<casacore::Complex> AgridOS(scimath::PaddingUtils::paddedShape(image.shape(),osfactor),0.);

    // this is how it's done in SynthesisParamsHelper::saveImageParameter():
    //imagePixels.resize(scimath::PaddingUtils::paddedShape(ip.shape(name),osfactor));
    //scimath::PaddingUtils::fftPad(ip.valueF(name),imagePixels);

    // destroy Agrid before resizing image. Small memory saving.
    {
        // Set up scratch arrays and lattices for changing resolution
        casacore::Array<casacore::Complex> Agrid(image.shape());
        casacore::ArrayLattice<casacore::Complex> Lgrid(Agrid);
 
        // copy image into a complex scratch space
        casacore::convertArray<casacore::Complex,float>(Agrid, image);
 
        // renormalise based on the imminent padding
        if (norm) {
            Agrid *= static_cast<float>(osfactor*osfactor);
        } 

        // fft to uv
        casacore::LatticeFFT::cfft2d(Lgrid, casacore::True);
 
        // "extract" original Fourier grid into the central portion of the new Fourier grid
        casacore::Array<casacore::Complex> subGrid = scimath::PaddingUtils::extract(AgridOS,osfactor);
        subGrid = Agrid;

        // ifft back to image and return the real part
        casacore::ArrayLattice<casacore::Complex> LgridOS(AgridOS);
        casacore::LatticeFFT::cfft2d(LgridOS, casacore::False);
    }

    image.resize(AgridOS.shape());
    image = real(AgridOS);

}

/// @brief remove Fourier zero-padding region to re-establish original resolution after cleaning
/// @param[in] osfactor extra oversampling factor
/// @param[in] image input oversampled image
/// @return downsampled image
/// @todo move osfactor to itsOsFactor to enforce consistency between oversample() & downsample()?
/// @todo use scimath::PaddingUtils::fftPad? Downsampling may not be supported at this stage
void ImageCleaningSolver::downsample(casacore::Array<float> &image, const float osfactor) const
{

    ASKAPCHECK(osfactor >= 1.0,
        "Oversampling factor in the solver is supposed to be greater than or equal to 1.0, you have "<<osfactor);
    if (osfactor == 1.) return;

    casacore::Array<casacore::Complex> Agrid;

    // destroy AgridOS before resizing image. Small memory saving.
    {
        // Set up scratch arrays and lattices for changing resolution
        casacore::Array<casacore::Complex> AgridOS(image.shape());
        casacore::ArrayLattice<casacore::Complex> LgridOS(AgridOS);
 
        // copy image into a complex scratch space
        casacore::convertArray<casacore::Complex,float>(AgridOS, image);
 
        // fft to uv. This is what we need to degrid from, so no need to renormalise
        casacore::LatticeFFT::cfft2d(LgridOS, casacore::True);
 
        // extract the central portion of the Fourier grid
        Agrid = scimath::PaddingUtils::extract(AgridOS,osfactor);
 
        // renormalise based on the imminent padding
        //Agrid /= static_cast<float>(osfactor*osfactor);
 
        // ifft back to image and return the real part
        casacore::ArrayLattice<casacore::Complex> Lgrid(Agrid);
        casacore::LatticeFFT::cfft2d(Lgrid, casacore::False);
    }

    image.resize(Agrid.shape());
    image = real(Agrid);

}


/// @brief configure basic parameters of the solver
/// @details This method encapsulates extraction of basic solver parameters from the parset.
/// @param[in] parset parset's subset (should have solver.Clean or solver.Dirty removed)
void ImageCleaningSolver::configure(const LOFAR::ParameterSet &parset)
{
  ImageSolver::configure(parset);
  setPaddingFactor(parset.getFloat("padding",1.));
  ASKAPLOG_INFO_STR(logger,"Solver padding of "<<paddingFactor()<<" will be used");
}

void ImageCleaningSolver::setInverseCouplingMatrix(casa::Matrix<casa::Double> &InverseMatrix) {
    ImageSolver::setInverseCouplingMatrix(InverseMatrix);
}

casa::Matrix<casa::Double> ImageCleaningSolver::getInverseCouplingMatrix() {
    return ImageSolver::getInverseCouplingMatrix();
}

} // namespace synthesis

} // namespace askap
