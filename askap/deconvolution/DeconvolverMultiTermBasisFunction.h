/// @file DeconvolverMultiTermBasisFunction.h
/// @brief Class for a BasisFunction-Clean-based deconvolver
/// @details This interface class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
/// This version can deal with multiple terms
/// D = B(0)*I(0) + B(1)*I(1) + B(2)*I(2)
/// The most common example is where the B's are the spectral dirty beams
/// and the I(0) are the Taylor series approximation to the frequency
/// dependent frequencies.
/// @ingroup Deconvolver
///
///
/// @copyright (c) 2007,2024 CSIRO
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
/// @author Mark Wieringa <Mark.Wieringa@csiro.au>
///

#ifndef ASKAP_SYNTHESIS_DECONVOLVERMULTITERMBASISFUNCTION_H
#define ASKAP_SYNTHESIS_DECONVOLVERMULTITERMBASISFUNCTION_H

#include <string>

#include <boost/shared_ptr.hpp>

#include <askap/deconvolution/DeconvolverBase.h>
#include <askap/deconvolution/BasisFunction.h>
#include <askap/utils/DeconvolveTimerUtils.h>


namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver using the BasisFunction Clean algorithm
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a dirty image, psf optionally using a mask and a weights image.
        /// This algorithm is similar to the MultiScale Clean (Cornwell 2009) with changes
        /// to improve performance and flexibility.
        ///
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverBasisFunction<double, DComplex>
        /// @ingroup Deconvolver
        template<class T, class FT>
        class DeconvolverMultiTermBasisFunction : public DeconvolverBase<T, FT> {

            public:

                typedef boost::shared_ptr<DeconvolverMultiTermBasisFunction<T, FT>> ShPtr;

                /// @brief Construct from dirty image and psf
                /// @details Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We use a vector
                /// here to allow multiple dirty images and PSFs for the
                /// same model (e.g. as in MFS)
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                /// @param[in] psf Point Spread Function containing 2*nTaylor-1 terms (array)
                DeconvolverMultiTermBasisFunction(Vector<Array<T>>& dirty,
                                                  Vector<Array<T>>& psf,
                                                  Vector<Array<T>>& psfLong);

                /// @brief Construct from dirty image and psf
                /// @details Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We keep this
                /// version for compatibility
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverMultiTermBasisFunction(Array<T>& dirty, Array<T>& psf);

                virtual ~DeconvolverMultiTermBasisFunction();

                /// @brief Set the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction.
                /// @param[in] bf Shared pointer to basisfunction instance
                void setBasisFunction(boost::shared_ptr<BasisFunction<T>> bf);

                /// @brief Return the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction
                boost::shared_ptr<BasisFunction<T>> basisFunction();

                /// @brief Set the type of solution used in finding the optimum component
                /// @details When trying to find the optimum component we can use MAXBASE
                /// to find the peak over term 0 bases or MAXCHISQ to also use the higher 
                /// order taylor terms in the optimization.
                /// @params[in] solutionType, specify either MAXBASE or MAXCHISQ
                void setSolutionType(const std::string& solutionType);

                /// @brief Perform the minor cycle clean iterations
                /// @details This is where the actual clean is implemented, uses OpenMP to speed it up
                void ManyIterations();

                /// @brief Perform the deconvolution
                /// @details This is the main deconvolution method.
                virtual bool deconvolve();

                /// @brief Initialize the deconvolution
                /// @details Initialise e.g., get residuals for each base and set weighted mask
                virtual void initialise();

                /// @brief Finalise the deconvolution
                /// @details Finalise the deconvolution by updating the residuals
                virtual void finalise();

                /// @brief configure basic parameters of the solver
                /// @details This method encapsulates extraction of basic solver parameters from the parset.
                /// @param[in] parset parset
                virtual void configure(const LOFAR::ParameterSet &parset);

                /// @brief Update only the dirty image
                /// @details Update an existing deconvolver for a changed dirty image
                /// @param[in] dirty Dirty image (array)
                /// @param[in] term term to update
                virtual void updateDirty(Array<T>& dirty, uInt term = 0);

                /// @brief Update only the dirty images
                /// @details Update an existing deconvolver for a changed dirty images.
                /// @param[in] dirty Dirty image (vector of arrays)
                virtual void updateDirty(Vector<Array<T>>& dirty);

                /// @brief export the scale mask
                /// @details Give access to the scale mask used in the deconvolution
                /// @return Matrix<T> scale mask (bitmask of scales for each pixel)
                const Matrix<T> scaleMask();

                /// @brief import initial scale mask
                /// @details Load an initial scale mask to use in the deconvolution
                /// @param[in]scaleMask a Matrix<T> with bitmask of scales for each pixel
                void setScaleMask(const Matrix<T>& scaleMask);

                /// @brief set noise map to spatially variable noise thresholds
                /// @details to allow spatially variant threshold we need a map of how
                /// the noise varies across the image. The noise map should be normalised
                /// to the average noise
                /// @param[in] noiseMap a Matrix<T> with normalised noise across the image
                /// @param[in] noiseBoxSize box size used to calculate the noise map
                void setNoiseMap(const Matrix<T>& noiseMap, uInt noiseBoxSize)
                { itsNoiseMap = noiseMap; itsNoiseBoxSize = noiseBoxSize;}

                /// @brief release working memory not needed between major cycles
                /// @details Deconvolvers can use a lot of memory, try to release as
                /// much as possible without losing essential state
                virtual void releaseMemory();

            private:

                /// @brief Initialise the PSFs
                /// @details only need to do this once per change in basis functions
                void initialisePSF();

                /// @brief Initialise the residual bases if required
                /// @details Recalculates the residual images for each basis
                void initialiseResidual();

                /// @brief Initialise the deep clean masks if required
                /// @details In practice masks are no longer used, instead a list of active pixels is kept for each base
                void initialiseMask();

                /// @brief Initialise the basisfunctions
                /// @details Calculate the basisfunction & PSFs at size needed
                /// @param[in] force, if true, force an update even if basisfunction hasn't changed
                virtual void initialiseForBasisFunction(bool force);

                /// @brief Find the next peak for clean
                /// @details Uses the specified algorithm to find the next optimal peak to subtract
                /// @param[out] optimumBase the base with the optimal peak
                /// @param[out] absPeakPos IPosition giving location of peak
                /// @param[out] absPeakVal Value of the peak
                /// @param[out] absPeakValScale Value of the peak, scaled with local noise noise estimate
                /// @param[in]  firstCycle bool value indicating we are in the first major cycle
                /// @param[in,out] highPixels list of high pixels that could be needed for each base
                /// @param[in,out] sectionTimer timer used to keep track of processing time in sections of the clean
                /// @param[in,out] maxPos only used as shared variable across threads
                /// @param[in,out] maxVal only used as shared variable across threads
                /// @param[in,out] maxValScaled only used as shared variable across threads
                /// @param[in] weights array with the data weights and/or mask, pixel values are multiplied by this before peak is determined,
                ///  allowed to be empty, in which case no weighting is done
                /// @param[in] neqchisq, work array for MAXCHISQ peak finding, can be empty
                /// @param[in] coefficients, work arrays for MAXCHISQ peak finding, can be empty
                void chooseComponent(uInt& optimumBase, IPosition& absPeakPos, T& absPeakVal, T& absPealValScaled, bool firstCycle,
                    const std::vector<std::vector<uInt>>&highPixels, askap::utils::SectionTimer& sectionTimer,
                    IPosition& maxPos, T& maxVal, T& maxValScaled, const Matrix<T>& weights, Matrix<T>& negchisq, Vector<Matrix<T>>& coefficients);

                /// @brief Fill the vector of high (or active) pixels for each base
                /// @details Rather than working with entire images, the peak search uses a list of high pixels that
                /// may be needed - with a cutoff a bit below the clean threshold
                /// @param[out] highPixels, list of high pixels to be considered in peak searching
                /// @param[in] weight, array with the data weights and/or mask
                void fillHighPixelList(std::vector<std::vector<uInt>>& highPixels, const Matrix<T>& weight);

                /// @brief Work out the slice of the residual and psf to use
                /// @details Based on the peak position, work out the size and start of the slice
                /// in the residual image and PSF we need
                /// @param[in] absPeakPos the peak position
                /// @param[out] shape the shape of the residual and psf slice to use
                /// @param[out] residualStart start of the residual slice
                /// @param[out] psfStart start of the PSF slice
                void getResidualAndPSFSlice(const IPosition& absPeakPos,
                    IPosition& shape, IPosition& residualStart, IPosition& psfStart);

                /// @brief Add the optimum basisfunction component to the model
                /// @details Updates the model by adding the current component multiplied by its corresponding basisfunction
                /// @param[in] peakValues the value at the peak location in each of the residual bases
                /// @param[in] shape the shape of the residual and psf slice to use
                /// @param[in] residualStart start of the residual slice
                /// @param[in] psfStart start of the PSF slice
                /// @param[in] optimumBase the base with the optimal peak
                /// @param[in,out] model only used as shared variable across threads
                void addComponentToModel(const Vector<T>& peakValues,
                        const IPosition& shape, const IPosition& resStart, const IPosition& psfStart,
                        const uInt optimumBase, Matrix<T>& model);

                /// @brief Subtract the basisfunction component from all residual bases and terms
                /// @details Update the residuals by subtracting the current component multiplied by all relevant PSF crossterms
                /// @param[in] peakValues the value at the peak location in each of the residual bases
                /// @param[in] shape the shape of the residual and psf slice to use
                /// @param[in] residualStart start of the residual slice
                /// @param[in] psfStart start of the PSF slice
                /// @param[in] optimumBase the base with the optimal peak
                void subtractPSF(const Vector<T>& peakValues,
                    const IPosition& shape, const IPosition& resStart, const IPosition& psfStart,
                    uInt optimumBase);

                /// @brief Subtract the basisfunction component from all residual bases and terms
                /// @details Update the residuals by subtracting the current component multiplied by all relevant PSF crossterms
                /// This version uses the list of active pixels instead of the full arrays
                /// @param[in] peakValues the value at the peak location in each of the residual bases
                /// @param[in] peakPos the peak position
                /// @param[in] optimumBase the base with the optimal peak
                /// @param[in] useHighPixels specify true if using the highPixels vectors of active pixels
                /// @param[in] highPixels, list of high pixels to be considered
                void subtractPSFPixels(const Vector<T>& peakValues, const IPosition& peakPos,
                    uInt optimumBase, bool useHighPixels, const std::vector<std::vector<uInt>>& highPixels);

                /// @brief Long vector of PSFs
                /// @details Vector of PSFs, with 2*nterm-1 Taylor terms
                Vector<Matrix<T>> itsPsfLongVec;

                /// @brief Residual images convolved with basis functions
                /// @details The shape of this is nbases x nterms x [nx,ny]
                Vector<Vector<Matrix<T>>> itsResidualBasis;

                /// @brief The maximum number of scale we can deal with
                /// @details Number of scales is limited to 24 as we write the scale mask out as a float image
                /// we could use std::numeric_limits<T>::digits, but for double
                /// the uInt used doesn't have enough range - anyway 24 scales seems plenty
                const uInt itsMaxScales = 24;

                /// @brief The set of pixels used for each scale/base
                /// @details We keep a list of pixels we've already cleaned for each scale
                /// The pixels are recorded as uInt offset from the start of the image.
                /// We would need to use ulong/size_t here instead of uInt for images > 64k^2
                std::vector<std::set<uInt>> itsScalePixels;

                /// @brief Point spread functions convolved with cross terms
                /// @details The shape of these is [nbases,nbases]x[nterms,nterms]x[nxsub,nysub], where [nxsub,nysub]
                /// is the PFS sub shape (set by psfwidth parameter)
                Matrix<Matrix<Matrix<T>>> itsPSFCrossTerms;

                /// @brief The coupling between different terms for each basis
                /// @details The shape of this is a vector of [nbases] with a matrices of shape [nterms,nterms]
                Vector<Matrix<double>> itsCouplingMatrix;

                /// @brief Inverse of the coupling matrix [nterms,nterms][nbases]
                /// @details The shape of this is a vector of [nbases] with a matrices of shape [nterms,nterms]
                Vector<Matrix<double>> itsInverseCouplingMatrix;

                /// @brief Pointer to the basis function used in the deconvolution
                /// @details This will point to a MultiScaleBasisFunction with a number of scales
                boost::shared_ptr<BasisFunction<T>> itsBasisFunction;

                /// @brief The flux subtracted on each term and scale 
                /// @details The shape of this is a vector of [nbases] with vectors of [nterms] values
                Vector< Vector<T>> itsTermBaseFlux;

                /// @brief Flag to indicate the dirty / residual images have been updated
                /// @details When the dirty image is updated we need to recalculate the residuals for each base
                bool itsDirtyChanged;

                /// @brief FLag to indicate the BasisFunction has been updated
                /// @details When the basisfunction changes, we need to recalculate residuals and psf cross terms
                bool itsBasisFunctionChanged;

                /// @brief Set the type of solution used in finding the optimum component
                /// @details When trying to find the optimum component we can use MAXBASE
                /// to find the peak over term 0 bases or MAXCHISQ to also use the higher 
                /// order taylor terms in the optimization.
                std::string itsSolutionType;

                /// @brief Use a list of high / active pixels
                /// @details If true, use lists of pixels for each base instead of images to do the minor cycles
                bool itsUsePixelLists;

                /// @brief pixel list tolerance
                /// @details Multiply the clean threshold by (1-tolerance) to decide which
                /// pixels go in the pixel list for a major cycle. Suggested value 0.1
                float itsPixelListTolerance;

                /// @brief pixel list n sigma limit
                /// @details Don't put pixels in the pixellist with amplitude < limit*noise
                /// The noise is determined separately for each scale / residual basis
                float itsPixelListNSigma;

                /// @brief pixellist range of number of pixels
                /// @details Avoid putting way too many pixels in the list
                /// Try to get the number between the first and second entry times
                /// the maximum number of iterations
                std::vector<float> itsPixelListNPixRange;

                /// @brief Read a pre-existing scale mask with the name given if not empty
                /// @details Load an initial scale mask to use in the deconvolution, the mask
                /// file will contain a bitmask specifying the active scales at each location
                std::string itsScaleMaskName;

                /// @brief image (map) of the noise level across the image
                /// @details to allow a spatially variant threshold we need a map of how
                /// the noise varies across the image. The noise map should be normalised
                /// to the average noise. Note that it may be smaller than the residual image
                /// by a factor of about itsNoiseBoxSize. Also note that this noise map is really
                /// more a map of poor dynamic range areas in the image.
                Matrix<T> itsNoiseMap;

                /// @brief size of the box used to calculate the noise map
                /// @details to allow spatially variant threshold we need a map of how
                /// the noise varies across the image. There are various ways such a noise map
                /// could be calculated, but setting a window or box size to define a local region
                /// is a common requirememt
                uInt itsNoiseBoxSize;

                /// @brief Use a pixel increment for larger scales
                /// @details If true, only consider every n'th pixel when doing peak searching for
                /// larger scales, where n = 1 for the first two scales, then doubles for each
                /// subsequent scale. When using pixellists, this avoids putting lots of closely
                /// spaced, correlated pixels in the list for large scales.
                bool itsUseIncrements;
        };

    } // namespace synthesis

} // namespace askap

#include <askap/deconvolution/DeconvolverMultiTermBasisFunction.tcc>

#endif
