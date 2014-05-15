/// @file DeconvolverBasisFunction.tcc
/// @brief Class for a deconvolver based on CLEANing with basis functions.
/// @details This concrete class defines a deconvolver used to estimate an
/// image from a residual image, psf optionally using a weights image.
/// @ingroup Deconvolver
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <scimath/Mathematics/MatrixMathLA.h>

// Local package includes
#include <measurementequation/SynthesisParamsHelper.h>
#include <deconvolution/DeconvolverBasisFunction.h>
#include <deconvolution/MultiScaleBasisFunction.h>

ASKAP_LOGGER(decbflogger, ".deconvolution.basisfunction");

namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver based on the BasisFunction Clean
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a residual image, psf optionally using a weights image.
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverBasisFunction<Double, DComplex>
        /// @ingroup Deconvolver
 
        template<class T, class FT>
        DeconvolverBasisFunction<T, FT>::DeconvolverBasisFunction(Vector<Array<T> >& dirty,
                                                                  Vector<Array<T> >& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf),
                itsUseCrossTerms(true), itsDecouple(true),
                itsDecouplingAlgorithm("diagonal")
        {
        };

        template<class T, class FT>
        DeconvolverBasisFunction<T, FT>::DeconvolverBasisFunction(Array<T>& dirty,
                                                                  Array<T>& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf),
                itsUseCrossTerms(true), itsDecouple(true),
                itsDecouplingAlgorithm("diagonal")
        {
        };

        template<class T, class FT>
        DeconvolverBasisFunction<T, FT>::~DeconvolverBasisFunction()
        {
        };

        template<class T, class FT>
        void DeconvolverBasisFunction<T, FT>::setBasisFunction(boost::shared_ptr<BasisFunction<T> > bf)
        {
            itsBasisFunction = bf;
        };

        template<class T, class FT>
        boost::shared_ptr<BasisFunction<T> > DeconvolverBasisFunction<T, FT>::basisFunction()
        {
            return itsBasisFunction;
        };

        template<class T, class FT>
        void DeconvolverBasisFunction<T, FT>::configure(const LOFAR::ParameterSet& parset)
        {
            DeconvolverBase<T, FT>::configure(parset);

            // Make the basis function
            {
                std::vector<float> defaultScales(3);
                defaultScales[0] = 0.0;
                defaultScales[1] = 10.0;
                defaultScales[2] = 30.0;
                const std::vector<float> scales = parset.getFloatVector("scales", defaultScales);

                ASKAPLOG_INFO_STR(decbflogger, "Constructing Multiscale basis function with scales " << scales);
                const Bool orthogonal = parset.getBool("orthogonal", false);

                if (orthogonal) {
                    ASKAPLOG_DEBUG_STR(decbflogger, "Multiscale basis functions will be orthogonalised");
                }

                itsBasisFunction = BasisFunction<Float>::ShPtr(new MultiScaleBasisFunction<Float>(scales,
                                   orthogonal));
            }
            itsUseCrossTerms = parset.getBool("usecrossterms", true);

            if (itsUseCrossTerms) {
                ASKAPLOG_DEBUG_STR(decbflogger, "Will use crossterms in subtraction");
            }

            itsDecouplingAlgorithm = parset.getString("decouplingalgorithm", "diagonal");
        }

        template<class T, class FT>
        void DeconvolverBasisFunction<T, FT>::finalise()
        {
            this->updateResiduals(this->itsModel);

            const Array<T> ones(this->itsL1image(0).shape(), 1.0);
            const T l0Norm(sum(ones(abs(this->itsL1image(0)) > T(0.0))));
            const T l1Norm(sum(abs(this->itsL1image(0))));
            ASKAPLOG_INFO_STR(decbflogger, "L0 norm = " << l0Norm << ", L1 norm   = " << l1Norm
                                  << ", Flux = " << sum(this->model()));

            for (uInt scale = 0; scale < itsScaleFlux.nelements(); scale++) {
                ASKAPLOG_INFO_STR(decbflogger, "   Scale " << scale << " Flux = " << itsScaleFlux(scale));
            }

        }

        template<class T, class FT>
        void DeconvolverBasisFunction<T, FT>::initialise()
        {
            DeconvolverBase<T, FT>::initialise();

            ASKAPLOG_INFO_STR(decbflogger, "Initialising Basis Function deconvolver");

            Int psfWidth = this->model().shape()(0);
            IPosition subPsfShape(2, 0, 0);

            // Only use the specified psfWidth if it makes sense
            if ((this->control()->psfWidth() > 0) && (this->control()->psfWidth() < psfWidth)) {
                psfWidth = this->control()->psfWidth();
                ASKAPLOG_INFO_STR(decbflogger, "Using subregion of Psf : size " << psfWidth
                                      << " pixels");
                subPsfShape = IPosition(2, psfWidth, psfWidth);
            } else {
                subPsfShape = IPosition(2, this->model().shape()(0), this->model().shape()(1));
            }

            this->itsBasisFunction->initialise(this->model().shape());
            initialiseResidual();
            this->itsBasisFunction->initialise(subPsfShape);
            initialisePSF();

            if (this->itsDecouplingAlgorithm == "basis") {
                // Decoupling using inverse coupling matrix generate orthogonal basis functions
                ASKAPLOG_INFO_STR(decbflogger, "Decoupling using inverse coupling matrix generate orthogonal basis functions");
                const Matrix<Double> inverseCouplingMatrix(this->itsInverseCouplingMatrix.copy());
                this->itsBasisFunction->initialise(this->model().shape());
                itsBasisFunction->multiplyArray(inverseCouplingMatrix);
                initialiseResidual();
                this->itsBasisFunction->initialise(subPsfShape);
                itsBasisFunction->multiplyArray(inverseCouplingMatrix);
                this->itsBasisFunction->multiplyArray(inverseCouplingMatrix);
                initialisePSF();
                //  SynthesisParamsHelper::saveAsCasaImage("BasisFunctionAfterInverseDecoupling.tab",
                //                         this->itsBasisFunction->basisFunction());
                //  SynthesisParamsHelper::saveAsCasaImage("ResidualsAfterInverseDecoupling.tab",
                //                         this->itsResidualBasisFunction);
            } else if (this->itsDecouplingAlgorithm == "residuals") {
                // Decoupling using inverse coupling matrix applied to basis and residuals
                ASKAPLOG_INFO_STR(decbflogger, "Decoupling using inverse coupling matrix applied to basis and residuals");
                const Array<T> invBF(applyInverse(this->itsInverseCouplingMatrix, this->itsBasisFunction->basisFunction()));
                this->itsBasisFunction->basisFunction() = invBF.copy();

                const Array<T> invRes(applyInverse(this->itsInverseCouplingMatrix, this->itsResidualBasisFunction));
                this->itsResidualBasisFunction = invRes.copy();

                if (itsUseCrossTerms) {
                    ASKAPLOG_DEBUG_STR(decbflogger, "Overriding usecrossterms since it makes no sense in this case");
                    this->itsUseCrossTerms = false;
                }

                //  SynthesisParamsHelper::saveAsCasaImage("BasisFunctionAfterResidualsDecoupling.tab",
                //                         this->itsBasisFunction->basisFunction());
                //  SynthesisParamsHelper::saveAsCasaImage("ResidualsAfterResidualsDecoupling.tab",
                //                         this->itsResidualBasisFunction);
            } else if (this->itsDecouplingAlgorithm == "inverse") {
                // Correcting coupling at subtraction phase with inverse coupling matrix
                ASKAPLOG_DEBUG_STR(decbflogger, "Correcting coupling at subtraction phase with inverse coupling matrix");
            } else if (this->itsDecouplingAlgorithm == "sqrtdiagonal") {
                // Correcting coupling at subtraction phase with inverse diag(coupling matrix)
                ASKAPLOG_DEBUG_STR(decbflogger, "Correcting coupling at subtraction phase with inverse sqrt(diag(coupling matrix))");
            } else if (this->itsDecouplingAlgorithm == "diagonal") {
                // Correcting coupling at subtraction phase with inverse diag(coupling matrix)
                ASKAPLOG_DEBUG_STR(decbflogger, "Correcting coupling at subtraction phase with inverse diag(coupling matrix)");
            } else if (this->itsDecouplingAlgorithm == "psfscales") {
                // Correcting coupling at subtraction phase with inverse psfscales
                ASKAPLOG_DEBUG_STR(decbflogger, "Correcting coupling at subtraction phase with inverse psfscales");
            } else if (this->itsDecouplingAlgorithm == "sqrtpsfscales") {
                // Correcting coupling at subtraction phase with inverse psfscales
                ASKAPLOG_DEBUG_STR(decbflogger, "Correcting coupling at subtraction phase with inverse sqrt(psfscales)");
            } else {
                // Correcting coupling at subtraction phase with inverse diag(coupling matrix)
                ASKAPLOG_DEBUG_STR(decbflogger, "Correcting coupling at subtraction phase with inverse diag(coupling matrix)");
            }

            const uInt nScales(this->itsBasisFunction->numberBases());
            const IPosition l1Shape(3, this->model().shape()(0), this->model().shape()(1), nScales);

            this->itsL1image.resize(this->itsNumberTerms);
            this->itsL1image(0).resize(l1Shape);
            this->itsL1image(0).set(0.0);
        }

        template<class T, class FT>
        void DeconvolverBasisFunction<T, FT>::initialiseResidual()
        {

            ASKAPCHECK(this->itsBasisFunction, "Basis function not initialised");

            this->state()->resetInitialObjectiveFunction();

            ASKAPLOG_DEBUG_STR(decbflogger, "Calculating cache of images");

            ASKAPLOG_DEBUG_STR(decbflogger, "Shape of basis functions "
                                   << this->itsBasisFunction->basisFunction().shape());

            const IPosition stackShape(this->itsBasisFunction->basisFunction().shape());

            itsResidualBasisFunction.resize(stackShape);

            Cube<FT> basisFunctionFFT(this->itsBasisFunction->basisFunction().shape());
            casa::setReal(basisFunctionFFT, this->itsBasisFunction->basisFunction());
            scimath::fft2d(basisFunctionFFT, true);

            Array<FT> residualFFT(this->dirty().shape().nonDegenerate());
            residualFFT.set(FT(0.0));
            casa::setReal(residualFFT, this->dirty().nonDegenerate());
            scimath::fft2d(residualFFT, true);

            Array<FT> work(this->model().nonDegenerate().shape());
            ASKAPLOG_DEBUG_STR(decbflogger,
                               "Calculating convolutions of residual image with basis functions");

            for (uInt term = 0; term < this->itsBasisFunction->numberBases(); term++) {

                ASKAPASSERT(basisFunctionFFT.xyPlane(term).nonDegenerate().shape().conform(residualFFT.nonDegenerate().shape()));
                work = conj(basisFunctionFFT.xyPlane(term).nonDegenerate()) * residualFFT.nonDegenerate();
                scimath::fft2d(work, false);

                // basis function * residual
                ASKAPLOG_DEBUG_STR(decbflogger, "Basis function(" << term
                                       << ") * Residual: max = " << max(real(work))
                                       << " min = " << min(real(work)));

                Cube<T>(itsResidualBasisFunction).xyPlane(term) = real(work);

            }
        }

        template<class T, class FT>
        void DeconvolverBasisFunction<T, FT>::initialisePSF()
        {
            // For the psf convolutions, we only need a small part of the
            // basis functions so we recalculate for that size
            Int psfWidth = this->model(0).shape()(0);

            // Only use the specified psfWidth if it makes sense
            if ((this->control()->psfWidth() > 0) && (this->control()->psfWidth() < psfWidth)) {
                psfWidth = this->control()->psfWidth();
                ASKAPLOG_DEBUG_STR(decbflogger, "Using subregion of Psf : size " << psfWidth
                                       << " pixels");
            }

            IPosition subPsfShape(2, psfWidth, psfWidth);

            Array<FT> work(subPsfShape);

            ASKAPLOG_DEBUG_STR(decbflogger, "Shape of basis functions "
                                   << this->itsBasisFunction->basisFunction().shape());

            const IPosition stackShape(this->itsBasisFunction->basisFunction().shape());

            // Now transform the basis functions
            Cube<FT> basisFunctionFFT(this->itsBasisFunction->basisFunction().shape());
            casa::setReal(basisFunctionFFT, this->itsBasisFunction->basisFunction());
            scimath::fft2d(basisFunctionFFT, true);

            this->itsPSFBasisFunction.resize(stackShape);

            this->itsScaleFlux.resize(stackShape(2));
            this->itsScaleFlux.set(T(0));

            // Calculate XFR for the subsection only
            Array<FT> subXFR(subPsfShape);

            const uInt nx(this->psf().shape()(0));
            const uInt ny(this->psf().shape()(1));

            const IPosition subPsfStart(2, nx / 2 - psfWidth / 2, ny / 2 - psfWidth / 2);
            const IPosition subPsfEnd(2, nx / 2 + psfWidth / 2 - 1, ny / 2 + psfWidth / 2 - 1);
            const IPosition subPsfStride(2, 1, 1);

            Slicer subPsfSlicer(subPsfStart, subPsfEnd, subPsfStride, Slicer::endIsLast);
            casa::IPosition minPos;
            casa::IPosition maxPos;
            T minVal, maxVal;
            casa::minMax(minVal, maxVal, minPos, maxPos, this->psf(0).nonDegenerate()(subPsfSlicer));
            ASKAPLOG_DEBUG_STR(decbflogger, "Maximum of PSF(0) = " << maxVal << " at " << maxPos);
            ASKAPLOG_DEBUG_STR(decbflogger, "Minimum of PSF(0) = " << minVal << " at " << minPos);
            this->itsPeakPSFVal = maxVal;
            this->itsPeakPSFPos(0) = maxPos(0);
            this->itsPeakPSFPos(1) = maxPos(1);

            const IPosition subPsfPeak(2, this->itsPeakPSFPos(0), this->itsPeakPSFPos(1));
            ASKAPLOG_DEBUG_STR(decbflogger, "Peak of PSF subsection at  " << subPsfPeak);
            ASKAPLOG_DEBUG_STR(decbflogger, "Shape of PSF subsection is " << subPsfShape);

            casa::setReal(subXFR, this->psf().nonDegenerate()(subPsfSlicer));
            scimath::fft2d(subXFR, true);

            // Now we have all the ingredients to calculate the convolutions
            // of basis function with psf's, etc.
            ASKAPLOG_DEBUG_STR(decbflogger, "Calculating convolutions of Psfs with basis functions");
            itsPSFScales.resize(this->itsBasisFunction->numberBases());

            for (uInt term = 0; term < this->itsBasisFunction->numberBases(); term++) {
                // basis function * psf
                ASKAPASSERT(basisFunctionFFT.xyPlane(term).nonDegenerate().shape().conform(subXFR.shape()));
                work = conj(basisFunctionFFT.xyPlane(term).nonDegenerate()) * subXFR;
                scimath::fft2d(work, false);
                Cube<T>(this->itsPSFBasisFunction).xyPlane(term) = real(work);

                ASKAPLOG_DEBUG_STR(decbflogger, "Basis function(" << term << ") * PSF: max = " << max(real(work)) << " min = " << min(real(work)));

                itsPSFScales(term) = max(real(work));
            }

            ASKAPLOG_DEBUG_STR(decbflogger, "Calculating double convolutions of PSF with basis functions");
            const IPosition crossTermsShape(4, psfWidth, psfWidth,
                                            this->itsBasisFunction->numberBases(),
                                            this->itsBasisFunction->numberBases());
            ASKAPLOG_DEBUG_STR(decbflogger, "Shape of cross terms " << crossTermsShape);
            itsPSFCrossTerms.resize(crossTermsShape);
            IPosition crossTermsStart(4, 0);
            IPosition crossTermsEnd(crossTermsShape - 1);
            IPosition crossTermsStride(4, 1);

            Array<FT> crossTermsPSFFFT(crossTermsShape);
            crossTermsPSFFFT.set(T(0));

            for (uInt term = 0; term < this->itsBasisFunction->numberBases(); term++) {
                crossTermsStart(2) = term;
                crossTermsEnd(2) = term;

                for (uInt term1 = 0; term1 < this->itsBasisFunction->numberBases(); term1++) {
                    crossTermsStart(3) = term1;
                    crossTermsEnd(3) = term1;
                    casa::Slicer crossTermsSlicer(crossTermsStart, crossTermsEnd, crossTermsStride, Slicer::endIsLast);
                    crossTermsPSFFFT(crossTermsSlicer).nonDegenerate() =
                        basisFunctionFFT.xyPlane(term).nonDegenerate() *
                        conj(basisFunctionFFT.xyPlane(term1)).nonDegenerate() * subXFR;
                }

            }

            this->itsCouplingMatrix.resize(itsBasisFunction->numberBases(), itsBasisFunction->numberBases());
            scimath::fft2d(crossTermsPSFFFT, true);
            this->itsPSFCrossTerms = real(crossTermsPSFFFT) / T(crossTermsShape(0) * crossTermsShape(1));

            for (uInt term = 0; term < this->itsBasisFunction->numberBases(); term++) {
                crossTermsStart(2) = term;
                crossTermsEnd(2) = term;

                for (uInt term1 = 0; term1 < this->itsBasisFunction->numberBases(); term1++) {
                    crossTermsStart(3) = term1;
                    crossTermsEnd(3) = term1;
                    casa::Slicer crossTermsSlicer(crossTermsStart, crossTermsEnd, crossTermsStride, Slicer::endIsLast);
                    casa::IPosition minPos;
                    casa::IPosition maxPos;
                    T minVal, maxVal;
                    casa::minMax(minVal, maxVal, minPos, maxPos, this->itsPSFCrossTerms(crossTermsSlicer));
                    this->itsCouplingMatrix(term, term1) = Double(maxVal);
                }

                this->itsCouplingMatrix(term, term) += Double(this->control()->lambda());
            }

            ASKAPLOG_DEBUG_STR(decbflogger, "Coupling matrix " << this->itsCouplingMatrix);
            this->itsInverseCouplingMatrix.resize(this->itsCouplingMatrix.shape());
            invertSymPosDef(this->itsInverseCouplingMatrix, this->itsDetCouplingMatrix, this->itsCouplingMatrix);
            ASKAPLOG_DEBUG_STR(decbflogger, "Coupling matrix determinant " << this->itsDetCouplingMatrix);
            ASKAPLOG_DEBUG_STR(decbflogger, "Inverse coupling matrix " << this->itsInverseCouplingMatrix);
            // Checked that the inverse really is an inverse.
            Matrix<T> identity(this->itsCouplingMatrix.shape(), 0.0);
            const uInt nRows(this->itsCouplingMatrix.nrow());
            const uInt nCols(this->itsCouplingMatrix.ncolumn());

            for (uInt row = 0; row < nRows; row++) {
                for (uInt col = 0; col < nCols; col++) {
                    identity(row, col) = sum(this->itsCouplingMatrix.row(row) * this->itsInverseCouplingMatrix.column(col));
                }
            }

            ASKAPLOG_DEBUG_STR(decbflogger, "Coupling matrix * inverse " << identity);


            // Now look at coupling between adjacent scales: this works well if the
            // scales are ordered.
            for (uInt term = 0; term < this->itsBasisFunction->numberBases() - 1; term++) {
                double det = this->itsCouplingMatrix(term, term) * this->itsCouplingMatrix(term + 1, term + 1) -
                             this->itsCouplingMatrix(term, term + 1) * this->itsCouplingMatrix(term + 1, term);
                ASKAPLOG_DEBUG_STR(decbflogger, "Independence between scales " << term << " and "
                                       << term + 1 << " = " << det);
            }
        }

        template<class T, class FT>
        bool DeconvolverBasisFunction<T, FT>::deconvolve()
        {
            this->initialise();

            ASKAPLOG_INFO_STR(decbflogger, "Performing BasisFunction CLEAN for "
                                  << this->control()->targetIter() << " iterations");

            do {
                this->oneIteration();
                this->monitor()->monitor(*(this->state()));
                this->state()->incIter();
            } while (!this->control()->terminate(*(this->state())));

            ASKAPLOG_INFO_STR(decbflogger, "Performed BasisFunction CLEAN for "
                                  << this->state()->currentIter() << " iterations");

            ASKAPLOG_INFO_STR(decbflogger, this->control()->terminationString());

            this->finalise();
            return True;
        }

        // This contains the heart of the BasisFunction Clean algorithm
        // The residual image and psfs are intrinsically two dimensional
        // but are expanded by projection onto the basis functions
        template<class T, class FT>
        bool DeconvolverBasisFunction<T, FT>::oneIteration()
        {
            // Find peak in residual image cube. This cube is full sized.
            casa::IPosition minPos;
            casa::IPosition maxPos;
            T minVal(0.0), maxVal(0.0);
            // Here the weights image is used as a weight in the determination
            // of the maximum i.e. it finds the max in weight . residual. The values
            // returned are without the weight
            minMaxMaskedScales(minVal, maxVal, minPos, maxPos, this->itsResidualBasisFunction,
                               this->weight(0));
            casa::IPosition absPeakPos;

            if (abs(minVal) < abs(maxVal)) {
                absPeakPos = maxPos;
            } else {
                absPeakPos = minPos;
            }

            // Find the peak values for each scale. Set the stopping criterion
            // to be the maximum of the maxima. Here we use the weighted
            // value since that's what we are interested in.
            const uInt nScales(this->itsBasisFunction->numberBases());
            Vector<T> peakValues(nScales);
            IPosition peakPos(absPeakPos);
            // If we are using residual decoupling, we need to
            // couple the peakvalues
            // Apply the inverse of the sqrt(diagonal values) to get the peak values
            Vector<T> coupledPeakValues(nScales);

            for (uInt scale = 0; scale < nScales; scale++) {
                peakPos(2) = scale;
                coupledPeakValues(scale) = this->itsResidualBasisFunction(peakPos);
            }

            if (itsDecouplingAlgorithm == "residuals") {
                // This is a special case - the residuals are already decoupled.
                peakValues = coupledPeakValues.copy();
            } else if (itsDecouplingAlgorithm == "inverse") {
                peakValues = apply(this->itsInverseCouplingMatrix, coupledPeakValues);
            } else if (itsDecouplingAlgorithm == "diagonal") {
                for (uInt scale = 0; scale < nScales; scale++) {
                    peakPos(2) = scale;
                    peakValues(scale) = coupledPeakValues(scale)
                                        / (this->itsCouplingMatrix(scale, scale));
                }
            } else if (itsDecouplingAlgorithm == "sqrtdiagonal") {
                for (uInt scale = 0; scale < nScales; scale++) {
                    peakPos(2) = scale;
                    peakValues(scale) = coupledPeakValues(scale)
                                        / sqrt(this->itsCouplingMatrix(scale, scale));
                }
            } else if (itsDecouplingAlgorithm == "psfscales") {
                for (uInt scale = 0; scale < nScales; scale++) {
                    peakPos(2) = scale;
                    peakValues(scale) = coupledPeakValues(scale)
                                        / this->itsPSFScales(scale);
                }
            } else if (itsDecouplingAlgorithm == "sqrtpsfscales") {
                for (uInt scale = 0; scale < nScales; scale++) {
                    peakPos(2) = scale;
                    peakValues(scale) = coupledPeakValues(scale)
                                        / sqrt(this->itsPSFScales(scale));
                }
            } else {
                ASKAPTHROW(AskapError, "Unknown decoupling algorithm " << itsDecouplingAlgorithm);
            }

            uInt optimumScale(0);
            T absPeakVal(0.0);

            for (uInt scale = 0; scale < nScales; scale++) {
                if (abs(peakValues(scale)) > abs(absPeakVal)) {
                    absPeakVal = peakValues(scale);
                    optimumScale = scale;
                }
            }

            // If we decoupled by residuals we need to recouple before
            // subtracting from the residuals
            if (this->itsDecouplingAlgorithm == "residuals") {
                peakValues = apply(this->itsCouplingMatrix, peakValues);
            } else {
                // Only the peak is useful
                for (uInt scale = 0; scale < nScales; scale++) {
                    if (scale != optimumScale) {
                        peakValues(scale) = T(0.0);
                    }
                }
            }

            if (this->state()->initialObjectiveFunction() == 0.0) {
                this->state()->setInitialObjectiveFunction(abs(absPeakVal));
            }

            this->state()->setPeakResidual(abs(absPeakVal));
            this->state()->setObjectiveFunction(abs(absPeakVal));
            this->state()->setTotalFlux(sum(this->model()));

            const casa::IPosition residualShape(this->itsResidualBasisFunction.shape());
            const casa::IPosition psfShape(this->itsPSFBasisFunction.shape());

            const casa::uInt ndim(this->itsResidualBasisFunction.shape().size());

            casa::IPosition residualStart(ndim, 0), residualEnd(ndim, 0), residualStride(ndim, 1);
            casa::IPosition psfStart(ndim, 0), psfEnd(ndim, 0), psfStride(ndim, 1);
            casa::IPosition psfCrossTermsStart(ndim + 1, 0), psfCrossTermsEnd(ndim + 1, 0), psfCrossTermsStride(ndim + 1, 1);

            const casa::IPosition modelShape(this->model().shape());
            const casa::uInt modelNdim(this->model().shape().size());
            casa::IPosition modelStart(modelNdim, 0), modelEnd(modelNdim, 0), modelStride(modelNdim, 1);

            // Wrangle the start, end, and shape into consistent form. It took me
            // quite a while to figure this out (slow brain day) so it may be
            // that there are some edge cases for which it fails.

            for (uInt dim = 0; dim < 2; dim++) {
                residualStart(dim) = max(0, Int(absPeakPos(dim) - psfShape(dim) / 2));
                residualEnd(dim) = min(Int(absPeakPos(dim) + psfShape(dim) / 2 - 1), Int(residualShape(dim) - 1));
                // Now we have to deal with the PSF. Here we want to use enough of the
                // PSF to clean the residual image.
                psfStart(dim) = max(0, Int(this->itsPeakPSFPos(dim) - (absPeakPos(dim) - residualStart(dim))));
                psfEnd(dim) = min(Int(this->itsPeakPSFPos(dim) - (absPeakPos(dim) - residualEnd(dim))),
                                  Int(psfShape(dim) - 1));

                psfCrossTermsStart(dim) = psfStart(dim);
                psfCrossTermsEnd(dim) = psfEnd(dim);

                modelStart(dim) = residualStart(dim);
                modelEnd(dim) = residualEnd(dim);
            }

            casa::Slicer modelSlicer(modelStart, modelEnd, modelStride, Slicer::endIsLast);

            // Add to model
            // Note that the model is only two dimensional. We could make it three dimensional
            // and keep the model layers separate
            // We loop over all terms and ignore those with no flux
            const casa::uInt nterms(this->itsResidualBasisFunction.shape()(2));

            for (uInt term = 0; term < nterms; term++) {
                if (abs(peakValues(term)) > 0.0) {
                    psfStart(2) = psfEnd(2) = term;
                    casa::Slicer psfSlicer(psfStart, psfEnd, psfStride, Slicer::endIsLast);
                    typename casa::Array<T> modelSlice = this->model()(modelSlicer).nonDegenerate();
                    modelSlice += this->control()->gain() * peakValues(term) *
                                  this->itsBasisFunction->basisFunction()(psfSlicer).nonDegenerate();
                }
            }

            // Keep track of strengths and locations of components
            for (uInt term = 0; term < nterms; term++) {
                if (abs(peakValues(term)) > 0.0) {
                    IPosition l1PeakPos(3, absPeakPos(0), absPeakPos(1), term);
                    casa::Slicer modelSlicer(modelStart, modelEnd, modelStride, Slicer::endIsLast);
                    this->itsL1image(0)(l1PeakPos) += this->control()->gain() * abs(peakValues(term));
                    this->itsScaleFlux(term) += this->control()->gain() * peakValues(term);
                }
            }

            // Subtract PSFs
            for (uInt term = 0; term < nterms; term++) {
                if (abs(peakValues(term)) > 0.0) {
                    psfStart(2) = psfEnd(2) = term;
                    casa::Slicer psfSlicer(psfStart, psfEnd, psfStride, Slicer::endIsLast);
                    residualStart(2) = residualEnd(2) = term;
                    casa::Slicer residualSlicer(residualStart, residualEnd, residualStride, Slicer::endIsLast);
                    typename casa::Array<T> residualBFSlice = this->itsResidualBasisFunction(residualSlicer).nonDegenerate();
                    residualBFSlice -= this->control()->gain() * peakValues(term) * this->itsPSFBasisFunction(psfSlicer).nonDegenerate();
                }
            }

            if (itsUseCrossTerms) {
                for (uInt term1 = 0; term1 < nterms; term1++) {
                    if (abs(peakValues(term1)) > 0.0) {
                        for (uInt term = 0; term < nterms; term++) {
                            if (term != term1) {
                                residualStart(2) = term;
                                residualEnd(2) = term;
                                casa::Slicer residualSlicer(residualStart, residualEnd, residualStride, Slicer::endIsLast);

                                psfCrossTermsStart(2) = term1;
                                psfCrossTermsEnd(2) = term1;
                                psfCrossTermsStart(3) = term;
                                psfCrossTermsEnd(3) = term;
                                casa::Slicer psfCrossTermsSlicer(psfCrossTermsStart, psfCrossTermsEnd, psfCrossTermsStride, Slicer::endIsLast);
                                typename casa::Array<T> residualBFSlice = this->itsResidualBasisFunction(residualSlicer).nonDegenerate();
                                residualBFSlice -= this->control()->gain() * peakValues(term1) *
                                                   this->itsPSFCrossTerms(psfCrossTermsSlicer).nonDegenerate();
                            }
                        }
                    }
                }
            }

            return True;
        }

        template<class T, class FT>
        void DeconvolverBasisFunction<T, FT>::minMaxMaskedScales(T& minVal, T& maxVal,
                IPosition& minPos, IPosition& maxPos,
                const Array<T>& dataArray,
                const Array<T>& weightArray)
        {
            const Cube<T> data(dataArray);
            bool isWeighted(weightArray.shape().nonDegenerate().conform(data.xyPlane(0).shape()));

            const uInt nScales = data.shape()(2);

            Vector<T> sMaxVal(nScales);
            Vector<T> sMinVal(nScales);
            Vector<IPosition> sMinPos(nScales);
            Vector<IPosition> sMaxPos(nScales);
            {
                if (isWeighted) {
                    for (uInt scale = 0; scale < nScales; scale++) {
                        casa::minMaxMasked(sMinVal(scale), sMaxVal(scale), sMinPos(scale), sMaxPos(scale),
                                           Cube<T>(dataArray).xyPlane(scale), weightArray.nonDegenerate());
                    }
                } else {
                    for (uInt scale = 0; scale < nScales; scale++) {
                        casa::minMax(sMinVal(scale), sMaxVal(scale), sMinPos(scale), sMaxPos(scale),
                                     Cube<T>(dataArray).xyPlane(scale));
                    }
                }
            }
            minPos = IPosition(3, sMinPos(0)(0), sMinPos(0)(1), 0);
            maxPos = IPosition(3, sMaxPos(0)(0), sMaxPos(0)(1), 0);
            minVal = sMinVal(0);
            maxVal = sMaxVal(0);

            for (uInt scale = 1; scale < nScales; scale++) {
                if (sMinVal(scale) <= minVal) {
                    minVal = sMinVal(scale);
                    minPos = IPosition(3, sMinPos(scale)(0), sMinPos(scale)(1), scale);
                }

                if (sMaxVal(scale) >= maxVal) {
                    maxVal = sMaxVal(scale);
                    maxPos = IPosition(3, sMaxPos(scale)(0), sMaxPos(scale)(1), scale);
                }
            }

            // If weighting (presumably with weights) was done we need to
            // look up the original values (without the weights).
            minVal = data.xyPlane(minPos(2))(minPos);
            maxVal = data.xyPlane(maxPos(2))(maxPos);
        }
        template<class T, class FT>
        Vector<T> DeconvolverBasisFunction<T, FT>::findCoefficients(const Matrix<Double>& invCoupling,
                const Vector<T>& peakValues)
        {
            const uInt nRows(invCoupling.nrow());
            const uInt nCols(invCoupling.ncolumn());
            Vector<T> coefficients(nRows);

            for (uInt row = 0; row < nRows; row++) {
                coefficients(row) = T(0.0);

                for (uInt col = 0; col < nCols; col++) {
                    coefficients(row) += T(invCoupling(row, col)) * peakValues(col);
                }
            }

            return coefficients;
        }

        template<class T, class FT>
        Array<T> DeconvolverBasisFunction<T, FT>::applyInverse(const Matrix<Double>& invCoupling,
                const Array<T> dataArray)
        {
            Array<T> invDataArray(dataArray.shape(), 0.0);

            const uInt nRows(invCoupling.nrow());
            const uInt nCols(invCoupling.ncolumn());
            const uInt nx = dataArray.shape()(0);
            const uInt ny = dataArray.shape()(1);

            for (uInt j = 0; j < ny; j++) {
                for (uInt i = 0; i < nx; i++) {

                    IPosition currentPosCol(3, i, j, 0);
                    IPosition currentPosRow(3, i, j, 0);

                    for (uInt row = 0; row < nRows; row++) {
                        currentPosRow(2) = row;

                        for (uInt col = 0; col < nCols; col++) {
                            currentPosCol(2) = col;
                            invDataArray(currentPosRow) += T(invCoupling(row, col)) * dataArray(currentPosCol);
                        }
                    }

                }
            }

            return invDataArray;
        }

        template<class T, class FT>
        Array<T> DeconvolverBasisFunction<T, FT>::apply(const Matrix<Double>& coupling,
                const Vector<T> dataVector)
        {
            Vector<T> vecDataVector(dataVector.shape(), 0.0);

            const uInt nRows(coupling.nrow());
            const uInt nCols(coupling.ncolumn());
            const uInt nx = dataVector.nelements();

            for (uInt i = 0; i < nx; i++) {
                for (uInt row = 0; row < nRows; row++) {
                    for (uInt col = 0; col < nCols; col++) {
                        vecDataVector(row) += T(coupling(row, col)) * dataVector(col);
                    }
                }
            }

            return vecDataVector;
        }

    }
}
// namespace synthesis
// namespace askap
