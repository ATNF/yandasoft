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
#include <askap/askap/AskapLogging.h>
#include <casacore/casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/MaskArrMath.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/scimath/Mathematics/MatrixMathLA.h>
#include <askap/scimath/fft/FFT2DWrapper.h>

// Local package includes
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/deconvolution/DeconvolverBasisFunction.h>
#include <askap/deconvolution/MultiScaleBasisFunction.h>

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
        void DeconvolverBasisFunction<T, FT>::setBasisFunction(const boost::shared_ptr<BasisFunction<T> >& bf)
        {
            itsBasisFunction = bf;
        };

        template<class T, class FT>
        const boost::shared_ptr<BasisFunction<T> >& DeconvolverBasisFunction<T, FT>::basisFunction() const
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

                // MV: a bit of technical debt highlighted by casacore's interface change. In principle, we could've
                // had scales as std::vector in the interface to avoid the explicit construction (in this particular case,
                // there is no benefit of using casacore::Vector)
                itsBasisFunction = BasisFunction<Float>::ShPtr(new MultiScaleBasisFunction<Float>(casacore::Vector<float>(scales),
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

            const Array<T> ones(this->itsL1image(0).shape(), static_cast<T>(1.0));
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
                this->itsBasisFunction->allBasisFunctions() = applyInverse(this->itsInverseCouplingMatrix, this->itsBasisFunction->allBasisFunctions());

                this->itsResidualBasisFunction = applyInverse(this->itsInverseCouplingMatrix, this->itsResidualBasisFunction);

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

            this->itsL1image.resize(this->nTerms());
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
                                   << this->itsBasisFunction->shape());

            const IPosition stackShape(this->itsBasisFunction->shape());

            itsResidualBasisFunction.resize(stackShape);

            Cube<FT> basisFunctionFFT(this->itsBasisFunction->shape(), FT(0.));

            // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
            // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
            // but one has to generate basis functions with the right order (or reorder them in a separate routine - leave it for the future)
            // Limit number of fft threads to 8 (more is slower for our fft sizes)
            scimath::FFT2DWrapper<FT> fft2d(true,8);

            // do explicit loop over basis functions here (the original code relied on iterator in old
            // fft2d and, therefore, low level representation of the basis function stack). This way
            // we have more control over the array structure and can transition to the more efficient order
            const casacore::uInt nBases = this->itsBasisFunction->numberBases();
            for (uInt base = 0; base < nBases; ++base) {
                 casacore::Matrix<FT> fftBuffer = basisFunctionFFT.xyPlane(base);
                 casacore::setReal(fftBuffer, this->itsBasisFunction->basisFunction(base));
                 fft2d(fftBuffer, true);
            }

            Matrix<FT> residualFFT(this->dirty().shape().nonDegenerate(2), static_cast<FT>(0.));
            casacore::setReal(residualFFT, this->dirty().nonDegenerate(2));
            fft2d(residualFFT, true);

            // there is no bypass of complex initialisation in casacore-3.4
            Matrix<FT> work(this->model().shape().nonDegenerate(2));
            ASKAPLOG_DEBUG_STR(decbflogger,
                               "Calculating convolutions of residual image with basis functions");

            for (uInt term = 0; term < nBases; ++term) {

                ASKAPASSERT(basisFunctionFFT.xyPlane(term).shape().conform(residualFFT.shape()));
                work = conj(basisFunctionFFT.xyPlane(term)) * residualFFT;
                fft2d(work, false);

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

            Matrix<FT> work(subPsfShape);

            ASKAPLOG_DEBUG_STR(decbflogger, "Shape of basis functions "
                                   << this->itsBasisFunction->shape());

            const IPosition stackShape(this->itsBasisFunction->shape());

            // Now transform the basis functions
            Cube<FT> basisFunctionFFT(this->itsBasisFunction->shape(), 0.);

            // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
            // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
            // but one has to generate basis functions with the right order (or reorder them in a separate routine - leave it for the future)
            // Limit number of fft threads to 8 (more is slower for our fft sizes)
            scimath::FFT2DWrapper<FT> fft2d(true,8);

            // do explicit loop over basis functions here (the original code relied on iterator in old
            // fft2d and, therefore, low level representation of the basis function stack). This way
            // we have more control over the array structure and can transition to the more efficient order
            const casacore::uInt nBases = this->itsBasisFunction->numberBases();
            for (uInt base = 0; base < nBases; ++base) {
                 casacore::Matrix<FT> fftBuffer = basisFunctionFFT.xyPlane(base);
                 casacore::setReal(fftBuffer, this->itsBasisFunction->basisFunction(base));
                 fft2d(fftBuffer, true);
            }

            this->itsPSFBasisFunction.resize(stackShape);

            this->itsScaleFlux.resize(stackShape(2));
            this->itsScaleFlux.set(T(0));

            // Calculate XFR for the subsection only
            Matrix<FT> subXFR(subPsfShape, static_cast<FT>(0.));

            const uInt nx(this->psf().shape()(0));
            const uInt ny(this->psf().shape()(1));

            const IPosition subPsfStart(2, nx / 2 - psfWidth / 2, ny / 2 - psfWidth / 2);
            const IPosition subPsfEnd(2, nx / 2 + psfWidth / 2 - 1, ny / 2 + psfWidth / 2 - 1);
            const IPosition subPsfStride(2, 1, 1);

            Slicer subPsfSlicer(subPsfStart, subPsfEnd, subPsfStride, Slicer::endIsLast);
            this->validatePSF(subPsfSlicer);

            const casacore::IPosition subPsfPeak = this->getPeakPSFPosition().getFirst(2);
            ASKAPLOG_DEBUG_STR(decbflogger, "Peak of PSF subsection at  " << subPsfPeak);
            ASKAPLOG_DEBUG_STR(decbflogger, "Shape of PSF subsection is " << subPsfShape);

            casacore::setReal(subXFR, this->psf().nonDegenerate(2)(subPsfSlicer));
            fft2d(subXFR, true);

            // Now we have all the ingredients to calculate the convolutions
            // of basis function with psf's, etc.
            ASKAPLOG_DEBUG_STR(decbflogger, "Calculating convolutions of Psfs with basis functions");
            itsPSFScales.resize(this->itsBasisFunction->numberBases());

            for (uInt term = 0; term < nBases; ++term) {
                // basis function * psf
                ASKAPASSERT(basisFunctionFFT.xyPlane(term).nonDegenerate().shape().conform(subXFR.shape()));
                work = conj(basisFunctionFFT.xyPlane(term).nonDegenerate()) * subXFR;
                fft2d(work, false);
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
                    casacore::Slicer crossTermsSlicer(crossTermsStart, crossTermsEnd, crossTermsStride, Slicer::endIsLast);
                    crossTermsPSFFFT(crossTermsSlicer).nonDegenerate(2) =
                        basisFunctionFFT.xyPlane(term) * conj(basisFunctionFFT.xyPlane(term1)) * subXFR;
                    // MV: now FFT has to be done explicitly per plane. In principle, it is possible to parallelise but leave it for the future
                    // use reference semantics of casacore arrays to get the right interface
                    casacore::Matrix<FT> tempMatrix(crossTermsPSFFFT(crossTermsSlicer).nonDegenerate(0));
                    fft2d(tempMatrix, true);
                }

            }

            this->itsCouplingMatrix.resize(itsBasisFunction->numberBases(), itsBasisFunction->numberBases());
            // the original wrapper is a method in scimath namespace, new wrapper is a local variable with the same name
            // it also could do explicit iteration over extra dimensions. Now it has to be done explicitly - do it inside the
            // nested for-loop above
            // It was probably a bug for a while that we passed the whole array here, the iteration over extra dimensions was made explicit
            // long time ago, but the problem wasn't triggered because this code is actually unused - need to clear this technical debt if we magically
            // have some spare time
            //scimath::fft2d(crossTermsPSFFFT, true);

            this->itsPSFCrossTerms = real(crossTermsPSFFFT) / T(crossTermsShape(0) * crossTermsShape(1));

            for (uInt term = 0; term < this->itsBasisFunction->numberBases(); term++) {
                crossTermsStart(2) = term;
                crossTermsEnd(2) = term;

                for (uInt term1 = 0; term1 < this->itsBasisFunction->numberBases(); term1++) {
                    crossTermsStart(3) = term1;
                    crossTermsEnd(3) = term1;
                    casacore::Slicer crossTermsSlicer(crossTermsStart, crossTermsEnd, crossTermsStride, Slicer::endIsLast);
                    casacore::IPosition minPos;
                    casacore::IPosition maxPos;
                    T minVal, maxVal;
                    casacore::minMax(minVal, maxVal, minPos, maxPos, this->itsPSFCrossTerms(crossTermsSlicer));
                    this->itsCouplingMatrix(term, term1) = Double(maxVal);
                }

                this->itsCouplingMatrix(term, term) += Double(this->control()->lambda());
            }

            ASKAPLOG_DEBUG_STR(decbflogger, "Coupling matrix " << this->itsCouplingMatrix);
            this->itsInverseCouplingMatrix.resize(this->itsCouplingMatrix.shape());
            invertSymPosDef(this->itsInverseCouplingMatrix, this->itsDetCouplingMatrix, this->itsCouplingMatrix);
            ASKAPLOG_DEBUG_STR(decbflogger, "Coupling matrix determinant " << this->itsDetCouplingMatrix);
            ASKAPLOG_DEBUG_STR(decbflogger, "Inverse coupling matrix " << this->itsInverseCouplingMatrix);

            // the following two methods are only reporting to the log with DEBUG severity,
            // there is no point doing this additional math in the production mode
            #ifdef ASKAP_DEBUG

            // double-check that the inverse really is an inverse (but only by writing the product to the log
            this->reportOnCouplingMatrix();

            // Now look at coupling between adjacent scales: this works well if the
            // scales are ordered.
            this->reportOnAdjacentScaleCoupling();
            #endif
        }

        template<typename T, typename FT>
        void DeconvolverBasisFunction<T, FT>::reportOnCouplingMatrix() const {
            Matrix<T> identity(this->itsCouplingMatrix.shape(), static_cast<T>(0.0));
            const uInt nRows(this->itsCouplingMatrix.nrow());
            const uInt nCols(this->itsCouplingMatrix.ncolumn());
            ASKAPDEBUGASSERT(this->itsCouplingMatrix.shape() == this->itsInverseCouplingMatrix.shape());

            for (uInt row = 0; row < nRows; row++) {
                for (uInt col = 0; col < nCols; col++) {
                    identity(row, col) = sum(this->itsCouplingMatrix.row(row) * this->itsInverseCouplingMatrix.column(col));
                }
            }

            ASKAPLOG_DEBUG_STR(decbflogger, "Coupling matrix * inverse " << identity);
        }

        template<typename T, typename FT>
        void DeconvolverBasisFunction<T, FT>::reportOnAdjacentScaleCoupling() const {
            ASKAPDEBUGASSERT(this->itsBasisFunction);
            // Look at coupling between adjacent scales: this works well if the scales are ordered.
            const casacore::Matrix<casacore::Double>& cm = this->itsCouplingMatrix;
            for (uInt term = 0; term < this->itsBasisFunction->numberBases() - 1; term++) {
                ASKAPDEBUGASSERT(term < cm.nrow() && term < cm.ncolumn());
                const double det = cm(term, term) * cm(term + 1, term + 1) - cm(term, term + 1) * cm(term + 1, term);
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
            casacore::IPosition minPos;
            casacore::IPosition maxPos;
            T minVal(0.0), maxVal(0.0);
            // Here the weights image is used as a weight in the determination
            // of the maximum i.e. it finds the max in weight . residual. The values
            // returned are without the weight
            minMaxMaskedScales(minVal, maxVal, minPos, maxPos, this->itsResidualBasisFunction,
                               this->weight(0));
            casacore::IPosition absPeakPos;

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

            const casacore::IPosition residualShape(this->itsResidualBasisFunction.shape());
            const casacore::IPosition psfShape(this->itsPSFBasisFunction.shape());

            const casacore::uInt ndim(residualShape.size());
            ASKAPDEBUGASSERT(ndim > 2);

            casacore::IPosition residualStart(ndim, 0), residualEnd(ndim, 0), residualStride(ndim, 1);
            casacore::IPosition psfStart(ndim, 0), psfEnd(ndim, 0), psfStride(ndim, 1);
            casacore::IPosition psfCrossTermsStart(ndim + 1, 0), psfCrossTermsEnd(ndim + 1, 0), psfCrossTermsStride(ndim + 1, 1);

            const casacore::IPosition modelShape(this->model().shape());
            const casacore::uInt modelNdim(this->model().shape().size());
            casacore::IPosition modelStart(modelNdim, 0), modelEnd(modelNdim, 0), modelStride(modelNdim, 1);

            // Wrangle the start, end, and shape into consistent form. It took me
            // quite a while to figure this out (slow brain day) so it may be
            // that there are some edge cases for which it fails.

            const casacore::IPosition peakPSFPos = this->getPeakPSFPosition();
            ASKAPDEBUGASSERT(peakPSFPos.nelements() >= 2);

            for (uInt dim = 0; dim < 2; dim++) {
                residualStart(dim) = max(0, Int(absPeakPos(dim) - psfShape(dim) / 2));
                residualEnd(dim) = min(Int(absPeakPos(dim) + psfShape(dim) / 2 - 1), Int(residualShape(dim) - 1));
                // Now we have to deal with the PSF. Here we want to use enough of the
                // PSF to clean the residual image.
                psfStart(dim) = max(0, Int(peakPSFPos(dim) - (absPeakPos(dim) - residualStart(dim))));
                psfEnd(dim) = min(Int(peakPSFPos(dim) - (absPeakPos(dim) - residualEnd(dim))),
                                  Int(psfShape(dim) - 1));

                psfCrossTermsStart(dim) = psfStart(dim);
                psfCrossTermsEnd(dim) = psfEnd(dim);

                modelStart(dim) = residualStart(dim);
                modelEnd(dim) = residualEnd(dim);
            }

            casacore::Slicer modelSlicer(modelStart, modelEnd, modelStride, Slicer::endIsLast);

            // Add to model
            // Note that the model is only two dimensional. We could make it three dimensional
            // and keep the model layers separate
            // We loop over all terms and ignore those with no flux
            const casacore::uInt nterms(this->itsResidualBasisFunction.shape()(2));
            ASKAPDEBUGASSERT(nterms == this->itsBasisFunction->numberBases());

            for (uInt term = 0; term < nterms; term++) {
                if (abs(peakValues(term)) > 0.0) {
                    //psfStart(2) = psfEnd(2) = term;
                    // slicer operates on a 2D image
                    casacore::Slicer psfSlicer(psfStart.getFirst(2), psfEnd.getFirst(2), psfStride.getFirst(2), Slicer::endIsLast);
                    typename casacore::Array<T> modelSlice = this->model()(modelSlicer).nonDegenerate();
                    modelSlice += this->control()->gain() * peakValues(term) *
                                  this->itsBasisFunction->basisFunction(term)(psfSlicer).nonDegenerate();
                }
            }

            // Keep track of strengths and locations of components
            for (uInt term = 0; term < nterms; term++) {
                if (abs(peakValues(term)) > 0.0) {
                    IPosition l1PeakPos(3, absPeakPos(0), absPeakPos(1), term);
                    casacore::Slicer modelSlicer(modelStart, modelEnd, modelStride, Slicer::endIsLast);
                    this->itsL1image(0)(l1PeakPos) += this->control()->gain() * abs(peakValues(term));
                    this->itsScaleFlux(term) += this->control()->gain() * peakValues(term);
                }
            }

            // Subtract PSFs
            for (uInt term = 0; term < nterms; term++) {
                if (abs(peakValues(term)) > 0.0) {
                    psfStart(2) = psfEnd(2) = term;
                    casacore::Slicer psfSlicer(psfStart, psfEnd, psfStride, Slicer::endIsLast);
                    residualStart(2) = residualEnd(2) = term;
                    casacore::Slicer residualSlicer(residualStart, residualEnd, residualStride, Slicer::endIsLast);
                    typename casacore::Array<T> residualBFSlice = this->itsResidualBasisFunction(residualSlicer).nonDegenerate();
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
                                casacore::Slicer residualSlicer(residualStart, residualEnd, residualStride, Slicer::endIsLast);

                                psfCrossTermsStart(2) = term1;
                                psfCrossTermsEnd(2) = term1;
                                psfCrossTermsStart(3) = term;
                                psfCrossTermsEnd(3) = term;
                                casacore::Slicer psfCrossTermsSlicer(psfCrossTermsStart, psfCrossTermsEnd, psfCrossTermsStride, Slicer::endIsLast);
                                typename casacore::Array<T> residualBFSlice = this->itsResidualBasisFunction(residualSlicer).nonDegenerate();
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
                        casacore::minMaxMasked(sMinVal(scale), sMaxVal(scale), sMinPos(scale), sMaxPos(scale),
                                           Cube<T>(dataArray).xyPlane(scale), weightArray.nonDegenerate());
                    }
                } else {
                    for (uInt scale = 0; scale < nScales; scale++) {
                        casacore::minMax(sMinVal(scale), sMaxVal(scale), sMinPos(scale), sMaxPos(scale),
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
            // MV: I don't understand the following - it seems like the dimensions mismatch - I'll replace it
            // (it probably worked fine because the last dimension goes beyond 3 dimensions of a cube and remains unchecked)
            //minVal = data.xyPlane(minPos(2))(minPos);
            //maxVal = data.xyPlane(maxPos(2))(maxPos);
            minVal = data(minPos);
            maxVal = data(maxPos);
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

        /// @brief basically matrix multiplication across the basis function domain
        /// @details This method applies inverse coupling matrix to every pixel of the data array
        /// @param[in] invCoupling inverse coupling matrix (i.e. the matrix multiplied by the data array)
        /// @param[in] dataArray data array (with basis function decomposition)
        /// @return new data array after the inverse coupling matrix is applied
        template<class T, class FT>
        Array<T> DeconvolverBasisFunction<T, FT>::applyInverse(const Matrix<Double>& invCoupling,
                const Array<T>& dataArray)
        {
            Array<T> invDataArray(dataArray.shape(), static_cast<T>(0.0));

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
            Vector<T> vecDataVector(dataVector.shape(), static_cast<T>(0.0));

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
