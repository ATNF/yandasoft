/// @file DeconvolverMultiTermBasisFunction.tcc
/// @brief Class for a deconvolver based on CLEANing with basis functions.
/// @details This concrete class defines a deconvolver used to estimate an
/// image from a residual image, psf optionally using a weights image.
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

#include <string>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(decmtbflogger, ".deconvolution.multitermbasisfunction");
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/profile/AskapProfiler.h>

#include <askap/deconvolution/DeconvolverMultiTermBasisFunction.h>
#include <askap/deconvolution/MultiScaleBasisFunction.h>
#include <askap/scimath/utils/OptimizedArrayMathUtils.h>
#include <askap/scimath/utils/OptimizedArrayMathUtils.h>
//#include <mpi.h>

#include <askap/scimath/fft/FFT2DWrapper.h>

namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver based on the BasisFunction Clean
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a residual image, psf optionally using a weights image.
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverMultiTermBasisFunction<double, DComplex>
        /// @ingroup Deconvolver

        template<class T> 
        void absMaxPos(T& maxVal, T& maxValScaled, IPosition& maxPos, const Matrix<T>& im,
            const Matrix<T>& mask, const std::vector<uInt>& pixels, const Matrix<T>& noise, uInt boxSize, uInt increment) 
        {
            // Set Shared Values
            maxVal = T(0.0);
            maxValScaled = T(0.0);
            // Set Private Values
            T maxVal_private(0), maxValScaled_private(0);
            uInt maxIndex_private = 0;
            ASKAPASSERT(im.contiguousStorage() && mask.contiguousStorage());
            const T* pIm = im.data();
            const T* pMask = mask.data();
            const bool useMask = mask.size() > 0;
            const bool usePixels = pixels.size() > 0;
            const bool useNoise = noise.size() > 0 && boxSize > 0;
            const uInt nrow = im.nrow();
            ASKAPDEBUGASSERT(!useNoise || (noise.nrow() >= (nrow - 1)/boxSize && noise.ncolumn() >= (im.ncolumn() - 1)/boxSize));
            const uInt n = (usePixels ? pixels.size() : im.size());
            ASKAPDEBUGASSERT(nrow > 0);

            #pragma omp for schedule(static)
            for (uInt j = 0; j < n; j++ ) {
                const uInt pixel = (usePixels ? pixels[j] : j);
                // skip pixels (for larger scales)
                if (increment == 1|| ((pixel % nrow)%increment == 0 &&
                         (pixel / nrow)%increment == 0)) {
                    T val = abs(pIm[pixel]);
                    if (useMask) val *= pMask[pixel];
                    T testVal = val;
                    if (useNoise) {
                        const uInt row = (pixel % nrow) / boxSize;
                        const uInt col = (pixel / nrow) / boxSize;
                        testVal /= noise(row, col);
                    } 
                    if (testVal > maxValScaled_private) {
                        maxValScaled_private = testVal;
                        maxVal_private = val;
                        maxIndex_private = pixel;
                    }
                }
            }
            #pragma omp critical
            {
                if (maxValScaled_private > maxValScaled) {
                    maxVal = maxVal_private;
                    maxValScaled = maxValScaled_private;
                    maxPos(0) = maxIndex_private % nrow;
                    maxPos(1) = maxIndex_private / nrow;
                }
            }
            #pragma omp barrier
        }

        template<class T, class FT>
        DeconvolverMultiTermBasisFunction<T, FT>::DeconvolverMultiTermBasisFunction(Vector<Array<T>>& dirty,
                Vector<Array<T>>& psf,
                Vector<Array<T>>& psfLong)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf), itsDirtyChanged(True), itsBasisFunctionChanged(True),
                itsSolutionType("MAXBASE"), itsUsePixelLists(true),
                itsPixelListTolerance(0.1), itsPixelListNSigma(4.0), itsPixelListNPixRange(std::vector<float>({2.0,10.0})),
                itsNoiseBoxSize(0)
        {
            ASKAPLOG_DEBUG_STR(decmtbflogger, "There are " << this->nTerms() << " terms to be solved");

            ASKAPCHECK(psfLong.size() == (2*this->nTerms() - 1), "Long PSF vector has incorrect length " << psfLong.size());
            itsPsfLongVec.resize(2*this->nTerms() - 1);

            for (uInt term = 0; term < (2*this->nTerms() - 1); ++term) {
                ASKAPCHECK(psfLong(term).nonDegenerate().shape().size() == 2, "PSF(" << term << ") has too many dimensions " << psfLong(term).shape());
                itsPsfLongVec(term).reference(psfLong(term).nonDegenerate());
            }
        }

        template<class T, class FT>
        DeconvolverMultiTermBasisFunction<T, FT>::DeconvolverMultiTermBasisFunction(Array<T>& dirty,
                Array<T>& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf), itsDirtyChanged(True), itsBasisFunctionChanged(True),
                itsSolutionType("MAXBASE"), itsUsePixelLists(true),
                itsPixelListTolerance(0.1), itsPixelListNSigma(4.0), itsPixelListNPixRange(std::vector<float>({2.0,10.0})),
                itsNoiseBoxSize(0)

        {
            ASKAPLOG_DEBUG_STR(decmtbflogger, "There is only one term to be solved");
            itsPsfLongVec.resize(1);
            itsPsfLongVec(0) = psf;
        }

        template<class T, class FT>
        DeconvolverMultiTermBasisFunction<T, FT>::~DeconvolverMultiTermBasisFunction()
        {
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::setSolutionType(const std::string& sol)
        {
            itsSolutionType = sol;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::setBasisFunction(boost::shared_ptr<BasisFunction<T>> bf)
        {
            itsBasisFunction = bf;
            itsBasisFunctionChanged = True;
        }

        template<class T, class FT>
        boost::shared_ptr<BasisFunction<T>> DeconvolverMultiTermBasisFunction<T, FT>::basisFunction()
        {
            return itsBasisFunction;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::updateDirty(Array<T>& dirty, uInt term)
        {
            DeconvolverBase<T, FT>::updateDirty(dirty, term);
            itsDirtyChanged = True;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::updateDirty(Vector<Array<T>>& dirtyVec)
        {
            DeconvolverBase<T, FT>::updateDirty(dirtyVec);
            itsDirtyChanged = True;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::configure(const LOFAR::ParameterSet& parset)
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::configure");
            // Configuring the base class causes the testspectral test to fail.
            // I suspect because ImageAMSMFSolver replaces the DeconvolverControl object.
            // Task cdeconvolver-mpi will set configurebase="true" because it needs this configuration done.
            if (parset.getString("configurebase","false")=="true") {
                DeconvolverBase<T, FT>::configure(parset); 
            }

            // Make the basis function
            std::vector<float> defaultScales({0.0,10.0,30.0});
            std::vector<float> scales = parset.getFloatVector("scales", defaultScales);
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Constructing Multiscale basis function with scales "
                                   << scales);
            bool orthogonal = parset.getBool("orthogonal", false);
            if (orthogonal) {
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Multiscale basis functions will be orthogonalised");
            }

            // MV: a bit of technical debt highlighted by casacore's interface change. In principle, we could've
            // had scales as std::vector in the interface to avoid the explicit construction (in this particular case,
            // there is no benefit of using Vector)
            itsBasisFunction = BasisFunction<float>::ShPtr(new MultiScaleBasisFunction<float>(Vector<float>(scales),
                               orthogonal));

            String solutionType = parset.getString("solutiontype", "MAXBASE");
            if(solutionType!="MAXCHISQ") {
               solutionType="MAXBASE";
            }

            if (solutionType == "MAXBASE") {
                itsSolutionType = solutionType;
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Component search to maximise over bases");
            } else {
                itsSolutionType = "MAXCHISQ";
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Component search to find maximum in chi-squared");
            }
            itsScaleMaskName = parset.getString("scalemask","");
            if (itsScaleMaskName != "") {
                ASKAPLOG_INFO_STR(decmtbflogger, "Read scale mask from image: "<<itsScaleMaskName);
                setScaleMask(SynthesisParamsHelper::imageHandler().read(itsScaleMaskName).nonDegenerate());
            }
            itsUsePixelLists = parset.getBool("usepixellists", true);
            // pixellists only implemented for MAXBASE
            if (itsUsePixelLists && itsSolutionType != "MAXBASE") {
                ASKAPLOG_WARN_STR(decmtbflogger,"Disabled usepixellists because of solutiontype "<<itsSolutionType);
                itsUsePixelLists = false;
            }
            if (itsUsePixelLists) {
                ASKAPLOG_INFO_STR(decmtbflogger, "Using pixel lists with active (high) pixels");
            }
            itsPixelListTolerance = parset.getFloat("usepixellists.tolerance",0.1);
            itsPixelListNSigma = parset.getFloat("usepixellists.nsigma",4.0);
            itsPixelListNPixRange = parset.getFloatVector("usepixellists.npixrange",std::vector<float>({2.0,10.0}));
            ASKAPCHECK(itsPixelListNPixRange.size()==2,"npixrange needs to have 2 values");
            ASKAPCHECK(itsPixelListNPixRange[0]<itsPixelListNPixRange[1],"first value of npixrange needs to be smaller than second");

            itsUseIncrements = parset.getBool("useincrements",false);
            if (itsUseIncrements) {
                ASKAPLOG_INFO_STR(decmtbflogger, "Using larger pixel increments for larger scales");
            }
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::finalise()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::finalise");
            this->updateResiduals(this->itsModel);

            // debug message
            for (uInt base = 0; base < itsTermBaseFlux.size(); base++) {
                for (uInt term = 0; term < itsTermBaseFlux(base).size(); term++) {
                    ASKAPLOG_DEBUG_STR(decmtbflogger, "   Term(" << term << "), Base(" << base
                                           << "): Flux = " << itsTermBaseFlux(base)(term));
                }
            }

            // info message
            for (uInt base = 0; base < itsTermBaseFlux.size(); base++) {
              ASKAPLOG_INFO_STR(decmtbflogger,"Total flux for scale "<<base<<" : "<<itsTermBaseFlux(base)(0));
            }
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialise()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialise");
            // This one does not appear to be needed
            //DeconvolverBase<T, FT>::initialise();

            // Initialise residuals
            initialiseResidual();

            // Initialise masks
            initialiseMask();

            // Force change in basis function
            initialiseForBasisFunction(true);
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialiseResidual()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialiseResidual");

            if (!itsDirtyChanged) {
                return;
            }

            // Initialise the basis function for residual calculations.
            ASKAPCHECK(itsBasisFunction, "Basis function not initialised");
            itsBasisFunction->initialise(this->dirty(0).shape());

            const uInt nBases(itsBasisFunction->numberBases());
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of basis functions "
                                   << itsBasisFunction->shape()<<" number of bases "<<nBases);

            itsResidualBasis.resize(nBases);
            for (uInt base = 0; base < nBases; base++) {
                itsResidualBasis(base).resize(this->nTerms());
            }

            // Calculate residuals convolved with bases [nx,ny][nterms][nbases]

            ASKAPLOG_INFO_STR(decmtbflogger,
                              "Calculating convolutions of residual images with basis functions");
            askap::utils::Timer timer;
            timer.start();

            // We can do bases in parallel & reduce #threads for FFTs - but it uses more memory,
            //  this is OK for continuum, but not for spectral line
            // const uInt nthreads = LOFAR::OpenMP::maxThreads();
            // const uInt nFFTthreads = min(8u,max(1u, nthreads / nBases));
            // scimath::FFT2DWrapper<FT> fft2d(true,8), nFFTthreads);
            // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
            // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
            // Limit number of fft threads to 8 (more is slower for our fft sizes)
            scimath::FFT2DWrapper<FT> fft2d(true,8);
            //#pragma omp parallel
            {
                //#pragma omp for
                for (uInt base = 0; base < nBases; base++) {
                     // Calculate transform of basis function [nx,ny,nbases]
                     const Matrix<T>& bfRef(itsBasisFunction->basisFunction(base));
                     Matrix<FT> basisFunctionFFT(bfRef.shape(), 0.);
                     setReal(basisFunctionFFT, bfRef);
                     fft2d(basisFunctionFFT, true);

                     for (uInt term = 0; term < this->nTerms(); term++) {

                        // Calculate transform of residual image
                        Matrix<FT> residualFFT(this->dirty(term).shape().nonDegenerate(), 0.);
                        setReal(residualFFT, this->dirty(term).nonDegenerate());
                        fft2d(residualFFT, true);

                        // Calculate product and transform back
                        ASKAPASSERT(basisFunctionFFT.shape().conform(residualFFT.shape()));

                        //the following line is equivalent to the optimised version called below
                        //residualFFT *= conj(basisFunctionFFT);
                        utility::multiplyByConjugate(residualFFT, basisFunctionFFT);

                        fft2d(residualFFT, false);

                        // temporary object is ok here because we do an assignment to uninitialised array later on
                        Matrix<T> work(real(residualFFT));
    #ifdef ASKAP_DEBUG
                        ASKAPLOG_DEBUG_STR(decmtbflogger, "Basis(" << base
                                                << ")*Residual(" << term << "): max = " << max(work)
                                                << " min = " << min(work));
    #endif
                        itsResidualBasis(base)(term).reference(work);
                    }
                }
            }
            timer.stop();
            ASKAPLOG_INFO_STR(decmtbflogger,
                              "Time to calculate residual images * basis functions: "<< timer.elapsedTime() << " sec");
        }
        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialiseMask()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialiseMask");
            ASKAPLOG_DEBUG_STR(decmtbflogger, "initialiseMask called");

            // check if we need the masks
            if (this->control()->targetObjectiveFunction2()==0) {
                return;
            }
            // check if we've already done this
            if (itsScalePixels.size()>0) {
                return;
            }
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Initialising deep clean masks");

            ASKAPCHECK(itsBasisFunction, "Basis function not initialised");

            uInt nBases(itsBasisFunction->numberBases());

            // Resize array that keeps track of pixels to clean at each scale
            itsScalePixels.resize(nBases);
            return;
        }
        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialiseForBasisFunction(bool force)
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialiseForBasisFunction");
            if (!force && !itsBasisFunctionChanged) {
                return;
            }

            ASKAPLOG_DEBUG_STR(decmtbflogger,
                               "Updating Multi-Term Basis Function deconvolver for change in basis function");

            IPosition subPsfShape(this->findSubPsfShape());

            // Use a smaller size for the psfs if specified.
            itsBasisFunction->initialise(subPsfShape);

            ASKAPLOG_DEBUG_STR(decmtbflogger, "Initialising for PSFs: shape = " << subPsfShape);
            initialisePSF();

            itsBasisFunctionChanged = False;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialisePSF()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialisePSF");

            if (!itsBasisFunctionChanged) {
                return;
            }

            ASKAPCHECK(itsBasisFunction, "Basis function not initialised");

            ASKAPLOG_DEBUG_STR(decmtbflogger,
                               "Updating Multi-Term Basis Function deconvolver for change in basis function");
            const IPosition subPsfShape(this->findSubPsfShape());

            Array<FT> work(subPsfShape);

            const uInt nBases(itsBasisFunction->numberBases());
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of basis functions "
                                   << itsBasisFunction->shape()<< " number of bases "<<nBases);

            const IPosition stackShape(itsBasisFunction->shape());

            // Now transform the basis functions. These may be a different size from
            // those in initialiseResidual so we don't keep either
            Cube<FT> basisFunctionFFT(stackShape,FT(0.));

            // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
            // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
            // Limit number of fft threads to 8 (more is slower for our fft sizes)
            scimath::FFT2DWrapper<FT> fft2d(true,8);

            // do explicit loop over basis functions here (the original code relied on iterator in
            // fft2d and, therefore, low level representation of the basis function stack). This way
            // we have more control over the array structure and can transition to the more efficient order
            for (uInt base = 0; base < nBases; ++base) {
                 // casacore arrays have reference semantics, no copying occurs in the following
                 Matrix<FT> fftBuffer = basisFunctionFFT.xyPlane(base);
                 setReal(fftBuffer, itsBasisFunction->basisFunction(base));
                 fft2d(fftBuffer, true);
            }

            itsTermBaseFlux.resize(nBases);
            for (uInt base = 0; base < nBases; base++) {
                itsTermBaseFlux(base).resize(this->nTerms());
                itsTermBaseFlux(base) = 0.0;
            }

            const uInt nx(this->psf(0).shape()(0));
            const uInt ny(this->psf(0).shape()(1));

            const IPosition subPsfStart(2, (nx - subPsfShape(0)) / 2, (ny - subPsfShape(1)) / 2);
            Slicer subPsfSlicer(subPsfStart, subPsfShape);
            // check just in case
            ASKAPCHECK(subPsfSlicer.length() == subPsfShape, "Slicer selected length of " <<
                subPsfSlicer.length() << " is different from requested shape " << subPsfShape);

            this->validatePSF(subPsfSlicer);

            const IPosition subPsfPeak=this->getPeakPSFPosition().getFirst(2);
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Peak of PSF subsection at  " << subPsfPeak);
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of PSF subsection is " << subPsfShape);

            // Calculate XFR for the subsection only. We need all PSF's up to
            // 2*nTerms-1
            ASKAPCHECK(itsPsfLongVec.size() == (2*this->nTerms() - 1),
                "PSF long vector has wrong length " << itsPsfLongVec.size());

            // Calculate all the transfer functions
            Vector<Array<FT>> subXFRVec(2*this->nTerms() - 1);
            for (uInt term1 = 0; term1 < subXFRVec.size(); ++term1) {
                subXFRVec(term1).resize(subPsfShape);
                // rely on reference semantics of casa arrays
                // MV: we can probably change subXFRVec to be a vector of matrices to reduce technical debt
                Matrix<FT> subXFRTerm1(subXFRVec(term1));
                subXFRTerm1.set(0.0);
                setReal(subXFRTerm1, itsPsfLongVec(term1).nonDegenerate()(subPsfSlicer));
                fft2d(subXFRTerm1, true);
                // we only need conjugated FT of subXFRVec (or real part of it, which doesn't change with conjugation),
                // it is better to compute conjugation in situ now and don't do it on the fly later
                utility::conjugateComplexArray(subXFRTerm1);
            }
            // Calculate residuals convolved with bases [nx,ny][nterms][nbases]

            // the following line is the original code which we do now in an optimised way
            //const T normPSF = sum(real(subXFRVec(0))) / subXFRVec(0).size();
            const T normPSF = utility::realPartMean(subXFRVec(0));
            ASKAPLOG_DEBUG_STR(decmtbflogger, "PSF effective volume = " << normPSF);

            itsPSFCrossTerms.resize(nBases, nBases);
            for (uInt base = 0; base < nBases; base++) {
                for (uInt base1 = 0; base1 < nBases; base1++) {
                    itsPSFCrossTerms(base, base1).resize(this->nTerms(), this->nTerms());
                }
            }

            itsCouplingMatrix.resize(nBases);
            for (uInt base1 = 0; base1 < nBases; base1++) {
                itsCouplingMatrix(base1).resize(this->nTerms(), this->nTerms());
                for (uInt base2 = base1; base2 < nBases; base2++) {
                    for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                        for (uInt term2 = term1; term2 < this->nTerms(); ++term2) {

                            // the following expression is what we had here originally. It is replaced by an optimised
                            // method allowing us to avoid creation of temporary objects (+ it is normally faster if OMP is used)
                            // note, the procedure doesn't have conj(subXFRVec(term1 + term2)) and for this we conjugated the whole
                            // subXFRVec(term1 + term2) above, when it is filled with values
                            //work = conj(basisFunctionFFT.xyPlane(base1)) * basisFunctionFFT.xyPlane(base2) *
                            //       conj(subXFRVec(term1 + term2)) / normPSF;
                            utility::calculateNormalisedProduct(work, basisFunctionFFT.xyPlane(base1), basisFunctionFFT.xyPlane(base2), subXFRVec(term1 + term2), normPSF);

                            //use reference semantics to get the right interface, we can probably change the interface to matrix to reduce technical debt
                            Matrix<FT> workMtr(work);
                            fft2d(workMtr, false);

                            ASKAPLOG_DEBUG_STR(decmtbflogger, "Base(" << base1 << ")*Base(" << base2
                                                   << ")*PSF(" << term1 + term2
                                                   << "): max = " << max(real(work))
                                                   << " min = " << min(real(work))
                                                   << " centre = " << real(work(subPsfPeak)));
                            // Remember that Array reuses the same memory where possible so this
                            // apparent redundancy does not cause any memory bloat
                            // I don't think that is true here: simple assigment does not share memory only the copy constructor does
                            // Need to use .reference() to get the behavior wanted
                            itsPSFCrossTerms(base1, base2)(term1, term2) = real(work);
                            itsPSFCrossTerms(base2, base1)(term1, term2).reference(itsPSFCrossTerms(base1, base2)(term1, term2));
                            itsPSFCrossTerms(base1, base2)(term2, term1).reference(itsPSFCrossTerms(base1, base2)(term1, term2));
                            itsPSFCrossTerms(base2, base1)(term2, term1).reference(itsPSFCrossTerms(base1, base2)(term1, term2));
                            if (base1 == base2) {
                                const T subPsfPeakValue =  real(work(subPsfPeak));
                                itsCouplingMatrix(base1)(term1, term2) = subPsfPeakValue;
                                itsCouplingMatrix(base1)(term2, term1) = subPsfPeakValue;
                            }
                        }
                    }
                }
            }

            ASKAPLOG_DEBUG_STR(decmtbflogger, "Calculating inverses of coupling matrices");

            // Invert the coupling matrices and check for correctness
            itsInverseCouplingMatrix.resize(nBases);
            Vector<double> detCouplingMatrix(nBases);

            for (uInt base = 0; base < nBases; base++) {
                itsInverseCouplingMatrix(base).resize(this->nTerms(), this->nTerms());
                ASKAPLOG_INFO_STR(decmtbflogger, "Coupling matrix(" << base << ")="
                                       << itsCouplingMatrix(base).row(0));
                for (uInt term = 1; term < this->nTerms(); ++term) {
                    ASKAPLOG_INFO_STR(decmtbflogger, "                   "
                                       << itsCouplingMatrix(base).row(term));
                }
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Calculating matrix inverse by Cholesky decomposition");
                invertSymPosDef(itsInverseCouplingMatrix(base),
                                detCouplingMatrix(base), itsCouplingMatrix(base));
                ASKAPLOG_INFO_STR(decmtbflogger, "Coupling matrix determinant(" << base << ") = "
                                       << detCouplingMatrix(base));
                ASKAPLOG_INFO_STR(decmtbflogger, "Inverse coupling matrix(" << base
                                       << ")=" << itsInverseCouplingMatrix(base).row(0));
                for (uInt term = 1; term < this->nTerms(); ++term) {
                    ASKAPLOG_INFO_STR(decmtbflogger, "                           "
                                       << itsInverseCouplingMatrix(base).row(term));
                }
            }
            itsBasisFunctionChanged = False;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::ManyIterations()
        {
            // Need to (re)set the subPsfShape and call validatePSF to get correct psf peak pixels as it gets unset by DeconvolverBase::initialise
            const IPosition subPsfShape(this->findSubPsfShape());
            const uInt nx(this->psf(0).shape()(0));
            const uInt ny(this->psf(0).shape()(1));
            const IPosition subPsfStart(2, (nx - subPsfShape(0)) / 2, (ny - subPsfShape(1)) / 2);
            this->validatePSF(Slicer(subPsfStart, subPsfShape));

            const uInt nBases(itsBasisFunction->numberBases());
            IPosition absPeakPos(2, 0);
            T absPeakVal(0.0), absPeakValScaled(0.0);
            uInt optimumBase(0);
            Vector<T> peakValues(this->nTerms());
            IPosition maxPos(2, 0);
            T maxVal(0.0), maxValScaled(0.0);
            Vector<Matrix<T>> coefficients(this->nTerms());
            Matrix<T> weights;
            Matrix<T> negchisq;
            IPosition residualShape(2);
            IPosition psfShape(2);
            const bool isWeighted((this->itsWeight.size() > 0) &&
                (this->weight(0).shape().nonDegenerate().conform(itsResidualBasis(0)(0).shape())));
            IPosition shape(2,0), resStart(2,0), psfStart(2,0);

            // temporary matrix reference
            Matrix<T> mat;

            // Timers for analysis
            const int no_timers = 9;
            askap::utils::SectionTimer sectionTimer(no_timers);

	      	// Termination
	      	int converged;

            bool fillHighPixels = itsUsePixelLists;
            bool listScalePixels = true;
            bool firstCycle = true;
            if (this->control()->targetIter() != 0) {

            std::vector<std::vector<uInt>> highPixels;
            if (itsUsePixelLists) {
                highPixels.resize(nBases);
            }

            askap::utils::Timer timer;
            timer.start();

            #pragma omp parallel
            {
                #pragma omp master
                {
                    const uInt nthreads = LOFAR::OpenMP::numThreads();
                    if (nthreads>1) ASKAPLOG_INFO_STR(decmtbflogger, "Cleaning using "<<nthreads<< " threads");
                }


                // =============== Set weights =======================

                // Section 0
                sectionTimer.start(0);

                if (isWeighted) {

                    #pragma omp single
                    weights = this->weight(0).nonDegenerate();

                    // Check weights for contiguity
                    if (!weights.contiguousStorage()) {
                        ASKAPLOG_WARN_STR(decmtbflogger, "weights (sec 0) is not contiguous\n");
                    }

                    if  (itsSolutionType == "MAXCHISQ") {
                        const uInt n = weights.size();
                        T* pWeights = weights.data();
                        // square weights for MAXCHISQ
                        #pragma omp for schedule(static)
                        for (size_t i = 0; i < n; i++ ) {
                            pWeights[i] *= pWeights[i];
                        }
                    }
                }

                sectionTimer.stop(0);

                // Commence cleaning iterations
                do {

                    // Reset peak pos
                    absPeakPos = 0;
                    // Reset peak Val
                    absPeakVal = 0.0;
                    absPeakValScaled = 0.0;
                    // Reset optimum base
                    optimumBase = 0;

                    // =============== Choose Component =======================
                    // Note we pass the section timer and use slots for section 1, 2 and 3
                    chooseComponent(optimumBase, absPeakPos, absPeakVal, absPeakValScaled, firstCycle, highPixels,
                        sectionTimer, maxPos, maxVal, maxValScaled, weights, negchisq, coefficients);

                    // Now that we know the location of the peak found using one of the
                    // above methods we can look up the values of the residuals. Remember
                    // that we have to decouple the answer

                    // Section 4
                    sectionTimer.start(4);

                    #pragma omp single
                    {
                        for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                            peakValues(term1) = 0.0;
                            for (uInt term2 = 0; term2 < this->nTerms(); ++term2) {
                                peakValues(term1) +=
                                    T(itsInverseCouplingMatrix(optimumBase)(term1, term2)) *
                                    itsResidualBasis(optimumBase)(term2)(absPeakPos);
                            }
                        }

                        // Record location of peak
                        ASKAPDEBUGASSERT(itsScalePixels.size()==0 || optimumBase < itsScalePixels.size());
                        if (itsScalePixels.size()) itsScalePixels[optimumBase].insert(absPeakPos[0]+absPeakPos[1]*nx);
                        // Take square root to get value comparable to peak residual
                        if (itsSolutionType == "MAXCHISQ") {
                            absPeakVal = sqrt(max(T(0.0), absPeakVal));
                            absPeakValScaled = sqrt(max(T(0.0), absPeakValScaled));
                        }
                    } // End of omp single section

                    // End of section 4
                    sectionTimer.stop(4);

                    // Section 5
                    sectionTimer.start(5);
                    #pragma omp single
                    {

                        if (this->state()->initialObjectiveFunction() == 0.0) {
                            this->state()->setInitialObjectiveFunction(abs(absPeakValScaled));
                        }
                        this->state()->setPeakResidual(abs(absPeakVal));
                        this->state()->setObjectiveFunction(abs(absPeakValScaled));
                    } // End of single

                    #pragma omp master
                    {
                        if (listScalePixels && this->control()->deepCleanMode()) {
                            for (uInt base = 0; base < nBases; base++) {
                                ASKAPLOG_DEBUG_STR(decmtbflogger,"Base "<<base<<" has "<<itsScalePixels[base].size()<<" active pixels for deep clean");
                            }
                            listScalePixels = false;
                        }
                    }

                    #pragma omp single
                    {
                        float sumFlux = 0.0;
                        for (uInt base=0; base < nBases; base++) {
                            sumFlux += itsTermBaseFlux(base)(0);
                        }
                        this->state()->setTotalFlux(sumFlux);
                        ASKAPLOG_DEBUG_STR(decmtbflogger,"Peak="<<absPeakVal<<", Pos="<< absPeakPos <<", Base="<<optimumBase<<", Total flux = "<<sumFlux);
                    }
                    // End of section 5
                    sectionTimer.stop(5);

                    // Section 6
                    sectionTimer.start(6);
                    getResidualAndPSFSlice(absPeakPos, shape, resStart, psfStart);
                    addComponentToModel(peakValues, shape, resStart, psfStart, optimumBase, mat);
                    // End of section 6
                    sectionTimer.stop(6);

                    // Section 7
                    sectionTimer.start(7);
                    const bool useScalePixels = this->control()->deepCleanMode();
                    const bool useHighPixels = itsUsePixelLists && !useScalePixels && !firstCycle ;
                    if (useScalePixels || useHighPixels) {
                        subtractPSFPixels(peakValues, absPeakPos, optimumBase, useHighPixels, highPixels);
                    } else {
                        subtractPSF(peakValues, shape, resStart, psfStart, optimumBase);
                    }
                    // End of section 7
                    sectionTimer.stop(7);

                    // Section 8
                    sectionTimer.start(8);
                    // if needed fill the list of high pixels we'll use for peak finding and cleaning residuals
                    if (fillHighPixels && itsUsePixelLists && !useScalePixels) {
                        fillHighPixelList(highPixels,weights);
                        #pragma omp single
                        fillHighPixels = false;
                    }

                    #pragma omp master
                    {
						this->monitor()->monitor(*(this->state()));
						this->state()->incIter();
                    }

                    //End of all iterations
                    #pragma omp barrier

					#pragma omp single
					{
						converged = this->control()->terminate(*(this->state()));
                        firstCycle = false;
					}
                    // End of section 8
                    sectionTimer.stop(8);

                } while (!converged);

            } // End of parallel section

            timer.stop();

            // Report Times
            ASKAPLOG_INFO_STR(decmtbflogger,
                              "Time for minor cycles: "<< timer.elapsedTime()<<" sec");
            sectionTimer.summary();

            ASKAPLOG_INFO_STR(decmtbflogger, "Performed Multi-Term BasisFunction CLEAN for "
                                  << this->state()->currentIter() << " iterations");
            ASKAPLOG_INFO_STR(decmtbflogger, this->control()->terminationString());

            } else {
              ASKAPLOG_INFO_STR(decmtbflogger,
                  "Bypassed Multi-Term BasisFunction CLEAN due to 0 iterations in the setup");
            }

        } // End of many iterations function

        // the maxPos, maxVal and maxValScaled arguments are shared state variables, the return values are not used
        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::chooseComponent(uInt& optimumBase, IPosition& absPeakPos, 
            T& absPeakVal, T& absPeakValScaled, bool firstCycle, const std::vector<std::vector<uInt>>&highPixels, 
            askap::utils::SectionTimer& sectionTimer, IPosition& maxPos, T& maxVal, T& maxValScaled, 
            const Matrix<T>& weights, Matrix<T>& negchisq, Vector<Matrix<T>>& coefficients)
        {
            const uInt nBases(itsBasisFunction->numberBases());

            for (uInt base = 0; base < nBases; base++) {

                #pragma omp single
                {
                    maxPos = 0;
                    maxVal = 0.0;
                    maxValScaled = 0.0;
                }
                bool haveMask = weights.size()>0;

                // We implement various approaches to finding the peak. The first is the cheapest
                // and evidently the best (according to Urvashi).

                // Look for the maximum in term=0 for this base
                if (itsSolutionType == "MAXBASE") {

                    // Section 1 Timer
                    sectionTimer.start(1);

                    const Matrix<T>& res =itsResidualBasis(base)(0);
                    // initialise list of pixels depending on mode we're in
                    const bool deepClean = this->control()->deepCleanMode();
                    const bool useHighPixels = !deepClean && itsUsePixelLists && !firstCycle;
                    const std::vector<uInt>& pixels (deepClean ? 
                        std::vector<uInt>(itsScalePixels[base].begin(),itsScalePixels[base].end()) :
                        ( useHighPixels ? highPixels[base] : std::vector<uInt>()));
                    const uInt increment = itsUseIncrements && base > 0 ? 1 << (base-1) : 1;
                    if (!(deepClean||useHighPixels) || pixels.size()>0) {
                        absMaxPos(maxVal,maxValScaled,maxPos,res,weights,pixels,itsNoiseMap,itsNoiseBoxSize,increment);
                    }
                    // In performing the search for the peak across bases, we want to take into account
                    // the SNR so we normalise out the coupling matrix for term=0 to term=0.
                    #pragma omp single
                    {
                        T couplingFactor = sqrt(itsCouplingMatrix(base)(0, 0));
                        ASKAPDEBUGASSERT(couplingFactor > 0);
                        maxVal /= couplingFactor;
                        maxValScaled /= couplingFactor;
                    }
                    sectionTimer.stop(1);

                } else if (itsSolutionType == "MAXCHISQ") {

                    // section 2
                    sectionTimer.start(2);

                    for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {

                        #pragma omp single
                        {
                            coefficients(term1).resize(this->dirty(0).shape().nonDegenerate());
                            coefficients(term1).set(T(0.0));
                            ASKAPDEBUGASSERT(coefficients(term1).contiguousStorage());
                        }

                        for (uInt term2 = 0; term2 < this->nTerms(); ++term2) {
                            T* coeff_pointer = coefficients(term1).data();
                            const T* res_pointer = itsResidualBasis(base)(term2).data();
                            #pragma omp for schedule(static)
                            for (size_t index = 0; index < coefficients(term1).size(); index++) {
                                coeff_pointer[index] += res_pointer[index] *
                                       T(itsInverseCouplingMatrix(base)(term1,term2));
                            }
                        }
                    } // End of for loop over terms

                    sectionTimer.stop(2);

                    sectionTimer.start(3);
                    #pragma omp single
                    {
                        negchisq.resize(this->dirty(0).shape().nonDegenerate());
                        negchisq.set(T(0.0));
                        ASKAPDEBUGASSERT(negchisq.contiguousStorage());
                    }

                    T* negchisq_pointer = negchisq.data();
                    for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                        const T* coeff_pointer = coefficients(term1).data();
                        const T* res_pointer = itsResidualBasis(base)(term1).data();
                        #pragma omp for schedule(static)
                        for (size_t index = 0; index < negchisq.size(); index++) {
                            negchisq_pointer[index] += coeff_pointer[index]*res_pointer[index];
                        }
                    }
                    // initialise list of pixels depending on mode we're in
                    const std::vector<uInt>& pixels (this->control()->deepCleanMode() ? 
                        std::vector<uInt>(itsScalePixels[base].begin(),itsScalePixels[base].end()) :
                        std::vector<uInt>());

                    absMaxPos(maxVal,maxValScaled,maxPos,negchisq,weights,pixels,itsNoiseMap,itsNoiseBoxSize,1);

                    // End of section 3
                    sectionTimer.stop(3);
                } // End of else decision

                #pragma omp single
                {
                    // We use the maxVal to find the optimum base
                    if (abs(maxValScaled) > absPeakValScaled) {
                            optimumBase = base;
                            absPeakVal = abs(maxVal);
                            absPeakValScaled = abs(maxValScaled);
                            absPeakPos = maxPos;
                    }
                }

            } // End of iteration over number of bases

        }


        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::fillHighPixelList(std::vector<std::vector<uInt>>&highPixels, const Matrix<T>& weights)
        {
            const T level = this->control()->level(*(this->state()),itsPixelListTolerance);
            const bool haveMask = weights.size()>0;
            const uInt nBases = highPixels.size();
            // Added code to limit the number of high pixels collected - above 1e5 things get slow
            // no more than 10 x nIter pixels to be collected per base
            const uInt upperLimit = itsPixelListNPixRange[1] * this->control()->targetIter();
            // but try to get at least 2 x nIter pixels
            const uInt lowerLimit = itsPixelListNPixRange[0] * this->control()->targetIter();
            #pragma omp for schedule(static)
            for (uInt base = 0; base < nBases; base++) {
                const uInt increment = itsUseIncrements && base > 0 ? 1 << (base-1) : 1;
                const Matrix<T>& res = itsResidualBasis(base)(0);
                // get a quick estimate of the rms using 1% of pixels
                ASKAPDEBUGASSERT(res.nrow()>10 && res.ncolumn()>10);
                const float sigma = 1.48f * madfm(res(Slice(0,res.nrow()/10,10),Slice(0,res.ncolumn()/10,10)));
                ASKAPDEBUGASSERT(res.contiguousStorage());
                ASKAPDEBUGASSERT(!haveMask || weights.contiguousStorage());
                const T* pRes = res.data();
                const T* pMask = weights.data();
                const uInt n = res.size();
                const uInt nRow = res.nrow();
                // check we don't overflow uInt
                ASKAPDEBUGASSERT(n==res.size());
                std::vector<uInt>& pixels = highPixels[base];
                // scale the level down for larger scales, but not below n sigma
                const float cutoff = max(itsPixelListNSigma*sigma,sqrt(itsCouplingMatrix(base)(0, 0)) * level);
                const float prev = this->state()->initialObjectiveFunction();
                float step = 1.0;
                float scale = 0.0;
                bool tooManyPixels = true;
                bool tooFewPixels = false;
                float trialCutoff = cutoff;
                while (tooManyPixels || tooFewPixels) {
                    pixels.clear();
                    float trialCutoff = cutoff + scale * (prev-cutoff);
                    if (itsScalePixels.size()>0) {
                        auto it = itsScalePixels[base].begin();
                        // add pixels to the list if they are above the cutoff or cleaned before
                        // uses the fact that itsScalePixels[base] is a sorted set
                        for (uInt j = 0; j < n; j++ ) {
                            // skip pixels: use the increment in x and y
                            const bool skip = increment!=1  && (((j % nRow) % increment) != 0 ||
                                ((j / nRow) % increment) != 0);
                            const T val = (skip ? 0 : (haveMask ? abs(pRes[j] * pMask[j]) : abs(pRes[j])));
                            const bool doInc = (it != itsScalePixels[base].end()) && (*it == j);
                            if (val > trialCutoff || doInc) {
                                pixels.push_back(j);
                            }
                            if (doInc) {
                                ++it;
                            }
                        }
                    } else {
                        // add pixels above the cutoff to the list
                        for (uInt j = 0; j < n; j++ ) {
                            // skip pixels: use the increment in x and y
                            const bool skip = increment!=1  && (((j % nRow) % increment) != 0 ||
                                ((j / nRow) % increment) != 0);
                            const T val = (skip ? 0 : (haveMask ? abs(pRes[j] * pMask[j]) : abs(pRes[j])));
                            if (val > trialCutoff) {
                                pixels.push_back(j);
                            }
                        }
                    }
                    // Check if we have way too many components (>>maxiter) and change the cutoff until we don't
                    // We do a binary search to land between tooFewPixels and tooManyPixels
                    // Except if we never have too many
                    int nPixels = pixels.size();
                    if (itsScalePixels.size()>0) nPixels -= itsScalePixels[base].size();
                    tooManyPixels = nPixels > upperLimit;
                    step /= 2;
                    if (tooManyPixels) {
                        scale += step;
                    }
                    tooFewPixels = nPixels < lowerLimit && scale > 0.;
                    if (tooFewPixels) {
                        scale -= step;
                    }
                }
                ASKAPLOG_DEBUG_STR(decmtbflogger,"Base "<<base<<" is using "<<highPixels[base].size()<<" pixels above "<<trialCutoff <<" sigma = "<<sigma);
            }
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::getResidualAndPSFSlice(const IPosition& absPeakPos,
            IPosition& shape, IPosition& residualStart, IPosition& psfStart)
        {
            IPosition residualShape(this->dirty(0).shape().getFirst(2));
            IPosition psfShape(itsBasisFunction->shape().getFirst(2));

            IPosition residualEnd(2, 0);

            const IPosition peakPSFPos = this->getPeakPSFPosition();
            ASKAPDEBUGASSERT(peakPSFPos.size() >= 2);
            // work out the array sections we need
            for (uInt dim = 0; dim < 2; dim++) {
                residualStart(dim) = max(0, Int(absPeakPos(dim) - psfShape(dim) / 2));
                residualEnd(dim) = min(Int(absPeakPos(dim) + psfShape(dim) / 2 - 1), Int(residualShape(dim) - 1));
                // Now we have to deal with the PSF. Here we want to use enough of the
                // PSF to clean the residual image.
                psfStart(dim) = max(0, Int(peakPSFPos(dim) - (absPeakPos(dim) - residualStart(dim))));
            }
            shape = residualEnd - residualStart + 1; // +1 added
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::addComponentToModel(const Vector<T>& peakValues,
            const IPosition& shape, const IPosition& resStart, const IPosition& psfStart,
            const uInt optimumBase, Matrix<T>& model)
        {
            // Add to model
            const Matrix<T>& basisFunc = itsBasisFunction->basisFunction(optimumBase);
            // We loop over all terms for the optimum base and ignore those terms with no flux
            for (uInt term = 0; term < this->nTerms(); ++term) {
                if (abs(peakValues(term)) > 0.0) {
                    const T amp = this->control()->gain() * peakValues(term);
                    #pragma omp single
                    {
                        model.reference(this->model(term).nonDegenerate());
                        itsTermBaseFlux(optimumBase)(term) += amp;
                    }
                    #pragma omp for schedule(static)
                    for (uInt j = 0; j < shape(1); j++ ) {
                        T* pModel = &model(resStart(0), resStart(1) + j);
                        const T* pBfn = &basisFunc(psfStart(0), psfStart(1) + j);
                        for (uInt i = 0; i < shape(0); i++ ) {
                            pModel[i] += amp * pBfn[i];
                        }
                    }
                }
            }
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::subtractPSF(const Vector<T>& peakValues,
            const IPosition& shape, const IPosition& resStart, const IPosition& psfStart,
            uInt optimumBase)
        {
            const IPosition psfShape(itsBasisFunction->shape().getFirst(2));
            const uInt nBases(itsBasisFunction->numberBases());

            // Subtract PSFs, including base-base crossterms
            for (uInt term1 = 0; term1 < this->nTerms(); term1++) {
                for (uInt term2 = 0; term2 < this->nTerms(); term2++) {
                    if (abs(peakValues(term2)) > 0.0) {
                        const T amp = this->control()->gain() * peakValues(term2);
                        for (uInt base = 0; base < nBases; base++) {
                            // optimise the following code
                            // itsResidualBasis(base)(term1)(residualSlicer) -=
                            //     control()->gain() * peakValues(term2) *
                            //     itsPSFCrossTerms(base, optimumBase)(term1, term2)(psfSlicer);
                            Matrix<T>& res = itsResidualBasis(base)(term1);
                            const Matrix<T>& psf = itsPSFCrossTerms(base, optimumBase)(term1, term2);
                            #pragma omp for schedule(static)
                            for (uInt j = 0; j < shape(1); j++) {
                                T* pRes = &res(resStart(0), resStart(1) + j);
                                const T* pPsf = &psf(psfStart(0), psfStart(1) + j);
                                for (uInt i = 0; i < shape(0); i++) {
                                    pRes[i] -= amp * pPsf[i];
                                }
                            }
                        }
                    }
                }
            }

        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::subtractPSFPixels(const Vector<T>& peakValues,
            const IPosition& peakPos, uInt optimumBase, bool useHighPixels,
            const std::vector<std::vector<uInt>>& highPixels)
        {
            const IPosition psfShape(itsBasisFunction->shape().getFirst(2));
            const uInt nBases(itsBasisFunction->numberBases());
            // Subtract PSFs, including base-base crossterms
            for (uInt term1 = 0; term1 < this->nTerms(); term1++) {
                for (uInt term2 = 0; term2 < this->nTerms(); term2++) {
                    if (abs(peakValues(term2)) > 0.0) {
                        const T amp = this->control()->gain() * peakValues(term2);
                        for (uInt base = 0; base < nBases; base++) {
                            Matrix<T>& res = itsResidualBasis(base)(term1);
                            const Matrix<T>& psf = itsPSFCrossTerms(base, optimumBase)(term1, term2);
                            // loop over pixels, work out if in range, if so, subtract psf
                            const std::vector<uInt>& pixels = (useHighPixels ? highPixels[base] :
                                std::vector<uInt>(itsScalePixels[base].begin(),itsScalePixels[base].end()));
                            const uInt n = pixels.size();
                            const uInt nrow = res.nrow();
                            #pragma omp for schedule(static)
                            for (uInt k = 0; k < n; k++) {
                                const uInt x = pixels[k] % nrow;
                                const uInt y = pixels[k] / nrow;
                                // need to use long here, xpsf & ypsf go negative
                                const long xpsf = x - peakPos(0) + psfShape(0)/2;
                                const long ypsf = y - peakPos(1) + psfShape(1)/2;
                                if (xpsf >= 0 && xpsf < psfShape(0) &&
                                    ypsf >= 0 && ypsf < psfShape(1)) {
                                    res(x,y) -= amp * psf(xpsf,ypsf);
                                }
                            }
                        }
                    }
                }
            }
        }

        template<class T, class FT>
        bool DeconvolverMultiTermBasisFunction<T, FT>::deconvolve()
        {
            // This is the parallel version of deconvolve using ManyIterations()
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::deconvolve");
            initialise();
            askap::utils::Timer timer;
            timer.start();
            ManyIterations();
            timer.stop();
            finalise();
            ASKAPLOG_INFO_STR(decmtbflogger, "Time Required: "<< timer.elapsedTime());
            // signal failure and finish the major cycles if we started to diverge
            return (this->control()->terminationCause() != DeconvolverControl<T>::DIVERGED);
        }

        /// @brief export the scale mask
        /// @detail Access the scale mask used during deconvolution, this is a bitmask
        /// where a bit is set if the corresponding scale was used for that pixel
        /// @param[in]scaleMask a Matrix<uInt> with bitmask of scales for each pixel
        template<class T, class FT>
        const Matrix<T> DeconvolverMultiTermBasisFunction<T, FT>::scaleMask()
        {
            Matrix<T> scaleMask(itsResidualBasis[0][0].shape().getFirst(2),T(0));
            for (int base=0; base<itsScalePixels.size(); base++) {
                const uInt scaleBit = 1<<base;
                for (const uInt& pixel : itsScalePixels[base]) {
                    scaleMask.data()[pixel] = static_cast<uInt>(scaleMask.data()[pixel]) | scaleBit;
                }
            }
            return scaleMask;
        }

        /// @brief import initial scale mask
        /// @detail Load an initial scale mask to use in the deconvolution. It is up to the user
        /// to make sure the number (<=24 for float) and size of the scales matches between deconvolution runs
        /// Clean will only look for components on a particular scale at pixels where the corresponding bit is set
        /// @param[in]scaleMask a Matrix<uInt> with bitmask of scales for each pixel
        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::setScaleMask(const Matrix<T>& scaleMask)
        {
            // Option to subset the input scaleMask to the size of the dirty image
            Matrix<T> mask;
            if (this->dirty(0).shape() == scaleMask.shape()) {
                mask = scaleMask;
            } else {
                ASKAPLOG_INFO_STR(decmtbflogger,"Shape mismatch of dirty image and scale mask, using central part of scale mask");
                IPosition blc = scaleMask.shape()/2 - this->dirty(0).shape()/2;
                IPosition trc = scaleMask.shape()/2 + this->dirty(0).shape()/2 - 1;
                mask = scaleMask(blc, trc);
            }
            ASKAPCHECK(itsBasisFunction, "Basis function not initialised");
            const uInt nBases(itsBasisFunction->numberBases());
            ASKAPCHECK(nBases <= itsMaxScales,"Scalemask only supports up to "<<itsMaxScales<<" scales");
            ASKAPASSERT(mask.contiguousStorage());
            itsScalePixels.resize(nBases);
            for (uInt i=0; i < mask.size(); i++) {
                const uInt val = static_cast<uInt>(mask.data()[i]);
                // Need to deal with multiple scales at same pixel
                if (val > 0) {
                    for (uInt scale = 0; scale < nBases; scale++) {
                        const uInt scaleBit = 1<<scale;
                        if (val & scaleBit) {
                            // this scale is present
                            itsScalePixels[scale].insert(i);
                        }
                        if (val == scaleBit) {
                            // no other scales to check
                            break;
                        }
                    }
                }
            }
            // we start deep cleaning straight away if scale mask is set
            this->control()->setDeepCleanMode();
            ASKAPLOG_INFO_STR(decmtbflogger, "Starting deep cleaning phase with provided scale mask");
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::releaseMemory()
        {
            DeconvolverBase<T, FT>::releaseMemory();
            uInt memory = 0;
            const uInt nBases(itsBasisFunction->numberBases());
            for (uInt base = 0; base < nBases; base++) {
                for (uInt term = 0; term < this->nTerms(); term++) {
                    memory += sizeof(imtype) * itsResidualBasis(base)(term).size();
                }
            }
            itsResidualBasis.resize();
            ASKAPLOG_DEBUG_STR(decbaselogger,"DeconvolverMultiTermBasisFunction released "<<memory/1024/1024<<" MB from residualBasis");
            ASKAPCHECK(itsBasisFunction, "Basis function not initialised");
            memory = sizeof(imtype) * itsBasisFunction->allBasisFunctions().size();
            itsBasisFunction->allBasisFunctions().resize();
            ASKAPLOG_DEBUG_STR(decbaselogger,"DeconvolverMultiTermBasisFunction released "<<memory/1024/1024<<" MB from basisfunctions");

        }

    }
}
// namespace synthesis
// namespace askap
