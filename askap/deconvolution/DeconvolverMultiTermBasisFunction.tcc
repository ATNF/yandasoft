/// @file DeconvolverMultiTermBasisFunction.tcc
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

#include <string>
#include <askap/askap/AskapLogging.h>
#include <casacore/casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MaskArrMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/scimath/Mathematics/MatrixMathLA.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/profile/AskapProfiler.h>
ASKAP_LOGGER(decmtbflogger, ".deconvolution.multitermbasisfunction");

#include <askap/deconvolution/DeconvolverMultiTermBasisFunction.h>
#include <askap/deconvolution/MultiScaleBasisFunction.h>
#include <askap/scimath/utils/OptimizedArrayMathUtils.h>
#include <mpi.h>

#include <askap/scimath/fft/FFT2DWrapper.h>

#ifdef USE_OPENACC
template<class T>
ACCManager<T>::ACCManager() {
    ASKAPLOG_INFO_STR(decmtbflogger,"In OPEN ACC mode instantiating manager");
}
template<class T>
ACCManager<T>::~ACCManager() {
    ASKAPLOG_INFO_STR(decmtbflogger,"Destructor FIXME delete the memory if required");
 // the residuals
    size_t nimages= nBases*nTerms;
    for (int i = 0; i< nimages ; i++) {
        T* tomove = (T *) residuals[i];
        #pragma acc exit data delete(tomove[0:npixels])
    }
    // the masks
    //
    for (int i = 0; i<nBases;i++) {
        T* tomove = (T *) masks[i];
        #pragma acc exit data delete(tomove[0:npixels])
    }

    T* tomove = (T *) maskToUse;
    #pragma acc exit data delete(tomove[0:npixels])

    tomove = (T *) weight;
    #pragma acc exit data delete(tomove[0:npixels])


}
template <class T>
void ACCManager<T>::CopyToDevice() {

    // the residuals
    size_t nimages= nBases*nTerms;
    for (int i = 0; i< nimages ; i++) {
        T* tomove = (T *) residuals[i];
        #pragma acc enter data copyin(tomove[0:npixels])
    }
    // the masks
    //
    for (int i = 0; i<nBases;i++) {
        T* tomove = (T *) masks[i];
        #pragma acc enter data copyin(tomove[0:npixels])
    }

    T* tomove = (T *) maskToUse;
    #pragma acc enter data copyin(tomove[0:npixels])

    tomove = (T *) weight;
    #pragma acc enter data copyin(tomove[0:npixels])

}

template <class T>
void ACCManager<T>::UpdateMask(int base) {
    T * basemask = (T *) masks[base];
    #pragma acc parallel loop present(maskToUse,basemask,weight)
    for (int i=0;i<npixels;i++) {
       maskToUse[i] = weight[i]*basemask[i];
    }
}
template <class T>
void ACCManager<T>::InitMask(int base) {

    #pragma acc parallel loop present(maskToUse,weight)
    for (int i=0;i<npixels;i++) {
       maskToUse[i] = weight[i];
    }
}
#endif

namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver based on the BasisFunction Clean
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a residual image, psf optionally using a weights image.
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverMultiTermBasisFunction<Double, DComplex>
        /// @ingroup Deconvolver


        template<class T>
        void absMaxPosOMP(T& maxVal, IPosition& maxPos, const Matrix<T>& im) {

            // Set Shared values
            maxVal = T(0.0);
            // Create and shared private values
            T maxVal_private(0.0);
            IPosition maxPos_private(2,0);
            const uInt ncol = im.ncolumn();
            const uInt nrow = im.nrow();
            #pragma omp for schedule(static)
            for (uInt j = 0; j < ncol; j++ ) {
                const T* pIm = &im(0,j);
                for (uInt i = 0; i < nrow; i++ ) {
                    T val = abs(*pIm++);
                    if (val > maxVal_private) {
                        maxVal_private = val;
                        maxPos_private(0) = i;
                        maxPos_private(1) = j;
                    }
                }
            }
            // Update shared max values and positions
            #pragma omp critical
            {
                if (maxVal_private > maxVal) {
                    maxVal = maxVal_private;
                    maxPos = maxPos_private;
                }
            }
            #pragma omp barrier
        }


        template<class T>
        void absMaxPosMaskedOMP(T& maxVal, IPosition& maxPos, const Matrix<T>& im, const Matrix<T>& mask) {

            // Set Shared Values
            maxVal = T(0.0);
            // Set Private Values
            T maxVal_private(0.0);
            IPosition maxPos_private(2,0);
            const uInt ncol = mask.ncolumn();
            const uInt nrow = mask.nrow();

            #pragma omp for schedule(static)
            for (uInt j = 0; j < ncol; j++ ) {
                const T* pIm = &im(0,j);
                const T* pMask = &mask(0,j);
                for (uInt i = 0; i < nrow; i++ ) {
                        T val = abs(*pIm++ * *pMask++);
                        if (val > maxVal_private) {
                            maxVal_private = val;
                            maxPos_private(0) = i;
                            maxPos_private(1) = j;
                        }
                }
            }
            #pragma omp critical
            {
                if (maxVal_private > maxVal) {
                    maxVal = maxVal_private;
                    maxPos = maxPos_private;
                }
            }
            #pragma omp barrier
        }

        template<class T, class FT>
        DeconvolverMultiTermBasisFunction<T, FT>::DeconvolverMultiTermBasisFunction(Vector<Array<T> >& dirty,
                Vector<Array<T> >& psf,
                Vector<Array<T> >& psfLong)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf), itsDirtyChanged(True), itsBasisFunctionChanged(True),
                itsSolutionType("MAXCHISQ"), itsDecoupled(false)
        {
            ASKAPLOG_DEBUG_STR(decmtbflogger, "There are " << this->nTerms() << " terms to be solved");

            ASKAPCHECK(psfLong.nelements() == (2*this->nTerms() - 1), "Long PSF vector has incorrect length " << psfLong.nelements());
            this->itsPsfLongVec.resize(2*this->nTerms() - 1);

            for (uInt term = 0; term < (2*this->nTerms() - 1); ++term) {
                ASKAPCHECK(psfLong(term).nonDegenerate().shape().nelements() == 2, "PSF(" << term << ") has too many dimensions " << psfLong(term).shape());
                this->itsPsfLongVec(term).reference(psfLong(term).nonDegenerate());
            }

        };

        template<class T, class FT>
        DeconvolverMultiTermBasisFunction<T, FT>::DeconvolverMultiTermBasisFunction(Array<T>& dirty,
                Array<T>& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf), itsDirtyChanged(True), itsBasisFunctionChanged(True),
                itsSolutionType("MAXCHISQ"), itsDecoupled(false)
        {
            ASKAPLOG_DEBUG_STR(decmtbflogger, "There is only one term to be solved");
            this->itsPsfLongVec.resize(1);
            this->itsPsfLongVec(0) = psf;
        };

        template<class T, class FT>
        DeconvolverMultiTermBasisFunction<T, FT>::~DeconvolverMultiTermBasisFunction()
        {
        };

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::setSolutionType(String sol)
        {
            itsSolutionType = sol;
        };

        template<class T, class FT>
        const String DeconvolverMultiTermBasisFunction<T, FT>::solutionType()
        {
            return itsSolutionType;
        };
        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::setDecoupled(Bool decoupled)
        {
            itsDecoupled=decoupled;
        };

        template<class T, class FT>
        const Bool DeconvolverMultiTermBasisFunction<T, FT>::decoupled()
        {
            return itsDecoupled;
        };

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::setBasisFunction(boost::shared_ptr<BasisFunction<T> > bf)
        {
            this->itsBasisFunction = bf;
            this->itsBasisFunctionChanged = True;
        };

        template<class T, class FT>
        boost::shared_ptr<BasisFunction<T> > DeconvolverMultiTermBasisFunction<T, FT>::basisFunction()
        {
            return itsBasisFunction;
        };

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::updateDirty(Array<T>& dirty, casacore::uInt term)
        {
            DeconvolverBase<T, FT>::updateDirty(dirty, term);
            this->itsDirtyChanged = True;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::updateDirty(Vector<Array<T> >& dirtyVec)
        {
            DeconvolverBase<T, FT>::updateDirty(dirtyVec);
            this->itsDirtyChanged = True;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::configure(const LOFAR::ParameterSet& parset)
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::configure");
            DeconvolverBase<T, FT>::configure(parset);

            // Make the basis function
            std::vector<float> defaultScales(3);
            defaultScales[0] = 0.0;
            defaultScales[1] = 10.0;
            defaultScales[2] = 30.0;
            std::vector<float> scales = parset.getFloatVector("scales", defaultScales);
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Constructing Multiscale basis function with scales "
                                   << scales);
            Bool orthogonal = parset.getBool("orthogonal", false);
            if (orthogonal) {
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Multiscale basis functions will be orthogonalised");
            }
            itsBasisFunction = BasisFunction<Float>::ShPtr(new MultiScaleBasisFunction<Float>(scales,
                               orthogonal));

            String solutionType = parset.getString("solutiontype", "MAXCHISQ");
            itsDecoupled = parset.getBool("decoupled", false);
            if (itsDecoupled) {
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Using decoupled residuals");
            }

            if (solutionType == "MAXBASE") {
                itsSolutionType = solutionType;
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Component search to maximise over bases");
            } else if (solutionType == "MAXTERM0") {
                itsSolutionType = solutionType;
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Component search to maximise Taylor term 0 over bases");
            } else {
                itsSolutionType = "MAXCHISQ";
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Component search to find maximum in chi-squared");
            }
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::finalise()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::finalise");
            this->updateResiduals(this->itsModel);

            // debug message
            for (uInt base = 0; base < itsTermBaseFlux.nelements(); base++) {
                for (uInt term = 0; term < itsTermBaseFlux(base).nelements(); term++) {
                    ASKAPLOG_DEBUG_STR(decmtbflogger, "   Term(" << term << "), Base(" << base
                                           << "): Flux = " << itsTermBaseFlux(base)(term));
                }
            }

            // info message
            for (uInt base = 0; base < itsTermBaseFlux.nelements(); base++) {
              ASKAPLOG_INFO_STR(decmtbflogger,"Total flux for scale "<<base<<" : "<<itsTermBaseFlux(base)(0));
            }
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialiseForBasisFunction(bool force)
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialiseForBasisFunction");
            if (!force && !this->itsBasisFunctionChanged) {
                return;
            }

            ASKAPLOG_DEBUG_STR(decmtbflogger,
                               "Updating Multi-Term Basis Function deconvolver for change in basis function");

            IPosition subPsfShape(this->findSubPsfShape());

            // Use a smaller size for the psfs if specified.
            this->itsBasisFunction->initialise(subPsfShape);

            ASKAPLOG_DEBUG_STR(decmtbflogger, "Initialising for PSFs: shape = " << subPsfShape);
            initialisePSF();

            itsBasisFunctionChanged = False;
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialise()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialise");
            DeconvolverBase<T, FT>::initialise();

            // Initialise residuals
            initialiseResidual();

            // Initialise masks
            initialiseMask();

#ifdef USE_OPENACC
            itsACCManager.CopyToDevice();
#endif
            // Force change in basis function
            initialiseForBasisFunction(true);

            //this->state()->resetInitialObjectiveFunction();
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialiseResidual()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialiseResidual");

            if (!this->itsDirtyChanged) {
                return;
            }

            // Initialise the basis function for residual calculations.
            ASKAPCHECK(this->itsBasisFunction, "Basis function not initialised");
            this->itsBasisFunction->initialise(this->dirty(0).shape());

            const uInt nBases(this->itsBasisFunction->numberBases());
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of basis functions "
                                   << this->itsBasisFunction->shape()<<" number of bases "<<nBases);

            itsResidualBasis.resize(nBases);
            for (uInt base = 0; base < nBases; base++) {
                itsResidualBasis(base).resize(this->nTerms());
            }

            // Calculate residuals convolved with bases [nx,ny][nterms][nbases]

            ASKAPLOG_INFO_STR(decmtbflogger,
                              "Calculating convolutions of residual images with basis functions");
            //const double start_time = MPI_Wtime();
            const time_t start_time = time(0);
            // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
            // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
            scimath::FFT2DWrapper<FT> fftWrapper(true);
            for (uInt base = 0; base < nBases; base++) {
                 // Calculate transform of basis function [nx,ny,nbases]
                 const Matrix<T> bfRef(this->itsBasisFunction->basisFunction(base));
                 Matrix<FT> basisFunctionFFT(bfRef.shape().nonDegenerate(2), 0.);
                 casacore::setReal(basisFunctionFFT, bfRef);
                 //scimath::fft2d(basisFunctionFFT, true);
                 fftWrapper(basisFunctionFFT, true);

                 for (uInt term = 0; term < this->nTerms(); term++) {

                    // Calculate transform of residual image
                    Matrix<FT> residualFFT(this->dirty(term).shape().nonDegenerate(), 0.);
                    casacore::setReal(residualFFT, this->dirty(term).nonDegenerate());
                    //scimath::fft2d(residualFFT, true);
                    fftWrapper(residualFFT, true);

                    // Calculate product and transform back
                    ASKAPASSERT(basisFunctionFFT.shape().conform(residualFFT.shape()));

                    //the following line is equivalent to the optimised version called below
                    //residualFFT *= conj(basisFunctionFFT);
                    utility::multiplyByConjugate(residualFFT, basisFunctionFFT);

                    //scimath::fft2d(residualFFT, false);
                    fftWrapper(residualFFT, false);

                    // temporary object is ok here because we do an assignment to uninitialised array later on
                    Matrix<T> work(real(residualFFT));

                    ASKAPLOG_DEBUG_STR(decmtbflogger, "Basis(" << base
                                           << ")*Residual(" << term << "): max = " << max(work)
                                           << " min = " << min(work));

                    this->itsResidualBasis(base)(term) = work;
                }
            }
            //const double end_time = MPI_Wtime();
            const time_t end_time = time(0);
            ASKAPLOG_INFO_STR(decmtbflogger,
                              "Time to calculate residual images * basis functions: "<<end_time-start_time<<" sec");
#ifdef USE_OPENACC
            itsACCManager.nBases = nBases;
            itsACCManager.nTerms = this->nTerms();
            itsACCManager.npixels = this->itsResidualBasis(0)(0).nelements();
            itsACCManager.nrows = this->itsResidualBasis(0)(0).shape()(1);
            itsACCManager.ncols = this->itsResidualBasis(0)(0).shape()(0);

            itsACCManager.residuals = new uInt64[itsACCManager.nBases*itsACCManager.nTerms];
            itsACCManager.deleteResiduals =  new casacore::Bool[itsACCManager.nBases*itsACCManager.nTerms];

            size_t idx = 0;
            for (uInt base = 0; base < nBases; base++) {
                // Calculate transform of residual images [nx,ny,nterms]
                for (uInt term = 0; term < this->nTerms(); ++term) {
                    Bool deleteIt;
                    itsACCManager.residuals[idx] = (uInt64) this->itsResidualBasis(base)(term).getStorage(deleteIt);
                    itsACCManager.deleteResiduals[idx] = deleteIt;
                    idx++;
                }
            }
#endif

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
            if (this->itsMask.nelements()>0) {
                return;
            }
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Initialising deep clean masks");

            ASKAPCHECK(this->itsBasisFunction, "Basis function not initialised");

            uInt nBases(this->itsBasisFunction->numberBases());

            this->itsMask.resize(nBases);

#ifdef USE_OPENACC

            Bool deleteIt;
            itsACCManager.tmpMask = this->itsWeight(0).nonDegenerate();
            itsACCManager.masks = new uInt64[nBases];
            itsACCManager.deleteMasks = new casacore::Bool[nBases];
            itsACCManager.weight = itsACCManager.tmpMask.getStorage(deleteIt);


#endif
            for (uInt base = 0; base < nBases; base++) {
                this->itsMask(base).resize(this->dirty(0).shape().nonDegenerate());
                this->itsMask(base).set(T(0.0));
#ifdef USE_OPENACC
                casacore::Bool deleteIt;
                itsACCManager.masks[base] = (uInt64) this->itsMask(base).getStorage(deleteIt);
                itsACCManager.deleteMasks[base] = deleteIt;
#endif

            }
#ifdef USE_OPENACC
            size_t npixels = this->itsMask(0).nelements();
            itsACCManager.maskToUse = new T[npixels];
#endif
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialisePSF()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialisePSF");

            if (!this->itsBasisFunctionChanged) {
                return;
            }

            ASKAPCHECK(this->itsBasisFunction, "Basis function not initialised");

            ASKAPLOG_DEBUG_STR(decmtbflogger,
                               "Updating Multi-Term Basis Function deconvolver for change in basis function");
            const IPosition subPsfShape(this->findSubPsfShape());

            Array<FT> work(subPsfShape);

            const uInt nBases(this->itsBasisFunction->numberBases());
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of basis functions "
                                   << this->itsBasisFunction->shape()<< " number of bases "<<nBases);

            const IPosition stackShape(this->itsBasisFunction->shape());

            // Now transform the basis functions. These may be a different size from
            // those in initialiseResidual so we don't keep either
            Cube<FT> basisFunctionFFT(stackShape,FT(0.));

            // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
            // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
            scimath::FFT2DWrapper<FT> fftWrapper(true);

            // do explicit loop over basis functions here (the original code relied on iterator in
            // fft2d and, therefore, low level representation of the basis function stack). This way
            // we have more control over the array structure and can transition to the more efficient order
            for (uInt base = 0; base < nBases; ++base) {
                 // casacore arrays have reference semantics, no copying occurs in the following
                 casacore::Matrix<FT> fftBuffer = basisFunctionFFT.xyPlane(base);
                 casacore::setReal(fftBuffer, this->itsBasisFunction->basisFunction(base));
                 //scimath::fft2d(fftBuffer, true);
                 fftWrapper(fftBuffer, true);
            }

            itsTermBaseFlux.resize(nBases);
            for (uInt base = 0; base < nBases; base++) {
                itsTermBaseFlux(base).resize(this->nTerms());
                itsTermBaseFlux(base) = 0.0;
            }

            const uInt nx(this->psf(0).shape()(0));
            const uInt ny(this->psf(0).shape()(1));

            const IPosition subPsfStart(2, (nx - subPsfShape(0)) / 2, (ny - subPsfShape(1)) / 2);
            //const IPosition subPsfEnd(2,(nx+subPsfShape(0))/2-1,(ny+subPsfShape(1))/2-1);
            //const IPosition subPsfStride(2,1,1);

            //Slicer subPsfSlicer(subPsfStart, subPsfEnd, subPsfStride, Slicer::endIsLast);
            Slicer subPsfSlicer(subPsfStart, subPsfShape);
            // check just in case
            ASKAPCHECK(subPsfSlicer.length() == subPsfShape, "Slicer selected length of " <<
                subPsfSlicer.length() << " is different from requested shape " << subPsfShape);

            this->validatePSF(subPsfSlicer);

            const casacore::IPosition subPsfPeak=this->getPeakPSFPosition().getFirst(2);
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Peak of PSF subsection at  " << subPsfPeak);
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of PSF subsection is " << subPsfShape);

            // Calculate XFR for the subsection only. We need all PSF's up to
            // 2*nTerms-1
            ASKAPCHECK(this->itsPsfLongVec.nelements() == (2*this->nTerms() - 1),
                "PSF long vector has wrong length " << this->itsPsfLongVec.nelements());

            // Calculate all the transfer functions
            Vector<Array<FT> > subXFRVec(2*this->nTerms() - 1);
            for (uInt term1 = 0; term1 < subXFRVec.nelements(); ++term1) {
                subXFRVec(term1).resize(subPsfShape);
                // rely on reference semantics of casa arrays 
                // MV: we can probably change subXFRVec to be a vector of matrices to reduce technical debt
                Matrix<FT> subXFRTerm1(subXFRVec(term1));
                subXFRTerm1.set(0.0);
                casacore::setReal(subXFRTerm1, this->itsPsfLongVec(term1).nonDegenerate()(subPsfSlicer));
                //scimath::fft2d(subXFRVec(term1), true);
                fftWrapper(subXFRTerm1, true);
                // we only need conjugated FT of subXFRVec (or real part of it, which doesn't change with conjugation), 
                // it is better to compute conjugation in situ now and don't do it on the fly later
                utility::conjugateComplexArray(subXFRTerm1);
            }
            // Calculate residuals convolved with bases [nx,ny][nterms][nbases]

            // the following line is the original code which we do now in an optimised way
            //const T normPSF = casacore::sum(casacore::real(subXFRVec(0))) / subXFRVec(0).nelements();
            const T normPSF = utility::realPartMean(subXFRVec(0));
            ASKAPLOG_DEBUG_STR(decmtbflogger, "PSF effective volume = " << normPSF);

            itsPSFCrossTerms.resize(nBases, nBases);
            for (uInt base = 0; base < nBases; base++) {
                for (uInt base1 = 0; base1 < nBases; base1++) {
                    itsPSFCrossTerms(base, base1).resize(this->nTerms(), this->nTerms());
                }
            }

            this->itsCouplingMatrix.resize(nBases);
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

                            //scimath::fft2d(work, false);
                            //use reference semantics to get the right interface, we can probably change the interface to matrix to reduce technical debt
                            Matrix<FT> workMtr(work);
                            fftWrapper(workMtr, false);

                            ASKAPLOG_DEBUG_STR(decmtbflogger, "Base(" << base1 << ")*Base(" << base2
                                                   << ")*PSF(" << term1 + term2
                                                   << "): max = " << max(real(work))
                                                   << " min = " << min(real(work))
                                                   << " centre = " << real(work(subPsfPeak)));
                            // Remember that casacore::Array reuses the same memory where possible so this
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
            this->itsInverseCouplingMatrix.resize(nBases);
            this->itsDetCouplingMatrix.resize(nBases);

            for (uInt base = 0; base < nBases; base++) {
                this->itsInverseCouplingMatrix(base).resize(this->nTerms(), this->nTerms());
                ASKAPLOG_INFO_STR(decmtbflogger, "Coupling matrix(" << base << ")="
                                       << this->itsCouplingMatrix(base).row(0));
                for (uInt term = 1; term < this->nTerms(); ++term) {
                    ASKAPLOG_INFO_STR(decmtbflogger, "                   "
                                       << this->itsCouplingMatrix(base).row(term));
                }
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Calculating matrix inverse by Cholesky decomposition");
                invertSymPosDef(this->itsInverseCouplingMatrix(base),
                                this->itsDetCouplingMatrix(base), this->itsCouplingMatrix(base));
                ASKAPLOG_INFO_STR(decmtbflogger, "Coupling matrix determinant(" << base << ") = "
                                       << this->itsDetCouplingMatrix(base));
                ASKAPLOG_INFO_STR(decmtbflogger, "Inverse coupling matrix(" << base
                                       << ")=" << this->itsInverseCouplingMatrix(base).row(0));
                for (uInt term = 1; term < this->nTerms(); ++term) {
                    ASKAPLOG_INFO_STR(decmtbflogger, "                           "
                                       << this->itsInverseCouplingMatrix(base).row(term));
                }
            }
            this->itsBasisFunctionChanged = False;
        }


        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::ManyIterations()
        {
            // Need to (re)set the subPsfShape and call validatePSF to get correct psf peak pixels as it gets unset by DeconvolverBase::initialise
            const IPosition subPsfShape(this->findSubPsfShape());
            const uInt nx(this->psf(0).shape()(0));
            const uInt ny(this->psf(0).shape()(1));
            const IPosition subPsfStart(2, (nx - subPsfShape(0)) / 2, (ny - subPsfShape(1)) / 2);
            const Slicer subPsfSlicer(subPsfStart, subPsfShape);
            this->validatePSF(subPsfSlicer);

            const uInt nBases(this->itsResidualBasis.nelements());
            IPosition absPeakPos(2, 0);
            T absPeakVal(0.0);
            float sumFlux;
            uInt optimumBase(0);
            Vector<T> peakValues(this->nTerms());
            Vector<T> maxValues(this->nTerms());
            Matrix<T> weights, mask, maskref;
            IPosition maxPos(2, 0);
            T maxVal(0.0);
            Bool haveMask;
            T norm;
            Vector<Array<T> > coefficients(this->nTerms());
            casa::Matrix<T> res, wt;
            Array<T> negchisq;
            casa::IPosition residualShape;
            casa::IPosition psfShape;
            bool isWeighted((this->itsWeight.nelements() > 0) &&
                (this->itsWeight(0).shape().nonDegenerate().conform(this->itsResidualBasis(0)(0).shape())));
            Vector<T> maxTermVals(this->nTerms());
            Vector<T> maxBaseVals(nBases);

            // Timers for analysis
            const int no_timers = 10;
            double TimerStart[no_timers], TimerStop[no_timers], Times[no_timers];
            memset(TimerStart, 0, sizeof(TimerStart));
            memset(TimerStop, 0, sizeof(TimerStop));
            memset(Times, 0, sizeof(Times));

		      	// Termination
		      	int converged;
            this->control()->maskNeedsResetting(true);

            if (this->control()->targetIter() != 0) {

            //const double start_time = MPI_Wtime();
            const time_t start_time = time(0);

            #pragma omp parallel
            {
                bool IsNotCont;

                // =============== Set weights =======================

                // Section 0
                #pragma omp single
                TimerStart[0] = MPI_Wtime();

                if (isWeighted) {

                    #pragma omp single
                    weights = this->itsWeight(0).nonDegenerate();

                    uInt ncol = weights.ncolumn();
                    uInt nrow = weights.nrow();
                    // Declare private versions of these
                    uInt i, j;

                    // Check weights for contiguity
                    if (!weights.contiguousStorage()) {
                        ASKAPLOG_WARN_STR(decmtbflogger, "weights (sec 0) is not contiguous\n");
                    }

                    if  (this->itsSolutionType == "MAXCHISQ") {
                        // square weights for MAXCHISQ
                        #pragma omp for schedule(static)
                        for (j = 0; j < ncol; j++ ) {
                            Vector<T> weightscol = weights.column(j);
                            T* pWeights = weightscol.getStorage(IsNotCont);
                            for (i = 0; i < nrow; i++ ) {
                                pWeights[i] = abs(*(pWeights+i) * (*(pWeights+i)));
                            }
                        }
                    }

                    if (nBases>1) {
                        // init scratch space for mask
                        #pragma omp single
                        mask = Matrix<T>(weights.shape(),0.0);
                        // Check mask base for contiguity
                        if (!mask.contiguousStorage()) {
                            ASKAPLOG_WARN_STR(decmtbflogger, "mask is not contiguous\n");
                        }
                    }

                }

                #pragma omp single
                { TimerStop[0] = MPI_Wtime(); Times[0] += (TimerStop[0]-TimerStart[0]); }

                // Commence cleaning iterations
                do {

                    // Reset peak pos
                    absPeakPos = 0;
                    // Reset peak Val
                    absPeakVal = 0.0;
                    // Reset optimum base
                    optimumBase = 0;

                    // =============== Set up deep cleaning mask =======================
                    // initialise a scratch space if mask needs to be reset each iter, otherwise just point at weights
                    if (isWeighted && this->control()->maskNeedsResetting()) {
                        uInt ncol = weights.ncolumn();
                        uInt nrow = weights.nrow();
                        if (this->control()->deepCleanMode()) {
                            if (nBases>1) {
                                #pragma omp single
                                maskref.reference(mask);
                            } else {
                                // only a single base, so multiple weights and mask now.
                                // reference the mask to the mask for base 0
                                #pragma omp single
                                maskref.reference(this->itsMask(0));
                                // multiply weights by the base 0 mask
                                #pragma omp for schedule(static)
                                for (uInt j = 0; j < ncol; j++ ) {
                                    Vector<T> weightscol = weights.column(j);
                                    T* pWeights = weightscol.getStorage(IsNotCont);
                                    Vector<T> maskcol = maskref.column(j);
                                    T* pMask = maskcol.getStorage(IsNotCont);
                                    for (uInt i = 0; i < nrow; i++ ) {
                                        T val = *(pWeights+i) * (*(pMask+i));
                                        pWeights[i] = val;
                                    }
                                }
                                // now reference the mask to the updated weigths
                                #pragma omp single
                                maskref.reference(weights);
                            }
                        } else {
                            // reference the mask to the weigths
                            #pragma omp single
                            maskref.reference(weights);
                        }
                        this->control()->maskNeedsResetting(false);
                    }

                    // =============== Choose Component =======================

                    for (uInt base = 0; base < nBases; base++) {

                        maxPos = 0;
                        maxVal = 0.0;

                        if (this->control()->deepCleanMode()) {

                            // Section 1 Timer
                            #pragma omp single
                            TimerStart[1] = MPI_Wtime();

                            if (isWeighted) {
                                // for a single base, the mask has already been set outside this loop

                                if (nBases>1) {
                                    uInt ncol = maskref.ncolumn();
                                    uInt nrow = maskref.nrow();
                                    Matrix<T> maskbase;
                                    maskbase.reference(this->itsMask(base));
                                    #pragma omp for schedule(static)
                                    for (uInt j = 0; j < ncol; j++ ) {
                                        Vector<T> weightscol = weights.column(j);
                                        T* pWeights = weightscol.getStorage(IsNotCont);
                                        Vector<T> maskbasecol = maskbase.column(j);
                                        T* pMaskBase = maskbasecol.getStorage(IsNotCont);
                                        Vector<T> maskcol = maskref.column(j);
                                        T* pMask = maskcol.getStorage(IsNotCont);
                                        for (uInt i = 0; i < nrow; i++ ) {
                                            pMask[i] = *(pWeights+i) * (*(pMaskBase+i));
                                        }
                                    }
                                }

                            } else {
                                #pragma omp single
                                maskref.reference(this->itsMask(base));
                            }

                            #pragma omp single
                            { TimerStop[1] = MPI_Wtime(); Times[1] += (TimerStop[1]-TimerStart[1]); }

                        }

                        #pragma omp single
                        haveMask = maskref.nelements()>0;

                        // We implement various approaches to finding the peak. The first is the cheapest
                        // and evidently the best (according to Urvashi).

                        // Look for the maximum in term=0 for this base
                        if (this->itsSolutionType == "MAXBASE") {

                            // Section 2 Timer
                            #pragma omp single
                            TimerStart[2] = MPI_Wtime();

                            #pragma omp single
                            res.reference(this->itsResidualBasis(base)(0));

                            if (haveMask) {
                                absMaxPosMaskedOMP(maxVal,maxPos,res,maskref);
                            } else {
                                absMaxPosOMP(maxVal,maxPos,res);
                            }

                            #pragma omp for schedule(static)
                            for (uInt term = 0; term < this->nTerms(); ++term) {
                                maxValues(term) = this->itsResidualBasis(base)(term)(maxPos);
                            }
                            // In performing the search for the peak across bases, we want to take into account
                            // the SNR so we normalise out the coupling matrix for term=0 to term=0.
                            #pragma omp single
                            {
                                norm = 1.0 / sqrt(this->itsCouplingMatrix(base)(0, 0));
                                maxVal *= norm;
                            }

                            #pragma omp single
                            { TimerStop[2] = MPI_Wtime(); Times[2] += (TimerStop[2]-TimerStart[2]); }

                        } else {  // Some other solver type than maxbase

                            // section 3
                            #pragma omp single
                            TimerStart[3] = MPI_Wtime();

                            for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {

                                #pragma omp single
                                {
                                    coefficients(term1).resize(this->dirty(0).shape().nonDegenerate());
                                    coefficients(term1).set(T(0.0));
                                }

                                for (uInt term2 = 0; term2 < this->nTerms(); ++term2) {
                                    T* coeff_pointer = coefficients(term1).getStorage(IsNotCont);
                                    T* res_pointer = (float*)this->itsResidualBasis(base)(term2).getStorage(IsNotCont);
                                    #pragma omp for schedule(static)
                                    for (int index = 0; index < coefficients(term1).nelements(); index++) {
                                        coeff_pointer[index] += res_pointer[index] *
                                               T(this->itsInverseCouplingMatrix(base)(term1,term2));
                                    }
                                }
                            } // End of for loop over terms

                            #pragma omp single
                            { TimerStop[3] = MPI_Wtime(); Times[3] += (TimerStop[3]-TimerStart[3]); }

                            if (this->itsSolutionType == "MAXTERM0") {

                                #pragma omp single
                                TimerStart[4] = MPI_Wtime();

                                #pragma omp single
                                res = coefficients(0);

                                if (haveMask) {
                                    absMaxPosMaskedOMP(maxVal, maxPos, res, maskref);
                                } else {
                                    absMaxPosOMP(maxVal, maxPos, res);
                                }

                                #pragma omp for schedule(static)
                                for (uInt term = 0; term < this->nTerms(); ++term) {
                                    maxValues(term) = coefficients(term)(maxPos);
                                }

                                #pragma omp single
                                { TimerStop[4] = MPI_Wtime(); Times[4] += (TimerStop[4]-TimerStart[4]); }

                            } else {
                                // MAXCHISQ
                                #pragma omp single
                                TimerStart[5] = MPI_Wtime();

                                #pragma omp single
                                {
                                    negchisq.resize(this->dirty(0).shape().nonDegenerate());
                                    negchisq.set(T(0.0));
                                }

                                T* negchisq_pointer = negchisq.getStorage(IsNotCont);
                                for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                                    T* coeff_pointer = coefficients(term1).getStorage(IsNotCont);
                                    T* res_pointer = this->itsResidualBasis(base)(term1).getStorage(IsNotCont);
                                    #pragma omp for schedule(static)
                                    for (int index = 0; index < negchisq.nelements(); index++) {
                                        negchisq_pointer[index] += coeff_pointer[index]*res_pointer[index];
                                    }
                                }

                                if (haveMask) {
                                    #pragma omp single
                                    res = negchisq;

                                    absMaxPosMaskedOMP(maxVal, maxPos, res, maskref);
                                } else {
                                    #pragma omp single
                                    res = negchisq;

                                    absMaxPosOMP(maxVal, maxPos, res);
                                }

                                // Small loop
                                #pragma omp for schedule(static)
                                for (uInt term = 0; term < this->nTerms(); ++term) {
                                            maxValues(term) = coefficients(term)(maxPos);
                                }

                                // End of section 5
                                #pragma omp single
                                { TimerStop[5] = MPI_Wtime(); Times[5] += (TimerStop[5]-TimerStart[5]); }

                            } // End of Maxterm0 or Maxchi solver decision
                        } // End of else decision

                        #pragma omp single
                        {
                            // We use the minVal and maxVal to find the optimum base
                            if (abs(maxVal) > absPeakVal) {
                                    optimumBase = base;
                                    absPeakVal = abs(maxVal);
                                    absPeakPos = maxPos;
                            }
                        }

                    } // End of iteration over number of bases

                    // Now that we know the location of the peak found using one of the
                    // above methods we can look up the values of the residuals. Remember
                    // that we have to decouple the answer

                    // Section 6
                    #pragma omp single
                    TimerStart[6] = MPI_Wtime();

                    #pragma omp single
                    {
                        for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                            peakValues(term1) = 0.0;
                            for (uInt term2 = 0; term2 < this->nTerms(); ++term2) {
                                peakValues(term1) +=
                                    T(this->itsInverseCouplingMatrix(optimumBase)(term1, term2)) *
                                    this->itsResidualBasis(optimumBase)(term2)(absPeakPos);
                            }
                        }

                        // Record location of peak in mask
                        if (this->itsMask.nelements()) this->itsMask(optimumBase)(absPeakPos)=T(1.0);

                        // Take square root to get value comparable to peak residual
                        if (this->itsSolutionType == "MAXCHISQ") {
                            absPeakVal = sqrt(max(T(0.0), absPeakVal));
                        }
                    } // End of omp single section

                    // End of section 6
                    #pragma omp single
                    { TimerStop[6] = MPI_Wtime(); Times[6] += (TimerStop[6]-TimerStart[6]); }

                    if (!this->control()->deepCleanMode() && !decoupled()) {

                        // **** Compute coupled residual ****

                        // Section 7
                        #pragma omp single
                        TimerStart[7] = MPI_Wtime();

                        for (uInt term = 0; term < this->nTerms(); term++) {
                            for (uInt base = 0; base < nBases; base++) {

                                maxPos(0) = 0; maxPos(1) = 0;
                                maxVal = 0.0;
                                #pragma omp single
                                res.reference(this->itsResidualBasis(base)(term));

                                if (isWeighted) {
                                    #pragma omp single
                                    wt.reference(this->itsWeight(0).nonDegenerate());

                                    absMaxPosMaskedOMP(maxVal, maxPos, res, wt);
                                } else {

                                    absMaxPosOMP(maxVal, maxPos, res);
                                }
                                // TODO: Do I need this barrier? No - absmax already has one
                                #pragma omp barrier

                                #pragma omp single
                                {
                                    maxBaseVals(base) = abs(this->itsResidualBasis(base)(term)(maxPos));
                                }
                            } // End of loop over bases

                            #pragma omp single
                            {
                                casa::IPosition minPos(1, 0);
                                casa::IPosition maxPos(1, 0);
                                T minVal(0.0), maxVal(0.0);
                                casa::minMax(minVal, maxVal, minPos, maxPos, maxBaseVals);
                                maxTermVals(term) = maxVal;
                            }
                        } // End of loop over terms

                        // End of Section 7
                        #pragma omp single
                        { TimerStop[7] = MPI_Wtime(); Times[7] += (TimerStop[7]-TimerStart[7]); }

                        #pragma omp single
                        {
                            casa::IPosition minPos(1, 0);
                            casa::IPosition maxPos(1, 0);
                            T minVal(0.0), maxVal(0.0);
                            casa::minMax(minVal, maxVal, minPos, maxPos, maxTermVals);
                            absPeakVal = maxVal;
                        }

                    } // End of coupled residual section

                    // TODO: Check this barrier?
                    #pragma omp barrier

                    // Section 8
                    #pragma omp single
                    TimerStart[8] = MPI_Wtime();

                    #pragma omp sections
                    {

                        #pragma omp section
                        {
                            if (this->state()->initialObjectiveFunction() == 0.0) {
                                this->state()->setInitialObjectiveFunction(abs(absPeakVal));
                            }
                        }

                        #pragma omp section
                        this->state()->setPeakResidual(abs(absPeakVal));

                        #pragma omp section
                        this->state()->setObjectiveFunction(abs(absPeakVal));

                        // Prepare the sum for the totalflux
                        #pragma omp section
                        sumFlux = 0.0;

                    } // End of single

                    float localsum = 0.0;
                    float* model_pointer = this->model(0).getStorage(IsNotCont);
                    #pragma omp for schedule(static) nowait
                    for (int index = 0; index < this->model(0).nelements(); index++) {
                        localsum += model_pointer[index];
                    }

                    #pragma omp critical
                    sumFlux += localsum;

                    // This barrier is required - no implicit barrier following criticals
                    #pragma omp barrier

                    // without OpenMP, this may be faster
                    //sumFlux = sum(this->model(0));

                    #pragma omp single
                    this->state()->setTotalFlux(sumFlux);

                    #pragma omp sections
                    {
                        // Now we adjust model and residual for this component
                        #pragma omp section
                        residualShape = this->dirty(0).shape().nonDegenerate();

                        #pragma omp section
                        {
                            psfShape(0) = this->itsBasisFunction->shape()(0),
                            psfShape(1) = this->itsBasisFunction->shape()(1);
                        }
                    }

                    casa::IPosition residualStart(2, 0), residualEnd(2, 0), residualStride(2, 1);
                    casa::IPosition psfStart(2, 0), psfEnd(2, 0), psfStride(2, 1);
                    const casa::IPosition modelShape(this->model(0).shape().nonDegenerate());
                    casa::IPosition modelStart(2, 0), modelEnd(2, 0), modelStride(2, 1);

                    // End of section 8
                    #pragma omp single
                    { TimerStop[8] = MPI_Wtime(); Times[8] += (TimerStop[8]-TimerStart[8]); }

                    // Section 9
                    #pragma omp single
                    TimerStart[9] = MPI_Wtime();

                    const casacore::IPosition peakPSFPos = this->getPeakPSFPosition();
                    ASKAPDEBUGASSERT(peakPSFPos.nelements() >= 2);
                    // that there are some edge cases for which it fails.
                    for (uInt dim = 0; dim < 2; dim++) {
                        residualStart(dim) = max(0, Int(absPeakPos(dim) - psfShape(dim) / 2));
                        residualEnd(dim) = min(Int(absPeakPos(dim) + psfShape(dim) / 2 - 1), Int(residualShape(dim) - 1));
                        // Now we have to deal with the PSF. Here we want to use enough of the
                        // PSF to clean the residual image.
                        psfStart(dim) = max(0, Int(peakPSFPos(dim) - (absPeakPos(dim) - residualStart(dim))));
                        psfEnd(dim) = min(Int(peakPSFPos(dim) - (absPeakPos(dim) - residualEnd(dim))),
                                        Int(psfShape(dim) - 1));
                        modelStart(dim) = residualStart(dim);
                        modelEnd(dim) = residualEnd(dim);
                    }

                    casa::Slicer psfSlicer(psfStart, psfEnd, psfStride, Slicer::endIsLast);
                    casa::Slicer residualSlicer(residualStart, residualEnd, residualStride, Slicer::endIsLast);
                    casa::Slicer modelSlicer(modelStart, modelEnd, modelStride, Slicer::endIsLast);

                    // Add to model
                    // We loop over all terms for the optimum base and ignore those terms with no flux
                    #pragma omp single
                    {
                        for (uInt term = 0; term < this->nTerms(); ++term) {
                            if (abs(peakValues(term)) > 0.0) {
                                casa::Array<float> slice = this->model(term).nonDegenerate()(modelSlicer);
                                slice += this->control()->gain() * peakValues(term) *
                                        this->itsBasisFunction->basisFunction(optimumBase).nonDegenerate()(psfSlicer);
                                this->itsTermBaseFlux(optimumBase)(term) += this->control()->gain() * peakValues(term);
                            }
                        }
                    }

/*
                    const uInt ni = residualEnd(0) - residualStart(0);
                    const uInt nj = residualEnd(1) - residualStart(1);
                    const uInt ri0 = residualStart(0);
                    const uInt rj0 = residualStart(1);
                    const uInt pi0 = psfStart(0);
                    const uInt pj0 = psfStart(1);

                    // Add to model
                    // We loop over all terms for the optimum base and ignore those terms with no flux
                    for (uInt term = 0; term < this->nTerms(); ++term) {
                        if (abs(peakValues(term)) > 0.0) {
                            const T amp = this->control()->gain() * peakValues(term);
                            casa::Matrix<T> mMdl, mBfn;
                            mMdl.reference(this->model(term).nonDegenerate()(modelSlicer));
                            mBfn.reference(Cube<T>(this->itsBasisFunction->basisFunction()).xyPlane(optimumBase).nonDegenerate()(psfSlicer));
                            #pragma omp for schedule(static)
                            for (uInt j = 0; j < nj; j++ ) {
                                Vector<T> mdlcol = mMdl.column(j);
                                T* pMdl = mdlcol.getStorage(IsNotCont);
                                Vector<T> bfncol = mBfn.column(j);
                                T* pBfn = bfncol.getStorage(IsNotCont);
                                for (uInt i = 0; i < ni; i++ ) {
                                    pMdl[i] += amp * (*(pBfn+i));
                                }
                            }

                            this->itsTermBaseFlux(optimumBase)(term) += this->control()->gain() * peakValues(term);

                        }
                    }
*/

                    #pragma omp single
                    {
                        // Subtract PSFs, including base-base crossterms
                        for (uInt term1 = 0; term1 < this->nTerms(); term1++) {
                            for (uInt term2 = 0; term2 < this->nTerms(); term2++) {
                                if (abs(peakValues(term2)) > 0.0) {
                                    for (uInt base = 0; base < nBases; base++) {
                                        // This can be done in parallel, but isnt worth it.
                                        this->itsResidualBasis(base)(term1)(residualSlicer) =
                                            this->itsResidualBasis(base)(term1)(residualSlicer)
                                            - this->control()->gain() * peakValues(term2) *
                                            this->itsPSFCrossTerms(base, optimumBase)(term1, term2)(psfSlicer);
                                    }
                                }
                            }
                        }
                    }

/*
                    // Subtract PSFs, including base-base crossterms
                    for (uInt term1 = 0; term1 < this->nTerms(); term1++) {
                        for (uInt term2 = 0; term2 < this->nTerms(); term2++) {
                            if (abs(peakValues(term2)) > 0.0) {
                                for (uInt base = 0; base < nBases; base++) {


                                    // This can be done in parallel, but isnt worth it.
                                    //this->itsResidualBasis(base)(term1)(residualSlicer) =
                                    //    this->itsResidualBasis(base)(term1)(residualSlicer)
                                    //    - this->control()->gain() * peakValues(term2) *
                                    //    this->itsPSFCrossTerms(base, optimumBase)(term1, term2)(psfSlicer);

                                    const T amp = this->control()->gain() * peakValues(term2);
                                    casa::Matrix<T> mRes, mPSF;
                                    mRes.reference(this->itsResidualBasis(base)(term1));
                                    mPSF.reference(this->itsPSFCrossTerms(base,optimumBase)(term1,term2));
                                    #pragma omp for schedule(static)
                                    for (uInt j = 0; j < nj; j++ ) {
                                        Vector<T> rescol = mRes.column(rj0+j);
                                        T* pRes = rescol.getStorage(IsNotCont);
                                        Vector<T> psfcol = mPSF.column(pj0+j);
                                        T* pPSF = psfcol.getStorage(IsNotCont);
                                        for (uInt i = 0; i < ni; i++ ) {
                                            pRes[ri0+i] -= amp * (*(pPSF+pi0+i));
                                        }
                                    }

                                }
                            }
                        }
                    }
*/

                    #pragma omp single
                    {
						this->monitor()->monitor(*(this->state()));
						this->state()->incIter();
                    }

                    // End of section 9
                    #pragma omp single
                    { TimerStop[9] = MPI_Wtime(); Times[9] += (TimerStop[9]-TimerStart[9]); }

                    //End of all iterations
                    #pragma omp barrier

					#pragma omp single
					{
						converged = this->control()->terminate(*(this->state()));
					}

                } while (!converged);

            } // End of parallel section

            //const double end_time = MPI_Wtime();
            const time_t end_time = time(0);
            ASKAPLOG_INFO_STR(decmtbflogger,
                              "Time for minor cycles: "<<end_time-start_time<<" sec");

            // Report Times
            double sum_time = 0.0;
            for (int i = 0; i < no_timers; i++) {
                ASKAPLOG_INFO_STR(decmtbflogger, "Section "<<i<<" Time: "<<Times[i]);
                sum_time += Times[i];
            }

            ASKAPLOG_INFO_STR(decmtbflogger, "Performed Multi-Term BasisFunction CLEAN for "
                                  << this->state()->currentIter() << " iterations");
            ASKAPLOG_INFO_STR(decmtbflogger, this->control()->terminationString());

            } else {
              ASKAPLOG_INFO_STR(decmtbflogger,
                  "Bypassed Multi-Term BasisFunction CLEAN due to 0 iterations in the setup");
            }

        } // End of many iterations function



        template<class T, class FT>
        bool DeconvolverMultiTermBasisFunction<T, FT>::deconvolve()
        {
            // This is the parallel version of deconvolve using ManyIterations()
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::deconvolve");
            this->initialise();
            //const double start_time = MPI_Wtime();
            const time_t start_time = time(0);
            this->ManyIterations();
            //const double end_time = MPI_Wtime();
            const time_t end_time = time(0);
            this->finalise();
            ASKAPLOG_INFO_STR(decmtbflogger, "Time Required: "<<end_time - start_time);
            // signal failure and finish the major cycles if we started to diverge
            return (this->control()->terminationCause() != DeconvolverControl<T>::DIVERGED);
        }

        // Helper function to replace minMaxMasked calls when we only need the abs maximum
        template<class T>
        void absMaxPosMasked(T& maxVal, IPosition&  maxPos,  const Matrix<T>& im, const Matrix<T>& mask)
        {
            maxVal = T(0);
            const uInt ncol = mask.ncolumn();
            const uInt nrow = mask.nrow();
            for (uInt j = 0; j < ncol; j++ ) {
                const T* pIm = &im(0,j);
                const T* pMask = &mask(0,j);
                for (uInt i = 0; i < nrow; i++ ) {
                    //T val = abs(mask(i,j) * im(i,j));
                    T val = abs(*pIm++ * *pMask++);
                    if (val > maxVal) {
                        maxVal = val;
                        maxPos(0) = i;
                        maxPos(1) = j;
                    }
                }
            }
        }

        template<class T>
        void absMaxPosMaskedACC(T& maxVal, int&  maxPos,  const T* im, const T* mask, uInt nele)
        {
            float maxValf = 0;
            int maxPosI = 0;
            //const uInt ncol = mask.ncolumn();
            //const uInt nrow = mask.nrow();
            // Parallel reduction with openacc
            #pragma acc parallel loop reduction(max:maxValf) present(mask, im)
            for (int i = 0; i < nele; i++ ) {
                T test = abs(im[i] * mask[i]);
                if (test > maxValf) {
                   maxValf = test;

                }
            }


            #pragma acc parallel loop present(mask,im)
            for (int i = 0; i < nele; i++ ) {
                if (abs(im[i] * mask[i]) == maxValf) {
                    maxPosI = i;
                }
            }

	        printf("MaxPosI = %d\n", maxPosI);

            maxVal = maxValf;
            maxPos = maxPosI;

           // Now we can change the shape to return (i,j) for the max value.
           // Do this later
           printf("DEBUG\tMS SUT Max value = %g, Location = %d\n", maxVal, maxPos);

        }


        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::getCoupledResidual(T& absPeakRes) {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction:::getCoupledResidual");
            const uInt nBases(this->itsResidualBasis.nelements());
            bool isWeighted((this->itsWeight.nelements() > 0) &&
                (this->itsWeight(0).shape().nonDegenerate().conform(this->itsResidualBasis(0)(0).shape())));

            Vector<T> maxTermVals(this->nTerms());
            Vector<T> maxBaseVals(nBases);

            for (uInt term = 0; term < this->nTerms(); term++) {
                for (uInt base = 0; base < nBases; base++) {
                    casacore::IPosition minPos(2, 0);
                    casacore::IPosition maxPos(2, 0);
                    T minVal(0.0), maxVal(0.0);
                    if (isWeighted) {
                        const casacore::Matrix<T> res = this->itsResidualBasis(base)(term);
                        const casacore::Matrix<T> wt = this->itsWeight(0).nonDegenerate();
                        absMaxPosMasked(maxVal, maxPos, res, wt);
                        //casacore::minMaxMasked(minVal, maxVal, minPos, maxPos, this->itsResidualBasis(base)(term),
                        //                   this->itsWeight(0).nonDegenerate());
                    } else {
                        casacore::minMax(minVal, maxVal, minPos, maxPos, this->itsResidualBasis(base)(term));
                    }
                    if (abs(minVal) > abs(maxVal)) {
                        maxBaseVals(base) = abs(this->itsResidualBasis(base)(term)(minPos));
                    }
                    else {
                        maxBaseVals(base) = abs(this->itsResidualBasis(base)(term)(maxPos));
                    }

                }
                casacore::IPosition minPos(1, 0);
                casacore::IPosition maxPos(1, 0);
                T minVal(0.0), maxVal(0.0);
                casacore::minMax(minVal, maxVal, minPos, maxPos,maxBaseVals);
                maxTermVals(term) = maxVal;
            }
            casacore::IPosition minPos(1, 0);
            casacore::IPosition maxPos(1, 0);
            T minVal(0.0), maxVal(0.0);
            casacore::minMax(minVal, maxVal, minPos, maxPos,maxTermVals);
            absPeakRes = maxVal;
        }
        // This contains the heart of the Multi-Term BasisFunction Clean algorithm
        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::chooseComponent(uInt& optimumBase,
                casacore::IPosition& absPeakPos, T& absPeakVal, Vector<T>& peakValues)
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction:::chooseComponent");

            const uInt nBases(this->itsResidualBasis.nelements());

	    // MS SUT Set up some manual timers
            const int NO_TIMERS = 12;
            static double MS_SUT_debug_timer[NO_TIMERS];
            static int MS_SUT_count = 0;  			 // Lazy iteration count
            double startTime[NO_TIMERS], endTime[NO_TIMERS]; 	// Wall times

            if (MS_SUT_count == 0){
                for (int timers = 0; timers < NO_TIMERS; timers++) {

                        MS_SUT_debug_timer[timers] = 0.0;
                }
            }
#ifdef USE_OPENACC
            int optimumIdx;
#endif
	    // MS SUT Set up some checksum variables
	        float testsum[] = {0.0, 0.0};


            absPeakVal = 0.0;

            ASKAPDEBUGASSERT(peakValues.nelements() <= this->nTerms());

            // Find the base having the peak value in term=0
            // Here the weights image is used as a weight in the determination
            // of the maximum i.e. it finds the max in weight . residual. The values
            // returned are without the weight
            bool isWeighted((this->itsWeight.nelements() > 0) &&
                (this->itsWeight(0).shape().nonDegenerate().conform(this->itsResidualBasis(0)(0).shape())));

            Vector<T> minValues(this->nTerms());
            Vector<T> maxValues(this->nTerms());

            // Set the mask - we need it for weighted search and deep clean
            Matrix<T> mask;
            if (isWeighted) {
                mask = this->itsWeight(0).nonDegenerate();
                if  (this->itsSolutionType == "MAXCHISQ") {
                    // square weights for MAXCHISQ
                    mask*=mask;
                }
            }

            for (uInt base = 0; base < nBases; base++) {

                // Find peak in residual image cube
                casacore::IPosition minPos(2, 0);
                casacore::IPosition maxPos(2, 0);
                int Idx;
                T minVal(0.0), maxVal(0.0);

                if (this->control()->deepCleanMode()) {
                    if (isWeighted) {
                        // recompute mask*weight for each new base
                        if (base>0) {
                            mask = this->itsWeight(0).nonDegenerate();
                            if  (this->itsSolutionType == "MAXCHISQ") {
                                // square weights for MAXCHISQ
                                mask*=mask;
                            }
                        }
#ifdef USE_OPENACC
                        itsACCManager.InitMask(base);
                        itsACCManager.UpdateMask(base);
#else
                        mask*=this->itsMask(base);
#endif

                    } else {
                        mask=this->itsMask(base);
                    }
                }

                else {
#ifdef USE_OPENACC
                itsACCManager.InitMask(base);
#endif
                }

                Bool haveMask=mask.nelements()>0;

                // We implement various approaches to finding the peak. The first is the cheapest
                // and evidently the best (according to Urvashi).

                // Look for the maximum in term=0 for this base
                if (this->itsSolutionType == "MAXBASE") {
                    if (haveMask) {
                        const casacore::Matrix<T> res = this->itsResidualBasis(base)(0);
#ifdef USE_OPENACC
                        bool deleteIm,deleteMa;
                        int nelements = res.nelements();
                        const T * im = (T *) itsACCManager.residuals[base*itsACCManager.nTerms];
                        const T * ma = (T *) itsACCManager.maskToUse;
                        printf("Check Array Locations im:%p ma:%p\n",im,ma);
                        absMaxPosMaskedACC(maxVal,Idx,im,ma,nelements);
                        const int y = Idx / itsACCManager.nrows;
                        const int x = Idx % itsACCManager.ncols;
                        maxPos(0) = x;
                        maxPos(1) = y;

                        printf("Check Max Locations (OpenACC): val=%f at: %d, %d, %d\n", maxVal, maxPos(0), maxPos(1), Idx);

                        //absMaxPosMasked(maxVal, maxPos, res, mask);
                        //printf("Check Max Locations (Serial): %d, %d\n", maxPos(0), maxPos(1));
#else
                        absMaxPosMasked(maxVal, maxPos, res, mask);
			           // printf("Check Max Locations (Serial): %d, %d\n", maxPos(0), maxPos(1));
//                      casacore::minMaxMasked(minVal, maxVal, minPos, maxPos, this->itsResidualBasis(base)(0),mask)

#endif
                    } else {
                        casacore::minMax(minVal, maxVal, minPos, maxPos, this->itsResidualBasis(base)(0));
                    }
                    for (uInt term = 0; term < this->nTerms(); ++term) {
                        minValues(term) = this->itsResidualBasis(base)(term)(minPos);
                        maxValues(term) = this->itsResidualBasis(base)(term)(maxPos);
                    }
                    // In performing the search for the peak across bases, we want to take into account
                    // the SNR so we normalise out the coupling matrix for term=0 to term=0.
                    const T norm(1 / sqrt(this->itsCouplingMatrix(base)(0, 0)));
                    maxVal *= norm;
                    minVal *= norm;
                } else {
                    // All these algorithms need the decoupled terms

                    // Decouple all terms using inverse coupling matrix
                    Vector<Array<T> > coefficients(this->nTerms());
                    for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                        coefficients(term1).resize(this->dirty(0).shape().nonDegenerate());
                        coefficients(term1).set(T(0.0));
                        for (uInt term2 = 0; term2 < this->nTerms(); ++term2) {
                            coefficients(term1) = coefficients(term1) +
                                                  T(this->itsInverseCouplingMatrix(base)(term1, term2)) *
                                                  this->itsResidualBasis(base)(term2);
                        }
                    }

                    if (this->itsSolutionType == "MAXTERM0") {
                        if (haveMask) {
                            casacore::minMaxMasked(minVal, maxVal, minPos, maxPos, coefficients(0),
                                               mask);
                        } else {
                            casacore::minMax(minVal, maxVal, minPos, maxPos, coefficients(0));
                        }
                        for (uInt term = 0; term < this->nTerms(); ++term) {
                            minValues(term) = coefficients(term)(minPos);
                            maxValues(term) = coefficients(term)(maxPos);
                        }
                    } else {
                        // MAXCHISQ
                        // Now form the criterion image and then search for the peak.
                        Array<T> negchisq(this->dirty(0).shape().nonDegenerate());
                        negchisq.set(T(0.0));
                        for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                            negchisq = negchisq + coefficients(term1) * this->itsResidualBasis(base)(term1);
                        }
                        // Need to take the square root to ensure that the SNR weighting is correct
                        //            ASKAPCHECK(min(negchisq)>0.0, "Negchisq has negative values");
                        //            negchisq=sqrt(negchisq);
                        //            SynthesisParamsHelper::saveAsCasaImage("negchisq.img",negchisq);
                        //            SynthesisParamsHelper::saveAsCasaImage("coefficients0.img",coefficients(0));
                        //            SynthesisParamsHelper::saveAsCasaImage("coefficients1.img",coefficients(1));
                        //            ASKAPTHROW(AskapError, "Written debug images");
                        // Remember that the weights must be squared.
                        if (haveMask) {
                            casacore::minMaxMasked(minVal, maxVal, minPos, maxPos, negchisq,
                                               mask);
                        } else {
                            casacore::minMax(minVal, maxVal, minPos, maxPos, negchisq);
                        }
                        for (uInt term = 0; term < this->nTerms(); ++term) {
                            minValues(term) = coefficients(term)(minPos);
                            maxValues(term) = coefficients(term)(maxPos);
                        }
                    }
                }

                // We use the minVal and maxVal to find the optimum base
                if (abs(minVal) > absPeakVal) {
                    optimumBase = base;
                    absPeakVal = abs(minVal);
                    absPeakPos = minPos;
#ifdef USE_OPENACC
                    optimumIdx = Idx;
#endif
                }
                if (abs(maxVal) > absPeakVal) {
                    optimumBase = base;
                    absPeakVal = abs(maxVal);
                    absPeakPos = maxPos;
#ifdef USE_OPENACC
                    optimumIdx = Idx;
#endif
                }
            }

            // Now that we know the location of the peak found using one of the
            // above methods we can look up the values of the residuals. Remember
            // that we have to decouple the answer
            for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                peakValues(term1) = 0.0;
                for (uInt term2 = 0; term2 < this->nTerms(); ++term2) {
                    peakValues(term1) +=
                        T(this->itsInverseCouplingMatrix(optimumBase)(term1, term2)) *
                        this->itsResidualBasis(optimumBase)(term2)(absPeakPos);
                }
            }

            // Record location of peak in mask
            if (this->itsMask.nelements()) this->itsMask(optimumBase)(absPeakPos)=T(1.0);
#ifdef USE_OPENACC
            T * tomove = (T *) itsACCManager.masks[optimumBase];
//          printf("device ptr %p offset: %d\n",tomove,optimumIdx);
            #pragma acc update device(tomove[optimumIdx:1])
#endif
            // Take square root to get value comparable to peak residual
            if (this->itsSolutionType == "MAXCHISQ") {
                absPeakVal = sqrt(max(T(0.0), absPeakVal));
            }

            // Not sure I agree with the this I think the absPeakVal should
            // be the absolute value of the peak residual
            // For deep cleaning we want to restrict the abspeakval to the mask
            // so we just use the value determined above
            if (!this->control()->deepCleanMode() && !decoupled()) getCoupledResidual(absPeakVal);


	    // Update collective time over all iterations
	    MS_SUT_debug_timer[0] += (endTime[0]-startTime[0]);
            MS_SUT_debug_timer[1] += (endTime[1]-startTime[1]);

	    // Report timings
	    // printf("DEBUG\t MS SUT Timings: Mask*Mask OpenACC: %g, Serial: %g\n", MS_SUT_debug_timer[0], MS_SUT_debug_timer[1]);

        }

        template<class T, class FT>
        bool DeconvolverMultiTermBasisFunction<T, FT>::oneIteration()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::oneIteration");

            // For the psf convolutions, we only need a small part of the
            // basis functions so we recalculate for that size
            IPosition subPsfShape(this->findSubPsfShape());

            const uInt nBases(this->itsResidualBasis.nelements());

            casacore::IPosition absPeakPos(2, 0);
            T absPeakVal(0.0);
            uInt optimumBase(0);
            Vector<T> peakValues(this->nTerms());
            chooseComponent(optimumBase, absPeakPos, absPeakVal, peakValues);

            // Report on progress
            // We want the worst case residual
            // T absPeakRes = max(abs(peakValues));



            //      ASKAPLOG_INFO_STR(decmtbflogger, "All terms: absolute max = " << absPeakRes << " at " << absPeakPos);
            //      ASKAPLOG_INFO_STR(decmtbflogger, "Optimum base = " << optimumBase);

            if (this->state()->initialObjectiveFunction() == 0.0) {
                this->state()->setInitialObjectiveFunction(abs(absPeakVal));
            }
            this->state()->setPeakResidual(abs(absPeakVal));
            this->state()->setObjectiveFunction(abs(absPeakVal));
            this->state()->setTotalFlux(sum(this->model(0)));

            //  Check if we should enter deep cleaning mode
            if (abs(absPeakVal) < this->control()->targetObjectiveFunction() &&
                this->control()->targetObjectiveFunction2()>0 &&
                abs(absPeakVal) > this->control()->targetObjectiveFunction2()) {
                if (!this->control()->deepCleanMode()) {
                    ASKAPLOG_INFO_STR(decmtbflogger, "Starting deep cleaning phase");
                }
                //setDeepCleanMode(True);
                this->control()->setDeepCleanMode();
            }

            // Now we adjust model and residual for this component
            const casacore::IPosition residualShape(this->dirty(0).shape().nonDegenerate());
            //IPosition subPsfStart(2, nx / 2 - subPsfShape(0) / 2, ny / 2 - subPsfShape(1) / 2);
            //IPosition subPsfEnd(2, nx / 2 + subPsfShape(0) / 2 - 1, ny / 2 + subPsfShape(1) / 2 - 1);
            //IPosition subPsfStride(2, 1, 1);

            //Slicer subPsfSlicer(subPsfStart, subPsfEnd, subPsfStride, Slicer::endIsLast);
            const casacore::IPosition psfShape = this->itsBasisFunction->shape().getFirst(2);

            casacore::IPosition residualStart(2, 0), residualEnd(2, 0), residualStride(2, 1);
            casacore::IPosition psfStart(2, 0), psfEnd(2, 0), psfStride(2, 1);

            const casacore::IPosition modelShape(this->model(0).shape().nonDegenerate());
            casacore::IPosition modelStart(2, 0), modelEnd(2, 0), modelStride(2, 1);

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

                modelStart(dim) = residualStart(dim);
                modelEnd(dim) = residualEnd(dim);
            }

            casacore::Slicer psfSlicer(psfStart, psfEnd, psfStride, Slicer::endIsLast);
            casacore::Slicer residualSlicer(residualStart, residualEnd, residualStride, Slicer::endIsLast);
            casacore::Slicer modelSlicer(modelStart, modelEnd, modelStride, Slicer::endIsLast);

            // Add to model
            // We loop over all terms for the optimum base and ignore those terms with no flux
            for (uInt term = 0; term < this->nTerms(); ++term) {
                if (abs(peakValues(term)) > 0.0) {
                    casacore::Array<float> slice = this->model(term).nonDegenerate()(modelSlicer);
                    slice += this->control()->gain() * peakValues(term) *
                             this->itsBasisFunction->basisFunction(optimumBase).nonDegenerate()(psfSlicer);
                    this->itsTermBaseFlux(optimumBase)(term) += this->control()->gain() * peakValues(term);
                }
            }

            // Subtract PSFs, including base-base crossterms
            for (uInt term1 = 0; term1 < this->nTerms(); ++term1) {
                for (uInt term2 = 0; term2 < this->nTerms(); ++term2) {
                    if (abs(peakValues(term2)) > 0.0) {
                        for (uInt base = 0; base < nBases; base++) {
                            this->itsResidualBasis(base)(term1)(residualSlicer) =
                                this->itsResidualBasis(base)(term1)(residualSlicer)
                                - this->control()->gain() * peakValues(term2) *
                                this->itsPSFCrossTerms(base, optimumBase)(term1, term2)(psfSlicer);

                        }
                    }
                }
            }

            return True;
        }

    }
}
// namespace synthesis
// namespace askap
