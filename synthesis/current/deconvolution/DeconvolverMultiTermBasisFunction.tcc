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
#include <askap/AskapLogging.h>
#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <scimath/Mathematics/MatrixMathLA.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <profile/AskapProfiler.h>
ASKAP_LOGGER(decmtbflogger, ".deconvolution.multitermbasisfunction");

#include <deconvolution/DeconvolverMultiTermBasisFunction.h>
#include <deconvolution/MultiScaleBasisFunction.h>

namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver based on the BasisFunction Clean
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a residual image, psf optionally using a weights image.
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverMultiTermBasisFunction<Double, DComplex>
        /// @ingroup Deconvolver

        template<class T, class FT>
        DeconvolverMultiTermBasisFunction<T, FT>::DeconvolverMultiTermBasisFunction(Vector<Array<T> >& dirty,
                Vector<Array<T> >& psf,
                Vector<Array<T> >& psfLong)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf), itsDirtyChanged(True), itsBasisFunctionChanged(True),
                itsSolutionType("MAXCHISQ")
        {
            ASKAPLOG_DEBUG_STR(decmtbflogger, "There are " << this->itsNumberTerms << " terms to be solved");

            ASKAPCHECK(psfLong.nelements() == (2*this->itsNumberTerms - 1), "Long PSF vector has incorrect length " << psfLong.nelements());
            this->itsPsfLongVec.resize(2*this->itsNumberTerms - 1);

            for (uInt term = 0; term < (2*this->itsNumberTerms - 1); term++) {
                ASKAPCHECK(psfLong(term).nonDegenerate().shape().nelements() == 2, "PSF(" << term << ") has too many dimensions " << psfLong(term).shape());
                this->itsPsfLongVec(term) = psfLong(term).nonDegenerate();
            }

        };

        template<class T, class FT>
        DeconvolverMultiTermBasisFunction<T, FT>::DeconvolverMultiTermBasisFunction(Array<T>& dirty,
                Array<T>& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf), itsDirtyChanged(True), itsBasisFunctionChanged(True),
                itsSolutionType("MAXCHISQ")
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
        void DeconvolverMultiTermBasisFunction<T, FT>::updateDirty(Array<T>& dirty, casa::uInt term)
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

            for (uInt base = 0; base < itsTermBaseFlux.nelements(); base++) {
                for (uInt term = 0; term < itsTermBaseFlux(base).nelements(); term++) {
                    ASKAPLOG_DEBUG_STR(decmtbflogger, "   Term(" << term << "), Base(" << base
                                           << "): Flux = " << itsTermBaseFlux(base)(term));
                }
            }

        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialiseForBasisFunction(bool force)
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialiseForBasisFunction");
            if (!force && !this->itsBasisFunctionChanged) return;

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

            // Force change in basis function
            initialiseForBasisFunction(true);

            this->state()->resetInitialObjectiveFunction();
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialiseResidual()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialiseResidual");

            if (!this->itsDirtyChanged) return;

            // Initialise the basis function for residual calculations.
            this->itsBasisFunction->initialise(this->dirty(0).shape());

            ASKAPCHECK(this->itsBasisFunction, "Basis function not initialised");

            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of basis functions "
                                   << this->itsBasisFunction->basisFunction().shape());

            uInt nBases(this->itsBasisFunction->numberBases());

            itsResidualBasis.resize(nBases);
            for (uInt base = 0; base < nBases; base++) {
                itsResidualBasis(base).resize(this->itsNumberTerms);
            }

            // Calculate residuals convolved with bases [nx,ny][nterms][nbases]
            // Calculate transform of PSF(0)
            Array<FT> xfrZero(this->psf(0).shape().nonDegenerate());
            xfrZero.set(FT(0.0));
            casa::setReal(xfrZero, this->psf(0).nonDegenerate());
            scimath::fft2d(xfrZero, true);
            T normPSF;
            normPSF = casa::sum(casa::real(xfrZero * conj(xfrZero))) / xfrZero.nelements();
            ASKAPLOG_DEBUG_STR(decmtbflogger, "PSF effective volume = " << normPSF);
            xfrZero = xfrZero / FT(normPSF);

            ASKAPLOG_DEBUG_STR(decmtbflogger,
                               "Calculating convolutions of residual images with basis functions");
            for (uInt base = 0; base < nBases; base++) {
                // Calculate transform of residual images [nx,ny,nterms]
                for (uInt term = 0; term < this->itsNumberTerms; term++) {

                    // Calculate transform of residual image
                    Array<FT> residualFFT(this->dirty(term).shape().nonDegenerate());
                    residualFFT.set(FT(0.0));
                    casa::setReal(residualFFT, this->dirty(term).nonDegenerate());
                    scimath::fft2d(residualFFT, true);

                    // Calculate transform of basis function [nx,ny,nbases]
                    Matrix<FT> basisFunctionFFT(this->dirty(term).shape().nonDegenerate());
                    basisFunctionFFT.set(FT(0.0));
                    casa::setReal(basisFunctionFFT, Cube<T>(this->itsBasisFunction->basisFunction()).xyPlane(base));
                    scimath::fft2d(basisFunctionFFT, true);

                    // Calculate product and transform back
                    Array<FT> work(this->dirty(term).nonDegenerate().shape());
                    ASKAPASSERT(basisFunctionFFT.shape().conform(residualFFT.shape()));
                    work = conj(basisFunctionFFT) * residualFFT * conj(xfrZero);
                    scimath::fft2d(work, false);

                    // basis function * psf
                    ASKAPLOG_DEBUG_STR(decmtbflogger, "Basis(" << base
                                           << ")*PSF(0)*Residual(" << term << "): max = " << max(real(work))
                                           << " min = " << min(real(work)));

                    this->itsResidualBasis(base)(term) = real(work);
                }
            }
        }

        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::initialisePSF()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::initialisePSF");

            if (!this->itsBasisFunctionChanged) return;

            ASKAPCHECK(this->itsBasisFunction, "Basis function not initialised");

            ASKAPLOG_DEBUG_STR(decmtbflogger,
                               "Updating Multi-Term Basis Function deconvolver for change in basis function");
            IPosition subPsfShape(this->findSubPsfShape());

            Array<FT> work(subPsfShape);

            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of basis functions "
                                   << this->itsBasisFunction->basisFunction().shape());

            IPosition stackShape(this->itsBasisFunction->basisFunction().shape());

            uInt nBases(this->itsBasisFunction->numberBases());

            // Now transform the basis functions. These may be a different size from
            // those in initialiseResidual so we don't keep either
            Cube<FT> basisFunctionFFT(this->itsBasisFunction->basisFunction().shape());
            basisFunctionFFT.set(FT(0.0));
            casa::setReal(basisFunctionFFT, this->itsBasisFunction->basisFunction());
            scimath::fft2d(basisFunctionFFT, true);

            itsTermBaseFlux.resize(nBases);
            for (uInt base = 0; base < nBases; base++) {
                itsTermBaseFlux(base).resize(this->itsNumberTerms);
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
            ASKAPCHECK(subPsfSlicer.length() == subPsfShape, "Slicer selected length of " << subPsfSlicer.length() <<
                       " is different from requested shape " << subPsfShape);

            casa::IPosition minPos;
            casa::IPosition maxPos;
            T minVal, maxVal;
            casa::minMax(minVal, maxVal, minPos, maxPos, this->psf(0).nonDegenerate()(subPsfSlicer));
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Maximum of PSF(0) = " << maxVal << " at " << maxPos);
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Minimum of PSF(0) = " << minVal << " at " << minPos);
            this->itsPeakPSFVal = maxVal;
            this->itsPeakPSFPos(0) = maxPos(0);
            this->itsPeakPSFPos(1) = maxPos(1);

            IPosition subPsfPeak(2, this->itsPeakPSFPos(0), this->itsPeakPSFPos(1));
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Peak of PSF subsection at  " << subPsfPeak);
            ASKAPLOG_DEBUG_STR(decmtbflogger, "Shape of PSF subsection is " << subPsfShape);

            // Calculate XFR for the subsection only. We need all PSF's up to
            // 2*nTerms-1
            ASKAPCHECK(this->itsPsfLongVec.nelements() == (2*this->itsNumberTerms - 1), "PSF long vector has wrong length " << this->itsPsfLongVec.nelements());

            // Calculate all the transfer functions
            Vector<Array<FT> > subXFRVec(2*this->itsNumberTerms - 1);
            for (uInt term1 = 0; term1 < (2*this->itsNumberTerms - 1); term1++) {
                subXFRVec(term1).resize(subPsfShape);
                subXFRVec(term1).set(0.0);
                casa::setReal(subXFRVec(term1), this->itsPsfLongVec(term1).nonDegenerate()(subPsfSlicer));
                scimath::fft2d(subXFRVec(term1), true);
            }
            // Calculate residuals convolved with bases [nx,ny][nterms][nbases]
            // Calculate transform of PSF(0)
            T normPSF;
            normPSF = casa::sum(casa::real(subXFRVec(0) * conj(subXFRVec(0)))) / subXFRVec(0).nelements();
            ASKAPLOG_DEBUG_STR(decmtbflogger, "PSF effective volume = " << normPSF);

            itsPSFCrossTerms.resize(nBases, nBases);
            for (uInt base = 0; base < nBases; base++) {
                for (uInt base1 = 0; base1 < nBases; base1++) {
                    itsPSFCrossTerms(base, base1).resize(this->itsNumberTerms, this->itsNumberTerms);
                    for (uInt term1 = 0; term1 < this->itsNumberTerms; term1++) {
                        for (uInt term2 = 0; term2 < this->itsNumberTerms; term2++) {
                            itsPSFCrossTerms(base, base1)(term1, term2).resize(subPsfShape);
                        }
                    }
                }
            }

            this->itsCouplingMatrix.resize(nBases);
            for (uInt base1 = 0; base1 < nBases; base1++) {
                itsCouplingMatrix(base1).resize(this->itsNumberTerms, this->itsNumberTerms);
                for (uInt base2 = base1; base2 < nBases; base2++) {
                    for (uInt term1 = 0; term1 < this->itsNumberTerms; term1++) {
                        for (uInt term2 = term1; term2 < this->itsNumberTerms; term2++) {
                            //          ASKAPLOG_DEBUG_STR(decmtbflogger, "Calculating convolutions of PSF("
                            //                << term1 << "+" << term2 << ") with basis functions");
                            work = conj(basisFunctionFFT.xyPlane(base1)) * basisFunctionFFT.xyPlane(base2) *
                                   subXFRVec(0) * conj(subXFRVec(term1 + term2)) / normPSF;
                            scimath::fft2d(work, false);
                            ASKAPLOG_DEBUG_STR(decmtbflogger, "Base(" << base1 << ")*Base(" << base2
                                                   << ")*PSF(" << term1 + term2
                                                   << ")*PSF(0): max = " << max(real(work))
                                                   << " min = " << min(real(work))
                                                   << " centre = " << real(work(subPsfPeak)));
                            // Remember that casa::Array reuses the same memory where possible so this
                            // apparent redundancy does not cause any memory bloat
                            itsPSFCrossTerms(base1, base2)(term1, term2) = real(work);
                            itsPSFCrossTerms(base2, base1)(term1, term2) = itsPSFCrossTerms(base1, base2)(term1, term2);
                            itsPSFCrossTerms(base1, base2)(term2, term1) = itsPSFCrossTerms(base1, base2)(term1, term2);
                            itsPSFCrossTerms(base2, base1)(term2, term1) = itsPSFCrossTerms(base1, base2)(term1, term2);
                            if (base1 == base2) {
                                itsCouplingMatrix(base1)(term1, term2) = real(work(subPsfPeak));
                                itsCouplingMatrix(base1)(term2, term1) = real(work(subPsfPeak));
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
                this->itsInverseCouplingMatrix(base).resize(this->itsNumberTerms, this->itsNumberTerms);
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Coupling matrix(" << base << ")="
                                       << this->itsCouplingMatrix(base));
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Calculating matrix inverse by Cholesky decomposition");
                invertSymPosDef(this->itsInverseCouplingMatrix(base),
                                this->itsDetCouplingMatrix(base), this->itsCouplingMatrix(base));
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Coupling matrix determinant(" << base << ") = "
                                       << this->itsDetCouplingMatrix(base));
                ASKAPLOG_DEBUG_STR(decmtbflogger, "Inverse coupling matrix(" << base
                                       << ")=" << this->itsInverseCouplingMatrix(base));
            }
            this->itsBasisFunctionChanged = False;
        }

        template<class T, class FT>
        bool DeconvolverMultiTermBasisFunction<T, FT>::deconvolve()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::deconvolve");
            this->initialise();

            if (this->control()->targetIter() != 0) {
                ASKAPLOG_INFO_STR(decmtbflogger, "Performing Multi-Term BasisFunction CLEAN for "
                                      << this->control()->targetIter() << " iterations");
                do {
                    this->oneIteration();
                    this->monitor()->monitor(*(this->state()));
                    this->state()->incIter();
                } while (!this->control()->terminate(*(this->state())));

                ASKAPLOG_INFO_STR(decmtbflogger, "Performed Multi-Term BasisFunction CLEAN for "
                                      << this->state()->currentIter() << " iterations");
                ASKAPLOG_INFO_STR(decmtbflogger, this->control()->terminationString());
            } else {
                ASKAPLOG_INFO_STR(decmtbflogger, "Bypassed Multi-Term BasisFunction CLEAN due to 0 iterations in the setup");
            }
            this->finalise();

            return True;
        }

        // This contains the heart of the Multi-Term BasisFunction Clean algorithm
        template<class T, class FT>
        void DeconvolverMultiTermBasisFunction<T, FT>::chooseComponent(uInt& optimumBase, casa::IPosition& absPeakPos,
                T& absPeakVal, Vector<T>& peakValues)
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction:::chooseComponent");
            const uInt nBases(this->itsResidualBasis.nelements());

            absPeakVal = 0.0;

            ASKAPDEBUGASSERT(peakValues.nelements() <= this->itsNumberTerms);

            // Find the base having the peak value in term=0
            // Here the weights image is used as a weight in the determination
            // of the maximum i.e. it finds the max in weight . residual. The values
            // returned are without the weight
            bool isWeighted((this->itsWeight.nelements() > 0) && (this->itsWeight(0).shape().nonDegenerate().conform(this->itsResidualBasis(0)(0).shape())));

            Vector<T> minValues(this->itsNumberTerms);
            Vector<T> maxValues(this->itsNumberTerms);

            for (uInt base = 0; base < nBases; base++) {

                // Find peak in residual image cube
                casa::IPosition minPos(2, 0);
                casa::IPosition maxPos(2, 0);
                T minVal(0.0), maxVal(0.0);

                // We implement various approaches to finding the peak. The first is the cheapest
                // and evidently the best (according to Urvashi).

                // Look for the maximum in term=0 for this base
                if (this->itsSolutionType == "MAXBASE") {
                    if (isWeighted) {
                        casa::minMaxMasked(minVal, maxVal, minPos, maxPos, this->itsResidualBasis(base)(0),
                                           this->itsWeight(0).nonDegenerate());
                    } else {
                        casa::minMax(minVal, maxVal, minPos, maxPos, this->itsResidualBasis(base)(0));
                    }
                    for (uInt term = 0; term < this->itsNumberTerms; term++) {
                        minValues(term) = this->itsResidualBasis(base)(term)(minPos);
                        maxValues(term) = this->itsResidualBasis(base)(term)(maxPos);
                    }
                    // In performing the search for the peak across bases, we want to take into account
                    // the SNR so we normalise out the coupling matrix for term=0 to term=0.
                    T norm(1 / sqrt(this->itsCouplingMatrix(base)(0, 0)));
                    maxVal *= norm;
                    minVal *= norm;
                } else {
                    // All these algorithms need the decoupled terms

                    // Decouple all terms using inverse coupling matrix
                    Vector<Array<T> > coefficients(this->itsNumberTerms);
                    for (uInt term1 = 0; term1 < this->itsNumberTerms; term1++) {
                        coefficients(term1).resize(this->dirty(0).shape().nonDegenerate());
                        coefficients(term1).set(T(0.0));
                        for (uInt term2 = 0; term2 < this->itsNumberTerms; term2++) {
                            coefficients(term1) = coefficients(term1)
                                                  + T(this->itsInverseCouplingMatrix(base)(term1, term2)) * this->itsResidualBasis(base)(term2);
                        }
                    }

                    if (this->itsSolutionType == "MAXTERM0") {
                        if (isWeighted) {
                            casa::minMaxMasked(minVal, maxVal, minPos, maxPos, coefficients(0),
                                               this->itsWeight(0).nonDegenerate());
                        } else {
                            casa::minMax(minVal, maxVal, minPos, maxPos, coefficients(0));
                        }
                        for (uInt term = 0; term < this->itsNumberTerms; term++) {
                            minValues(term) = coefficients(term)(minPos);
                            maxValues(term) = coefficients(term)(maxPos);
                        }
                    } else {
                        // MAXCHISQ
                        // Now form the criterion image and then search for the peak.
                        Array<T> negchisq(this->dirty(0).shape().nonDegenerate());
                        negchisq.set(T(0.0));
                        for (uInt term1 = 0; term1 < this->itsNumberTerms; term1++) {
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
                        if (isWeighted) {
                            casa::minMaxMasked(minVal, maxVal, minPos, maxPos, negchisq,
                                               this->itsWeight(0).nonDegenerate()*this->itsWeight(0).nonDegenerate());
                        } else {
                            casa::minMax(minVal, maxVal, minPos, maxPos, negchisq);
                        }
                        for (uInt term = 0; term < this->itsNumberTerms; term++) {
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
                }
                if (abs(maxVal) > absPeakVal) {
                    optimumBase = base;
                    absPeakVal = abs(maxVal);
                    absPeakPos = maxPos;
                }
            }

            // Now that we know the location of the peak found using one of the
            // above methods we can look up the values of the residuals. Remember
            // that we have to decouple the answer
            for (uInt term1 = 0; term1 < this->itsNumberTerms; ++term1) {
                peakValues(term1) = 0.0;
                for (uInt term2 = 0; term2 < this->itsNumberTerms; ++term2) {
                    peakValues(term1) +=
                        T(this->itsInverseCouplingMatrix(optimumBase)(term1, term2)) * this->itsResidualBasis(optimumBase)(term2)(absPeakPos);
                }
            }

            // Take square root to get value comparable to peak residual
            if (this->itsSolutionType == "MAXCHISQ") {
                absPeakVal = sqrt(max(T(0.0), absPeakVal));
            }
        }

        template<class T, class FT>
        bool DeconvolverMultiTermBasisFunction<T, FT>::oneIteration()
        {
            ASKAPTRACE("DeconvolverMultiTermBasisFunction::oneIteration");

            // For the psf convolutions, we only need a small part of the
            // basis functions so we recalculate for that size
            IPosition subPsfShape(this->findSubPsfShape());

            const uInt nBases(this->itsResidualBasis.nelements());

            casa::IPosition absPeakPos(2, 0);
            T absPeakVal(0.0);
            uInt optimumBase(0);
            Vector<T> peakValues(this->itsNumberTerms);
            chooseComponent(optimumBase, absPeakPos, absPeakVal, peakValues);

            // Report on progress
            // We want the worst case residual
            T absPeakRes = max(abs(peakValues));

            //      ASKAPLOG_INFO_STR(decmtbflogger, "All terms: absolute max = " << absPeakRes << " at " << absPeakPos);
            //      ASKAPLOG_INFO_STR(decmtbflogger, "Optimum base = " << optimumBase);

            if (this->state()->initialObjectiveFunction() == 0.0) {
                this->state()->setInitialObjectiveFunction(abs(absPeakVal));
            }
            this->state()->setPeakResidual(abs(absPeakRes));
            this->state()->setObjectiveFunction(abs(absPeakVal));
            this->state()->setTotalFlux(sum(this->model(0)));

            uInt nx(this->psf(0).shape()(0));
            uInt ny(this->psf(0).shape()(1));

            // Now we adjust model and residual for this component
            const casa::IPosition residualShape(this->dirty(0).shape().nonDegenerate());
            IPosition subPsfStart(2, nx / 2 - subPsfShape(0) / 2, ny / 2 - subPsfShape(1) / 2);
            IPosition subPsfEnd(2, nx / 2 + subPsfShape(0) / 2 - 1, ny / 2 + subPsfShape(1) / 2 - 1);
            IPosition subPsfStride(2, 1, 1);

            Slicer subPsfSlicer(subPsfStart, subPsfEnd, subPsfStride, Slicer::endIsLast);

            const casa::IPosition psfShape(2, this->itsBasisFunction->basisFunction().shape()(0),
                                           this->itsBasisFunction->basisFunction().shape()(1));

            casa::IPosition residualStart(2, 0), residualEnd(2, 0), residualStride(2, 1);
            casa::IPosition psfStart(2, 0), psfEnd(2, 0), psfStride(2, 1);

            const casa::IPosition modelShape(this->model(0).shape().nonDegenerate());
            casa::IPosition modelStart(2, 0), modelEnd(2, 0), modelStride(2, 1);

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

                modelStart(dim) = residualStart(dim);
                modelEnd(dim) = residualEnd(dim);
            }

            casa::Slicer psfSlicer(psfStart, psfEnd, psfStride, Slicer::endIsLast);
            casa::Slicer residualSlicer(residualStart, residualEnd, residualStride, Slicer::endIsLast);
            casa::Slicer modelSlicer(modelStart, modelEnd, modelStride, Slicer::endIsLast);

            // Add to model
            // We loop over all terms for the optimum base and ignore those terms with no flux
            for (uInt term = 0; term < this->itsNumberTerms; ++term) {
                if (abs(peakValues(term)) > 0.0) {
                    casa::Array<float> slice = this->model(term).nonDegenerate()(modelSlicer);
                    slice += this->control()->gain() * peakValues(term) *
                             Cube<T>(this->itsBasisFunction->basisFunction()).xyPlane(optimumBase).nonDegenerate()(psfSlicer);
                    this->itsTermBaseFlux(optimumBase)(term) += this->control()->gain() * peakValues(term);
                }
            }

            // Subtract PSFs, including base-base crossterms
            for (uInt term1 = 0; term1 < this->itsNumberTerms; term1++) {
                for (uInt term2 = 0; term2 < this->itsNumberTerms; term2++) {
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
