/// @file DeconvolverFiesta.tcc
/// @brief Class for a deconvolver based on the Fista
/// @details This concrete class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
/// @ingroup Deconvolver
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

#include <casacore/casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/STLIO.h>

#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(decfistalogger, ".deconvolution.fista");

#include <askap/deconvolution/DeconvolverFista.h>
#include <askap/deconvolution/MultiScaleBasisFunction.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/scimath/utils/ImageUtils.h>

namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver based on the Fista Clean
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a dirty image, psf optionally using a mask and a weights image.
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverFista<Double, DComplex>
        /// @ingroup Deconvolver

        template<class T, class FT>
        DeconvolverFista<T, FT>::~DeconvolverFista()
        {
        };

        template<class T, class FT>
        DeconvolverFista<T, FT>::DeconvolverFista(Vector<Array<T> >& dirty, Vector<Array<T> >& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf)
        {
            if (this->itsNumberDirtyTerms > 1) {
                throw(AskapError("FISTA deconvolver cannot perform multi-term deconvolutions"));
            }
        };

        template<class T, class FT>
        DeconvolverFista<T, FT>::DeconvolverFista(Array<T>& dirty, Array<T>& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf)
        {
        };

        template<class T, class FT>
        void DeconvolverFista<T, FT>::configure(const LOFAR::ParameterSet& parset)
        {
            DeconvolverBase<T, FT>::configure(parset);

            // Make the basis function
            {
                std::vector<float> defaultScales(3);
                defaultScales[0] = 0.0;
                defaultScales[1] = 10.0;
                defaultScales[2] = 30.0;
                std::vector<float> scales = parset.getFloatVector("scales", defaultScales);

                ASKAPLOG_INFO_STR(decfistalogger, "Constructing Multiscale basis function with scales " << scales);
                const Bool orthogonal = parset.getBool("orthogonal", false);
                if (orthogonal) {
                    ASKAPLOG_DEBUG_STR(decfistalogger, "Multiscale basis functions will be orthogonalised");
                }

                // MV: a bit of technical debt highlighted by casacore's interface change. In principle, we could've
                // had scales as std::vector in the interface to avoid the explicit construction (in this particular case,
                // there is no benefit of using casacore::Vector) + similar code repeats a few times, we could have a better reuse
                itsBasisFunction = BasisFunction<Float>::ShPtr(new MultiScaleBasisFunction<Float>(casacore::Vector<casacore::Float>(scales),
                                   orthogonal));
            }
        }

        template<class T, class FT>
        void DeconvolverFista<T, FT>::initialise()
        {
            DeconvolverBase<T, FT>::initialise();

            if (itsBasisFunction) {
                this->itsBasisFunction->initialise(this->model().shape());
                itsBasisFunctionTransform.resize(itsBasisFunction->shape());
                itsBasisFunctionTransform.set(0.);

                // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
                // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
                // but one has to generate basis functions with the right order (or reorder them in a separate routine - leave it for the future)
                // Limit number of fft threads to 8 (more is slower for our fft sizes)
                scimath::FFT2DWrapper<FT> fft2d(true,8);

                // do explicit loop over basis functions here (the original code relied on iterator in old
                // fft2d and, therefore, low level representation of the basis function stack). This way
                // we have more control over the array structure and can transition to the more efficient order
                const casacore::uInt nBases = itsBasisFunction->numberBases();
                for (uInt base = 0; base < nBases; ++base) {
                     casacore::Matrix<FT> fftBuffer = itsBasisFunctionTransform.xyPlane(base);
                     casacore::setReal(fftBuffer, this->itsBasisFunction->basisFunction(base));
                     fft2d(fftBuffer, true);
                }
            }

            ASKAPLOG_INFO_STR(decfistalogger, "Initialised FISTA solver");
        }

        template<class T, class FT>
        bool DeconvolverFista<T, FT>::deconvolve()
        {
            this->initialise();

            bool isMasked(this->itsWeight.nelements());
            if (!this->itsWeight(0).shape().conform(this->dirty().shape())) isMasked = false;

            Array<T> X, X_old, X_temp;

            X_temp.resize(this->model().shape());
            X_temp.set(T(0.0));

            X_old.resize(this->model().shape());
            X_old.set(T(0.0));

            X.resize(this->model().shape());
            X = this->model().copy();

            T absPeakVal;

            ASKAPLOG_INFO_STR(decfistalogger, "Performing Fista for " << this->control()->targetIter() << " iterations");

            this->updateResiduals(X);

            X_temp = X.copy();

            absPeakVal = max(abs(this->dirty()));

            //      T effectiveLambda(absPeakVal*this->control()->fractionalThreshold()+this->control()->lambda());
            T lambda(absPeakVal*(1.0 - this->control()->gain()));
            ASKAPLOG_INFO_STR(decfistalogger, "Effective lambda = " << lambda);

            T t_new = 1;

            T lipschitz(10.0);

            do {
                X_old = X_temp.copy();
                T t_old = t_new;

                this->updateResiduals(X);
                scimath::saveAsCasaImage("residuals.tab", this->dirty());

                X = X + this->dirty() / lipschitz;

                // Transform to other (e.g. multiscale) space
                Array<T> WX;
                scimath::saveAsCasaImage("X.tab", X);
                this->W(WX, X);
                scimath::saveAsCasaImage("W.tab", WX);

                // Now shrink the coefficients towards zero and clip those below
                // lambda/lipschitz.
                Array<T> shrink(WX.shape());

                Array<T> truncated(abs(WX) - lambda / lipschitz);
                shrink = truncated(truncated > T(0.0));
                shrink = sign(WX) * shrink;
                shrink(truncated < T(0.0)) = T(0.0);

                scimath::saveAsCasaImage("shrink.tab", WX);
                // Transform back from other (e.g. wavelet) space here

                this->WT(X_temp, shrink);
                scimath::saveAsCasaImage("WT.tab", X_temp);

                t_new = (T(1.0) + sqrt(T(1.0) + T(4.0) * square(t_old))) / T(2.0);
                X = X_temp + ((t_old - T(1.0)) / t_new) * (X_temp - X_old);
                {
                    casacore::IPosition minPos;
                    casacore::IPosition maxPos;
                    T minVal(0.0), maxVal(0.0);
                    if (isMasked) {
                        casacore::minMaxMasked(minVal, maxVal, minPos, maxPos,
                                           this->dirty(),
                                           this->itsWeight(0));
                    } else {
                        casacore::minMax(minVal, maxVal, minPos, maxPos, this->dirty());
                    }

                    ASKAPLOG_INFO_STR(decfistalogger, "   Maximum = " << maxVal << " at location " << maxPos);
                    ASKAPLOG_INFO_STR(decfistalogger, "   Minimum = " << minVal << " at location " << minPos);
                    if (abs(minVal) < abs(maxVal)) {
                        absPeakVal = abs(maxVal);
                    } else {
                        absPeakVal = abs(minVal);
                    }
                }

                T l1Norm = sum(abs(X_temp));
                T fit = casacore::sum(this->dirty() * this->dirty());
                T objectiveFunction(fit + lambda*l1Norm);
                this->state()->setPeakResidual(absPeakVal);
                this->state()->setObjectiveFunction(objectiveFunction);
                this->state()->setTotalFlux(sum(X_temp));

                if (absPeakVal < lambda) {
                    lambda *= 1.0 - this->control()->gain();
                    ASKAPLOG_INFO_STR(decfistalogger, "Setting Lagrange multiplier lambda = " << lambda);
                }
                this->monitor()->monitor(*(this->state()));
                this->state()->incIter();
            } while (!this->control()->terminate(*(this->state())));
            this->model() = X_temp.copy();
            this->updateResiduals(this->model());

            ASKAPLOG_INFO_STR(decfistalogger, "Performed Fista for " << this->state()->currentIter() << " iterations");

            ASKAPLOG_INFO_STR(decfistalogger, this->control()->terminationString());

            this->finalise();

            absPeakVal = casacore::max(casacore::abs(this->dirty()));

            this->state()->setPeakResidual(absPeakVal);
            this->state()->setObjectiveFunction(absPeakVal);

            return True;
        }

        template<class T, class FT>
        void DeconvolverFista<T, FT>::setBasisFunction(boost::shared_ptr<BasisFunction<T> > bf)
        {
            itsBasisFunction = bf;
        };

        template<class T, class FT>
        boost::shared_ptr<BasisFunction<T> > DeconvolverFista<T, FT>::basisFunction()
        {
            return itsBasisFunction;
        };

        // Apply the convolve operation - this is undecimated and redundant. The 2D image
        // will be expanded along the third axis. Since the basis functions are normalised
        // to unit volume, the values here should be roughly maintained.
        template <class T, class FT>
        void DeconvolverFista<T, FT>::W(Array<T>& out, const Array<T>& in)
        {
            if (itsBasisFunction) {
                // MV: the original code was quite untidy and with technical debt regarding degenerate dimensions and
                // taking shape. I have to change this as the new FFT wrapper interface gets the Matrix for 2D FFT rather than
                // an array with built in iteration over the third or any further dimensions. As we don't use this deconvolver now,
                // the change won't be well tested (+ I'm not sure it ever worked satisfactory, this is largely an experimental code).
                // See change sets associated with AXA-1031 for further details.
                casacore::Matrix<FT> inTransform(in.shape().nonDegenerate(2), static_cast<FT>(0.));
                casacore::Matrix<FT> outPlaneTransform(inTransform.shape());
                out.resize(itsBasisFunction->shape());
                casacore::Cube<T> outCube(out);
                casacore::setReal(inTransform, in.nonDegenerate(2));
                // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
                // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
                // but one has to generate basis functions with the right order (or reorder them in a separate routine - leave it for the future)
                // Limit number of fft threads to 8 (more is slower for our fft sizes)
                scimath::FFT2DWrapper<FT> fft2d(true,8);

                fft2d(inTransform, true);
                const uInt nPlanes(itsBasisFunction->shape()(2));
                for (uInt plane = 0; plane < nPlanes; plane++) {
                    outPlaneTransform = inTransform * itsBasisFunctionTransform.xyPlane(plane);
                    fft2d(outPlaneTransform, false);
                    outCube.xyPlane(plane) = real(outPlaneTransform);
                }
            } else {
                out = in.copy();
            }
        }

        // Apply the transpose of the W operation - this is undecimated and redundant so we need
        // to sum over the planes
        template <class T, class FT>
        void DeconvolverFista<T, FT>::WT(Array<T>& out, const Array<T>& in)
        {
            if (itsBasisFunction) {
                const Cube<T> inCube(in);
                // MV: the original code was quite untidy and with technical debt regarding degenerate dimensions and
                // taking shape. I have to change this as the new FFT wrapper interface gets the Matrix for 2D FFT rather than
                // an array with built in iteration over the third or any further dimensions. As we don't use this deconvolver now,
                // the change won't be well tested (+ I'm not sure it ever worked satisfactory, this is largely an experimental code).
                // See change sets associated with AXA-1031 for further details.
                casacore::Matrix<FT> inPlaneTransform(out.shape().nonDegenerate(2),static_cast<FT>(0.));
                casacore::Matrix<FT> outTransform(inPlaneTransform.shape());
                outTransform.set(FT(0.0));

                // To reconstruct, we filter out each basis from the cumulative sum
                // and then add the corresponding term from the in array.
                const uInt nPlanes(itsBasisFunction->shape()(2));

                casacore::setReal(inPlaneTransform, inCube.xyPlane(nPlanes - 1));
                // Do harmonic reorder as with the original wrapper (hence, pass true to the wrapper), it may be possible to
                // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
                // but one has to generate basis functions with the right order (or reorder them in a separate routine - leave it for the future)
                // Limit number of fft threads to 8 (more is slower for our fft sizes)
                scimath::FFT2DWrapper<FT> fft2d(true,8);

                fft2d(inPlaneTransform, true);
                outTransform = itsBasisFunctionTransform.xyPlane(nPlanes - 1) * (inPlaneTransform);

                for (uInt plane = 1; plane < nPlanes; plane++) {
                    inPlaneTransform.set(0.);
                    casacore::setReal(inPlaneTransform, inCube.xyPlane(nPlanes - 1 - plane));
                    fft2d(inPlaneTransform, true);
                    outTransform = outTransform
                                   + itsBasisFunctionTransform.xyPlane(nPlanes - 1 - plane) * (inPlaneTransform - outTransform);
                }
                fft2d(outTransform, false);
                out.nonDegenerate(2) = real(outTransform);
            } else {
                out = in.copy();
            }
        }

    } // namespace synthesis

} // namespace askap
