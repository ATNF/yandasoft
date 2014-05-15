/// @file EntropyBase.tcc
/// @brief Entropy operations as needed for Cornwell-Evans algorithm
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
#include <boost/shared_ptr.hpp>
#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/BasicMath/Math.h>
ASKAP_LOGGER(decentropybaselogger, ".deconvolution.entropy");

using namespace casa;

namespace askap {

    namespace synthesis {

        template<class T>
        EntropyBase<T>::EntropyBase() : itsAlpha(T(0.0)), itsBeta(T(0.0)), itsQ(T(40.0)), itsScale(1.0),
                itsTolerance(0.3), itsUseFluxConstraint(false)
        {
        };

        template<class T>
        EntropyBase<T>::~EntropyBase()
        {
        };

        template<class T>
        T EntropyBase<T>::entropy(const Array<T>& model)
        {
            throw(AskapError("Called base class entropy"));
            ASKAPCHECK(model.shape().nelements(), "Model has no elements");
            return 0.0;
        };

        template<class T>
        T EntropyBase<T>::entropy(const Array<T>& model, const Array<T>& mask)
        {
            throw(AskapError("Called base class entropy"));
            ASKAPCHECK(model.shape().conform(mask.shape()), "Model and mask images have different shapes");
            return 0.0;
        };

        template<class T>
        void EntropyBase<T>::gradEntropy(Array<T>& gradH, Array<T>& rHess, const Array<T>& model)
        {
            throw(AskapError("Called base class gradEntropy"));
            // Silence the compiler by doing meaningless comparisons
            ASKAPCHECK(gradH.shape().conform(rHess.shape()), "Gradient and Hessian images have different shapes");
            ASKAPCHECK(gradH.shape().conform(model.shape()), "Gradient and model images have different shapes");
        }

        template<class T>
        void EntropyBase<T>::gradEntropy(Array<T>& gradH, Array<T>& rHess, const Array<T>& model,
                                         const Array<T>& mask)
        {
            throw(AskapError("Called base class gradEntropy"));
            // Silence the compiler by doing meaningless comparisons
            ASKAPCHECK(gradH.shape().conform(rHess.shape()), "Gradient and Hessian images have different shapes");
            ASKAPCHECK(gradH.shape().conform(model.shape()), "Gradient and model images have different shapes");
            ASKAPCHECK(gradH.shape().conform(mask.shape()), "Gradient and mask images have different shapes");
        }

        template<class T>
        void EntropyBase<T>::setScale(const T scale)
        {
            this->itsScale = scale;
        };

        template<class T>
        T EntropyBase<T>::formLength(const Matrix<T>& GDG)
        {
            return GDG(H, H) +  square(this->itsAlpha) * GDG(C, C) + square(this->itsBeta)*GDG(F, F);
        }

        template<class T>
        Matrix<T> EntropyBase<T>::formGDG(const Array<T>& model, const Array<T>& residual)
        {
            ASKAPCHECK(model.shape().conform(residual.shape()), "Model and residual images have different shapes");

            Matrix<T> GDG(4, 4);
            GDG.set(0.0);

            Array<T> rHess(model.shape());
            Array<T> gradH(model.shape());
            gradEntropy(gradH, rHess, model);
            Array<T> gradC = -T(2.0) * residual;
            GDG(H, H) = sum(gradH * rHess * gradH);
            GDG(H, C) = sum(gradH * rHess * gradC);
            GDG(H, F) = sum(gradH * rHess);
            GDG(C, C) = sum(gradC * rHess * gradC);
            GDG(C, F) = sum(gradC * rHess);
            GDG(F, F) = sum(rHess);
            GDG(H, J) = GDG(H, H) -  this->itsAlpha * GDG(H, C) - this->itsBeta * GDG(H, F);
            GDG(C, J) = GDG(H, C) -  this->itsAlpha * GDG(C, C) - this->itsBeta * GDG(C, F);
            GDG(F, J) = GDG(H, F) -  this->itsAlpha * GDG(C, F) - this->itsBeta * GDG(F, F);
            GDG(J, J) = GDG(H, H) +  square(this->itsAlpha) * GDG(C, C)
                        + square(this->itsBeta) * GDG(F, F)  + 2 * this->itsAlpha * this->itsBeta * GDG(C, F)
                        - 2 * this->itsAlpha * GDG(H, C) - 2 * this->itsBeta * GDG(H, F);
            return GDG;
        }

        template<class T>
        Matrix<T> EntropyBase<T>::formGDGStep(const Array<T>& model, const Array<T>& residual, Array<T>& step)
        {
            ASKAPCHECK(model.shape().conform(residual.shape()), "Model and residual images have different shapes");

            Matrix<T> GDG(4, 4);
            GDG.set(0.0);

            // Probably need to write as iterators soon
            Array<T> rHess(model.shape());
            Array<T> gradH(model.shape());
            gradEntropy(gradH, rHess, model);
            Array<T> gradC = - T(2.0) * residual;
            Array<T> gradJ = gradH - this->itsAlpha * gradC - this->itsBeta;
            step = rHess * gradJ;
            GDG(H, H) = sum(gradH * rHess * gradH);
            GDG(H, C) = sum(gradH * rHess * gradC);
            GDG(H, F) = sum(gradH * rHess);
            GDG(C, C) = sum(gradC * rHess * gradC);
            GDG(C, F) = sum(gradC * rHess);
            GDG(F, F) = sum(rHess);
            GDG(H, J) = GDG(H, H) -  this->itsAlpha * GDG(H, C) - this->itsBeta * GDG(H, F);
            GDG(C, J) = GDG(H, C) -  this->itsAlpha * GDG(C, C) - this->itsBeta * GDG(C, F);
            GDG(F, J) = GDG(H, F) -  this->itsAlpha * GDG(C, F) - this->itsBeta * GDG(F, F);
            GDG(J, J) = GDG(H, H) +  square(this->itsAlpha) * GDG(C, C)
                        + square(this->itsBeta) * GDG(F, F)  + 2 * this->itsAlpha * this->itsBeta * GDG(C, F)
                        - 2 * this->itsAlpha * GDG(H, C) - 2 * this->itsBeta * GDG(H, F);

            ASKAPCHECK(model.shape().conform(step.shape()), "Model and step images have different shapes");

            return GDG;
        }

        template<class T>
        T EntropyBase<T>::formGDS(const Array<T>& model, const Array<T>& residual, const Array<T>& step)
        {
            ASKAPCHECK(model.shape().conform(step.shape()), "Model and step images have different shapes");
            ASKAPCHECK(model.shape().conform(residual.shape()), "Model and residual images have different shapes");

            // Probably need to write as iterators soon
            Array<T> gradH(model.shape());
            Array<T> rHess(model.shape());
            gradEntropy(gradH, rHess, model);
            return - sum(step *(gradH + T(2.0) * this->itsAlpha * residual - this->itsBeta));
        };

        template<class T>
        void EntropyBase<T>::changeAlphaBeta(const Matrix<T>& GDG, const T targetChisq, const T chisq,
                                             const T targetFlux, const T flux)
        {
            T length = GDG(H, H) + square(this->itsAlpha) * GDG(C, C) + square(this->itsBeta) * GDG(F, F);
            if (this->itsAlpha == 0.0 && this->itsBeta == 0.0) {
                length = GDG(F, F);
            }
            T normGrad = GDG(J, J) / length;
            if (this->itsAlpha == 0.0) {
                normGrad = 0.0;
            }
            if (normGrad < this->itsTolerance) {
                this->updateAlphaBeta(GDG, targetChisq, chisq, targetFlux, flux);
            } else {
                this->initialiseAlphaBeta(GDG);
            }
        };

        template<class T>
        Bool EntropyBase<T>::initialiseAlphaBeta(const Matrix<T>& GDG)
        {
            if (this->itsAlpha > T(0.0)) {
                return false;
            }
            if (!this->itsUseFluxConstraint) {
                this->itsAlpha = max(0.0,  abs(GDG(H, C) / GDG(C, C)));
                this->itsBeta = 0.0;
            } else {
                Double det = GDG(C, C) * GDG(F, F) - GDG(C, F) * GDG(C, F);
                this->itsAlpha = (GDG(F, F) * GDG(H, C) - GDG(C, F) * GDG(H, F)) / det;
                this->itsBeta = (GDG(C, C) * GDG(H, F) - GDG(C, F) * GDG(H, C)) / det;
            }
            return true;
        };

        template<class T>
        void EntropyBase<T>::updateAlphaBeta(const Matrix<T>& GDG, const T targetChisq, const T chisq,
                                             const T targetFlux, const T flux)
        {

            // code stolen from SDE's mem.f
            T length(this->formLength(GDG));

            Double a = GDG(C, J) / GDG(C, C);
            Double b = square(a) - (GDG(J, J) - this->itsTolerance * length) / GDG(C, C);

            Double dAMax;
            Double dAMin;
            Double dBMax;
            Double dBMin;

            if (b > 0.0) {
                b = sqrt(b);
                dAMax = a + b;
                dAMin = a - b;
            } else {
                dAMax = 0.0;
                dAMin = 0.0;
            }

            Double dChisq;
            Double dAlpha;
            Double dBeta;
            Double dFlux;

            if (! this->itsUseFluxConstraint) {
                dChisq =  chisq -  targetChisq + GDG(C, J);
                dAlpha = dChisq / GDG(C, C);

                dAlpha = max(dAMin, min(dAMax, dAlpha));
                this->itsAlpha = max(0.0, this->itsAlpha + dAlpha);
            } else {
                a = GDG(F, J) / GDG(F, F);
                b = square(a) - (GDG(J, J) - this->itsTolerance * length) / GDG(F, F);
                if (b > 0.0) {
                    b = sqrt(b);
                    dBMax = a + b;
                    dBMin = a - b;
                } else {
                    dBMax = 0.0;
                    dBMin = 0.0;
                }

                dChisq = chisq - targetChisq + GDG(C, J);
                dFlux  = flux  - targetFlux  + GDG(F, J);
                T det = GDG(C, C) * GDG(F, F) - square(GDG(F, C));
                dAlpha = (GDG(F, F) * dChisq - GDG(C, F) * dFlux) / det;
                dBeta = (GDG(C, C) * dFlux  - GDG(C, F) * dChisq) / det;

                dAlpha = max(dAMin, min(dAMax, dAlpha));
                this->itsAlpha = max(0.0, this->itsAlpha + dAlpha);
                dBeta    = max(dBMin, min(dBMax, dBeta));
                this->itsBeta  = itsBeta + dBeta;
            }
        };

    } // namespace synthesis

} // namespace askap
