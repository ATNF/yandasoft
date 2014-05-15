/// @file EntrophyBase.h
/// @brief Base class for a deconvolver
/// @details This interface class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
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

#ifndef ASKAP_SYNTHESIS_ENTROPYBASE_H
#define ASKAP_SYNTHESIS_ENTROPYBASE_H

#include <string>

#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Array.h>
#include <Common/ParameterSet.h>


namespace askap {

    namespace synthesis {

        // Base class
        template<class T>
        class EntropyBase {
            public:

                typedef boost::shared_ptr<EntropyBase<T> > ShPtr;

                enum GRADTYPE {H = 0, C, F, J };

                // Construct the basis class
                EntropyBase();

                // A virtual destructor may be necessary for use in derived classes.
                virtual ~EntropyBase();

                // Set the current scaling
                void setScale(T scale);

                // calculate the entropy for the whole image
                virtual T entropy(const casa::Array<T>& model);

                // calculate the entropy for the whole image
                virtual T entropy(const casa::Array<T>& model, const casa::Array<T>& mask);

                // calculate the gradient entropy for the whole image
                virtual void gradEntropy(casa::Array<T>& gradH, casa::Array<T>& rHess,
                                         const casa::Array<T>& model);

                // calculate the gradient entropy for the whole image
                virtual void gradEntropy(casa::Array<T>& gradH, casa::Array<T>& rHess,
                                         const casa::Array<T>& model, const casa::Array<T>& mask);

                // Form length
                T formLength(const casa::Matrix<T>& GDG);

                // calculate the Gradient dot Gradient matrix
                casa::Matrix<T> formGDG(const casa::Array<T>& model, const casa::Array<T>& residual);

                // calculate the Gradient dot Gradient matrix, calculate Step
                casa::Matrix<T> formGDGStep(const casa::Array<T>& model, const casa::Array<T>& residual,
                                      casa::Array<T>& step);

                // calculate Gradient dot Step
                T formGDS(const casa::Array<T>& model, const casa::Array<T>& residual, const casa::Array<T>& step);

                void setAlpha(const T alpha) {itsAlpha = alpha;};

                void setBeta(const T beta) {itsBeta = beta;};

                T alpha() {return itsAlpha;};

                T beta() {return itsBeta;};

                void setQ(const T Q) {itsQ = Q;};

                void setTolerance(const T tolerance) {itsTolerance = tolerance;};

                void setFluxConstraint(const casa::Bool useFluxConstraint) {itsUseFluxConstraint = useFluxConstraint;};

                void setMask(const casa::Array<T>& mask) {itsMask = mask.copy();};

                void setPrior(const casa::Array<T>& prior) {itsPrior = prior.copy();};

                void changeAlphaBeta(const casa::Matrix<T>& GDG, const T targetChisq, const T chisq,
                                     const T targetFlux, const T flux);

                void updateAlphaBeta(const casa::Matrix<T>& GDG, const T targetChisq, const T chisq,
                                     const T targetFlux, const T flux);

                // Do we need to initialise alpha and beta? If so, do so.
                casa::Bool initialiseAlphaBeta(const casa::Matrix<T>& GDG);

            protected:

                T itsAlpha;

                T itsBeta;

                T itsQ;

                T itsScale;

                T itsTolerance;

                casa::Bool itsUseFluxConstraint;

                casa::Array<T> itsMask;

                casa::Array<T> itsPrior;

        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/EntropyBase.tcc>

#endif
