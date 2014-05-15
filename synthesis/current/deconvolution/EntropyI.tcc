/// @file EntropyI.tcc
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

#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/BasicMath/Math.h>
#include <askap/AskapLogging.h>

ASKAP_LOGGER(decentropyilogger, ".deconvolution.entropyi");

using namespace casa;

namespace askap {

    namespace synthesis {

        template<class T>
        T EntropyI<T>::entropy(const Array<T>& model)
        {
            const T flux = sum(model);
            if (this->itsPrior.conform(model)) {
                return - sum(model*log(model)) / flux + log(model.shape().product());
            } else {
                return - sum(model*log(model)) / flux + log(model.shape().product());
            }

        }

        template<class T>
        T EntropyI<T>::entropy(const Array<T>& model, const Array<T>& mask)
        {
            const T flux = sum(mask * model);
            if (this->itsPrior.conform(model)) {
                return - sum(mask*model*log(model)) / flux + log(model.shape().product());
            } else {
                return - sum(mask*model*log(model)) / flux + log(model.shape().product());
            }

        }

        template<class T>
        void EntropyI<T>::gradEntropy(Array<T>& gradH, Array<T>& rHess, const Array<T>& model)
        {
            const T ggc = 2 * this->itsAlpha * this->itsQ;
            gradH.resize(model.shape());
            rHess.resize(model.shape());
            if (this->itsPrior.conform(model)) {
                gradH = -log(model / this->itsPrior);
            } else {
                gradH = -log(model);
            }
            rHess = model / (T(1.0) + ggc * model);
        }

        template<class T>
        void EntropyI<T>::gradEntropy(Array<T>& gradH, Array<T>& rHess, const Array<T>& model,
                                      const Array<T>& mask)
        {
            const T ggc = 2 * this->itsAlpha * this->itsQ;
            gradH.resize(model.shape());
            rHess.resize(model.shape());
            if (this->itsPrior.conform(model)) {
                gradH = -mask * log(model / this->itsPrior);
            } else {
                gradH = -mask * log(model);
            }
            rHess = mask * model / (T(1.0) + ggc * model);
        }

    } // namespace synthesis

} // namespace askap
