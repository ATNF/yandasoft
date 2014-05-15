/// @file Emptiness.tcc
/// @brief Emptiness operations as needed for Cornwell-Evans algorithm
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
#include <casa/Arrays/ArrayMath.h>
#include <casa/BasicMath/Math.h>

ASKAP_LOGGER(decemptinesslogger, ".deconvolution.emptiness");

#include <deconvolution/EntropyBase.h>

namespace askap {

    namespace synthesis {

        template<class T>
        T Emptiness<T>::entropy(const Array<T>& model)
        {
            ASKAPCHECK(this->itsScale > 0.0, "Scaling in Emptiness is invalid");
            return this->itsScale * sum((log(cosh((model - this->itsPrior) / this->itsScale)))) ;
        };

        template<class T>
        T Emptiness<T>::entropy(const Array<T>& model, const Array<T>& mask)
        {
            ASKAPCHECK(this->itsScale > 0.0, "Scaling in Emptiness is invalid");
            return this->itsScale * sum(mask *(log(cosh((model - this->itsPrior) / this->itsScale)))) ;
        };

        template<class T>
        void Emptiness<T>::gradEntropy(Array<T>& gradH, Array<T>& rHess, const Array<T>& model)
        {
            ASKAPCHECK(this->itsScale > 0.0, "Scaling in Emptiness is invalid");
            const T ggc = 2 * this->itsAlpha * this->itsQ;
            gradH.resize(model.shape());
            rHess.resize(model.shape());
            if (this->itsPrior.conform(model)) {
                gradH = -tanh((model - this->itsPrior) / this->itsScale);
            } else {
                gradH = -tanh(model / this->itsScale);
            }
            rHess = T(1.0) / ((T(1.0) - square(gradH)) / this->itsScale + ggc) ;
        }

        template<class T>
        void Emptiness<T>::gradEntropy(Array<T>& gradH, Array<T>& rHess, const Array<T>& model,
                                       const Array<T>& mask)
        {
            ASKAPCHECK(this->itsScale > 0.0, "Scaling in Emptiness is invalid");
            const T ggc = 2 * this->itsAlpha * this->itsQ;
            if (this->itsPrior.conform(model)) {
                gradH = -mask * tanh((model - this->itsPrior) / this->itsScale);
            } else {
                gradH = -mask * tanh(model / this->itsScale);
            }
            rHess = mask / ((T(1.0) - square(gradH)) / this->itsScale + ggc) ;
        }

    } // namespace synthesis

} // namespace askap
