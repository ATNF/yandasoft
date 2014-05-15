/// @file
/// @brief Basis Function for a point source
/// @details Holds basis function for a point source
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

#include <askap_synthesis.h>

#include <askap/AskapLogging.h>
#include <casa/aips.h>
ASKAP_LOGGER(decpointbaselogger, ".deconvolution.pointbasis");

#include <deconvolution/PointBasisFunction.h>

using namespace casa;

namespace askap {

    namespace synthesis {

        template<class T>
        PointBasisFunction<T>::PointBasisFunction() :
                BasisFunction<T>::BasisFunction()
        {
        };

        template<class T>
        PointBasisFunction<T>::PointBasisFunction(const IPosition shape) :
                BasisFunction<T>::BasisFunction()
        {
            initialise(shape);
        };

        template<class T>
        void PointBasisFunction<T>::initialise(const IPosition shape)
        {
            const IPosition centre(3, shape[0] / 2, shape[1] / 2, 0);
            BasisFunction<T>::itsBasisFunction.resize(IPosition(3, shape[0], shape[1], 1));
            BasisFunction<T>::itsBasisFunction(centre) = T(1.0);
        };

    } // namespace synthesis

} // namespace askap
