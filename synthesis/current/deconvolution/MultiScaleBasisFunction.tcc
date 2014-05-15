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
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>
#include <fft/FFTWrapper.h>

ASKAP_LOGGER(decmsbaselogger, ".deconvolution.multiscalebasisfunction");

#include <deconvolution/MultiScaleBasisFunction.h>

using namespace casa;

namespace askap {

    namespace synthesis {

        template<class T>
        MultiScaleBasisFunction<T>::MultiScaleBasisFunction(const Vector<Float>& scales,
                const Bool orthogonal) :
                BasisFunction<T>::BasisFunction(), itsScales(scales)
        {
            this->itsOrthogonal = orthogonal;
            this->itsNumberBases = scales.nelements();
            ASKAPLOG_INFO_STR(decmsbaselogger, "Initialising multiscale basis function with "
                                  << this->itsNumberBases << " scales: " << scales);
        }

        template<class T>
        MultiScaleBasisFunction<T>::MultiScaleBasisFunction(const IPosition shape,
                const Vector<Float>& scales,
                const Bool orthogonal) :
                BasisFunction<T>::BasisFunction(), itsScales(scales)
        {
            this->itsOrthogonal = orthogonal;
            this->itsNumberBases = scales.nelements();
            initialise(shape);
        }

        template<class T>
        void MultiScaleBasisFunction<T>::initialise(const IPosition shape)
        {
            BasisFunction<T>::initialise(shape);

            ASKAPLOG_INFO_STR(decmsbaselogger, "Calculating multiscale basis functions");
            // Make a cube and use referencing in Array
            casa::Cube<T> scaleCube(this->itsBasisFunction);

            for (uInt scale = 0; scale < itsScales.nelements(); scale++) {
                const Float scaleSize = itsScales(scale);
                ASKAPCHECK(scaleSize >= 0.0, "Scale size " << scale << " is not positive " << scaleSize);

                if (scaleSize < 1e-6) {
                    scaleCube(shape[0] / 2, shape[1] / 2, 0) = T(1.0);
                } else {
                    const Int nx = shape[0];
                    const Int ny = shape[1];
                    const Int mini = max(0, (Int)(nx / 2 - scaleSize));
                    const Int maxi = min(nx - 1, (Int)(nx / 2 + scaleSize));
                    const Int minj = max(0, (Int)(ny / 2 - scaleSize));
                    const Int maxj = min(ny - 1, (Int)(ny / 2 + scaleSize));

                    Float ypart = 0.0;
                    Float volume = 0.0;
                    Float rad2 = 0.0;
                    Float rad = 0.0;

                    for (Int j = minj; j <= maxj; j++) {
                        ypart = square((ny / 2 - (Double)(j)) / scaleSize);
                        for (Int i = mini; i <= maxi; i++) {
                            rad2 =  ypart + square((nx / 2 - (Double)(i)) / scaleSize);
                            if (rad2 < 1.0) {
                                if (rad2 <= 0.0) {
                                    rad = 0.0;
                                } else {
                                    rad = sqrt(rad2);
                                }
                                T value;
                                value = (1.0 - rad2) * spheroidal(rad);
                                scaleCube(i, j, scale) = value;
                                volume += value;
                            }
                        }
                    }
                    scaleCube.xyPlane(scale) = scaleCube.xyPlane(scale) / volume;
                }
            }
            if (this->itsOrthogonal) {
                ASKAPLOG_INFO_STR(decmsbaselogger, "Orthogonalising multi-scale basis function using Gram-Schmidt algorithm");
                this->gramSchmidt(this->itsBasisFunction);
            }
        }

        // Calculate the spheroidal function
        template<class T>
        T MultiScaleBasisFunction<T>::spheroidal(T nu)
        {

            if (nu <= 0) {
                return T(1.0);
            } else if (nu >= 1.0) {
                return T(0.0);
            } else {
                uInt np = 5;
                uInt nq = 3;
                Matrix<float> p(np, 2);
                Matrix<float> q(nq, 2);
                p(0, 0) = 8.203343e-2;
                p(1, 0) = -3.644705e-1;
                p(2, 0) =  6.278660e-1;
                p(3, 0) = -5.335581e-1;
                p(4, 0) =  2.312756e-1;
                p(0, 1) =  4.028559e-3;
                p(1, 1) = -3.697768e-2;
                p(2, 1) = 1.021332e-1;
                p(3, 1) = -1.201436e-1;
                p(4, 1) = 6.412774e-2;
                q(0, 0) = 1.0000000e0;
                q(1, 0) = 8.212018e-1;
                q(2, 0) = 2.078043e-1;
                q(0, 1) = 1.0000000e0;
                q(1, 1) = 9.599102e-1;
                q(2, 1) = 2.918724e-1;
                uInt part = 0;
                Float nuend = 0.0;
                if (nu >= 0.0 && nu < 0.75) {
                    part = 0;
                    nuend = 0.75;
                } else if (nu >= 0.75 && nu <= 1.00) {
                    part = 1;
                    nuend = 1.0;
                }

                Float top = p(0, part);
                Float delnusq = pow(nu, 2.0) - pow(nuend, 2.0);
                uInt k;
                for (k = 1; k < np; k++) {
                    top += p(k, part) * pow(delnusq, (Float)k);
                }
                Float bot = q(0, part);
                for (k = 1; k < nq; k++) {
                    bot += q(k, part) * pow(delnusq, (Float)k);
                }

                if (bot != 0.0) {
                    return T(top / bot);
                } else {
                    return T(0.0);
                }
            }
        }

    } // namespace synthesis

} // namespace askap
