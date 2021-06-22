/// @file BasisFunction.tcc
/// @brief Basis Function for a  source
/// @details Holds basis function for a  source
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

#include <askap/askap_synthesis.h>

#include <askap/askap/AskapLogging.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Cube.h>

using namespace casa;

namespace askap {

    namespace synthesis {

        template<class T>
        BasisFunction<T>::BasisFunction() : itsNumberBases(1)
        {
        };

        template<class T>
        BasisFunction<T>::BasisFunction(const IPosition& shape) : itsNumberBases(1), itsOrthogonal(false)
        {
            initialise(shape);
        };

        template<class T>
        void BasisFunction<T>::initialise(const IPosition& shape)
        {
            ASKAPASSERT(itsNumberBases);
            ASKAPASSERT(shape.nelements() >= 2);
            const IPosition bfShape = shape.getFirst(2).concatenate(IPosition(1, itsNumberBases));
            itsBasisFunction.resize(bfShape);
            itsBasisFunction.set(T(0.0));
            // we probably need to patch the code to remove expectations of the 3rd dimension to be passed here
            itsShape = bfShape;
        };

        template<class T>
        void BasisFunction<T>::multiplyArray(const Matrix<Double>& A)
        {
            Array<T> mDataArray(this->shape());
            mDataArray.set(T(0.0));

            const uInt nRows(A.nrow());
            const uInt nCols(A.ncolumn());
            ASKAPDEBUGASSERT(this->itsBasisFunction.shape() >= 2);
            const uInt nx = this->itsBasisFunction.shape()(0);
            const uInt ny = this->itsBasisFunction.shape()(1);

            for (uInt j = 0; j < ny; j++) {
                for (uInt i = 0; i < nx; i++) {
                    for (uInt row = 0; row < nRows; row++) {
                        const IPosition currentPosRow(3, i, j, row);
                        for (uInt col = 0; col < nCols; col++) {
                            const IPosition currentPosCol(3, i, j, col);
                            mDataArray(currentPosRow) += T(A(row, col)) * this->itsBasisFunction(currentPosCol);
                        }
                    }
                }
            }
            this->itsBasisFunction = mDataArray.copy();
            Cube<T> BF(this->itsBasisFunction);
            for (uInt i = 0; i < this->itsNumberBases; i++) {
                const T sumBF(sum(BF.xyPlane(i)));
                if (abs(sumBF) > 0.0) {
                    BF.xyPlane(i) = BF.xyPlane(i) / sumBF;
                }
            }
        }

        /// @brief return requested basis function
        /// @details index index of the basis function to return [0..numberBases()-1]
        /// @return basis function of interest
        /// @note due to reference semantics of casa arrays it is possible to change the 
        /// basis function despite the fact that this method is const.
        template<typename T>
        casacore::Matrix<T> BasisFunction<T>::basisFunction(casacore::uInt index) const
        {
           ASKAPASSERT(index < itsNumberBases);
           casacore::Cube<T> tmpCube(itsBasisFunction);
           return tmpCube.xyPlane(index);
        }

        template<class T>
        void BasisFunction<T>::gramSchmidt(Array<T>& bf)
        {
            // Use reference semantics of Casa Arrays.
            Cube<T> v(bf.nonDegenerate());

            const uInt nScales = bf.shape()(2);

            for (uInt j = 0; j < nScales; j++) {
                for (uInt i = 0; i < j; i++) {
                    // remove component in direction vi
                    const T proj(sum(v.xyPlane(j)*v.xyPlane(i)) / sum(v.xyPlane(i)*v.xyPlane(i)));
                    v.xyPlane(j) = v.xyPlane(j) - proj * v.xyPlane(i);
                }
                // Normalise to rms==1.0
                v.xyPlane(j) = v.xyPlane(j) / sum(v.xyPlane(j) * v.xyPlane(j));
            }
            // Now normalise each plane to unit sum. The basis functions will remain
            // orthogonal but not orthonormal
            for (uInt j = 0; j < nScales; j++) {
                const T rNorm(T(1.0) / sum(v.xyPlane(j)));
                v.xyPlane(j) = v.xyPlane(j) * rNorm;
            }
        }

    } // namespace synthesis
} // namespace askap
