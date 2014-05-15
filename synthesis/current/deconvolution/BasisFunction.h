/// @file BasisFunction.h
/// @brief Base class for Basis functions
/// @details Holder for basis functions used in deconvolution
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

#ifndef ASKAP_SYNTHESIS_BASISFUNCTION_H
#define ASKAP_SYNTHESIS_BASISFUNCTION_H

#include <string>

#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casa/Arrays/Array.h>

namespace askap {

    namespace synthesis {

        /// @brief Base class for basis functions used in deconvolution
        /// @details Basis functions can be used during deconvolution
        /// to represent the emission. Typical examples are Points (for Hogbom)
        /// and MultiScale (for BasisFunctionClean).
        /// @ingroup Deconvolver

        template<typename T> class BasisFunction {

            public:
                typedef boost::shared_ptr<BasisFunction<T> > ShPtr;

                virtual ~BasisFunction() {};

                /// @brief Construct empty basis function
                BasisFunction();

                /// @brief Construct from a specified shape
                /// param[in] shape Shape of desired basis function
                BasisFunction(const casa::IPosition shape);

                /// @brief Initialise from a specified shape (actually only first two axes)
                /// @detail Set the shape of the basis function and fill in the actual
                /// function.
                /// param[in] shape Shape of desired basis function on the first two axes.
                virtual void initialise(const casa::IPosition shape);

                /// @brief Return the number of bases in the basis function
                casa::uInt numberBases() const {return itsNumberBases;};

                /// @brief Return the basis function as an array
                /// @details The basis function is returned as an array
                /// of shape (nx, ny, nbases) where nx, ny are the
                /// lengths of the first two axes.
                virtual casa::Array<T>& basisFunction() {return itsBasisFunction;};

                /// @brief Multiply by a matrix on the third dimension
                virtual void multiplyArray(const casa::Matrix<casa::Double>& arr);

                /// @brief Orthogonalise using Gram Schmidt algorithm
                void gramSchmidt(casa::Array<T>& bf);

            protected:
                casa::IPosition itsShape;
                casa::Array<T> itsBasisFunction;
                casa::uInt itsNumberBases;
                casa::Bool itsOrthogonal;
        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/BasisFunction.tcc>

#endif
