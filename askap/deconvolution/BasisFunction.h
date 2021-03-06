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

#include <casacore/casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <casacore/casa/Arrays/Array.h>

namespace askap {

    namespace synthesis {

        /// @brief Base class for basis functions used in deconvolution
        /// @details Basis functions can be used during deconvolution
        /// to represent the emission. Typical examples are Points (for Hogbom)
        /// and MultiScale (for BasisFunctionClean).
        /// @ingroup Deconvolver

        template<typename T> class BasisFunction : public boost::noncopyable {

            public:
                typedef boost::shared_ptr<BasisFunction<T> > ShPtr;

                virtual ~BasisFunction() {};

                /// @brief Construct empty basis function
                BasisFunction();

                /// @brief Construct from a specified shape
                /// param[in] shape Shape of desired basis function
                explicit BasisFunction(const casacore::IPosition& shape);

                /// @brief Initialise from a specified shape (actually only first two axes)
                /// @detail Set the shape of the basis function and fill in the actual
                /// function.
                /// param[in] shape Shape of desired basis function on the first two axes.
                virtual void initialise(const casacore::IPosition& shape);

                /// @brief Return the number of bases in the basis function
                casacore::uInt numberBases() const {return itsNumberBases;};

                /// @brief shape of the basis
                /// @return shape of the basis function itself
                casacore::IPosition shape() const {return itsShape;}

                /// @brief Return the basis function as an array
                /// @details The basis function is returned as an array
                /// of shape (nx, ny, nbases) where nx, ny are the
                /// lengths of the first two axes.
                casacore::Array<T>& allBasisFunctions() {return itsBasisFunction;};

                /// @brief return requested basis function
                /// @details index index of the basis function to return [0..numberBases()-1]
                /// @return basis function of interest
                /// @note due to reference semantics of casa arrays it is possible to change the 
                /// basis function despite the fact that this method is const.
                casacore::Matrix<T> basisFunction(casacore::uInt index) const;

                /// @brief Multiply by a matrix on the third dimension
                void multiplyArray(const casacore::Matrix<casacore::Double>& arr);

            protected:
                /// @brief Orthogonalise using Gram Schmidt algorithm
                /// @param[in] bf cube to work with, 3rd dimension is the basis function index
                static void gramSchmidt(casacore::Array<T>& bf);

            private:
                casacore::IPosition itsShape;
            protected:
                casacore::Array<T> itsBasisFunction;
                casacore::uInt itsNumberBases;
                casacore::Bool itsOrthogonal;
        };

    } // namespace synthesis

} // namespace askap

#include <askap/deconvolution/BasisFunction.tcc>

#endif
