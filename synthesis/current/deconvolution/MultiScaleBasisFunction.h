/// @file 
/// @brief Holder fo multiscale functions
/// @details Holder for multiscale functions used in deconvolution
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

#ifndef ASKAP_SYNTHESIS_MULTISCALEBASISFUNCTION_H
#define ASKAP_SYNTHESIS_MULTISCALEBASISFUNCTION_H

#include <string>

#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casa/Arrays/Array.h>

#include <deconvolution/BasisFunction.h>

namespace askap {

    namespace synthesis {

        /// @brief Holder for multiscale basis functions used in MSClean
        /// @details The basis functions used here are those used in
        /// MSClean - invert parabolas in radius with a prolate
        /// spheroidal wave function multipled in. This looks
        /// Gaussian-ish but has the virtue of limited support.
        /// If you want to change the shape of the blobs, make
        /// another class to hold it.
        /// @ingroup Deconvolver

        template<typename T>
        class MultiScaleBasisFunction : public BasisFunction<T> {

            public:
                typedef boost::shared_ptr<MultiScaleBasisFunction<T> > ShPtr;

                /// @brief Construct
                /// @details Set up internals - shape not set yet
                /// @param[in] scales Scale (in pixels) of each blob
                /// @param[in] orthogonal Orthogonalise using Gram-Schmidt
                MultiScaleBasisFunction(const casa::Vector<casa::Float>& scales,
                                        const casa::Bool orthogonal = false);

                /// @brief Construct with specified shape
                /// @details Set up internals for specified shape
                /// @param[in] shape Shape of first two axes
                /// @param[in] scales Scale (in pixels) of each blob
                /// @param[in] orthogonal Orthogonalise using Gram-Schmidt
                MultiScaleBasisFunction(const casa::IPosition shape,
                                        const casa::Vector<casa::Float>& scales,
                                        const casa::Bool orthogonal = false);

                /// @brief Specify shape and fill in array
                /// @details The first two axes are set from shape, and the
                /// array is then filled in with the calculated values.
                /// @param[in] shape Shape of first two axes
                virtual void initialise(const casa::IPosition shape);

            private:
                /// Vector of scales (in pixels)
                casa::Vector<casa::Float> itsScales;

                /// Ancient routine (originally from F. Schwab) to calculate the PSWF.
                T spheroidal(T nu);
        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/MultiScaleBasisFunction.tcc>

#endif
