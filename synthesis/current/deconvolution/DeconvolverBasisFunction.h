/// @file DeconvolverBasisFunction.h
/// @brief Class for a BasisFunction-Clean-based deconvolver
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

#ifndef ASKAP_SYNTHESIS_DECONVOLVERBASISFUNCTION_H
#define ASKAP_SYNTHESIS_DECONVOLVERBASISFUNCTION_H

#include <string>

#include <boost/shared_ptr.hpp>
#include <casa/aips.h>
#include <casa/Arrays/Array.h>

#include <deconvolution/DeconvolverBase.h>
#include <deconvolution/DeconvolverState.h>
#include <deconvolution/DeconvolverControl.h>
#include <deconvolution/DeconvolverMonitor.h>
#include <deconvolution/BasisFunction.h>

namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver using the BasisFunction Clean algorithm
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a dirty image, psf optionally using a mask and a weights image.
        /// This algorithm is similar to the MultiScale Clean (Cornwell 2009) with changes
        /// to improve performance and flexibility.
        ///
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverBasisFunction<Double, DComplex>
        /// @ingroup Deconvolver
        template<class T, class FT>
        class DeconvolverBasisFunction : public DeconvolverBase<T, FT> {

            public:
                typedef boost::shared_ptr<DeconvolverBasisFunction<T, FT> > ShPtr;

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We use a vector
                /// here to allow multiple dirty images and PSFs for the
                /// same model (e.g. as in MFS)
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverBasisFunction(casa::Vector<casa::Array<T> >& dirty,
                                         casa::Vector<casa::Array<T> >& psf);

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We keep this
                /// version for compatibility
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverBasisFunction(casa::Array<T>& dirty, casa::Array<T>& psf);

                virtual ~DeconvolverBasisFunction();

                /// @brief Set the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction.
                /// @param[in] bf Shared pointer to basisfunction instance
                void setBasisFunction(boost::shared_ptr<BasisFunction<T> > bf);

                /// @brief Return the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction
                boost::shared_ptr<BasisFunction<T> > basisFunction();

                /// @brief Perform the deconvolution
                /// @detail This is the main deconvolution method.
                virtual bool deconvolve();

                /// @brief Initialize the deconvolution
                /// @detail Initialise e.g. set weighted mask
                virtual void initialise();

                /// @brief Finalise the deconvolution
                /// @detail Finalise the deconvolution
                virtual void finalise();

                /// @brief configure basic parameters of the solver
                /// @details This method encapsulates extraction of basic solver parameters from the parset.
                /// @param[in] parset parset
                virtual void configure(const LOFAR::ParameterSet &parset);

            private:

                /// @brief Perform the deconvolution
                /// @detail This is the main deconvolution method.
                bool oneIteration();

                void initialisePSF();

                void initialiseResidual();

                void minMaxMaskedScales(T& minVal, T& maxVal,
                                        casa::IPosition& minPos, casa::IPosition& maxPos,
                                        const casa::Array<T>& dataArray,
                                        const casa::Array<T>& maskArray);

                // Find the coefficients for each scale by applying the
                // inverse of the coupling matrix
                casa::Vector<T> findCoefficients(const casa::Matrix<casa::Double>& invCoupling,
                                                 const casa::Vector<T>& peakValues);

                casa::Array<T> applyInverse(const casa::Matrix<casa::Double>& invCoupling,
                                            const casa::Array<T> dataArray);

                casa::Array<T> apply(const casa::Matrix<casa::Double>& invCoupling,
                                     const casa::Vector<T> dataVector);

                void gramSchmidt(casa::Array<T>& bf);

                /// Residual images convolved with basis functions
                casa::Array<T> itsResidualBasisFunction;

                /// Point spread functions convolved with basis functions
                casa::Array<T> itsPSFBasisFunction;

                /// Use cross terms in the source removal step?
                casa::Bool itsUseCrossTerms;

                /// The coupling between different scales.
                casa::Matrix<casa::Double> itsCouplingMatrix;

                /// Inverse of the coupling matrix
                casa::Matrix<casa::Double> itsInverseCouplingMatrix;

                /// Determinant of the coupling Matrix
                casa::Double itsDetCouplingMatrix;

                /// Point spread functions convolved with cross terms of basis functions
                casa::Array<T> itsPSFCrossTerms;

                casa::Bool itsDecouple;

                casa::String itsDecouplingAlgorithm;

                /// Basis function used in the deconvolution
                boost::shared_ptr<BasisFunction<T> > itsBasisFunction;

                /// We keep track of the strength and location of each component
                /// identified. This allows calculation of the L1 norm of the
                /// model. We use clean to minimise the L1 so this is a good
                /// check to make. Ideally for stokes I, this should be equal to the flux.
                casa::Vector<casa::Array<T> > itsL1image;

                /// The flux subtracted on each scale
                casa::Vector<T> itsScaleFlux;

                /// The peak of the convolved PSF as a function of scale
                casa::Vector<T> itsPSFScales;
        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/DeconvolverBasisFunction.tcc>

#endif  // #ifndef I_DECONVOLVERBASISFUNCTION_H
