/// @file DeconvolverMultiTermBasisFunction.h
/// @brief Class for a BasisFunction-Clean-based deconvolver
/// @details This interface class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
/// This version can deal with multiple terms
/// D = B(0)*I(0) + B(1)*I(1) + B(2)*I(2)
/// The most common example is where the B's are the spectral dirty beams
/// and the I(0) are the Taylor series approximation to the frequency
/// dependent frequencies.
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

#ifndef ASKAP_SYNTHESIS_DECONVOLVERMULTITERMBASISFUNCTION_H
#define ASKAP_SYNTHESIS_DECONVOLVERMULTITERMBASISFUNCTION_H

#include <string>

#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
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
        class DeconvolverMultiTermBasisFunction : public DeconvolverBase<T, FT> {

            public:

                typedef boost::shared_ptr<DeconvolverMultiTermBasisFunction<T, FT> > ShPtr;


                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We use a vector
                /// here to allow multiple dirty images and PSFs for the
                /// same model (e.g. as in MFS)
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                /// @param[in] psf Point Spread Function containing 2*nTaylor-1 terms (array)
                DeconvolverMultiTermBasisFunction(casa::Vector<casa::Array<T> >& dirty,
                                                  casa::Vector<casa::Array<T> >& psf,
                                                  casa::Vector<casa::Array<T> >& psfLong);

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We keep this
                /// version for compatibility
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverMultiTermBasisFunction(casa::Array<T>& dirty, casa::Array<T>& psf);

                virtual ~DeconvolverMultiTermBasisFunction();

                /// @brief Set the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction.
                /// @param[in] bf Shared pointer to basisfunction instance
                void setBasisFunction(boost::shared_ptr<BasisFunction<T> > bf);

                /// @brief Return the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction
                boost::shared_ptr<BasisFunction<T> > basisFunction();

                /// @brief Set the type of solution used in finding the optimum component
                void setSolutionType(casa::String solutionType);

                /// @brief Get the type of solution used in finding the optimum component
                const casa::String solutionType();

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

                /// @brief Update only the dirty image
                /// @detail Update an existing deconvolver for a changed dirty image
                /// @param[in] dirty Dirty image (array)
                /// @param[in] term term to update
                virtual void updateDirty(casa::Array<T>& dirty, casa::uInt term = 0);

                /// @brief Update only the dirty images
                /// @detail Update an existing deconvolver for a changed dirty images.
                /// @param[in] dirty Dirty image (vector of arrays)
                virtual void updateDirty(casa::Vector<casa::Array<T> >& dirty);

            private:

                // Perform one iteration
                bool oneIteration();

                // Initialise the PSFs - only need to do this once per change in basis functions
                void initialisePSF();

                void initialiseResidual();

                // Initialise the PSFs - only need to do this once per change in basis functions
                virtual void initialiseForBasisFunction(bool force);

                void chooseComponent(uInt& optimumBase, casa::IPosition& absPeakPos, T& absPeakVal, Vector<T>& peakValues);

                // Long vector of PSFs
                casa::Vector<casa::Array<T> > itsPsfLongVec;

                /// Residual images convolved with basis functions, [nx,ny][nterms][nbases]
                casa::Vector<casa::Vector<casa::Array<T> > > itsResidualBasis;

                /// Point spread functions convolved with cross terms
                // [nx,ny][nterms,nterms][nbases,nbases]
                casa::Matrix<casa::Matrix<casa::Array<T> > > itsPSFCrossTerms;

                /// The coupling between different terms for each basis [nterms,nterms][nbases]
                casa::Vector<casa::Matrix<casa::Double> > itsCouplingMatrix;

                /// Inverse of the coupling matrix [nterms,nterms][nbases]
                casa::Vector<casa::Matrix<casa::Double> > itsInverseCouplingMatrix;

                /// Determinants of the coupling Matrix for each base [nbases]
                casa::Vector<casa::Double> itsDetCouplingMatrix;

                /// Basis function used in the deconvolution
                boost::shared_ptr<BasisFunction<T> > itsBasisFunction;

                /// The flux subtracted on each term and scale [nterms][nbases]
                casa::Vector< casa::Vector<T> > itsTermBaseFlux;

                casa::Bool itsDirtyChanged;

                casa::Bool itsBasisFunctionChanged;

                casa::String itsSolutionType;
        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/DeconvolverMultiTermBasisFunction.tcc>

#endif
