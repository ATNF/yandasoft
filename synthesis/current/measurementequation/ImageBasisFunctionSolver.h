/// @file
///
/// ImageBasisFunctionSolver: This solver calculates the dirty image (or equivalent)
/// for all parameters called image*
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
#ifndef SYNIMAGEBASISFUNCTIONSOLVER_H_
#define SYNIMAGEBASISFUNCTIONSOLVER_H_

#include <measurementequation/ImageCleaningSolver.h>

#include <lattices/Lattices/ArrayLattice.h>
#include <deconvolution/DeconvolverBasisFunction.h>

namespace askap {
    namespace synthesis {
        /// @brief BasisFunction solver for images.
        ///
        /// @ingroup measurementequation
        class ImageBasisFunctionSolver : public ImageCleaningSolver {
            public:

                /// @brief default constructor
                ImageBasisFunctionSolver();

                /// @brief default constructor
     	        ImageBasisFunctionSolver(casa::Vector<float>& scales);

                /// @brief Initialize this solver
                virtual void init();
 
                /// @brief Solve for parameters, updating the values kept internally
                /// The solution is constructed from the normal equations. The parameters named 
                /// image* are interpreted as images and solved for.
                /// @param[in] ip current model (to be updated)
                /// @param[in] quality Solution quality information
                virtual bool solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality);

                /// @brief Clone this object
                /// @return a shared pointer to the clone
                virtual askap::scimath::Solver::ShPtr clone() const;

	  /// @brief configure basic parameters of the solver
	  /// @details This method encapsulates extraction of basic solver parameters from the parset.
	  /// @param[in] parset parset's subset (should have solver.Clean or solver.Dirty removed)
	  virtual void configure(const LOFAR::ParameterSet &parset); 

	  virtual void setBasisFunction(BasisFunction<Float>::ShPtr bf);

	  BasisFunction<Float>::ShPtr basisFunction();

            protected:
                /// @brief Precondition the PSF and the dirty image
                /// @param[in] psf point spread function to precondition (in/out)
                /// @param[in] dirty dirty image to precondition (in/out)
                void preconditionNE(casa::ArrayLattice<float>& psf, casa::ArrayLattice<float>& dirty);

	  boost::shared_ptr<DeconvolverControl<Float> > itsControl;
      
	  boost::shared_ptr<DeconvolverMonitor<Float> > itsMonitor;

	  Bool itsUseCrossTerms;

	  BasisFunction<Float>::ShPtr itsBasisFunction;

            private:

        };

    }
}
#endif
