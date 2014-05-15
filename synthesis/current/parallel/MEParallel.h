/// @file
///
/// Provides generic methods for parallel algorithms using the measurement equation
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
#ifndef ASKAP_SYNTHESIS_MEPARALLEL_H_
#define ASKAP_SYNTHESIS_MEPARALLEL_H_

// System includes
#include <string>

// ASKAPsoft includes
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <fitting/INormalEquations.h>
#include <fitting/Equation.h>
#include <fitting/Solver.h>

// Loacl package includes
#include <parallel/SynParallel.h>

namespace askap
{
	namespace synthesis
	{
		/// @brief Support for parallel algorithms using the measurement equation
		///
		/// @details Support for parallel applications using the measurement equation
		/// classes. An application is derived from this abstract base. The model used is that the
		/// application has many prediffers and one solver, running in separate MPI processes
		/// or in one single thread. The solver is the master so the number of processes
		/// is one more than the number of prediffers. Each prediffer is currently given
		/// a separate data set.
		///
		/// The steps are:
		/// (a) define an initial model and distribute to all prediffers
		/// (b) calculate the normal equations for each data set (this part is
		/// distributed across the prediffers)
		/// (c) send all normal equations to the solver for merging
		/// (d) solve the merged normal equations
		/// (e) distribute the model to all prediffers and return to (b)
		///
		/// The caller is responsible for ensuring that the model is transferred correctly 
		/// before a CalcNE and after a SolveNE. For example:
		/// @code
		///		for (int cycle=0;cycle<nCycles;cycle++)
		///		{
		///			imager.os() << "*** Starting major cycle " << cycle << " ***" << std::endl;
		///			if(cycle>0) imager.receiveModel();
		///			imager.calcNE();
		///			imager.solveNE();
		///			// Broadcast the model
		///			if (cycle<(nCycles-1)) imager.broadcastModel();
		///		}
		/// @endcode
		/// The normal equations are transferred automatically between the calcNE and solveNE
		/// steps so the called does not need to be concerned about that.
		///
		/// If the number of nodes is 1 then everything occurs in the same process with
		/// no overall for transmission of model or normal equations.
		///
		/// @ingroup parallel
		class MEParallel : public SynParallel
		{
			public:

				/// @brief Constructor 
				/// @details The first parameter is needed solely for MPI, the second
                /// is the parset to be used in derived classes
				/// @param[in] comms communication object
				/// @param[in] parset parameter set      				
				MEParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset);

				virtual ~MEParallel();

				/// Calculate the normalequations (runs only in the workers)
				virtual void calcNE() = 0;

				/// Solve the normal equations (runs only in the master)
				virtual void solveNE() = 0;

				/// Write the model (runs only in the master)
				/// @param[in] postfix a string to be added to the file name
				virtual void writeModel(const std::string &postfix = std::string()) = 0;

				/// @brief Send the normal equations from this worker to the master
				void sendNE();

                /// @brief Receive the normal equations from all workers into this master
                void receiveNE();

                /// @brief Perform a reduction for normal equations from all
                /// workers to the master.
                void reduceNE(askap::scimath::INormalEquations::ShPtr ne);

			protected:
		
                // Point-to-point send normal equations
                // @param[in] ne    pointer to normal equations to send
                // @param[in] dest  rank of process to send normal equations to    
                void sendNormalEquations(const askap::scimath::INormalEquations::ShPtr ne, int dest);

                // Point-to-point receive normal equations
                // @param[in] source    rank of the process from which normal
                // equations will be received
                // @return a shared pointer, pointing to the received normal equations
                askap::scimath::INormalEquations::ShPtr receiveNormalEquations(int source);

				/// Holder for the normal equations
				askap::scimath::INormalEquations::ShPtr itsNe;

				/// Holder for the solver
				askap::scimath::Solver::ShPtr itsSolver;
				
				/// Holder for the equation
				askap::scimath::Equation::ShPtr itsEquation;
		};

	}
}
#endif
