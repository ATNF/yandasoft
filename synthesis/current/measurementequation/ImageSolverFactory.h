/// @file
///
/// ImageSolverFactory: Factory class for image solvers
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
#ifndef ASKAP_SYNTHESIS_IMAGESOLVERFACTORY_H_
#define ASKAP_SYNTHESIS_IMAGESOLVERFACTORY_H_

#include <fitting/Solver.h>
#include <fitting/Params.h>

#include <Common/ParameterSet.h>
#include <measurementequation/ImageSolver.h>

#include <boost/shared_ptr.hpp>

namespace askap
{
  namespace synthesis
  {
    /// @brief Construct image solvers according to parameters
    /// @ingroup measurementequation
    class ImageSolverFactory
    {
      public:
        ImageSolverFactory();
        
        ~ImageSolverFactory();

        /// @brief Make a shared pointer for an image solver
        /// @param ip Params for the solver
        /// @param parset ParameterSet containing description of
        /// solver to be constructed
        /// @return shared pointer to the solver
        static askap::scimath::Solver::ShPtr make(const LOFAR::ParameterSet& parset); 
        
        /// @brief Helper method to configure preconditioners
        /// @details This method parses reconditioner.Names and sets up
        /// given preconditioners. Parameters necessary for this setup are also
        /// extracted from parset
        /// @param[in] parset parameter set to extract the input from
        /// @param[in] solver shared pointer to the solver to be configured
        /// @note This method is made public (temporary) because the same
        /// functionality is required to setup the restore solver. In future, we
        /// expect to set up the restore solver through the same factory
        static void configurePreconditioners(const LOFAR::ParameterSet &parset,
                     const boost::shared_ptr<ImageSolver> &solver);
        
      protected:
        /// @brief Helper method to configure minor cycle threshold(s)
        /// @details This method parses threshold.minorcycle parameter
        /// of the parset file. The parameter can be either a single string
        /// or a vector of two strings. A number without units is interpreted
        /// as a fractional stopping threshold (w.r.t the peak residual), 
        /// so does the number with the percentage sign. An absolute flux given
        /// in Jy or related units is interpreted as an absolute threshold.
        /// Either one or both of these thresholds can be given in the same
        /// time.
        /// @param[in] parset parameter set to extract the input from
        /// @param[in] solver shared pointer to the solver to be configured
        static void configureThresholds(const LOFAR::ParameterSet &parset,
                     const boost::shared_ptr<ImageSolver> &solver);
                             
    };

  }
}
#endif
