/// @file
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

#include <casa/aips.h>
#include <casa/Quanta.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.meparsetinterface");

#include <askap/AskapError.h>

#include <measurementequation/MEParsetInterface.h>
#include <measurementequation/ImageSolver.h>
#include <measurementequation/ImageMultiScaleSolver.h>

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;
using namespace askap;

namespace askap
{
	namespace synthesis
	{

		void operator<<(askap::scimath::Solver::ShPtr& solver, 
				const LOFAR::ParameterSet &parset)
		{
			askap::scimath::Params params;
			if(parset.getString("solver")=="Clean")
			{
				std::vector<float> defaultScales(3);
				defaultScales[0]=0.0;
				defaultScales[1]=10.0;
				defaultScales[2]=30.0;
				std::vector<float> scales=parset.getFloatVector("solver.Clean.scales", defaultScales);
				solver = askap::scimath::Solver::ShPtr(new ImageMultiScaleSolver(casa::Vector<float>(scales)));
				ASKAPLOG_INFO_STR(logger, "Constructed image multiscale solver" );
				solver->setGain(parset.getFloat("solver.Clean.gain", 0.7));
				solver->setAlgorithm(parset.getString("solver.Clean.algorithm", "MultiScale"));
				solver->setVerbose(parset.getBool("solver.Clean.verbose", true));
				solver->setNiter(parset.getInt32("solver.Clean.niter", 100));
				casa::Quantity threshold;
				casa::Quantity::read(threshold, parset.getString("solver.Clean.threshold", "0Jy"));
				solver->setThreshold(threshold);
			}
			else
			{
				solver = askap::scimath::Solver::ShPtr(new ImageSolver);
				casa::Quantity threshold;
				casa::Quantity::read(threshold, parset.getString("solver.Dirty.threshold", "0Jy"));
				solver->setThreshold(threshold);
				ASKAPLOG_INFO_STR(logger, "Constructed dirty image solver" );
			}
		}
	}
}
