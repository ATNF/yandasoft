/// @file
/// 
/// @brief An abstract measurement equation.
/// @details To be able to use common code regardless on the time of the
/// measurement equaiton used (i.e. ComponentEquation, ImageFFTEquation, etc)
/// we need a common ancestor of the measurement equation classes.
/// askap::scimath::Equation is not specialised enough for this purpose.
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>


#include <measurementequation/IMeasurementEquation.h>

using namespace askap;
using namespace askap::synthesis;

/// @brief empty virtual descrtuctor to make the compiler happy
IMeasurementEquation::~IMeasurementEquation() {}

