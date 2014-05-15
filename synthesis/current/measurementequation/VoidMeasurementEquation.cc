/// @file
/// 
/// @brief A measurement equation, which does nothing.
/// @details The current calibration class requires a perfect measurement
/// equation. This class has been written to be able to use the same code 
/// for both applying a calibration and solving for parameters. It is 
/// a void measurement equation in the sense that it does nothing to the 
/// data or normal equations given to it.
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


#include <measurementequation/VoidMeasurementEquation.h>

using namespace askap;
using namespace askap::synthesis;

/// @brief Predict model visibilities for one accessor (chunk).
/// @details This prediction is done for single chunk of data only. 
/// It seems that all measurement equations should work with accessors 
/// rather than iterators (i.e. the iteration over chunks should be 
/// moved to the higher level, outside this class). 
void VoidMeasurementEquation::predict(accessors::IDataAccessor &) const {}

/// @brief Calculate the normal equation for one accessor (chunk).
/// @details This calculation is done for a single chunk of
/// data only (one iteration).It seems that all measurement
/// equations should work with accessors rather than iterators
/// (i.e. the iteration over chunks should be moved to the higher
/// level, outside this class). 
void VoidMeasurementEquation::calcEquations(const accessors::IConstDataAccessor &,
                          askap::scimath::INormalEquations&) const {}
