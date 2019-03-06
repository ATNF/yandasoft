/// @file
/// 
/// @brief A structural type joining together ImagingEquation and MultiChunkEquation
/// @details Because we deal with the double inheritance here, we need to overload
/// explicitly predict and calcImagingEquation methods, otherwise the methods
/// in the scimath::Equation tree are left pure abstract. Theoretically, we can
/// get rid of the double inheritance here and the need of this class would 
/// disappear. However, MultiChunkEquation is envisaged to be a temporary
/// class before we start to work with accessors only. Therefore, the inheritance
/// is left double to highlight the future interface 
/// (derived from IMeasurementEquation).
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


// own includes
#include <measurementequation/ImagingMultiChunkEquation.h>
#include <askap/AskapError.h>
#include <measurementequation/NormalEquationsTypeError.h>

// std includes
#include <stdexcept>

namespace askap {

namespace synthesis {

/// @brief Standard constructor, which remembers data iterator.
/// @param[in] idi data iterator
ImagingMultiChunkEquation::ImagingMultiChunkEquation(const accessors::IDataSharedIter& idi) :
         MultiChunkEquation(idi) {}


/// @brief Calculate the normal equations for the iterator
/// @details This version iterates through all chunks of data and
/// calls an abstract method declared in IMeasurementEquation for each 
/// individual accessor (each iteration of the iterator). 
/// @param[in] ne Normal equations
void ImagingMultiChunkEquation::calcImagingEquations(askap::scimath::ImagingNormalEquations& ne) const
{
  MultiChunkEquation::calcEquations(ne);
}

/// @brief Calculate the normal equation for one accessor (chunk).
/// @details This calculation is done for a single chunk of
/// data only (one iteration).It seems that all measurement
/// equations should work with accessors rather than iterators
/// (i.e. the iteration over chunks should be moved to the higher
/// level, outside this class). This method overrides an abstract method
/// of MultiChunkEquation. It calls calcImagingEquations(chunk,ne)  with 
/// ne converted to ImagingNormalEquations
/// @param[in] chunk a read-write accessor to work with
/// @param[in] ne Normal equations
void ImagingMultiChunkEquation::calcEquations(const accessors::IConstDataAccessor &chunk,
                   askap::scimath::INormalEquations& ne) const
{
  try {
     calcImagingEquations(chunk, 
             dynamic_cast<askap::scimath::ImagingNormalEquations&>(ne));
  }
  catch (const std::bad_cast &bc) {
     ASKAPTHROW(NormalEquationsTypeError, "An attempt to use incompatible type of "
                 "the normal equations class with a derivative from "
                 "ImagingMultiChunkEquation. It accepts only ImagingNormalEquations "
                 "and derivatives. This exception probably indicates a logic error");    
  }
}

/// @brief Predict model visibility for the iterator.
/// @details This version of the predict method iterates
/// over all chunks of data and calls an abstract method declared
/// in IMeasurementEquation for each accessor. 
void ImagingMultiChunkEquation::predict() const      
{
  MultiChunkEquation::predict();
}

/// @brief Calculate the normal equations for the iterator
/// @details This version iterates through all chunks of data and
/// calls an abstract method declared in IMeasurementEquation for each 
/// individual accessor (each iteration of the iterator)
/// @param[in] ne Normal equations
void ImagingMultiChunkEquation::calcEquations(askap::scimath::INormalEquations& ne) const
{
  MultiChunkEquation::calcEquations(ne);
}
 

} // namespace synthesis

} // namespace askap
