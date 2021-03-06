/// @file
/// 
/// @brief A structural type joining together GenericEquation and MultiChunkEquation
/// @details Because we deal here with double inheritance, we need to overload
/// explicitly predict and calcGenericEquation methods, otherwise the methods
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

#ifndef GENERIC_MULTI_CHUNK_EQUATION_H
#define GENERIC_MULTI_CHUNK_EQUATION_H

// own include
#include <askap/measurementequation/MultiChunkEquation.h>
#include <askap/scimath/fitting/GenericEquation.h>
#include <askap/scimath/fitting/GenericNormalEquations.h>


namespace askap {

namespace synthesis {


/// @brief A structural type joining together GenericEquation and MultiChunkEquation
/// @details Because we deal here with double inheritance, we need to overload
/// explicitly predict and calcGenericEquation methods, otherwise the methods
/// in the scimath::Equation tree are left pure abstract. Theoretically, we can
/// get rid of the double inheritance here and the need of this class would 
/// disappear. However, MultiChunkEquation is envisaged to be a temporary
/// class before we start to work with accessors only. Therefore, the inheritance
/// is left double to highlight the future interface 
/// (derived from IMeasurementEquation).
/// @ingroup measurementequation
struct GenericMultiChunkEquation : virtual public MultiChunkEquation,
                              virtual public askap::scimath::GenericEquation    
{  
  /// @brief Standard constructor, which remembers data iterator.
  /// @param[in] idi data iterator
  GenericMultiChunkEquation(const accessors::IDataSharedIter& idi);

  /// @brief Calculate the normal equations for the iterator
  /// @details This version iterates through all chunks of data and
  /// calls an abstract method declared in IMeasurementEquation for each 
  /// individual accessor (each iteration of the iterator). This method is
  /// overriden in this class to do a proper type conversion.
  /// @param[in] ne Normal equations
  virtual void calcGenericEquations(askap::scimath::GenericNormalEquations& ne) const;

  /// @brief Calculate the normal equation for one accessor (chunk).
  /// @details This calculation is done for a single chunk of
  /// data only (one iteration).It seems that all measurement
  /// equations should work with accessors rather than iterators
  /// (i.e. the iteration over chunks should be moved to the higher
  /// level, outside this class). This method overrides an abstract method
  /// of MultiChunkEquation. It calls calcGenericEquation  with ne converted
  /// to GenericNormalEquations
  /// @param[in] chunk a read-write accessor to work with
  /// @param[in] ne Normal equations
  virtual void calcEquations(const accessors::IConstDataAccessor &chunk,
                          askap::scimath::INormalEquations& ne) const;
  
  /// @brief Calculate the normal equation for one accessor (chunk).
  /// @details This calculation is done for a single chunk of
  /// data only (one iteration).It seems that all measurement
  /// equations should work with accessors rather than iterators
  /// (i.e. the iteration over chunks should be moved to the higher
  /// level, outside this class). 
  /// This class translates calls to calcEquations(chunk, ne) to 
  /// this method and does type conversion.
  /// @param[in] chunk a read-write accessor to work with
  /// @param[in] ne Normal equations
  virtual void calcGenericEquations(const accessors::IConstDataAccessor &chunk,
                          askap::scimath::GenericNormalEquations& ne) const = 0;

  /// @brief Predict model visibility for the iterator.
  /// @details This version of the predict method iterates
  /// over all chunks of data and calls an abstract method declared
  /// in IMeasurementEquation for each accessor. 
  virtual void predict() const;  
 
  /// @brief Calculate the normal equations for the iterator
  /// @details This version iterates through all chunks of data and
  /// calls an abstract method declared in IMeasurementEquation for each 
  /// individual accessor (each iteration of the iterator)
  /// @param[in] ne Normal equations
  virtual void calcEquations(askap::scimath::INormalEquations& ne) const;
  
  using MultiChunkEquation::predict;
};


} // namespace synthesis

} // namespace askap

#endif // #ifndef GENERIC_MULTI_CHUNK_EQUATION_H
