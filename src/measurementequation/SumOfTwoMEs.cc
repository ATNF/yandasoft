/// @file
/// 
/// @brief A measurement equation, which is a sum of two measurement
/// equations.
/// @details For simulation it is necessary to be able to add noise 
/// to the simulated visibilities. One way of doing this is to write
/// a special measurement equation which predict noise and use a 
/// composite equation when a prediction must be made. Such an equation
/// can't be solved with a regular solver (due to a stochastic nature
/// of the problem statistical estimators are needed), but prediction
/// would work. Another application of this class is a composite imaging
/// equation where the model is composed from an image and a list of 
/// components. If there are many other additive effects to be implemented
/// and/or solution for parameters is required, the measurement equation
/// corresponding to the random visibility noise generator can be reorganized
/// into a template + individual effects in a similar way to that how 
/// CalibrationME template is written. 
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


#include <measurementequation/SumOfTwoMEs.h>
#include <dataaccess/MemBufferDataAccessor.h>
#include <measurementequation/NormalEquationsTypeError.h>
#include <askap/AskapError.h>

#include <stdexcept>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

/// @brief Constructor   
/// @details Creates a new composite measurement equation equivalent
/// to a sum of the given equations. Equations passed as parameters 
/// are not changed.
/// @param[in] first a shared pointer to the first equation
/// @param[in] second a shared pointer to the second equation
SumOfTwoMEs::SumOfTwoMEs(
              const boost::shared_ptr<IMeasurementEquation const> &first,
              const boost::shared_ptr<IMeasurementEquation const> &second,
              const IDataSharedIter &it) : MultiChunkEquation(it),
                itsFirstME(first), itsSecondME(second) {}

/// @brief Predict model visibilities for one accessor (chunk).
/// @details This prediction is done for single chunk of data only. 
/// It seems that all measurement equations should work with accessors 
/// rather than iterators (i.e. the iteration over chunks should be 
/// moved to the higher level, outside this class).
/// @param[in] chunk a read-write accessor to work with 
void SumOfTwoMEs::predict(IDataAccessor &chunk) const 
{
  MemBufferDataAccessor secondResult(chunk);
  ASKAPDEBUGASSERT(itsFirstME);
  ASKAPDEBUGASSERT(itsSecondME);
  itsSecondME->predict(secondResult);
  itsFirstME->predict(chunk);
  casa::Cube<casa::Complex> &rwVis = chunk.rwVisibility();
  rwVis += secondResult.visibility();
}

/// @brief Calculate the normal equation for one accessor (chunk).
/// @details This calculation is done for a single chunk of
/// data only (one iteration).It seems that all measurement
/// equations should work with accessors rather than iterators
/// (i.e. the iteration over chunks should be moved to the higher
/// level, outside this class). 
/// @note This method will work correctly only if two parts of the
/// equation are completely independent. If there is a common 
/// parameter for both parts, normal equations on that parameter will be
/// wrong because the cross terms are omitted. This class is currently
/// seen to be used for simulations (where only predict method is used),
/// therefore it is not an issue. However, if a proper functionality is 
/// required, the only way to achieve it is to use a similar approach
/// to CalibrationME template and plug in effects.
/// @param[in] chunk a read-write accessor to work with
/// @param[in] ne Normal equations
void SumOfTwoMEs::calcEquations(const IConstDataAccessor &chunk,
                          askap::scimath::INormalEquations &ne) const 
{
  ASKAPDEBUGASSERT(itsFirstME);
  ASKAPDEBUGASSERT(itsSecondME);
  
  // each calcEquations assumes that the target visibilities are given in
  // the chunk supplied as the input parameter. However, in a composite
  // equation the chunk passed as input corresponds to the whole equation.
  // Therefore, it can not be simply passed to the individual parts.
  // The value produced by another part has to be subtracted from the 
  // chunk's visibility cube before calling calcEquations for the 
  // individual parts.
  MemBufferDataAccessor resultOfCalculatedPart(chunk);
  MemBufferDataAccessor resultOfOtherPart(chunk);
  itsSecondME->predict(resultOfOtherPart);
  resultOfCalculatedPart.rwVisibility() = 
                 chunk.visibility() - resultOfOtherPart.visibility();
  try {
     itsFirstME->calcEquations(resultOfCalculatedPart,ne);
  }
  // ignore the error, if type of the measurement equation is incompatible with
  // the given normal equations. This just means that there is no dependency
  // of this part of the composite equation on the parameters represented by
  // the normal equations.
  catch (const NormalEquationsTypeError &) {}
 
  itsFirstME->predict(resultOfOtherPart);
  resultOfCalculatedPart.rwVisibility() = 
                 chunk.visibility() - resultOfOtherPart.visibility();
  try {
     itsSecondME->calcEquations(resultOfCalculatedPart, ne);
  }
  // ignore the error, if type of the measurement equation is incompatible with
  // the given normal equations. This just means that there is no dependency
  // of this part of the composite equation on the parameters represented by
  // the normal equations.
  catch (const NormalEquationsTypeError &) {}
}


/// Clone this into a shared pointer
/// @return shared pointer to a copy
scimath::Equation::ShPtr SumOfTwoMEs::clone() const 
{
  return scimath::Equation::ShPtr(new SumOfTwoMEs(*this));
}

/// Predict the data from the parameters.
void SumOfTwoMEs::predict() const
{
  // probably we could get away with just "using" keyword
  MultiChunkEquation::predict();
}

/// Calculate the normal equations for the given data and parameters
/// @param ne Normal equations to be filled
void SumOfTwoMEs::calcEquations(scimath::INormalEquations& ne) const
{ 
  // probably we could get away with just "using" keyword
  MultiChunkEquation::calcEquations(ne);
}
