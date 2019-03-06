/// @file
/// 
/// @brief An adapter to make imaging equation a derivative from IMeasurementEquation
/// @details The current imaging code works with iterators, rather than accessors.
/// Although ImagingMultiChunkEquation allows to take this iterator dependency
/// in stages, it is still a lot of work to convert calcEquations and predict
/// methods of a typical imaging measurement equation to be able to derive it 
/// from this class. This adapter allows to translate calls to the virtual
/// methods of IMeasurementEquation to the appropriate call of the 
/// iterator-dependent measurement equation. I hope this adapter is temporary.
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

#include <measurementequation/ImagingEquationAdapter.h>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

#include <dataaccess/FakeSingleStepIterator.h>
#include <measurementequation/MultiChunkEquation.h>
#include <measurementequation/ImageFFTEquation.h>
#include <askap/AskapError.h>

/// @brief constructor
/// @details This constructor initializes a fake iterator. Actual measurement
/// equation is set up via a call to assign method (templated)
ImagingEquationAdapter::ImagingEquationAdapter() : 
        itsIterAdapter(new FakeSingleStepIterator) {}

/// @brief access to parameters
/// @details This call is translated to itsActualEquation. We need to override
/// this method of scimath::Equation, as otherwise an empty parameter class
/// initialized in the default constructor would always be returned.
/// @return const reference to paramters of this equation
const scimath::Params& ImagingEquationAdapter::parameters() const
{
  ASKAPCHECK(itsActualEquation, 
     "assign method should be called before first usage of ImagingEquationAdapter");
  return itsActualEquation->parameters();
}
   
/// @brief set parameters
/// @details This call is translated to itsActualEquation.
/// @param[in] par new parameters
void ImagingEquationAdapter::setParameters(const scimath::Params &par)
{
  ASKAPCHECK(itsActualEquation, 
     "assign method should be called before first usage of ImagingEquationAdapter");
  itsActualEquation->setParameters(par);  
}
   
/// @brief predict visibilities
/// @details This call is translated to itsActualEquation.
void ImagingEquationAdapter::predict() const
{
  ASKAPCHECK(itsActualEquation, 
     "assign method should be called before first usage of ImagingEquationAdapter");
  // there will be an exception if this class is initialized with the type,
  // which works with the iterator directly and bypasses accessor-based method.
  itsActualEquation->predict();
}
   
/// @brief calculate normal equations
/// @details This call is translated to itsActualEquation.
/// @param[in] ne normal equations to be updated
void ImagingEquationAdapter::calcEquations(scimath::INormalEquations &ne) const
{ 
  ASKAPCHECK(itsActualEquation, 
     "assign method should be called before first usage of ImagingEquationAdapter");
  
  // there will be an exception if this class is initialized with the type,
  // which works with the iterator directly and bypasses accessor-based method.
  itsActualEquation->calcEquations(ne);
}
   
/// @brief clone this "composite" equation
/// @details The operations performed by this method are more complex 
/// than just copy constructor, because we store shared pointers to 
/// the iterator adapter and underlying measurement equation. They both
/// have to be cloned properly.
/// @return a shared pointer with the cloned version
scimath::Equation::ShPtr ImagingEquationAdapter::clone() const
{
  boost::shared_ptr<ImagingEquationAdapter> result(new ImagingEquationAdapter(*this));
  try {
     if (itsIterAdapter) {
         // we need an explicit type conversion, because shared iter would
         // dereference the iterator according to its interface
         const boost::shared_ptr<IDataIterator> basicIt = itsIterAdapter;
         ASKAPDEBUGASSERT(basicIt);
     
         const FakeSingleStepIterator &it = dynamic_cast<const FakeSingleStepIterator&>(*basicIt);
         result->itsIterAdapter = IDataSharedIter(new FakeSingleStepIterator(it));
     }
     if (itsActualEquation) {
         result->itsActualEquation = itsActualEquation->clone();
     }
  }
  catch (const std::bad_cast &bc) {
     ASKAPTHROW(AskapError, "Bad cast inside ImagingEquationAdapter::clone, most likely this means "
                 "there is a logical error");
  }
  return result;
}
   
/// @brief accessor-based version of predict
/// @details This version of predict is implemented via iterator-based
/// version of itsIterAdapter.
/// @param[in] chunk a chunk to be filled with predicted data
void ImagingEquationAdapter::predict(IDataAccessor &chunk) const
{     
   boost::shared_ptr<FakeSingleStepIterator> it = 
             itsIterAdapter.dynamicCast<FakeSingleStepIterator>();
   if (!it) {
       ASKAPTHROW(AskapError, "Bad cast inside ImagingEquationAdapter::predict, most likely this means "
              "there is a logical error"); 
   }
   it->assignDataAccessor(chunk);  
   predict();
   it->detachAccessor();      
}
   
/// @brief accessor-based version of calcEquations
/// @details This version of calcEquations is implemented via iterator-based
/// version of itsIterAdapter.
/// @param[in] chunk a chunk of data to work with
/// @param[in] ne normal equations to update
void ImagingEquationAdapter::calcEquations(const IConstDataAccessor &chunk,
             scimath::INormalEquations &ne) const
{
  boost::shared_ptr<FakeSingleStepIterator> it = 
            itsIterAdapter.dynamicCast<FakeSingleStepIterator>();
  if (!it) {
      ASKAPTHROW(AskapError, "Bad cast inside ImagingEquationAdapter::calcEquations, most likely this means "
                 "there is a logical error");
  }
  it->assignConstDataAccessor(chunk);  
  calcEquations(ne);
  it->detachAccessor();      
}       
 
 
/// @brief assign the actual measurement equation to an adapter
/// @details This method assigns a measurement set equation, which has 
/// already been constructed. Templated methods construct a new object.
/// @param[in] me a shared pointer to a measurement set
void ImagingEquationAdapter::assign(const scimath::Equation::ShPtr &me)
{ 
  ASKAPDEBUGASSERT(me);
  itsActualEquation = me; 
  boost::shared_ptr<MultiChunkEquation> multiChunkME = 
                 boost::dynamic_pointer_cast<MultiChunkEquation>(me);
  if (multiChunkME) {
      // multi-chunk MEs hold own copies of the iterator, which has to be
      // substituted to a fake iterator here
      multiChunkME->setIterator(itsIterAdapter);
  } else {
     // a bit ugly solution, but it should go away when we stop using
     // iterator-based measurement equations
     boost::shared_ptr<ImageFFTEquation> imageFFTME = 
                 boost::dynamic_pointer_cast<ImageFFTEquation>(me);
     // an instance of ImageFFTEquation holds own copy of the iterator, which
     // has to be substituted to a fake iterator here
     imageFFTME->setIterator(itsIterAdapter);        
  }
}
 
