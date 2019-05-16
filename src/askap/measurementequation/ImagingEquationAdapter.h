/// @file
/// 
/// @brief An adapter to make imaging equation a derivative from IMeasurementEquation
/// @details The current imaging code works with iterators, rather than accessors.
/// Although ImagingMultiChunkEquation allows to take this iterator dependency out
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

#ifndef IMAGING_EQUATION_ADAPTER_H
#define IMAGING_EQUATION_ADAPTER_H

#include <fitting/Equation.h>
#include <fitting/Params.h>
#include <fitting/INormalEquations.h>
#include <measurementequation/IMeasurementEquation.h>
#include <gridding/IVisGridder.h>
#include <dataaccess/SharedIter.h>
#include <dataaccess/IDataIterator.h>
#include <dataaccess/IDataAccessor.h>

namespace askap {

namespace synthesis {

/// @brief An adapter to make imaging equation a derivative from IMeasurementEquation
/// @details The current imaging code works with iterators, rather than accessors.
/// Although ImagingMultiChunkEquation allows to take this iterator dependency out
/// in stages, it is still a lot of work to convert calcEquations and predict
/// methods of a typical imaging measurement equation to be able to derive it 
/// from this class. This adapter allows to translate calls to the virtual
/// methods of IMeasurementEquation to the appropriate call of the 
/// iterator-dependent measurement equation. I hope this adapter is temporary.
/// @ingroup measurementequation
struct ImagingEquationAdapter : virtual public IMeasurementEquation,
                      virtual public askap::scimath::Equation
{
   /// @brief constructor
   /// @details This constructor initializes a fake iterator. Actual measurement
   /// equation is set up via a call to assign method (templated)
   ImagingEquationAdapter();
   
   /// @brief assign actual measurement equation to an adapter
   /// @details This templated method constructs the actual measurement equation
   /// of an appropriate type and sets up itsActualEquation. itsIterAdapter is
   /// passed as an iterator. This version accepts just the parameters.
   /// @param[in] par input parameters
   template<typename ME>
   void assign(const scimath::Params &par) 
   { itsActualEquation.reset(new ME(par, itsIterAdapter)); }
   
   /// @brief assign actual measurement equation to an adapter
   /// @details This templated method constructs the actual measurement equation
   /// of an appropriate type and sets up itsActualEquation. itsIterAdapter is
   /// passed as an iterator. This version accepts parameters and a gridder.
   /// @param[in] par input parameters
   /// @param[in] gridder input gridder (passed as shared pointer)
   template<typename ME>
   void assign(const scimath::Params &par, const IVisGridder::ShPtr &gridder) 
   { itsActualEquation.reset(new ME(par, itsIterAdapter, gridder)); }
   
   /// @brief assign the actual measurement equation to an adapter
   /// @details This method assigns a measurement set equation, which has 
   /// already been constructed. Templated methods construct a new object.
   /// @param[in] me a shared pointer to a measurement set
   void assign(const scimath::Equation::ShPtr &me);
   
   /// @brief access to parameters
   /// @details This call is translated to itsActualEquation. We need to override
   /// this method of scimath::Equation, as otherwise an empty parameter class
   /// initialized in the default constructor would always be returned.
   /// @return const reference to paramters of this equation
   virtual const scimath::Params& parameters() const;
   
   /// @brief set parameters
   /// @details This call is translated to itsActualEquation.
   /// @param[in] par new parameters
   virtual void setParameters(const scimath::Params &par);
   
   /// @brief predict visibilities
   /// @details This call is translated to itsActualEquation.
   virtual void predict() const;
   
   /// @brief calculate normal equations
   /// @details This call is translated to itsActualEquation.
   /// @param[in] ne normal equations to be updated
   virtual void calcEquations(scimath::INormalEquations &ne) const;
   
   /// @brief clone this "composite" equation
   /// @details The operations performed by this method are more complex 
   /// than just copy constructor, because we store shared pointers to 
   /// the iterator adapter and underlying measurement equation. They both
   /// have to be cloned properly.
   /// @return a shared pointer with the cloned version
   scimath::Equation::ShPtr clone() const;        
   
   /// @brief accessor-based version of predict
   /// @details This version of predict is implemented via iterator-based
   /// version of itsIterAdapter.
   /// @param[in] chunk a chunk to be filled with predicted data
   virtual void predict(accessors::IDataAccessor &chunk) const;
   
   /// @brief accessor-based version of calcEquations
   /// @details This version of calcEquations is implemented via iterator-based
   /// version of itsIterAdapter.
   /// @param[in] chunk a chunk of data to work with
   /// @param[in] ne normal equations to update
   virtual void calcEquations(const accessors::IConstDataAccessor &chunk,
             scimath::INormalEquations &ne) const;
   
private:
   /// @brief iterator adapter
   accessors::IDataSharedIter itsIterAdapter;
   /// @brief actual measurement equation 
   scimath::Equation::ShPtr itsActualEquation;
};


} // namespace synthesis

} // namespace askap

#endif // #ifndef IMAGING_EQUATION_ADAPTER_H

