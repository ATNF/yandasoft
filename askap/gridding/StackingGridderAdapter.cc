/// @file StackingGridderAdapter.cpp
/// @brief An adapter that buffers the input data.
/// @details Essentially when accessing data, instead of gridding every time
/// step, the input is buffered so it can be transposed.
//
///  @author Stephen Ord (c) (CSIRO)  2021
//

#include <stdio.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>


#include "StackingGridderAdapter.h"

ASKAP_LOGGER(logger, ".gridding.stackinggridderadapter");

using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

/// @brief initialise the adapter
/// @details
/// @param[in] gridder a shared pointer to the gridder to be wrapped by this adapter
/// @param[in] nsteps how many time-steps to buffer
StackingGridderAdapter::StackingGridderAdapter(const boost::shared_ptr<IVisGridder> &gridder) {
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief copy constructor
/// @details We need this because the gridder doing actual work is held by a shared pointer,
/// which is a non-trivial type
/// @param[in] other an object to copy from
StackingGridderAdapter::StackingGridderAdapter(const StackingGridderAdapter &other) {
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief clone a copy of this gridder
/// @return shared pointer to the clone
boost::shared_ptr<IVisGridder> StackingGridderAdapter::clone() {
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief initialise the gridding
/// @details
/// @param[in] axes axes specifications
/// @param[in] shape Shape of output image: cube: u,v,pol,chan
/// @param[in] dopsf Make the psf?

void StackingGridderAdapter::initialiseGrid(const scimath::Axes& axes,
             const casacore::IPosition& shape, const bool dopsf,
                            const bool)
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief grid the visibility data.
/// @param[in] acc const data accessor to work with
void StackingGridderAdapter::grid(accessors::IConstDataAccessor& acc)
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief form the final output image
/// @param[in] out output double precision image or PSF
void StackingGridderAdapter::finaliseGrid(casacore::Array<imtype>& out)
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief finalise weights
/// @details Form the sum of the convolution function squared, multiplied by the weights for each
/// different convolution function. This is used in the evaluation of the second derivative.
/// @param[in] out output double precision sum of weights images
void StackingGridderAdapter::finaliseWeights(casacore::Array<imtype>& out)
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief initialise the degridding
/// @param[in] axes axes specifications
/// @param[in] image input image cube: u,v,pol,chan
void StackingGridderAdapter::initialiseDegrid(const scimath::Axes& axes,
       const casacore::Array<imtype>& image)
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief make context-dependant changes to the gridder behaviour
/// @param[in] context context description
void StackingGridderAdapter::customiseForContext(const std::string &context)
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief set visibility weights
/// @param[in] viswt shared pointer to visibility weights
void StackingGridderAdapter::initVisWeights(const IVisWeights::ShPtr &viswt)
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief degrid the visibility data.
/// @param[in] acc non-const data accessor to work with
void StackingGridderAdapter::degrid(accessors::IDataAccessor& acc)
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief finalise degridding
void StackingGridderAdapter::finaliseDegrid()
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @brief check whether the model is empty
/// @details A simple check allows us to bypass heavy calculations if the input model
/// is empty (all pixels are zero). This makes sense for degridding only.
/// @brief true, if the model is empty
bool StackingGridderAdapter::isModelEmpty() const
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}

/// @ brief check whether the adapter can stack
bool StackingGridderAdapter::isStackingGridder() const
{
  ASKAPTHROW(AskapError, "This StackingGridderAdapter method is not yet implemented");
}


