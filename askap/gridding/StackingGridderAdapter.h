/// @file
///
/// @brief Gridder adapter perform buffering of input data
/// @details it is becoming clear that a number of schemes to 
/// segment the input dataset will need to be prototyped and tested.
/// The current accessors are not set up to be too flexible in this regard
/// And a change to implement say an input TIME chunck of all the data 
/// may not be easy to asbstract. So this class is simply here to provide
/// a buffer of the whole input allocation - without changing the selectors
/// and accessors.
///
/// @copyright (c) 2021 CSIRO
/// @author Stephen Ord <stephen.ord@csiro.au>


#ifndef STACKING_GRIDDER_ADAPTER_H
#define STACKING_GRIDDER_ADAPTER_H

#include <askap/gridding/IVisGridder.h>
#include <askap/dataaccess/MemBufferDataAccessor.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/IDataIterator.h>

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief Gridder adapter to perform buffering of input data
/// @details Intended to allow further segmentation of input visset
/// perhaps by derived classes
/// @ingroup gridding
class StackingGridderAdapter : virtual public IVisGridder
{
public:
  
   typedef boost::shared_ptr<StackingGridderAdapter> ShPtr;
   
   /// @brief initialise the adapter
   /// @details
   /// @param[in] gridder a shared pointer to the gridder to be wrapped by this adapter
   /// @param[in] nsteps how many time-steps to buffer
   StackingGridderAdapter(const boost::shared_ptr<IVisGridder> &gridder);

   /// @brief copy constructor
   /// @details We need this because the gridder doing actual work is held by a shared pointer,
   /// which is a non-trivial type
   /// @param[in] other an object to copy from
   StackingGridderAdapter(const StackingGridderAdapter &other);

   /// @brief clone a copy of this gridder
   /// @return shared pointer to the clone
   virtual boost::shared_ptr<IVisGridder> clone();

   /// @brief initialise the gridding
   /// @details
   /// @param[in] axes axes specifications
   /// @param[in] shape Shape of output image: cube: u,v,pol,chan
   /// @param[in] dopsf Make the psf?
   virtual void initialiseGrid(const scimath::Axes& axes,
                const casacore::IPosition& shape, const bool dopsf = true,
                const bool dopcf=false);

   /// @brief grid the visibility data.
   /// @param[in] acc const data accessor to work with
   virtual void grid(accessors::IConstDataAccessor& acc);

   /// @brief form the final output image
   /// @param[in] out output double precision image or PSF
   virtual void finaliseGrid(casacore::Array<imtype>& out);

   /// @brief finalise weights
   /// @details Form the sum of the convolution function squared, multiplied by the weights for each
   /// different convolution function. This is used in the evaluation of the second derivative.
   /// @param[in] out output double precision sum of weights images
   virtual void finaliseWeights(casacore::Array<imtype>& out);

   /// @brief initialise the degridding
   /// @param[in] axes axes specifications
   /// @param[in] image input image cube: u,v,pol,chan
   virtual void initialiseDegrid(const scimath::Axes& axes,
					const casacore::Array<imtype>& image);

   /// @brief make context-dependant changes to the gridder behaviour
   /// @param[in] context context description
   virtual void customiseForContext(const std::string &context);

   /// @brief set visibility weights
   /// @param[in] viswt shared pointer to visibility weights
   virtual void initVisWeights(const IVisWeights::ShPtr &viswt);

   /// @brief degrid the visibility data.
   /// @param[in] acc non-const data accessor to work with
   virtual void degrid(accessors::IDataAccessor& acc);

   /// @brief finalise degridding
   virtual void finaliseDegrid();

   /// @brief check whether the model is empty
   /// @details A simple check allows us to bypass heavy calculations if the input model
   /// is empty (all pixels are zero). This makes sense for degridding only.
   /// @brief true, if the model is empty
   virtual bool isModelEmpty() const;

   /// @ brief check whether the adapter can stack
   virtual bool isStackingGridder() const;

private:

   /// @brief wrapped gridder doing actual job
   boost::shared_ptr<IVisGridder> itsGridder;

   /// @brief flag that the model is empty for degridding
   /// @details It allows to bypass expensive image regridding
   bool itsModelIsEmpty;
   
   /// @brief buffer for the input visibilites
   /// @details this buffer should be reused

   std::vector<accessors::MemBufferDataAccessor> itsVisBuffer;
 
   
   };
} // namespace synthesis

} // namespace askap

#endif // #ifndef STACKING_GRIDDER_ADAPTER_H
