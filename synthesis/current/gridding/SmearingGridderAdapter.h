/// @file
///
/// @brief Gridder adapter to simulate smearing effects 
/// @details When in the degridding mode, this adapter runs the wrapped 
/// gridder at finer time and frequency resolution and averages the result
/// to simulate bandwidth at time-average smearing. This can be done for both
/// simulation and the forward prediction step of the ordinary imaging.
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
///

#ifndef SMEARING_GRIDDER_ADAPTER_H
#define SMEARING_GRIDDER_ADAPTER_H

#include <gridding/IVisGridder.h>
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief Gridder adapter to simulate smearing effects 
/// @details When in the degridding mode, this adapter runs the wrapped 
/// gridder at finer time and frequency resolution and averages the result
/// to simulate bandwidth at time-average smearing. This can be done for both
/// simulation and the forward prediction step of the ordinary imaging.
/// @ingroup gridding
class SmearingGridderAdapter : virtual public IVisGridder 
{
public:
   /// @brief initialise the adapter
   /// @details
   /// @param[in] gridder a shared pointer to the gridder to be wrapped by this adapter
   /// @param[in] bandwidth effective bandwidth of a single spectral channel (we don't have access to it and
   /// in general it may be different from spectral resolution), same units as used in the other part of 
   /// the gridding package (Hz)
   /// @param[in] nFreqSteps number of frequency points in the simulation (default is 1, i.e. no simulation)
   /// @note time-average smearing is not yet implemented
   SmearingGridderAdapter(const boost::shared_ptr<IVisGridder> &gridder,
           const double bandwidth,
           const casa::uInt nFreqSteps = 1);


   /// @brief copy constructor
   /// @details We need this because the gridder doing actual work is held by a shared pointer,
   /// which is a non-trivial type
   /// @param[in] other an object to copy from
   SmearingGridderAdapter(const SmearingGridderAdapter &other);
   
   /// @brief clone a copy of this gridder
   /// @return shared pointer to the clone
   virtual boost::shared_ptr<IVisGridder> clone();

   /// @brief initialise the gridding
   /// @details
   /// @param[in] axes axes specifications
   /// @param[in] shape Shape of output image: cube: u,v,pol,chan
   /// @param[in] dopsf Make the psf?
   virtual void initialiseGrid(const scimath::Axes& axes,
                const casa::IPosition& shape, const bool dopsf = true);

   /// @brief grid the visibility data.
   /// @param[in] acc const data accessor to work with
   virtual void grid(accessors::IConstDataAccessor& acc);

   /// @brief form the final output image
   /// @param[in] out output double precision image or PSF
   virtual void finaliseGrid(casa::Array<double>& out);

   /// @brief finalise weights
   /// @details Form the sum of the convolution function squared, multiplied by the weights for each
   /// different convolution function. This is used in the evaluation of the second derivative.
   /// @param[in] out output double precision sum of weights images
   virtual void finaliseWeights(casa::Array<double>& out);

   /// @brief initialise the degridding
   /// @param[in] axes axes specifications
   /// @param[in] image input image cube: u,v,pol,chan
   virtual void initialiseDegrid(const scimath::Axes& axes,
					const casa::Array<double>& image);

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
      
private:

   /// @brief wrapped gridder doing actual job
   boost::shared_ptr<IVisGridder> itsGridder;
   
   /// @brief bandwidth of a single channel (same units as everywhere else in gridding -> Hz)
   const double itsBandwidth;
   
   /// @brief number of frequency steps to simulate (1 means a void operation, pass accessor as it is)
   const casa::uInt itsNFreqSteps;
                  
   /// @brief flag that the model is empty for degridding
   /// @details It allows to bypass expensive image regridding
   bool itsModelIsEmpty;
   };


} // namespace synthesis

} // namespace askap

#endif // #ifndef SMEARING_GRIDDER_ADAPTER_H

