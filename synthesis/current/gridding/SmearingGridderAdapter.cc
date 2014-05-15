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

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.snapshotimaginggridderadapter");

#include <gridding/SmearingGridderAdapter.h>
#include <askap/AskapError.h>
#include <dataaccess/SmearingAccessorAdapter.h>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

/// @brief initialise the adapter
/// @details
/// @param[in] gridder a shared pointer to the gridder to be wrapped by this adapter
/// @param[in] bandwidth effective bandwidth of a single spectral channel (we don't have access to it and
/// in general it may be different from spectral resolution), same units as used in the other part of 
/// the gridding package (Hz)
/// @param[in] nFreqSteps number of frequency points in the simulation (default is 1, i.e. no simulation)
/// @note time-average smearing is not yet implemented
SmearingGridderAdapter::SmearingGridderAdapter(const boost::shared_ptr<IVisGridder> &gridder,
           const double bandwidth, const casa::uInt nFreqSteps) :
           itsBandwidth(bandwidth), itsNFreqSteps(nFreqSteps), itsModelIsEmpty(true)
{
  ASKAPCHECK(gridder, "SmearingGridderAdapter should only be initialised with a valid gridder");
  itsGridder = gridder->clone();
  ASKAPCHECK(itsNFreqSteps > 0, "SmearingGridderAdapter should only be initialised with a positive number of frequency integration steps");
  ASKAPCHECK(itsBandwidth > 0, "SmearingGridderAdapter should only be initialised with a positive bandwidth");
}           


/// @brief copy constructor
/// @details We need this because the gridder doing actual work is held by a shared pointer,
/// which is a non-trivial type
/// @param[in] other an object to copy from
SmearingGridderAdapter::SmearingGridderAdapter(const SmearingGridderAdapter &other) :
    itsBandwidth(other.itsBandwidth), itsNFreqSteps(other.itsNFreqSteps), itsModelIsEmpty(other.itsModelIsEmpty)   
{   
  ASKAPCHECK(other.itsGridder, 
      "copy constructor of SmearingGridderAdapter got an object somehow set up with an empty gridder");
  itsGridder = other.itsGridder->clone();  
}  
   
/// @brief clone a copy of this gridder
/// @return shared pointer to the clone
boost::shared_ptr<IVisGridder> SmearingGridderAdapter::clone() 
{  
  boost::shared_ptr<SmearingGridderAdapter> newOne(new SmearingGridderAdapter(*this));
  return newOne;
}

/// @brief initialise the gridding
/// @details
/// @param[in] axes axes specifications
/// @param[in] shape Shape of output image: cube: u,v,pol,chan
/// @param[in] dopsf Make the psf?
void SmearingGridderAdapter::initialiseGrid(const scimath::Axes& axes,
                const casa::IPosition& shape, const bool dopsf)
{
   ASKAPDEBUGASSERT(itsGridder);
   itsGridder->initialiseGrid(axes,shape,dopsf);
}

/// @brief grid the visibility data.
/// @param[in] acc const data accessor to work with
void SmearingGridderAdapter::grid(accessors::IConstDataAccessor& acc)
{
   ASKAPDEBUGASSERT(itsGridder);
   itsGridder->grid(acc);
}

/// @brief form the final output image
/// @param[in] out output double precision image or PSF
void SmearingGridderAdapter::finaliseGrid(casa::Array<double>& out)
{
   ASKAPDEBUGASSERT(itsGridder);
   itsGridder->finaliseGrid(out);
}

/// @brief finalise weights
/// @details Form the sum of the convolution function squared, multiplied by the weights for each
/// different convolution function. This is used in the evaluation of the second derivative.
/// @param[in] out output double precision sum of weights images
void SmearingGridderAdapter::finaliseWeights(casa::Array<double>& out)
{
   ASKAPDEBUGASSERT(itsGridder);
   itsGridder->finaliseWeights(out);
}

/// @brief initialise the degridding
/// @param[in] axes axes specifications
/// @param[in] image input image cube: u,v,pol,chan
void SmearingGridderAdapter::initialiseDegrid(const scimath::Axes& axes,
					const casa::Array<double>& image)
{
   ASKAPDEBUGASSERT(itsGridder);
   itsGridder->initialiseDegrid(axes,image);
   itsModelIsEmpty = itsGridder->isModelEmpty();
}

/// @brief make context-dependant changes to the gridder behaviour
/// @param[in] context context description
void SmearingGridderAdapter::customiseForContext(const std::string &context)
{
   ASKAPDEBUGASSERT(itsGridder);
   itsGridder->customiseForContext(context);
}
			
/// @brief set visibility weights
/// @param[in] viswt shared pointer to visibility weights
void SmearingGridderAdapter::initVisWeights(const IVisWeights::ShPtr &viswt)
{
   ASKAPDEBUGASSERT(itsGridder);
   itsGridder->initVisWeights(viswt);
}

/// @brief degrid the visibility data.
/// @param[in] acc non-const data accessor to work with  
void SmearingGridderAdapter::degrid(accessors::IDataAccessor& acc)
{
   ASKAPDEBUGASSERT(itsGridder);
   if (itsNFreqSteps < 2) {
       // void operation, as we have a single integration step only
       itsGridder->degrid(acc); 
   } else {
       // we need a new adapter here allowing to change frequency information
       // this part is to be written
       SmearingAccessorAdapter accBuffer(acc);
       accBuffer.useFrequencyBuffer();
       const casa::uInt centralStep = itsNFreqSteps / 2; 
       ASKAPDEBUGASSERT(centralStep > 0);
       const double freqInc = itsBandwidth / double(itsNFreqSteps - 1);
       // for an odd number of integration steps we get one point exactly at the centre of the channel,
       // for an even number - equidistant on both sides, hence the offset
       const double freqOff = itsNFreqSteps % 2 == 0 ? freqInc / 2. : 0.;
       const casa::uInt nChan = acc.nChannel();
       ASKAPDEBUGASSERT(accBuffer.rwFrequency().nelements() == nChan);
       accBuffer.rwVisibility().set(0.0);
       for (casa::uInt it = 0; it < itsNFreqSteps; ++it) {
            // do the integration via Trapezium method
            // first deal with the first and the last point because we scale down the result later
            const casa::uInt step = (it == 0 ? 0 : (it == 1 ? itsNFreqSteps - 1 : it - 1));
            ASKAPDEBUGASSERT(step < itsNFreqSteps); 
            // fill the new frequency vector
            for (casa::uInt chan = 0; chan<nChan; ++chan) {
                 accBuffer.rwFrequency()[chan] = acc.frequency()[chan] + freqOff + 
                          double(casa::Int(step) - casa::Int(centralStep)) * freqInc;
            }
            itsGridder->degrid(accBuffer);
            if (it == 1) {
                // first and last have a weight of 0.5, other elements have a weight of 1.
                // this would automatically revert to the rectangle method if only two points are available
                accBuffer.rwVisibility() *= float(0.5);
            }
       }
       acc.rwVisibility() = accBuffer.visibility() / float(itsNFreqSteps - 1);       
   }
}

/// @brief finalise degridding
void SmearingGridderAdapter::finaliseDegrid()
{
   ASKAPDEBUGASSERT(itsGridder);
   itsGridder->finaliseDegrid();
}

/// @brief check whether the model is empty
/// @details A simple check allows us to bypass heavy calculations if the input model
/// is empty (all pixels are zero). This makes sense for degridding only.
/// @brief true, if the model is empty
bool SmearingGridderAdapter::isModelEmpty() const
{
  return itsModelIsEmpty;
}


