/// @file
///
/// Provides generic methods for parallel algorithms
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///
#ifndef ASKAP_SYNTHESIS_SYNPARALLEL_H_
#define ASKAP_SYNTHESIS_SYNPARALLEL_H_

#include <fitting/Equation.h>
#include <fitting/Solver.h>
#include <fitting/INormalEquations.h>
#include <fitting/Params.h>
#include <Common/ParameterSet.h>
#include <gridding/IVisGridder.h>


#include <askapparallel/AskapParallel.h>
#include <measures/Measures/MFrequency.h>

namespace askap
{
  namespace synthesis
  {
    /// @brief Support for parallel algorithms 
    ///
    /// @details Support for parallel applications in the area.
    /// An application is derived from this abstract base. The model used is that the
    /// application has many workers and one master, running in separate MPI processes
    /// or in one single thread. The master is the master so the number of processes
    /// is one more than the number of workers. 
    ///
    /// If the number of nodes is 1 then everything occurs in the same process with
    /// no overall for transmission of model.
    ///
    /// @ingroup parallel
    class SynParallel
    {
  public:

      /// @brief Constructor 
      /// @details The first parameter is needed solely for MPI, the second
      /// is the parset to be used in derived classes
      /// @param[in] comms communications object
      /// @param[in] parset parameter set      
      SynParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset);

      virtual ~SynParallel();

      /// Return the model
      askap::scimath::Params::ShPtr& params();

      /// @brief Broadcast the model to all workers or a group of workers
      void broadcastModel();

      /// @brief Receive the model from the master
      void receiveModel();

      /// Substitute %w by worker number, and %n by number of workers (one less than the number
      // of nodes). This allows workers to do different work! This just calls
      // through to the AskapParallel version of substitute()
      std::string substitute(const std::string& s) const;

  protected:
      /// @brief helper method to indentify model parameters to broadcast
      /// @details We use itsModel to buffer some derived images like psf, weights, etc
      /// which are not required for prediffers. It just wastes memory and CPU time if
      /// we broadcast them. At the same time, some auxilliary parameters like peak
      /// residual value need to be broadcast (so the major cycle can terminate in workers).
      /// This method returns the vector with all parameters to be broadcast. By default
      /// it returns all parameter names. This method is supposed to be overridden in
      /// derived classes (e.g. ImagerParallel) where a different behavior is needed.
      /// @return a vector with parameters to broadcast
      virtual std::vector<std::string> parametersToBroadcast() const;
  
      /// @brief actual implementation of the model broadcast
      /// @details This method is only supposed to be called from the master.
      /// @param[in] model the model to send
      void broadcastModelImpl(const scimath::Params &model);

      /// @brief actual implementation of the model receive
      /// @details This method is only supposed to be called from workers. 
      /// There should be one to one match between the number of calls to 
      /// broadcastModelImpl and receiveModelImpl.
      /// @param[in] model the model to fill
      void receiveModelImpl(scimath::Params &model);
      
      
      /// @brief obtain parameter set
      /// @details to be used in derived classes
      /// @return reference to the parameter set object
      inline const LOFAR::ParameterSet& parset() const { return itsParset;}

      /// @brief read the models from parset file to the given params object
      /// @details The model can be composed from both images and components. This
      /// method populates Params object by adding model data read from the parset file.
      /// The model is given by shared pointer because the same method can be used for both
      /// simulations and calibration (the former populates itsModel, the latter populates
      /// itsPerfectModel) 
      /// @param[in] pModel shared pointer to the params object (must exist)
      void readModels(const scimath::Params::ShPtr &pModel) const;

      /// The model
      askap::scimath::Params::ShPtr itsModel;

      /// Class for communications
      askap::askapparallel::AskapParallel& itsComms;

      /// obtain frequency reference frame
      inline casa::MFrequency::Ref getFreqRefFrame() const { return itsFreqRefFrame;}

      /// @brief helper method to create and configure gridder
      /// @details It is expected to be called from the constructor of derived classes
      /// @param[in] comms communications object
      /// @param[in] parset parameter set      
      static IVisGridder::ShPtr createGridder(const askap::askapparallel::AskapParallel& comms, 
                           const LOFAR::ParameterSet& parset);
  private:
      /// @brief parameter set to get the parameters from
      LOFAR::ParameterSet itsParset;
 
      /// @brief reference frame for frequency
      /// @details We may want to simulate/image in different reference frames.
      /// This field contains the reference frame selected in the parset.
      casa::MFrequency::Ref itsFreqRefFrame;    
    };

  }
}
#endif
