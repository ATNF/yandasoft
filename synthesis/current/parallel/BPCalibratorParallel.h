/// @file
///
/// BPCalibratorParallel: part of the specialised tool to do optimised bandpass calibration with
/// limited functionality. Unlike CalibratorParallel, this class
///
///      * solves for bandpass only
///      * works only with preaveraging calibration approach
///      * does not support multiple chunks in time (i.e. only one solution is made for the whole dataset)
///      * does not support data distribution except per beam 
///      * does not support a distributed model (e.h. with individual workers dealing with individual Taylor terms)
///      * does not require exact match between number of workers and number of channel chunks, data are dealt with
///        serially by each worker with multiple iterations over data, if required.
///      * solves normal equations at the worker level in the parallel case
///
/// This specialised tool matches closely BETA needs and will be used for BETA initially (at least until we converge
/// on the best approach to do bandpass calibration). The lifetime of this tool is uncertain at present. In many
/// instances the code is quick and dirty, just to suit our immediate needs.  
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
#ifndef ASKAP_SYNTHESIS_BP_CALIBRATOR_PARALLEL_H
#define ASKAP_SYNTHESIS_BP_CALIBRATOR_PARALLEL_H

// own includes
#include <askapparallel/AskapParallel.h>
#include <parallel/MEParallelApp.h>
#include <Common/ParameterSet.h>
#include <gridding/IVisGridder.h>
#include <measurementequation/IMeasurementEquation.h>
#include <dataaccess/SharedIter.h>
#include <calibaccess/ICalSolutionSource.h>
#include <utils/MultiDimPosIter.h>


// std includes
#include <utility>

// boost includes
#include <boost/shared_ptr.hpp>

namespace askap
{
  namespace synthesis
  {
    /// @brief Support for parallel algorithms implementing bandpass calibration
    /// @details
    /// BPCalibratorParallel: part of the specialised tool to do optimised bandpass calibration with
    /// limited functionality. Unlike CalibratorParallel, this class
    ///
    ///      * solves for bandpass only
    ///      * works only with preaveraging calibration approach
    ///      * does not support multiple chunks in time (i.e. only one solution is made for the whole dataset)
    ///      * does not support data distribution except per beam 
    ///      * does not support a distributed model (e.h. with individual workers dealing with individual Taylor terms)
    ///      * does not require exact match between number of workers and number of channel chunks, data are dealt with
    ///        serially by each worker with multiple iterations over data, if required.
    ///      * solves normal equations at the worker level in the parallel case
    ///
    /// This specialised tool matches closely BETA needs and will be used for BETA initially (at least until we converge
    /// on the best approach to do bandpass calibration). The lifetime of this tool is uncertain at present. In many
    /// instances the code is quick and dirty, just to suit our immediate needs.
    ///  
    /// @ingroup parallel
    class BPCalibratorParallel : public MEParallelApp
    {
      public:

      /// @brief Constructor from ParameterSet
      /// @details The parset is used to construct the internal state. We could
      /// also support construction from a python dictionary (for example).
      /// @param[in] comms communication object
      /// @param[in] parset ParameterSet for inputs
      BPCalibratorParallel(askap::askapparallel::AskapParallel& comms,
          const LOFAR::ParameterSet& parset);

      /// @brief method which does the main job
      /// @details it iterates over all channels/beams and writes the result.
      /// In the parallel mode each worker iterates over their own portion of work and
      /// then sends the result to master for writing.
      void run(); 

      protected:      

      // virtual methods of the abstract base, define them as protected because they
      // are no longer supposed to be called directly from the application level
    
      /// @brief Calculate the normal equations (runs in workers)
      /// @details Model, either image-based or component-based, is used in conjunction with 
      /// CalibrationME to calculate the generic normal equations. 
      virtual void calcNE();

      /// @brief Solve the normal equations (runs in workers)
      /// @details Parameters of the calibration problem are solved for here
      virtual void solveNE();
      
      /// @brief Write the results (runs in master)
      /// @details The solution (calibration parameters) is reported via solution accessor
      /// @param[in] postfix a string to be added to the file name (unused in this class)
	  virtual void writeModel(const std::string &postfix = std::string());
                    
      /// @brief create measurement equation
      /// @details This method initialises itsEquation with shared pointer to a proper type.
      /// It uses internal flags to create a correct type (i.e. polarisation calibration or
      /// just antenna-based gains). Parameters are passed directly to the constructor of 
      /// CalibrationME template.
      /// @param[in] dsi data shared iterator 
      /// @param[in] perfectME uncorrupted measurement equation
      void createCalibrationME(const accessors::IDataSharedIter &dsi, 
                const boost::shared_ptr<IMeasurementEquation const> &perfectME);
  
      /// @brief helper method to rotate all phases
      /// @details This method rotates the phases of all gains in itsModel
      /// to have the phase of itsRefGain exactly 0. This operation does
      /// not seem to be necessary for SVD solvers, however it simplifies
      /// "human eye" analysis of the results (otherwise the phase degeneracy
      /// would make the solution different from the simulated gains).
      /// @note The method throws exception if itsRefGain is not among
      /// the parameters of itsModel
      void rotatePhases();
      
      /// @brief helper method to extract solution time from NE.
      /// @details To be able to time tag the calibration solutions we add
      /// start and stop times extracted from the dataset as metadata to normal
      /// equations. It allows us to send these times to the master, which
      /// ultimately writes the calibration solution. Otherwise, these times 
      /// could only be obtained in workers who deal with the actual data.
      /// @return solution time (seconds since 0 MJD)
      /// @note if no start/stop time metadata are present in the normal equations
      /// this method returns 0.
      double solutionTime() const;
            
  private:
      /// @brief helper method to remove the dataset name from a parset
      /// @details We deal with multiple measurement sets in a dit different
      /// way from the other synthesis applications (they are not per worker
      /// here). This method allows to remove the string with measurement sets
      /// in the parset passed to base classes and replace it by empty string
      /// @param[in] parset input parset 
      /// @return a copy without the dataset keyword
      static LOFAR::ParameterSet emptyDatasetKeyword(const LOFAR::ParameterSet &parset);


      /// @brief read the model from parset file and populate itsPerfectModel
      /// @details This method is common between several classes and probably
      /// should be pushed up in the class hierarchy
      inline void readModels() const { SynParallel::readModels(itsPerfectModel); }
      
      /// @brief number of antennas to solve for
      /// @return number of antennas to solve for
      inline casa::uInt nAnt() const { return parset().getInt32("nAnt", 36); }
      
      /// @brief number of beams to solve for
      /// @return number of beams to solve for
      inline casa::uInt nBeam() const { return parset().getInt32("nBeam", 1); }
      
      /// @brief number of channels to solve for
      /// @return number of channels to solve for
      inline casa::uInt nChan() const { return parset().getInt32("nChan", 304); }
      
      /// @brief extract current beam/channel pair from the iterator
      /// @details This method encapsulates interpretation of the output of itsWorkUnitIterator.cursor() for workers and
      /// in the serial mode. However, it extracts the current beam and channel info out of the model for the master
      /// in the parallel case. This is done because calibration data are sent to the master asynchronously and there is no
      /// way of knowing what iteration in the worker they correspond to without looking at the data.
      /// @return pair of beam (first) and channel (second) indices
      std::pair<casa::uInt, casa::uInt> currentBeamAndChannel() const; 

      /// @brief helper method to invalidate curremt solution
      void invalidateSolution(); 

      /// @brief verify that the current solution is valid
      /// @details We use a special keywork 'invalid' in the model to 
      /// signal that a particular solution failed. for whatever reason. 
      /// This flag is checked to avoid writing the solution (which would 
      /// automatically set validity flag
      /// @return true, if the current solution is valid
      bool validSolution() const;
 
      /// Calculate normal equations for one data set, channel and beam
      /// @param[in] ms Name of data set
      /// @param[in] chan channel to work with
      /// @param[in] beam beam to work with
      void calcOne(const std::string& ms, const casa::uInt chan, const casa::uInt beam);
      
      /// @brief send current model to the master
      /// @details This method is supposed to be called from workers in the parallel mode and
      /// sends the current results to the master rank
      void sendModelToMaster() const;
      
      /// @brief asynchronously receive model from one of the workers
      /// @details This method is supposed to be used in the master rank in the parallel mode. It
      /// waits until the result becomes available from any of the workers and then stores it 
      /// in itsModel. 
      void receiveModelFromWorker();
      
      /// uncorrupted model
      askap::scimath::Params::ShPtr itsPerfectModel;
      
      /// @brief reference antenna (index)
      /// @details Negative number means no referencing required
      int itsRefAntenna;
      
      /// @brief name of the parameter taken as a reference
      /// @details empty string means no referencing is required
      //wasim was here 
      /*
      std::string itsRefGain;
      */
      std::string itsRefGainXX;
      std::string itsRefGainYY;
                                   
      /// @brief solution source to store the result
      /// @details This object is initialised by the master. It stores the solution
      /// in parset file, casa table or a database.
      boost::shared_ptr<accessors::ICalSolutionSource> itsSolutionSource;

      /// @brief solution accessor used to stage the results in memory
      /// @details This object is initialised by the master. It provides a way to store
      /// the solutions in memory, until we write out at the end.
      boost::shared_ptr<accessors::ICalSolutionAccessor> itsSolAcc;
        
      /// @brief shared pointer to measurement equation correspondent to the perfect model
      /// @details It is handy to store the perfect measurement equation, so it is not
      /// recreated every time for each solution interval. 
      boost::shared_ptr<IMeasurementEquation const> itsPerfectME;
      
      /// @brief iterator over channels and beams
      /// @details This class allows us to split work domain between a number of workers (=iteration chunks)
      scimath::MultiDimPosIter itsWorkUnitIterator;
      
      /// @brief solution ID to work with
      /// @details This field should only be used if itsSolutionIDValid is true
      long itsSolutionID;
      
      /// @brief solution ID validity flag
      bool itsSolutionIDValid;
    };

  }
}
#endif // #ifndef ASKAP_SYNTHESIS_BP_CALIBRATOR_PARALLEL_H

