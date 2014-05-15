/// @file
///
/// Support for parallel statistics accululation to advise on imaging parameters
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

#include <parallel/AdviseParallel.h>
#include <askap/AskapError.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>
#include <dataaccess/SharedIter.h>
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

#include <fitting/INormalEquations.h>

#include <casa/aips.h>
#include <casa/OS/Timer.h>


#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief a helper adapter to reuse the existing MW framework.
/// @details We could've moved it to a separate file, but it is used
/// only in this particular cc file at the moment.
struct EstimatorAdapter : public scimath::INormalEquations {

  /// @brief constructor
  /// @details
  /// @param[in] estimator statistics estimator to work with (reference semantics)
  explicit EstimatorAdapter(const boost::shared_ptr<VisMetaDataStats> &estimator) :
           itsEstimator(estimator) {ASKAPDEBUGASSERT(itsEstimator);}

  /// @brief Clone this into a shared pointer
  /// @details "Virtual constructor" - creates a copy of this object. Derived
  /// classes must override this method to instantiate the object of a proper 
  /// type.
  virtual INormalEquations::ShPtr clone() const
  {
    boost::shared_ptr<VisMetaDataStats> newEstimator(new VisMetaDataStats(*itsEstimator));
    boost::shared_ptr<EstimatorAdapter> result(new EstimatorAdapter(newEstimator));
    return result;
  }

  /// @brief reset the normal equation object
  /// @details After a call to this method the object has the same pristine
  /// state as immediately after creation with the default constructor
  virtual void reset()
      { itsEstimator->reset(); }
  
  /// @brief Merge these normal equations with another
  /// @details Combining two normal equations depends on the actual class type
  /// (different work is required for a full matrix and for an approximation).
  /// This method must be overriden in the derived classes for correct 
  /// implementation. 
  /// This means that we just add
  /// @param[in] src an object to get the normal equations from
  virtual void merge(const INormalEquations& src)  {
    try {
       const EstimatorAdapter& ea = dynamic_cast<const EstimatorAdapter&>(src);
       itsEstimator->merge(*(ea.itsEstimator));
    }
    catch (const std::bad_cast &bc) {
       ASKAPTHROW(AskapError, "Unsupported type of normal equations used with the estimator adapter: "<<bc.what());
    }
  }
  
  /// @brief write the object to a blob stream
  /// @param[in] os the output stream
  virtual void writeToBlob(LOFAR::BlobOStream& os) const 
    { itsEstimator->writeToBlob(os); }

  /// @brief read the object from a blob stream
  /// @param[in] is the input stream
  /// @note Not sure whether the parameter should be made const or not 
  virtual void readFromBlob(LOFAR::BlobIStream& is) 
    { itsEstimator->readFromBlob(is); }
      
  
  /// @brief obtain shared pointer
  /// @return shared pointer to the estimator
  inline boost::shared_ptr<VisMetaDataStats> get() const { return itsEstimator;}
  
  
  /// @brief stubbed method for this class
  /// @return nothing, throws an exception
  virtual const casa::Matrix<double>& normalMatrix(const std::string &, 
                        const std::string &) const 
      { ASKAPTHROW(AskapError, "Method is not supported"); }
      
  /// @brief stubbed method for this class
  /// @return nothing, throws an exception
  virtual const casa::Vector<double>& dataVector(const std::string &) const
      { ASKAPTHROW(AskapError, "Method is not supported"); }
                        
  /// @brief stubbed method for this class
  /// @return nothing, throws an exception
  virtual std::vector<std::string> unknowns() const 
      { ASKAPTHROW(AskapError, "Method is not supported"); }
      
private:
  /// @brief estimator to work with, reference semantics
  boost::shared_ptr<VisMetaDataStats> itsEstimator;      
};

// actual AdviseParallel implementation 

/// @brief Constructor from ParameterSet
/// @details The parset is used to construct the internal state. We could
/// also support construction from a python dictionary (for example).
/// The command line inputs are needed solely for MPI - currently no
/// application specific information is passed on the command line.
/// @param comms communication object 
/// @param parset ParameterSet for inputs
AdviseParallel::AdviseParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset) :
    MEParallelApp(comms, addMissingFields(parset)), itsTangentDefined(false)
{
   itsWTolerance = parset.getDouble("wtolerance",-1.);
   if (parset.isDefined("tangent")) {
       const std::vector<std::string> direction = parset.getStringVector("tangent");
       ASKAPCHECK(direction.size() == 3, "Direction should have exactly 3 parameters, you have "<<direction.size());
       ASKAPCHECK(direction[2] == "J2000", "Only J2000 is implemented at the moment, you have requested "<<direction[2]);
      
       const double ra = SynthesisParamsHelper::convertQuantity(direction[0],"rad");
       const double dec = SynthesisParamsHelper::convertQuantity(direction[1],"rad");
       itsTangent = casa::MVDirection(ra,dec);
       itsTangentDefined = true;
   }
   itsNe.reset();
}    

/// @brief make the estimate
/// @details This method iterates over one or more datasets, accumulates and aggregates statistics. If
/// tangent point is not defined, two iterations are performed. The first one is to estimate the tangent
/// point and the second to obtain  
void AdviseParallel::estimate()
{
   if (itsTangentDefined) {
       ASKAPLOG_INFO_STR(logger, "Using explicitly defined tangent point "<<printDirection(itsTangent)<<" (J2000)");
       itsEstimator.reset(new VisMetaDataStats(itsTangent, itsWTolerance));
       // we only need one iteration here
   } else {
       // the first iteration is just to get an estimate for the tangent point
       itsEstimator.reset(new VisMetaDataStats);
   }
   calcNE();
   if (!itsTangentDefined) {
       // second iteration is necessary
       if (itsComms.isMaster()) {
           ASKAPDEBUGASSERT(params());
           params()->add("tangent",itsEstimator->centre().get());
       }
       broadcastModel();
       receiveModel();
       ASKAPCHECK(params()->has("tangent"), "tangent is not defined. There is likely to be a problem with model broadcast/receive");
       const casa::Vector<casa::Double> tangent = params()->value("tangent");
       ASKAPCHECK(tangent.nelements() == 2, "Expect a 2-element vector for tangent, you have "<<tangent);       
       itsTangent = casa::MVDirection(tangent);
       itsTangentDefined = true;
       ASKAPLOG_INFO_STR(logger, "Using tangent "<<printDirection(itsTangent)<<" (estimated most central direction)");
       // now all ranks should have the same value of itsTangent & itsTangentDefined, ready for the second iteration
       itsEstimator.reset(new VisMetaDataStats(itsTangent, itsWTolerance));
       calcNE();        
   }
}
   
/// @brief perform the accumulation for the given dataset
/// @details This method iterates over the given dataset, predicts visibilities according to the
/// model and subtracts these model visibilities from the original visibilities in the dataset.
/// This is the core operation of the doSubtraction method, which manages the parallel aspect of it.
/// All actual calculations are done inside this helper method.
/// @param[in] ms measurement set name
void AdviseParallel::calcOne(const std::string &ms)
{
   casa::Timer timer;
   timer.mark();
   ASKAPLOG_INFO_STR(logger, "Performing iteration to accumulate metadata statistics for " << ms );
   ASKAPDEBUGASSERT(itsEstimator);
   
   accessors::TableDataSource ds(ms, accessors::TableDataSource::MEMORY_BUFFERS, dataColumn());
   ds.configureUVWMachineCache(uvwMachineCacheSize(),uvwMachineCacheTolerance());      
   accessors::IDataSelectorPtr sel=ds.createSelector();
   sel << parset();
   accessors::IDataConverterPtr conv=ds.createConverter();
   conv->setFrequencyFrame(getFreqRefFrame(), "Hz");
   conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
   conv->setEpochFrame(); // time since 0 MJD
   accessors::IDataSharedIter it=ds.createIterator(sel, conv);
   for (; it.hasMore(); it.next()) {
        // iteration over the dataset
        itsEstimator->process(*it);
   }
   
   ASKAPLOG_INFO_STR(logger, "Finished iteration for "<< ms << " in "<< timer.real()
                   << " seconds ");    
}
      
/// @brief calculate "normal equations", i.e. statistics for this dataset
void AdviseParallel::calcNE()
{
   ASKAPDEBUGASSERT(itsEstimator);
   itsNe.reset(new EstimatorAdapter(itsEstimator));
   if (itsComms.isWorker()) {
       ASKAPCHECK(itsNe, "Statistics estimator (stored as NormalEquations) is not defined");
       if (itsComms.isParallel()) {
           calcOne(measurementSets()[itsComms.rank()-1]);
           sendNE();
       } else {
          for (size_t iMs=0; iMs<measurementSets().size(); ++iMs) {
               calcOne(measurementSets()[iMs]);
          }
       }       
   }
   if (itsComms.isMaster()) {
       receiveNE();
   }
   // after reduction the adapter may have a different object (may create copies and reset them)
   // therefore, we update the shared pointer which would take care of the reference counting for us
   // and destroy the original object, if necessary.
   ASKAPDEBUGASSERT(itsNe);
   boost::shared_ptr<EstimatorAdapter> ea = boost::dynamic_pointer_cast<EstimatorAdapter>(itsNe);
   ASKAPASSERT(ea);
   itsEstimator = ea->get();
}

/// @brief helper method to get statistics estimator
/// @return const reference to the statistics estimator
/// @note An exception is thrown if the estimator is not defined or the method is called from
/// worker process.
const VisMetaDataStats& AdviseParallel::estimator() const
{
   ASKAPCHECK(itsComms.isMaster(), "AdviseParallel::estimator() is supposed to be called from the master process only");
   ASKAPCHECK(itsEstimator, "estimator is not defined!");
   return *itsEstimator;
}

/// @brief summarise stats into log
/// @details This method just summarises all stats accumulated in the call to estimate() method
/// into log. Nothing is done for worker process.
void AdviseParallel::summary() const
{
   if (itsComms.isMaster()) {
       const VisMetaDataStats &stats = estimator();
       ASKAPLOG_INFO_STR(logger, "AdviseParallel::summary - statistics for the visibility dataset:");
       ASKAPLOG_INFO_STR(logger, "  Total number of visibilities (ignoring polarisation): "<<stats.nVis());
       ASKAPLOG_INFO_STR(logger, "  Largest U: "<<stats.maxU()<<" wavelengths");
       ASKAPLOG_INFO_STR(logger, "  Largest V: "<<stats.maxV()<<" wavelengths");
       ASKAPLOG_INFO_STR(logger, "  Largest W: "<<stats.maxW()<<" wavelengths");
       if (itsWTolerance >= 0.) {
           ASKAPLOG_INFO_STR(logger, "  Largest residual W: "<<stats.maxResidualW()<<" wavelengths");
       } else {
           ASKAPLOG_INFO_STR(logger, "  Largest residual W - not determined");
       }
       // the following assumes the standard accessor units selection
       ASKAPLOG_INFO_STR(logger, "  Frequency from "<<stats.minFreq()/1e6<<" to "<<stats.maxFreq()/1e6<<" MHz");
       ASKAPLOG_INFO_STR(logger, "  Number of antennas: "<<stats.nAntennas());
       ASKAPLOG_INFO_STR(logger, "  Number of beams: "<<stats.nBeams());
       if (itsTangentDefined) {
           ASKAPLOG_INFO_STR(logger, "  Assumed tangent point: "<<printDirection(itsTangent)<<" (J2000)");
       }
       ASKAPLOG_INFO_STR(logger, "  Most central pointing direction of the field: "<<printDirection(stats.centre())<<" (J2000)");
       const std::pair<double,double> offsets = stats.maxOffsets();
       ASKAPLOG_INFO_STR(logger, "  Largest beam offsets from the central direction: "<<offsets.first/casa::C::pi*180.<<" , "<<offsets.second/casa::C::pi*180.<<" deg");
       const double cellSize = stats.squareCellSize(); // in arcsec
       ASKAPLOG_INFO_STR(logger, "  Estimated square cell size: "<<cellSize<<" arcsec");

       if (itsTangentDefined) {
           ASKAPLOG_INFO_STR(logger, "  Distance of the 'average' pointing direction from the tangent point: "<<itsTangent.separation(stats.centre())*180./casa::C::pi<<" deg");
       }
       for (int stage = 0; stage<2; ++stage) {
            const std::string what = stage == 0 ? "'average' pointing" : "tangent point";
            const double fieldSize = stats.squareFieldSize(stage == 1); // in deg
            ASKAPLOG_INFO_STR(logger, "  Estimated square field size (about the "<<what<<"): "<<fieldSize<<" deg");
            const long imgSize = long(fieldSize * 3600 / cellSize) + 1;
            ASKAPLOG_INFO_STR(logger, "  --- or the minimum image size: "<<imgSize<<" x "<<imgSize<<" pixels");
       }
   }
} 

/// @brief a hopefully temporary method to define missing fields in parset
/// @details We reuse some code for general synthesis application, but it requires some
/// parameters (like gridder) to be defined. This method fills the parset with stubbed fields.
/// Hopefully, it is a temporary approach.
/// @param parset ParameterSet for inputs
/// @return new parset 
LOFAR::ParameterSet AdviseParallel::addMissingFields(const LOFAR::ParameterSet& parset)
{
  LOFAR::ParameterSet result(parset);
  if (!result.isDefined("gridder")) {
      result.add("gridder","SphFunc");
  }
  return result;
}


} // namespace synthesis

} // namespace askap
