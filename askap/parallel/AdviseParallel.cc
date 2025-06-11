/// @file
///
/// Support for parallel statistics accumulation to advise on imaging parameters
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

#include <askap/parallel/AdviseParallel.h>
#include <askap/askap/AskapError.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");
#include <askap/profile/AskapProfiler.h>


#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/Solver.h>

#include <askap/scimath/utils/EstimatorAdapter.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/OS/Timer.h>

#include <Blob/BlobString.h>
#include <Blob/BlobIBufString.h>
#include <Blob/BlobOBufString.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobOStream.h>


#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

namespace askap {

// put specialisation of EstimatorAdapterMonitor in the scimath namespace where all related definitions are
namespace scimath {

/// @brief specialisation of default monitor class for the type used here
/// @details This monitor just prints appropriate messages to the log. If this specialisation is
/// removed (e.g. commented out), the default template will be used resulting in void operation which may
/// be even optimised by the compiler
template<>
struct EstimatorAdapterMonitor<EstimatorAdapter<synthesis::VisMetaDataStats> > {
   inline static void aboutToMerge(const EstimatorAdapter<synthesis::VisMetaDataStats>& thisEstimator,
                                   const EstimatorAdapter<synthesis::VisMetaDataStats>& otherEstimator) {
       ASKAPDEBUGASSERT(thisEstimator.get());
       ASKAPLOG_DEBUG_STR(logger, "Merge This Largest residual W: "<<thisEstimator.get()->maxW()<<" wavelengths, percentile = "<<thisEstimator.get()->wPercentile());
       ASKAPDEBUGASSERT(otherEstimator.get());
       ASKAPLOG_DEBUG_STR(logger, "With other Largest residual W: "<<otherEstimator.get()->maxW()<<" wavelengths, percentile = "<<otherEstimator.get()->wPercentile());
   };
}; // struct EstimatorAdapterMonitor

} // namespace scimath

namespace synthesis {

// actual AdviseParallel implementation

/// @brief Constructor from ParameterSet
/// @details The parset is used to construct the internal state. We could
/// also support construction from a python dictionary (for example).
/// The command line inputs are needed solely for MPI - currently no
/// application specific information is passed on the command line.
/// @param comms communication object
/// @param parset ParameterSet for inputs
AdviseParallel::AdviseParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset) :
    MEParallelApp(comms, addMissingFields(parset),true), itsTangentDefined(false)
{
   init(parset);
}

void AdviseParallel::init(const LOFAR::ParameterSet& parset)
{
   // update the parset in the base class
   setParset(addMissingFields(parset));
   itsWTolerance = parset.getDouble("wtolerance",-1.);
   if (parset.isDefined("tangent")) {
       const std::vector<std::string> direction = parset.getStringVector("tangent");
       ASKAPCHECK(direction.size() == 3, "Direction should have exactly 3 parameters, you have "<<direction.size());
       ASKAPCHECK(direction[2] == "J2000", "Only J2000 is implemented at the moment, you have requested "<<direction[2]);
       itsTangent = asMDirection(direction).getValue();
       itsTangentDefined = true;
   }
   // keep the historical default for now, taking the flags into account can be a lot slower
   itsIncludeFlaggedData = parset.getBool("includeflagged",true);
   ASKAPLOG_INFO_STR(logger, (itsIncludeFlaggedData ? "I" : "Not i")<< "ncluding flagged data in meta data statistics ");
   // set w percentile if present
   itsWPercentile = parset.getDouble("wpercentile",99.9)/100.0;
   ASKAPCHECK(itsWPercentile>0 && itsWPercentile<1,"wpercentile value needs to be between 0 and 100");
   itsNe.reset();
}

/// @brief make the estimate
/// @details This method iterates over one or more datasets, accumulates and aggregates statistics. If
/// tangent point is not defined, two iterations are performed. The first one is to estimate the tangent
/// point and the second to obtain the stats
void AdviseParallel::estimate()
{
    ASKAPTRACE("AdviseParallel::estimate");

   if (itsTangentDefined) {
       ASKAPLOG_INFO_STR(logger, "Using explicitly defined tangent point "<<printDirection(itsTangent)<<" (J2000)");
       itsEstimator.reset(new VisMetaDataStats(itsTangent, itsWTolerance, itsIncludeFlaggedData,itsWPercentile));
       // we only need one iteration here
   } else {
       // the first iteration is just to get an estimate for the tangent point
       itsEstimator.reset(new VisMetaDataStats(true,itsWPercentile));
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
       const casacore::Vector<imtype> tangent = params()->valueT("tangent");
       ASKAPCHECK(tangent.nelements() == 2, "Expect a 2-element vector for tangent, you have "<<tangent);
       itsTangent = casacore::MVDirection(tangent(0),tangent(1));
       itsTangentDefined = true;
       ASKAPLOG_INFO_STR(logger, "Using tangent "<<printDirection(itsTangent)<<" (estimated most central direction)");
       // now all ranks should have the same value of itsTangent & itsTangentDefined, ready for the second iteration
       itsEstimator.reset(new VisMetaDataStats(itsTangent, itsWTolerance, itsIncludeFlaggedData,itsWPercentile));
       calcNE();
   }
   // to synchronise stats held by all workers (wouldn't be necessary if the only thing we wanted was to print summary)
   broadcastStatistics();
}

/// @brief set estimator
/// @details Normally, the estimator is set inside estimate() call with the right parameters
/// (i.e. depending on whether the tangent point is known or not). For a greater reuse of the code it is
/// handy to be able to set a custom estimator and then use calcOne method directly for accumulation.
/// This method sets the estimator (to be reset if estimate() method is called).
/// @param[in] estimator shared pointer to custom estimator to work with
void AdviseParallel::setCustomEstimator(const boost::shared_ptr<VisMetaDataStats> &estimator)
{
   itsEstimator=estimator;
}

/// @brief helper method to broadcast statistics to all workers
/// @details It seems better conceptually, if all ranks hold the same statistics at the end of
/// the estimate() call. This allows workers to update their own parsets in the parallel way.
/// This helper method encapsulates all actions required to broadcast statistics estimator from
/// the master to all workers.
void AdviseParallel::broadcastStatistics()
{
  ASKAPDEBUGTRACE("AdviseParallel::broadcastStatistics");
  // only need to do something in the parallel case
  if (itsComms.isParallel()) {
      LOFAR::BlobString bs;
      bs.resize(0);
      ASKAPASSERT(itsEstimator);
      const int currentVersion = 1;
      if (itsComms.isMaster()) {
          // master has up to date stats, it sends them to all workers
          LOFAR::BlobOBufString bob(bs);
          LOFAR::BlobOStream out(bob);
          out.putStart("statistics", currentVersion);
          itsEstimator->writeToBlob(out);
          out.putEnd();
          itsComms.broadcastBlob(bs ,0);
      }
      if (itsComms.isWorker()) {
          // workers receive stats from the master
          itsComms.broadcastBlob(bs, 0);
          LOFAR::BlobIBufString bib(bs);
          LOFAR::BlobIStream in(bib);
          const int version=in.getStart("statistics");
          ASKAPASSERT(version == currentVersion);
          itsEstimator->readFromBlob(in);
          in.getEnd();
      }
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
   ASKAPTRACE("AdviseParallel::calcOne");

   casacore::Timer timer;
   timer.mark();
   ASKAPLOG_INFO_STR(logger, "Performing iteration to accumulate metadata statistics for " << ms);
   ASKAPDEBUGASSERT(itsEstimator);

   accessors::TableDataSource ds(ms, accessors::TableDataSource::MEMORY_BUFFERS, dataColumn());

   ASKAPLOG_DEBUG_STR(logger, "Initialised accessor");
   ds.configureUVWMachineCache(uvwMachineCacheSize(),uvwMachineCacheTolerance());

   ASKAPLOG_DEBUG_STR(logger, "Initialised UVW cache");
   accessors::IDataSelectorPtr sel=ds.createSelector();
   ASKAPDEBUGASSERT(sel);
   ASKAPLOG_DEBUG_STR(logger, "Initialised data selector");

   sel << parset();
   ASKAPLOG_DEBUG_STR(logger, "Filled selector\n" << parset());
   accessors::IDataConverterPtr conv=ds.createConverter();
   ASKAPDEBUGASSERT(conv);
   ASKAPLOG_DEBUG_STR(logger, "Initialised converter");
   conv->setFrequencyFrame(getFreqRefFrame(), "Hz");
   conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));
   conv->setEpochFrame(); // time since 0 MJD

   ASKAPLOG_DEBUG_STR(logger, "Initialise iterator");
   for (accessors::IConstDataSharedIter it=ds.createConstIterator(sel, conv); it.hasMore(); it.next()) {
        // iteration over the dataset
        itsEstimator->process(*it);

   }

   ASKAPLOG_INFO_STR(logger, "Finished iteration for "<< ms << " in "<< timer.real()
                   << " seconds");
}

/// @brief calculate "normal equations", i.e. statistics for this dataset
void AdviseParallel::calcNE()
{
   ASKAPDEBUGASSERT(itsEstimator);
   ASKAPTRACE("AdviseParallel::calcNE");

   itsNe.reset(new scimath::EstimatorAdapter<VisMetaDataStats>(itsEstimator));
   // the following call avoids merging of normal equations inside the solver (which we don't use)
   // this is an artefact of reusing Master/Worker framework
   itsSolver.reset(new scimath::Solver);
   if (itsComms.isWorker()) {
       ASKAPCHECK(itsNe, "Statistics estimator (stored as NormalEquations) is not defined");
       if (itsComms.isParallel()) {
           ASKAPLOG_DEBUG_STR(logger,"calcNE parallel mode  - rank = "<<itsComms.rank()<<" #ms = "<<measurementSets().size());
           // if we are in imager's joint imaging mode we do the MS's in serial
           // this is yet another patch for imager's broken advise design
           if (parset().getBool("updatedirection",false)) {
               for (size_t iMs=0; iMs<measurementSets().size(); ++iMs) {
                    calcOne(measurementSets()[iMs]);
               }
           } else {
               calcOne(measurementSets()[itsComms.rank()-1]);
           }
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
   boost::shared_ptr<scimath::EstimatorAdapter<VisMetaDataStats> > ea = boost::dynamic_pointer_cast<scimath::EstimatorAdapter<VisMetaDataStats> >(itsNe);
   ASKAPASSERT(ea);
   itsEstimator = ea->get();
   // we don't need normal equation adapter any more
   itsNe.reset();
}

/// @brief helper method to get statistics estimator
/// @return const reference to the statistics estimator
/// @note An exception is thrown if the estimator is not defined or the method is called from
/// worker process.
const VisMetaDataStats& AdviseParallel::estimator() const
{
   //broadcastStatistics() at the end of estimate() means that this can be called on workers as well. Needed for advise in Cimager.
   //ASKAPCHECK(itsComms.isMaster(), "AdviseParallel::estimator() is supposed to be called from the master process only");
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
       ASKAPLOG_INFO_STR(logger, "  Largest beam offsets from the central direction: "<<offsets.first/casacore::C::pi*180.<<" , "<<offsets.second/casacore::C::pi*180.<<" deg");
       const double cellSize = stats.squareCellSize(); // in arcsec
       ASKAPLOG_INFO_STR(logger, "  Estimated square cell size: "<<cellSize<<" arcsec");

       if (itsTangentDefined) {
           ASKAPLOG_INFO_STR(logger, "  Distance of the 'average' pointing direction from the tangent point: "<<itsTangent.separation(stats.centre())*180./casacore::C::pi<<" deg");
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
