/// @file
///
/// Support for parallel statistics accumulation to advise on imaging parameters
///
/// @copyright (c) 2016 CSIRO
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
/// @author Stephen Ord <stephen.ord@csiro.au>
///

#include <askap/distributedimager/AdviseDI.h>
#include <askap/parallel/ImagerParallel.h>
#include <askap/askap/AskapError.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapUtil.h>

ASKAP_LOGGER(logger, ".adviseDI");

#include <askap/profile/AskapProfiler.h>


#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/Solver.h>

#include <casacore/casa/BasicSL.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>
#include <casacore/ms/MSOper/MSReader.h>
#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/iostream.h>
#include <casacore/casa/namespace.h>
#include <casacore/casa/Quanta/MVTime.h>

#include <Blob/BlobString.h>
#include <Blob/BlobIBufString.h>
#include <Blob/BlobOBufString.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobOStream.h>


#include <vector>
#include <string>
#include <float.h>
#include <math.h>
#include <boost/shared_ptr.hpp>


namespace askap {

namespace synthesis {

using utility::toString;

float frequency_tolerance = 0.0;

bool compare_tol (const casacore::MFrequency& X, const casacore::MFrequency& Y) {

    if (abs(X.getValue() - Y.getValue()) <= frequency_tolerance) {
      return true;
    }
    else {
      return false;
    }
}

bool compare (const casacore::MFrequency& X, const casacore::MFrequency& Y) {
    return (X.getValue() == Y.getValue());
}

bool custom_lessthan (const casacore::MFrequency& X, const casacore::MFrequency& Y) {
    return (X.getValue() < Y.getValue());
}
bool in_range(double end1, double end2, double test) {


    if (test >= end1 && test <= end2) {
        ASKAPLOG_DEBUG_STR(logger,"Test frequency " << test \
        << " between " << end1 << " and " << end2);
        return true;
    }
    if (test <= end1 && test >= end2) {
        ASKAPLOG_DEBUG_STR(logger,"Test frequency " << test \
        << " between " << end1 << " and " << end2);
        return true;
    }
    ASKAPLOG_DEBUG_STR(logger,"Test frequency " << test \
    << " NOT between " << end1 << " and " << end2);
    return false;
}
// actual AdviseDI implementation

/// @brief Constructor from ParameterSet
/// @details The parset is used to construct the internal state. We could
/// also support construction from a python dictionary (for example).
/// This is needed becuase the default AdviseParallel assumes a master/worker
/// distribution that may not be the case.

/// The command line inputs are needed solely for MPI - currently no
/// application specific information is passed on the command line.
/// @param comms communication object
/// @param parset ParameterSet for inputs
AdviseDI::AdviseDI(askap::cp::CubeComms& comms, LOFAR::ParameterSet& parset) :
    AdviseParallel(comms,parset),itsParset(parset)
{
    isPrepared = false;
    itsFreqRefFrame = casacore::MFrequency::Ref(casacore::MFrequency::TOPO);
    itsWorkUnitCount=0;
}

/// @brief obtain metadata stats for one dataset
/// @details Hopefully this is a temporary method until this class (or the whole imager) can be
/// redesigned, at least to avoid accessing measurement set from first principles. Currently, this
/// method adds to uglyness of the code - it uses the functionality of the base class (i.e. the originally
/// written estimator) to get the basic info for a single measurement set without relying on the correct
/// parallel distribution. The returned statistics is valid until the next call to this method or until
/// the full parallel advise is done.
/// @param[in] ms measurement set to work with
/// @return const reference to the object populated with resulted metadata statistics
const VisMetaDataStats& AdviseDI::computeVisMetaDataStats(const std::string &ms)
{
   // a very basic estimator with all defaults - we currently use it only to get tangent point
   boost::shared_ptr<VisMetaDataStats> stats(new VisMetaDataStats());
   // work with this estimator until the next call to this method or until the proper parallel procedure is performed
   // via estimate() - note, it seems to be tied down to the original design of cimager and doesn't map on what the
   // current imager is doing / its distribution pattern
   setCustomEstimator(stats);
   // now analyse given measurement set
   calcOne(ms);
   // in principle we could've returned stats as it would be the same, but it is cleaner this way, shared pointer is kept in the base class
   return estimator();
}

void AdviseDI::prepare() {
    // this assumes only a sinlge spectral window - must generalise
    ASKAPLOG_INFO_STR(logger,"Running prepare");
    // Read from the configruation the list of datasets to process
    const vector<string> ms = getDatasets();
    ASKAPLOG_INFO_STR(logger,"Data set list is " << ms);
    const unsigned int nWorkers = itsComms.nProcs() - 1;
    ASKAPLOG_DEBUG_STR(logger, "nWorkers " << nWorkers);
    ASKAPCHECK(nWorkers > 0, "This code is intended to be executed in parallel with at least two ranks available");

    const unsigned int nGroups = itsComms.nGroups();
    ASKAPLOG_DEBUG_STR(logger, "nGroups " << nGroups);
    ASKAPCHECK(nGroups > 0, "Expect at least one rank group to be available");

    const unsigned int nWorkersPerGroup = nWorkers/nGroups;
    ASKAPLOG_DEBUG_STR(logger, "nWorkersPerGroup " << nWorkersPerGroup);
    ASKAPCHECK(nWorkersPerGroup > 0, "The number of workers per group is expected to be positive, i.e. you get less workers than groups");

    int nchanpercore = itsParset.getInt32("nchanpercore", 1);
    ASKAPLOG_DEBUG_STR(logger,"nchanpercore " << nchanpercore);

    const int nwriters = itsParset.getInt32("nwriters", 1);
    ASKAPLOG_DEBUG_STR(logger,"nwriters " << nwriters);
    frequency_tolerance = itsParset.getDouble("channeltolerance",0.0);

    ASKAPCHECK(nwriters > 0 ,"Number of writers must be greater than zero");
    ASKAPCHECK(nwriters <= nWorkers ,"Number of writers must less than or equal to number of workers");

    /// Get the channel range
    /// The imager ususally uses the Channels keyword in the parset to
    /// defined the work unit that it will attempt. That does not work in this



    casacore::uInt srow = 0;
    chanFreq.resize(ms.size());
    chanWidth.resize(ms.size());
    effectiveBW.resize(ms.size());
    resolution.resize(ms.size());
    centre.resize(ms.size());

    // Not really sure what to do for multiple ms or what that means in this
    // context... but i'm doing it any way - probably laying a trap for myself

    // Iterate over all measurement sets
    // combine all the channels into a list...
    // these measurement sets may now be from different epochs they should not
    // have different channel ranges - but it is possible that the channel range
    // may have been broken up into chunks.

    // need to calculate the allocations too.

    casacore::uInt totChanIn = 0;
    ASKAPLOG_INFO_STR(logger,"Testing for Channels");
    vector<int> itsChannels = getChannels();
    ASKAPLOG_INFO_STR(logger,"Testing for Frequencies");
    vector<double> itsFrequencies = getFrequencies();

    bool user_defined_channels = false;
    bool user_defined_frequencies = false;

    if (itsChannels[2] > 0) {
        user_defined_channels = true;
        ASKAPLOG_INFO_STR(logger,"User specified channels");
    }
    else {
        ASKAPLOG_INFO_STR(logger,"User has not specified channels");
    }
    if (itsFrequencies[2] != 0) {
        user_defined_frequencies = true;
        ASKAPLOG_INFO_STR(logger,"User specified frequencies explicitly");
    }
    else {
        ASKAPLOG_INFO_STR(logger,"User has not specified frequencies explicitly");
    }
    if (user_defined_channels && user_defined_frequencies) {
        ASKAPLOG_WARN_STR(logger,
            "User has specified both Channels AND Frequencies - Frequency range will take preference");
        user_defined_channels = false;
    }

    for (unsigned int n = 0; n < ms.size(); ++n) {
        chanFreq[n].resize(0);
        chanWidth[n].resize(0);
        effectiveBW[n].resize(0);
        resolution[n].resize(0);
        centre[n].resize(0);

        // MV: the design of this class is very ugly from the C++ point of view, it needs to be redesigned to get more
        // structure if we want to extend it further (or even debug - I suspect the issues I am working with now is
        // just a tip of an iceberg). Also, it needs to access MSs through standard interfaces and mimic the same access patterns
        // as the actual imager (or do appropriate estimates). To address the immediate problem of too simplistic interpretation of
        // MSs, I simply added another iteration through data accessor interface provided by base class (and btw, the association rather than
        // inheritance may be more appropriate here, with possible refactoring). This essentially would introduce another iteration over
        // metadata. But, in this simplistic form, when only tangent point is required, data won't be touched, so hopefully no huge performance
        // penalty. Besides, this additional code is well encapsulated keeping (additional) technical debt to a minimum.
        ASKAPLOG_DEBUG_STR(logger, "Assessing " << ms[n] << " via the standard interfaces");
        const VisMetaDataStats &mdStats = computeVisMetaDataStats(ms[n]);
        //

        // Open the input measurement set
        ASKAPLOG_DEBUG_STR(logger, "Opening " << ms[n] << " filecount " << n );
        const casacore::MeasurementSet in(ms[n]);
        const casacore::ROMSColumns srcCols(in);
        const casacore::ROMSSpWindowColumns& sc = srcCols.spectralWindow();
        const casacore::ROMSFieldColumns& fc = srcCols.field();
        const casacore::ROMSObservationColumns& oc = srcCols.observation();
        const casacore::ROMSAntennaColumns& ac = srcCols.antenna();
        const casacore::ROArrayColumn<casacore::Double> times = casacore::ROArrayColumn<casacore::Double>(oc.timeRange());
        const casacore::ROArrayColumn<casacore::Double> ants = casacore::ROArrayColumn<casacore::Double>(ac.position());
        const casacore::uInt thisRef = casacore::ROScalarColumn<casacore::Int>(in.spectralWindow(),"MEAS_FREQ_REF")(0);

        casacore::uInt thisChanIn = 0;
        srow = sc.nrow()-1;

        ASKAPCHECK(srow==0,"More than one spectral window not currently supported in adviseDI");

        // get the channel selection - first append all the input channels

        size_t chanStart=0;
        size_t chanStop=0;
        size_t chanStep=0;


        chanStart = 0;
        chanStop = casacore::ROScalarColumn<casacore::Int>(in.spectralWindow(),"NUM_CHAN")(0);
        chanStep = 1;
        thisChanIn = 0;

        for (uint i = chanStart; i < chanStop; i = i + chanStep) {
            chanFreq[n].push_back(sc.chanFreq()(srow)(casacore::IPosition(1, i)));
            chanWidth[n].push_back(sc.chanWidth()(srow)(casacore::IPosition(1, i)));
            effectiveBW[n].push_back(sc.effectiveBW()(srow)(casacore::IPosition(1, i)));
            resolution[n].push_back(sc.resolution()(srow)(casacore::IPosition(1, i)));
            thisChanIn++;
        }

        totChanIn = totChanIn + thisChanIn;

        
        // MV: see comments above, this is a somewhat ugly approach to get the correct phase centre
        itsTangent.push_back(mdStats.centre());
        itsDirVec.push_back(casa::Vector<casacore::MDirection>(1,casacore::MDirection(itsTangent[n], casacore::MDirection::J2000)));
        const casacore::Vector<casacore::MDirection> oldDirVec(fc.phaseDirMeasCol()(0));
        ASKAPLOG_DEBUG_STR(logger, "Tangent point for "<<ms[n]<<" : "<<printDirection(itsTangent[n])<<
                           " (J2000), old way: "<<printDirection(oldDirVec(0).getValue())<<" (J2000)");
        

        // the old way to get phase centre/tangent:
        //itsDirVec.push_back(fc.phaseDirMeasCol()(0));
        //itsTangent.push_back(itsDirVec[n](0).getValue());

        // Read the position on Antenna 0
        Array<casacore::Double> posval;
        ants.get(0,posval,true);
        vector<double> pval = posval.tovector();

        MVPosition mvobs(Quantity(pval[0], "m").getBaseValue(),
        Quantity(pval[1], "m").getBaseValue(),
        Quantity(pval[2], "m").getBaseValue());

        itsPosition.push_back(MPosition(mvobs,casacore::MPosition::ITRF));

        // Get the Epoch
        Array<casacore::Double> tval;
        vector<double> tvals;

        times.get(0,tval,true);
        tvals = tval.tovector();
        double mjd = tvals[0]/(86400.);
        casacore::MVTime dat(mjd);

        itsEpoch.push_back(MVEpoch(dat.day()));

        itsRef = thisRef;

        ASKAPLOG_INFO_STR(logger, "Completed filecount " << n);
    }


    // ASKAPLOG_INFO_STR(logger, "Assuming tangent point shared: "<<printDirection(itsTangent[0])<<" (J2000)");

    // the frequencies
    itsFFrameFrequencies.resize(0);
    itsInputFrequencies.resize(0);
    itsRequestedFrequencies.resize(0);
    // the work allocations
    itsAllocatedFrequencies.resize(nWorkersPerGroup);
    itsAllocatedWork.resize(nWorkers);

    // setup frequency frame
    const Bool bc = itsParset.getBool("barycentre",false);

    std::string freqFrame = itsParset.getString("freqframe","topo");

    if (bc == true) { // legacy
        freqFrame = "bary";
    }

    if (freqFrame == "topo") {
        ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as topocentric");
        itsFreqRefFrame = casacore::MFrequency::Ref(casacore::MFrequency::TOPO);
        itsFreqType = casacore::MFrequency::TOPO;
    } else if (freqFrame == "lsrk") {
        ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as lsrk");
        itsFreqRefFrame = casacore::MFrequency::Ref(casacore::MFrequency::LSRK);
        itsFreqType = casacore::MFrequency::LSRK;
    } else if (freqFrame == "bary") {
        ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as barycentric");
        itsFreqRefFrame = casacore::MFrequency::Ref(casacore::MFrequency::BARY);
        itsFreqType = casacore::MFrequency::BARY;
    } else {
        ASKAPTHROW(AskapError, "Unsupported frequency frame "<<freqFrame);
    }


    // At this point we now have each channel from each MS
    // in a unique array.
    // It's freq frame is itsRef - so in principle we support any input reference frame ....
    // first we need to sort and uniqify the list
    // then resize the list to get the channel range.
    // This is required becuase we are trying to form a unique
    // reference channel list from the input measurement sets

    // Then we have to use the user preferences to get the correct list of desired frequencies




    for (unsigned int n = 0; n < ms.size(); ++n) {



        // builds a list of all the channels

        for (unsigned int ch = 0; ch < chanFreq[n].size(); ++ch) {

            MeasFrame itsFrame(MEpoch(itsEpoch[n]),itsPosition[n],itsDirVec[n][0]);
            MFrequency::Ref refin(MFrequency::castType(itsRef),itsFrame); // the frame of the input channels



            itsInputFrequencies.push_back(MFrequency(MVFrequency(chanFreq[n][ch]),refin));

            /// The original scheme attempted to convert the input into the output frame
            /// and only keep those output (FFRAME) channels that matched.
            /// This was not deemed useful or helpful enough though. But I need to keep that
            /// mode as a default.

        }
    }

    /// uniquifying the lists

    bool (*custom_compare)(const casacore::MFrequency& , const casacore::MFrequency& ) = NULL;

    if (frequency_tolerance > 0) {
      custom_compare = compare_tol;
      ASKAPLOG_WARN_STR(logger,"Comparing frequencies with floating point tolerance of " << frequency_tolerance);
    }
    else {
      custom_compare = compare;
      ASKAPLOG_INFO_STR(logger,"Using standard compare (zero tolerance) for frequency allocations");
    }


    std::sort(itsInputFrequencies.begin(),itsInputFrequencies.end(), custom_lessthan);
    std::vector<casacore::MFrequency>::iterator topo_it;
    topo_it = std::unique(itsInputFrequencies.begin(),itsInputFrequencies.end(),custom_compare);
    itsInputFrequencies.resize(std::distance(itsInputFrequencies.begin(),topo_it));
    ASKAPLOG_DEBUG_STR(logger," Unique sizes Input " << itsInputFrequencies.size() << " Output " << itsFFrameFrequencies.size());


    // Now they are unique we need to get a list of desired output freqencies that meet
    // the requirements specified in the parset

    if (user_defined_channels) {
        // the user has specified nchan and a start channel
        // this is probably based upon the input frame as the user has
        // <probably> not done the maths to work out the output channel mapping so.
        // The channel width is unchanged - I may allow the channel width to be used
        // but that is probably another ticket.

        size_t n = itsChannels[0];
        size_t st = itsChannels[1];

        if (n > itsInputFrequencies.size()){
            ASKAPLOG_WARN_STR(logger, "Requested nchan > available channels; truncating");
        }

        for (unsigned int ch = 0; ch < itsInputFrequencies.size(); ++ch) {
            if (ch >= st) {
                if (itsRequestedFrequencies.size() < n)
                    itsRequestedFrequencies.push_back(itsInputFrequencies[ch]);
            }
        }
    }
    else if (user_defined_frequencies) {
        // This time the user has specified frequencies in the freqFrame frame
        // So we now fill the desired frequency array with that in mind.
        // Easy
        ASKAPLOG_WARN_STR(logger, "User requested frequency range is being used");
        size_t n=itsFrequencies[0];
        double width = itsFrequencies[2];
        double st = itsFrequencies[1] + width/2.0;
        ASKAPLOG_WARN_STR(logger, "Starting at " << st << " width " << width);
        for (unsigned int ch = 0; ch < n ; ch++) {
            itsRequestedFrequencies.push_back(MFrequency(Quantity(st+ch*width,"Hz"),itsFreqRefFrame));
        }

    }
    else {
        ASKAPLOG_WARN_STR(logger, "Full channel range is being used");
        for (unsigned int ch = 0; ch < itsInputFrequencies.size(); ++ch) {
           itsRequestedFrequencies.push_back(itsInputFrequencies[ch]);
        }
    }
    // Now we have a list of requested frequencies lets allocate them to nodes - some maybe empty.
    ASKAPLOG_INFO_STR(logger,
    " User requests " << itsRequestedFrequencies.size() << " cube " << " starting at " << itsRequestedFrequencies[0].getValue());

    if (itsRequestedFrequencies.size()/nWorkersPerGroup != nchanpercore) {
        ASKAPLOG_WARN_STR(logger,"Mismatch: #requested frequencies = "<<itsRequestedFrequencies.size()<<
        ", #workers = "<<nWorkersPerGroup<<", #chan per core = "<<nchanpercore);
        nchanpercore = itsRequestedFrequencies.size()/nWorkersPerGroup;
        if (nchanpercore == 0) {
            nchanpercore = 1;
        }
        ASKAPLOG_WARN_STR(logger,"Trying nchanpercore = "<<nchanpercore);
        const string param = "nchanpercore";
        const string pstr = toString(nchanpercore);
        ASKAPLOG_INFO_STR(logger, "  updating parameter " << param << ": " << pstr);
        if (itsParset.isDefined(param)) {
            itsParset.replace(param, pstr);
        } else {
            itsParset.add(param, pstr);
        }
    }
    ASKAPCHECK(itsRequestedFrequencies.size()/nWorkersPerGroup == nchanpercore,"Miss-match nchanpercore is incorrect");

    for (unsigned int ch = 0; ch < itsRequestedFrequencies.size(); ++ch) {

        ASKAPLOG_DEBUG_STR(logger,"Requested Channel " << ch << ":" << itsRequestedFrequencies[ch]);
        unsigned int allocation_index = floor(ch / nchanpercore);

        ASKAPLOG_DEBUG_STR(logger,"Allocating frequency "<< itsRequestedFrequencies[ch].getValue() \
        << " to worker " << allocation_index+1);

        itsAllocatedFrequencies[allocation_index].push_back(itsRequestedFrequencies[ch].getValue());
    }


    // Now for each allocated workunit we need to fill in the rest of the workunit
    // we now have a workUnit for each channel in the allocation - but not
    // for each Epoch.

    // We have to match the desired frequencies to those present in the data set.

    int globalChannel = 0;
    vector<int> itsBeams = getBeams();

    // initially lets just use the first beam in the list

    int myBeam = itsBeams[0];

    //
    for (unsigned int work = 0; work < itsAllocatedFrequencies.size(); ++work) {
        ASKAPLOG_DEBUG_STR(logger,"Allocating frequency channels for worker " << work);
        // loop over the measurement sets and find the local channel number
        // associated with the barycentric channel

        vector<double>& thisAllocation = itsAllocatedFrequencies[work];

        for (unsigned int frequency=0;frequency < thisAllocation.size();++frequency) {
            /// Global channels are in order so I can get the global channel count from here

            // need to allocate the measurement sets for this channel to this allocation
            // this may require appending new work units.
            ASKAPLOG_DEBUG_STR(logger,"Allocating " << thisAllocation[frequency] \
            << "Global channel " << globalChannel);

            bool allocated = false;
            for (unsigned int set=0;set < ms.size();++set){


                MeasFrame itsFrame(MEpoch(itsEpoch[set]),itsPosition[set],itsDirVec[set][0]);
                MFrequency::Ref refin(MFrequency::castType(itsRef),itsFrame); // the frame of the input channels
                MFrequency::Ref refout(itsFreqType,itsFrame); // the frame desired
                MFrequency::Convert forw(refin,refout); // from input to desired
                MFrequency::Convert backw(refout,refin); // from desired to input

                vector<int> lc;
                // try and match the converted frequency in the input data
                // now this returns all channels in a range.

                // Determine the range of the channels in the default case
                // we need to know the chanWidth
                // in the frequency case this is from itsFrequencies
                MVFrequency oneEdge;
                MVFrequency otherEdge;
                if (user_defined_frequencies) {
                    oneEdge = thisAllocation[frequency] - itsFrequencies[2]/2.;
                    otherEdge = thisAllocation[frequency] + itsFrequencies[2]/2.;
                }
                else {
                    oneEdge = thisAllocation[frequency] - chanWidth[set][0]/2.0;
                    otherEdge = thisAllocation[frequency] + chanWidth[set][0]/2.0;
                }

                // try and find the requested channels in the input dataset
                lc = matchall(set,backw(oneEdge).getValue(),backw(otherEdge).getValue());

                if (lc.size() > 0) {
                    size_t limit = lc.size();
                    if (lc.size() > 1) {
                        ASKAPLOG_WARN_STR(logger,"More that one channel has been found in this range. This is currently unsupported and will be ignored");
                        limit = 1;
                    }

                    // there is at least one channel of this frequency in the measurement set
                    for (size_t lc_part=0; lc_part < limit; lc_part++) {
                        cp::ContinuumWorkUnit wu;

                        wu.set_payloadType(cp::ContinuumWorkUnit::WORK);
                        wu.set_channelFrequency(thisAllocation[frequency]);
                        wu.set_beam(myBeam);

                        if (itsRequestedFrequencies.size() > 1)
                            wu.set_channelWidth(fabs(itsInputFrequencies[1].getValue() - itsInputFrequencies[0].getValue()));
                        else
                            wu.set_channelWidth(fabs(chanWidth[0][0]));

                        wu.set_localChannel(lc[lc_part]);
                        wu.set_globalChannel(globalChannel);
                        wu.set_dataset(ms[set]);
                        ASKAPDEBUGASSERT(work < itsAllocatedWork.size());
                        itsAllocatedWork[work].push_back(wu);
                        itsWorkUnitCount++;
                        ASKAPLOG_DEBUG_STR(logger,"MATCH Found desired freq " << thisAllocation[frequency] \
                        << " in local channel number " << lc << " ( " << chanFreq[set][lc[lc_part]] << " ) of width " << wu.get_channelWidth()  \
                        << " in set: " << ms[set] <<  " to rank " << work+1 << " this rank has " \
                        << itsAllocatedWork[work].size() << " of a total count " << itsWorkUnitCount \
                        << " the global channel is " << globalChannel);

                        allocated = true;
                    }
                }

            }
            globalChannel++;
            if (allocated == false)
            {
                ASKAPLOG_WARN_STR(logger,"Allocating FAIL Cannot match " << thisAllocation[frequency] << " for rank " << work+1 \
                << " in any set: ");
                // warn it does not match ....
                // have to increment the workcount for the cleanup.
                cp::ContinuumWorkUnit wu;
                wu.set_payloadType(cp::ContinuumWorkUnit::NA);
                wu.set_channelFrequency(thisAllocation[frequency]);
                wu.set_beam(myBeam);

                if (itsRequestedFrequencies.size() > 1)
                    wu.set_channelWidth(fabs(itsInputFrequencies[1].getValue() - itsInputFrequencies[0].getValue()));
                else
                    wu.set_channelWidth(fabs(chanWidth[0][0]));

                wu.set_localChannel(-1);
                wu.set_globalChannel(-1);

                ASKAPDEBUGASSERT(work < itsAllocatedWork.size());
                itsAllocatedWork[work].push_back(wu);
                itsWorkUnitCount++;

            }

        }

    }
        // expand the channels by the number of groups - this is cheap on memory and
        // allows easier indexing
        // But this is only really needed by the master
        /// Now if required we need to allocate the writers for a parallel writers
        /// The writers do not need to be dedicated cores - they can write in addition
        /// to their other duties.

    // This loop is trying to find a writer with work
    unsigned int nWorkersPerWriter = floor(itsAllocatedWork.size() / nwriters);
    int mywriter = 0;
    for (int wrk = 0; wrk < itsAllocatedWork.size(); wrk++) {
        if (nWorkersPerWriter>0) {
            mywriter = floor(wrk/nWorkersPerWriter)*nWorkersPerWriter;
        }
        bool has_work = false;
        while (has_work == false) {
            for (int unit = 0; unit < itsAllocatedWork[mywriter].size(); unit++) {
                if (itsAllocatedWork[mywriter][unit].get_payloadType() == cp::ContinuumWorkUnit::WORK) {
                    has_work = true;
                    break;
                }
            }
            if (has_work == true) {
                break;
            }
            else {
                mywriter++;
                ASKAPLOG_WARN_STR(logger,"Ran out of eligible writers will write myself");
                mywriter = wrk;
                break;

            }
        }

        for (int unit = 0; unit < itsAllocatedWork[wrk].size(); unit++) {
                    itsAllocatedWork[wrk][unit].set_writer(mywriter+1); // plus 1 for rank
                    ASKAPLOG_DEBUG_STR(logger,"Set rank " << wrk+1 << " writer to be rank " << mywriter+1);
        }

    }
    for (int grp = 1; grp < itsComms.nGroups(); grp++) {
        for (int wrk = 0; wrk < nWorkersPerGroup; wrk++) {
            ASKAPDEBUGASSERT(grp*nWorkersPerGroup + wrk < itsAllocatedWork.size());
            itsAllocatedWork[grp*nWorkersPerGroup+wrk] = itsAllocatedWork[wrk];

            itsWorkUnitCount=itsWorkUnitCount + itsAllocatedWork[wrk].size();

            ASKAPLOG_DEBUG_STR(logger,"Allocating rank " << grp*nWorkersPerGroup+wrk+1 \
            << " the same units as rank " << wrk+1 << "(" << itsAllocatedWork[wrk].size() << ")"<< " Count " << itsWorkUnitCount);
        }
    }

    isPrepared = true;
    ASKAPLOG_DEBUG_STR(logger, "Prepared the advice");
}

cp::ContinuumWorkUnit AdviseDI::getAllocation(int id) {
    ASKAPCHECK(id >= 0 && id < itsAllocatedWork.size(), "Requested work id = "<<id<<" is outside of bounds, itsAllocatedWork.size() = "<<itsAllocatedWork.size());
    cp::ContinuumWorkUnit rtn;
    if (itsAllocatedWork[id].empty() == true) {
        ASKAPLOG_WARN_STR(logger, "Stack is empty for " << id+1);
        rtn.set_payloadType(cp::ContinuumWorkUnit::DONE);
        return rtn;
    }
    else {
        rtn = itsAllocatedWork[id].back();
        itsAllocatedWork[id].pop_back();

    }
    int count=0;
    for (int alloca = 0 ; alloca < itsAllocatedWork.size() ; alloca++) {
        count = count + itsAllocatedWork[alloca].size();

    }
    itsWorkUnitCount = count;
    return rtn;
}
vector<int> AdviseDI::matchall(int ms_number,
casacore::MVFrequency oneEdge, casacore::MVFrequency otherEdge) const {
    /// return all the input channels in the range
    vector<int> matches;
    for (int ch=0 ; ch < chanFreq[ms_number].size(); ++ch) {
            ASKAPLOG_DEBUG_STR(logger, "looking for " << chanFreq[ms_number][ch] << "Hz");
            if (in_range(oneEdge.getValue(),otherEdge.getValue(),chanFreq[ms_number][ch])) {
                ASKAPLOG_DEBUG_STR(logger, "Found");
                matches.push_back(ch);
            }
    }

    return matches;
}

void AdviseDI::addMissingParameters() {
    this->addMissingParameters(this->itsParset,True);
}

void AdviseDI::updateDirectionFromWorkUnit(LOFAR::ParameterSet& parset, askap::cp::ContinuumWorkUnit& wu) {

  string wu_dataset = wu.get_dataset();
  std::vector<std::string> ms = getDatasets();
  const vector<string> imageNames = parset.getStringVector("Images.Names", false);
  ASKAPLOG_DEBUG_STR(logger,"Image names " << imageNames);

  for (unsigned int n = 0; n < ms.size(); ++n) {

    if (wu_dataset.compare(ms[n]) == 0) {
      // MATCH
      string param = "Images.direction";
      string pstr = "["+printLon(itsTangent[n])+","+printLat(itsTangent[n])+", J2000]";
      // Only J2000 is implemented at the moment.
      ASKAPLOG_INFO_STR(logger, "  updating parameter " << param << ": " << pstr);
      parset.replace(param, pstr);

      for (size_t img = 0; img < imageNames.size(); ++img) {
        param ="Images."+imageNames[img]+".direction";
        ASKAPLOG_INFO_STR(logger, "  updating parameter " << param << ": " << pstr);
        parset.replace(param, pstr);

      }
    }
  }
}

void AdviseDI::addMissingParameters(LOFAR::ParameterSet& parset, bool extra)
{

    ASKAPLOG_DEBUG_STR(logger,"Adding missing params ");

    if (isPrepared == true) {
        ASKAPLOG_DEBUG_STR(logger,"Prepared therefore can add frequency label for the output image");
        std::vector<casacore::MFrequency>::const_iterator begin_it;
        std::vector<casacore::MFrequency>::const_iterator end_it;

        ASKAPCHECK(!itsRequestedFrequencies.empty(),"No frequencies requested for ouput");
        begin_it = itsRequestedFrequencies.begin();
        end_it = itsRequestedFrequencies.end()-1;


        this->minFrequency = (*begin_it).getValue();
        this->maxFrequency = (*end_it).getValue();
        ASKAPLOG_DEBUG_STR(logger,"Min:Max frequency -- " << this->minFrequency << ":" << this->maxFrequency);


    }

   // test for missing image-specific parameters:

   // these parameters can be set globally or individually
   bool cellsizeNeeded = false;
   bool shapeNeeded = false;
   int nTerms = 1;

   string param;


   const vector<string> imageNames = parset.getStringVector("Images.Names", false);

   param = "Images.direction";

   if ( !parset.isDefined(param) ) {
       // Only J2000 is implemented at the moment.
       string pstr = "["+printLon(itsTangent[0])+","+printLat(itsTangent[0])+", J2000]";
       ASKAPLOG_INFO_STR(logger, "  Advising on parameter (getting from the tangent point on the first measurement set) " << param << ": " << pstr);
       itsParset.add(param, pstr);
   }
   param = "Images.restFrequency";

   if ( !parset.isDefined(param) ) {
       string pstr = "HI";
       ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << ": " << pstr);
       parset.add(param, pstr);
   }

   for (size_t img = 0; img < imageNames.size(); ++img) {

     if (!parset.isDefined("Images."+imageNames[img]+".cellsize")) {
         cellsizeNeeded = true;
     }
     if (!parset.isDefined("Images."+imageNames[img]+".shape")) {
         shapeNeeded = true;
     }

     param = "Images."+imageNames[img]+".frequency";
     if ( !parset.isDefined(param) && isPrepared == true) {
       const string key="Images."+imageNames[img]+".frequency";
       // changing this to match adviseParallel
       const double aveFreq = 0.5*(minFrequency+maxFrequency);
       const string val = "["+toString(aveFreq)+","+toString(aveFreq)+"]";
       ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << val);
       parset.add(key,val);
     }

     param ="Images."+imageNames[img]+".direction";
     if ( !parset.isDefined(param) ) {

       if (parset.isDefined("Images.direction") ) {
         const std::vector<std::string> direction = parset.getStringVector("Images.direction");
         const double ra = SynthesisParamsHelper::convertQuantity(direction[0],"rad");
         const double dec = SynthesisParamsHelper::convertQuantity(direction[1],"rad");
         const casacore::MVDirection itsDirection = casacore::MVDirection(ra,dec);

         // Only J2000 is implemented at the moment.
         string pstr = "["+printLon(itsDirection)+","+printLat(itsDirection)+", J2000]";
         ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << " (obtained from Images.direction) : " << pstr);
         itsParset.add(param, pstr);
       } else {
         // Only J2000 is implemented at the moment.
         string pstr = "["+printLon(itsTangent[0])+","+printLat(itsTangent[0])+", J2000]";
         ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << "(obtained from tangent point of first MS): " << pstr);
         parset.add(param, pstr);
       }
     }
     param = "Images."+imageNames[img]+".nterms"; // if nterms is set, store it for later
     if (parset.isDefined(param)) {
       if ((nTerms>1) && (nTerms!=itsParset.getInt(param))) {
         ASKAPLOG_WARN_STR(logger, "  Imaging with different nterms may not work");
       }
       nTerms = itsParset.getInt(param);
     }

     param = "Images."+imageNames[img]+".nchan";
     if ( !parset.isDefined("Images."+imageNames[img]+".nchan") ) {
         string pstr = "1";  // Default nchan is 1
         ASKAPLOG_INFO_STR(logger, "  Advising on parameter nchan : " << pstr);
         parset.add(param, pstr);

     }
   }

   if (nTerms > 1) { // check required MFS parameters
       param = "visweights"; // set to "MFS" if unset and nTerms > 1
       if (!parset.isDefined(param)) {
           string pstr="MFS";
           ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<" (obtained by default - we know no other weighting) : " << pstr);
           parset.add(param, pstr);
       }

       param = "visweights.MFS.reffreq"; // set to average frequency if unset and nTerms > 1
       if ((parset.getString("visweights")=="MFS")) {
           if (!parset.isDefined(param)) {
               const double aveFreq = 0.5*(minFrequency+maxFrequency);
               string val = toString(aveFreq);
               ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<" (using average frequency):  " << val);
               parset.add(param,val);
           }

       }
   }

   if (!extra) {
       return;
   }
   // test for general missing parameters:
   // if any of these is missing we run the stats
   const string wMaxGridder = ImagerParallel::wMaxAdviceNeeded(parset);
   bool missing = !parset.isDefined("nUVWMachines");
   missing |= (cellsizeNeeded && !parset.isDefined("Images.cellsize"));
   missing |= (shapeNeeded && !parset.isDefined("Images.shape"));
   missing |= (wMaxGridder!="");

   if (missing) {
       ASKAPLOG_INFO_STR(logger, "Remaining missing parameters require a pass through the data");

       if (!parset.isDefined("tangent")) {
           string pstr = "["+printLon(itsTangent[0])+","+printLat(itsTangent[0])+", J2000]";
           ASKAPLOG_INFO_STR(logger, "  Adding tangent for advise: " << pstr);
           parset.add("tangent", pstr);
           parset.add("tangent.advised", pstr);
       }

       // Use the snapshotimaging wtolerance as the advise wtolerance. But only if wmax is needed.
       bool snapShot = parset.getBool("gridder.snapshotimaging",false);
       param = "gridder.snapshotimaging.wtolerance";
       if (snapShot && parset.isDefined(param) && !parset.isDefined("wtolerance") && (wMaxGridder!="")) {
           string wtolerance = parset.getString(param);
           ASKAPLOG_INFO_STR(logger, "  Adding wtolerance for advise: " << wtolerance);
           parset.add("wtolerance", wtolerance);
           parset.add("wtolerance.advised", wtolerance);
       }
       if (wMaxGridder!="") {
           param = "gridder."+wMaxGridder+".wpercentile";
           if (parset.isDefined(param)) {
               string wpercentile = parset.get(param);
               parset.add("wpercentile", wpercentile);
               parset.add("wpercentile.advised", wpercentile);
           }
       }

       estimate(parset);
       const VisMetaDataStats &advice = estimator();

       param = "nUVWMachines"; // if the number of uvw machines is undefined, set it to the number of beams.
       if (!parset.isDefined(param)) {
           const string pstr = toString(advice.nBeams());
           ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr);
           parset.add(param, pstr);
       }

       param = "Images.cellsize";
       std::vector<double> cellSize(2, advice.squareCellSize());
       if (parset.isDefined(param)) {
           cellSize = SynthesisParamsHelper::convertQuantity(parset.getStringVector(param),"arcsec");
       } else if (cellsizeNeeded) {
           string pstr = "["+toString(cellSize[0])+"arcsec,"+toString(cellSize[1])+"arcsec]";
           ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr);
           parset.add(param, pstr);
       }

       param = "Images.shape"; // if image shape is undefined, use the advice.
       if (shapeNeeded && !parset.isDefined(param)) {
           const double fieldSize = advice.squareFieldSize(1); // in deg
           const long lSize = long(fieldSize * 3600 / cellSize[0]) + 1;
           const long mSize = long(fieldSize * 3600 / cellSize[1]) + 1;
           string pstr = "["+toString(lSize)+","+toString(mSize)+"]";
           ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr);
           parset.add(param, pstr);
       }

       const string gridder = ImagerParallel::wMaxAdviceNeeded(parset); // returns empty string if wmax is not required.
       if (gridder!="") {
           param = "gridder."+gridder+".wmax"; // if wmax is undefined but needed, use the advice.
           if (!parset.isDefined(param)) {
               // both should be set if either is, but include the latter to ensure that maxResidualW is generated.
               if (parset.isDefined("gridder.snapshotimaging.wtolerance") && parset.isDefined("wtolerance") ) {
                   const string pstr = toString(advice.maxResidualW()); // could use parset.getString("gridder.snapshotimaging.wtolerance");
                   ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr);
                   parset.add(param, pstr);
               } else if (parset.isDefined("gridder."+gridder+".wpercentile")){
                   const string pstr = toString(advice.wPercentile());
                   ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr << " with specified percentile value");
                   parset.add(param, pstr);
                   // if we're setting from percentile, turn on clipping
                   param = "gridder."+gridder+".wmaxclip";
                   if (!parset.isDefined(param)) {
                       parset.add(param, "true");
                   }
               } else {
                   const string pstr = toString(advice.maxW());
                   ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr);
                   parset.add(param, pstr);
               }
           }
       }

       ASKAPCHECK(!parset.isDefined("Images.extraoversampling"), "Images.extraoversampling cannot be set by user");
       if (parset.getBool("Images.nyquistgridding",false) || parset.isDefined("Images.griddingcellsize")) {

           // need to make sure that extraOsFactor results in an integer number of pixels,
           // which could get complicated for rectangular grids and pixels.
           ASKAPCHECK(cellSize[0]==cellSize[1], "nyquistgridding only set up for square pixels");

           const std::vector<int> imSize = parset.getInt32Vector("Images.shape");
           std::vector<double> gCellSize(2);
           if (parset.isDefined("Images.griddingcellsize")) {
               gCellSize =
                   SynthesisParamsHelper::convertQuantity(parset.getStringVector("Images.griddingcellsize"),"arcsec");
               ASKAPCHECK(gCellSize[0]==gCellSize[1], "nyquistgridding only set up for square pixels");
               ASKAPCHECK(gCellSize[0]>=cellSize[0], "griddingcellsize must not be less than cellsize");
           }
           else {
               const double uv_max = max(advice.maxU(), advice.maxV());
               const double fov = cellSize[0] * imSize[0] / 3600. * M_PI / 180.0;
               const double wk_max = 6/fov + advice.maxW()*fov;
               ASKAPASSERT(uv_max > 0);
               // calculate the resolution in arcsec corresponding to the smallest grid that will fit all the data
               //  - the reciprical of twice the longest baseline plus the largest w support
               gCellSize[0] = 0.5 / (uv_max + wk_max) * 3600.0 * 180.0 / M_PI;
               gCellSize[1] = gCellSize[0];
           }

           // nominal ratio between gridding resolution and cleaning resolution
           double extraOsFactor = gCellSize[0]/cellSize[0];
           // now tweak the ratio to result in an integer number of pixels and reset the gridding cell size
           double nPix = ceil(double(imSize[0])/extraOsFactor);
           // also ensure that it is even
           nPix += int(nPix) % 2;
           // reset the extra multiplicative factor
           extraOsFactor = double(imSize[0]) / nPix;

           gCellSize[0] = cellSize[0] * extraOsFactor;
           gCellSize[1] = gCellSize[0];

           ASKAPLOG_INFO_STR(logger, "  Adding new parameter extraoversampling = "<<extraOsFactor);
           ASKAPLOG_INFO_STR(logger, "  Changing cellsize from "<<parset.getStringVector("Images.cellsize")<<
                                     " to "<<"["<<gCellSize[0]<<"arcsec,"<<gCellSize[1]<<"arcsec]");
           ASKAPLOG_INFO_STR(logger, "  Changing shape from "<<parset.getInt32Vector("Images.shape")<<
                                     " to "<<"["<<long(nPix)<<","<<long(nPix)<<"]");

           {
               std::ostringstream pstr;
               pstr<<extraOsFactor;
               parset.add("Images.extraoversampling", pstr.str().c_str());
           }
           {
               std::ostringstream pstr;
               pstr<<"["<<gCellSize[0]<<"arcsec,"<<gCellSize[1]<<"arcsec]";
               parset.replace("Images.cellsize", pstr.str().c_str());
           }
           {
               std::ostringstream pstr;
               pstr<<"["<<long(nPix)<<","<<long(nPix)<<"]";
               parset.replace("Images.shape", pstr.str().c_str());
           }
       }

   }
   ASKAPLOG_DEBUG_STR(logger,"Done adding missing params ");

   if (parset.isDefined("tangent.advised")) {
       parset.remove("tangent");
       parset.remove("tangent.advised");
   }

   if (parset.isDefined("wtolerance.advised")) {
       parset.remove("wtolerance");
       parset.remove("wtolerance.advised");
   }

   if (parset.isDefined("wpercentile.advised")) {
       parset.remove("wpercentile");
       parset.remove("wpercentile.advised");
   }
}


// Utility function to get dataset names from parset.
std::vector<std::string> AdviseDI::getDatasets() const
{
    if (itsParset.isDefined("dataset") && itsParset.isDefined("dataset0")) {
        ASKAPTHROW(std::runtime_error,
            "Both dataset and dataset0 are specified in the parset");
    }

    // First look for "dataset" and if that does not exist try "dataset0"
    vector<string> ms;
    if (itsParset.isDefined("dataset")) {
        ms = itsParset.getStringVector("dataset", true);
    } else {
        string key = "dataset0";   // First key to look for
        long idx = 0;
        while (itsParset.isDefined(key)) {
            const string value = itsParset.getString(key);
            ms.push_back(value);
            key = "dataset"+toString(idx+1);
            ++idx;
        }
    }

    return ms;
}

/// the adviseDI should be smart enough to tell the difference between a straight list and
/// actual frequencies .... or is that too different.
std::vector<double> AdviseDI::getFrequencies() const {
    std::vector<double> f(3,0);

    if (!itsParset.isDefined("Frequencies")) {
    ASKAPLOG_WARN_STR(logger,
        "Frequencies keyword is not defined");

    }
    else {
        const std::vector<std::string> fstr = itsParset.getStringVector("Frequencies",true);
        ASKAPCHECK(fstr.size()>=3, "Require at least 3 elements for frequency");

        f[0] = utility::fromString<float>(fstr[0]);
        f[1] = utility::fromString<float>(fstr[1]);
        f[2] = utility::fromString<float>(fstr[2]);
         // just a start and stop - assume no averaging
        if (f[2] == 0) {
           ASKAPLOG_WARN_STR(logger,
           "Channel width not specified this setting will be ignored");
        }

    }
    return f;
}

std::vector<int> AdviseDI::getChannels() const {

    // channels should now behave more like the historical version
    // <nchan> <start> ...
    //

    std::vector<int> c(3,0);

    if (!itsParset.isDefined("Channels")) {
    ASKAPLOG_WARN_STR(logger,
        "Channels keyword is not defined");
        c[2] = -1;
    }
    else {
        const std::vector<std::string> cstr = itsParset.getStringVector("Channels",true);
        ASKAPCHECK(cstr.size() >= 2, "Require at least two elements for Channels keyword");
        c[0] = utility::fromString<float>(cstr[0]);
        const string wild = "%w";

        if (cstr[1].compare(wild) == 0) {
            ASKAPLOG_WARN_STR(logger,
        "Wild card in the Channel list");
            c[1] = 0;
            c[2] = -1; // this now images all the channels
        }
        else {
            c[1] = utility::fromString<int>(cstr[1]);
            c[2] = 1; // should maybe allow averaging
        }


    }

    // just a start and stop - assume no averaging
    if (c[2] == 0) {
        c[2] = 1;
    }
    // otherwise averaging needs to be performed -- see later

    return c;

}

void AdviseDI::updateComms() {

    cp::CubeComms& itsCubeComms = dynamic_cast< cp::CubeComms& >(itsComms);

    /// Go through the work allocations and set the writers
    for (int worker=0; worker < itsAllocatedWork.size() ; worker++) {
        itsCubeComms.addWorker(worker+1);
        double last_freq = itsAllocatedWork[worker][0].get_channelFrequency();
        for (int alloc=0; alloc < itsAllocatedWork[worker].size() ; alloc++) {
            double current_freq = itsAllocatedWork[worker][alloc].get_channelFrequency();
            // need to test whether this is a distinct channel or a different epoch for
            // the same channel
            if (itsAllocatedWork[worker][alloc].get_payloadType() == cp::ContinuumWorkUnit::NA) {
                ASKAPLOG_WARN_STR(logger,
                "NA Payload in updateComms");
            }
            else if (current_freq != last_freq || alloc == 0) {

                itsCubeComms.addWriter(itsAllocatedWork[worker][alloc].get_writer());

                itsCubeComms.addChannelToWriter(itsAllocatedWork[worker][alloc].get_writer(),worker+1);
                itsCubeComms.addChannelToWorker(worker+1);

                last_freq = current_freq;
            }
        }
    }
    /// First lets set up the cube
    const bool singleSink = itsParset.getBool("singleoutputfile",false);
    if (singleSink){
        itsCubeComms.setSingleSink();
    }
    else {
        itsCubeComms.setMultiSink();
    }


}
std::vector<int> AdviseDI::getBeams() const
{
    std::vector<int> bs;

    if (itsParset.isDefined("beams")) {
        bs = itsParset.getInt32Vector("beams",bs);

    }
    else {
        bs.push_back(0);
    }
    return bs;
}
double AdviseDI::getBaseFrequencyAllocation(int workerNumber) const {
    ASKAPDEBUGASSERT(workerNumber < itsAllocatedFrequencies.size());
    return itsAllocatedFrequencies[workerNumber][0];
}
} // namespace synthesis

} // namespace askap
