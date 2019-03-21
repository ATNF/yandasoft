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

#include <distributedimager/AdviseDI.h>
#include <askap/AskapError.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>
#include <dataaccess/SharedIter.h>
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>


ASKAP_LOGGER(logger, ".adviseDI");

#include <profile/AskapProfiler.h>


#include <fitting/INormalEquations.h>
#include <fitting/Solver.h>

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
#include <boost/shared_ptr.hpp>


namespace askap {

namespace synthesis {

float frequency_tolerance = 0.0;

bool compare_tol (const casa::MFrequency& X, const casa::MFrequency& Y) {

    if (abs(X.getValue() - Y.getValue()) <= frequency_tolerance) {
      return true;
    }
    else {
      return false;
    }
}

bool compare (const casa::MFrequency& X, const casa::MFrequency& Y) {
    return (X.getValue() == Y.getValue());
}

bool custom_lessthan (const casa::MFrequency& X, const casa::MFrequency& Y) {
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
    barycentre = false;
    itsFreqRefFrame = casa::MFrequency::Ref(casa::MFrequency::TOPO);
    itsWorkUnitCount=0;
}

void AdviseDI::prepare() {
    // this assumes only a sinlge spectral window - must generalise
    ASKAPLOG_INFO_STR(logger,"Running prepare");
    // Read from the configruation the list of datasets to process
    const vector<string> ms = getDatasets();
    ASKAPLOG_INFO_STR(logger,"Data set list is " << ms);
    unsigned int nWorkers = itsComms.nProcs() - 1;
    ASKAPLOG_DEBUG_STR(logger, "nWorkers " << nWorkers);

    unsigned int nWorkersPerGroup = nWorkers/itsComms.nGroups();
    ASKAPLOG_DEBUG_STR(logger, "nWorkersPerGroup " << nWorkersPerGroup);
    unsigned int nGroups = itsComms.nGroups();
    ASKAPLOG_DEBUG_STR(logger, "nGroups " << nGroups);
    const int nchanpercore = itsParset.getInt32("nchanpercore", 1);
    ASKAPLOG_DEBUG_STR(logger,"nchanpercore " << nchanpercore);
    const int nwriters = itsParset.getInt32("nwriters", 1);
    ASKAPLOG_DEBUG_STR(logger,"nwriters " << nwriters);
    frequency_tolerance = itsParset.getDouble("channeltolerance",0.0);

    ASKAPCHECK(nwriters > 0 ,"Number of writers must be greater than zero");

    /// Get the channel range
    /// The imager ususally uses the Channels keyword in the parset to
    /// defined the work unit that it will attempt. That does not work in this



    casa::uInt srow = 0;
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

    casa::uInt totChanIn = 0;



    for (unsigned int n = 0; n < ms.size(); ++n) {
        chanFreq[n].resize(0);
        chanWidth[n].resize(0);
        effectiveBW[n].resize(0);
        resolution[n].resize(0);
        centre[n].resize(0);
    // Open the input measurement set
        ASKAPLOG_DEBUG_STR(logger, "Opening " << ms[n] << " filecount " << n );
        const casa::MeasurementSet in(ms[n]);
        const casa::ROMSColumns srcCols(in);
        const casa::ROMSSpWindowColumns& sc = srcCols.spectralWindow();
        const casa::ROMSFieldColumns& fc = srcCols.field();
        const casa::ROMSObservationColumns& oc = srcCols.observation();
        const casa::ROMSAntennaColumns& ac = srcCols.antenna();
        const casa::ROArrayColumn<casa::Double> times = casa::ROArrayColumn<casa::Double>(oc.timeRange());
        const casa::ROArrayColumn<casa::Double> ants = casa::ROArrayColumn<casa::Double>(ac.position());
        const casa::uInt thisRef = casa::ROScalarColumn<casa::Int>(in.spectralWindow(),"MEAS_FREQ_REF")(0);

        casa::uInt thisChanIn = 0;
        srow = sc.nrow()-1;

        ASKAPCHECK(srow==0,"More than one spectral window not currently supported in adviseDI");

        // get the channel selection

        vector<int> itsChannels = getChannels();

        size_t chanStart=0;
        size_t chanStop=0;
        size_t chanStep=0;


        if (itsChannels[2] > 0) {
            // this also picks up whether the Channels keyword wad defined
            // we should be averaging
            // not sure how yet so just step for the moment
             // number of channels
            chanStart = itsChannels[1]; // either number of channels or start freq.
            chanStop = itsChannels[1] + itsChannels[0];
            chanStep = itsChannels[2];
            
            if (chanStop > casa::ROScalarColumn<casa::Int>(in.spectralWindow(),"NUM_CHAN")(0)) {
                chanStop = casa::ROScalarColumn<casa::Int>(in.spectralWindow(),"NUM_CHAN")(0);
                
            }
            thisChanIn = 0;
        }
        else {
            chanStart = 0;
            chanStop = casa::ROScalarColumn<casa::Int>(in.spectralWindow(),"NUM_CHAN")(0);
            chanStep = 1;

        }
        vector<double> itsFrequencies = getFrequencies();
        if (itsFrequencies[2] != 0) {
            // this means something has been set
            // The logic is the same as Channels <numchan> <start> <width>
            double freq_start = itsFrequencies[1];
            double delta = itsFrequencies[2];
            double freq_stop = freq_start + itsFrequencies[0]*delta;
            ASKAPLOG_INFO_STR(logger,"User has specified start:" << freq_start << " stop: " << freq_stop << " and width: " << delta);
        
        }


        ASKAPLOG_INFO_STR(logger, "Chan start: " << chanStart << " stop: " << chanStop << " step: " << chanStep);

        for (uint i = chanStart; i < chanStop; i = i + chanStep) {
          chanFreq[n].push_back(sc.chanFreq()(srow)(casa::IPosition(1, i)));
          chanWidth[n].push_back(sc.chanWidth()(srow)(casa::IPosition(1, i)));
          effectiveBW[n].push_back(sc.effectiveBW()(srow)(casa::IPosition(1, i)));
          resolution[n].push_back(sc.resolution()(srow)(casa::IPosition(1, i)));
          thisChanIn++;
        }

        totChanIn = totChanIn + thisChanIn;

        itsDirVec.push_back(fc.phaseDirMeasCol()(0));
        itsTangent.push_back(itsDirVec[n](0).getValue());

        // Read the position on Antenna 0
        Array<casa::Double> posval;
        ants.get(0,posval,true);
        vector<double> pval = posval.tovector();

        MVPosition mvobs(Quantity(pval[0], "m").getBaseValue(),
        Quantity(pval[1], "m").getBaseValue(),
        Quantity(pval[2], "m").getBaseValue());

        itsPosition.push_back(MPosition(mvobs,casa::MPosition::ITRF));

        // Get the Epoch
        Array<casa::Double> tval;
        vector<double> tvals;

        times.get(0,tval,true);
        tvals = tval.tovector();
        double mjd = tvals[0]/(86400.);
        casa::MVTime dat(mjd);

        itsEpoch.push_back(MVEpoch(dat.day()));

        itsRef = thisRef;

        ASKAPLOG_INFO_STR(logger, "Completed filecount " << n);
    }


    // ASKAPLOG_INFO_STR(logger, "Assuming tangent point shared: "<<printDirection(itsTangent[0])<<" (J2000)");


    itsFFrameFrequencies.resize(0);
    itsTopoFrequencies.resize(0);
    itsRequestedFrequencies.resize(0);
    
    // setup frequency frame
    const std::string freqFrame = itsParset.getString("freqframe","topo");
    if (freqFrame == "topo") {
        ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as topocentric");
        itsFreqRefFrame = casa::MFrequency::Ref(casa::MFrequency::TOPO);
        itsFreqType = casa::MFrequency::TOPO;
    } else if (freqFrame == "lsrk") {
        ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as lsrk");
        itsFreqRefFrame = casa::MFrequency::Ref(casa::MFrequency::LSRK);
        itsFreqType = casa::MFrequency::LSRK;
    } else if (freqFrame == "bary") {
        ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as barycentric");
        itsFreqRefFrame = casa::MFrequency::Ref(casa::MFrequency::BARY);
        itsFreqType = casa::MFrequency::BARY;
    } else {
        ASKAPTHROW(AskapError, "Unsupported frequency frame "<<freqFrame);
    }
    
   
    // At this point we now have each topocentric channel from each MS
    // in a unique array.
    // first we need to sort and uniqify the list
    // then resize the list to get the channel range.
    // This is required becuase we are trying to form a unique
    // reference channel list from the input measurement sets

    // This first loop just appends all the frequencies into 2 single arrays
    // the list of TOPO and FFRAME frequencies.


    itsAllocatedFrequencies.resize(nWorkersPerGroup);
    itsAllocatedWork.resize(nWorkers);

    for (unsigned int n = 0; n < ms.size(); ++n) {

        MeasFrame itsFrame(MEpoch(itsEpoch[n]),itsPosition[n],itsDirVec[n][0]);
        MFrequency::Ref refin(MFrequency::castType(itsRef),itsFrame);
        MFrequency::Ref refout(itsFreqType,itsFrame);
        MFrequency::Convert forw(refin,refout);
        MFrequency::Convert backw(refout,refin);

        // builds a list of all the FFRAME channels

        for (unsigned int ch = 0; ch < chanFreq[n].size(); ++ch) {

            ASKAPLOG_DEBUG_STR(logger, "CHECK --- File " << n << " Chan " << ch << " Freq " << chanFreq[n][ch]);
            /// possible output (desired frame) frequencies
            
            itsFFrameFrequencies.push_back(forw(chanFreq[n][ch]).getValue());
            /// actual input (topo) frequencies
            itsTopoFrequencies.push_back(MFrequency(MVFrequency(chanFreq[n][ch]),refin));

            /// The original scheme attempted to convert the input into the output frame
            /// and only keep those output (FFRAME) channels that matched.
            /// This was not deemed useful or helpful enough though. But I need to keep that
            /// mode as a default.
           
        }
    }

    /// uniquifying the lists

    bool (*custom_compare)(const casa::MFrequency& , const casa::MFrequency& ) = NULL;

    if (frequency_tolerance > 0) {
      custom_compare = compare_tol;
      ASKAPLOG_WARN_STR(logger,"Comparing frequencies with floating point tolerance of " << frequency_tolerance);
    }
    else {
      custom_compare = compare;
      ASKAPLOG_INFO_STR(logger,"Using standard compare for (zero tolerance) for freuqnecy allocations");
    }

    std::sort(itsFFrameFrequencies.begin(),itsFFrameFrequencies.end(), custom_lessthan);
    std::vector<casa::MFrequency>::iterator fframe_it;
    fframe_it = std::unique(itsFFrameFrequencies.begin(),itsFFrameFrequencies.end(),custom_compare);
    itsFFrameFrequencies.resize(std::distance(itsFFrameFrequencies.begin(),fframe_it));

    std::sort(itsTopoFrequencies.begin(),itsTopoFrequencies.end(), custom_lessthan);
    std::vector<casa::MFrequency>::iterator topo_it;
    topo_it = std::unique(itsTopoFrequencies.begin(),itsTopoFrequencies.end(),custom_compare);
    itsTopoFrequencies.resize(std::distance(itsTopoFrequencies.begin(),topo_it));
    ASKAPLOG_DEBUG_STR(logger," Unique sizes Topo " << itsTopoFrequencies.size() << " Bary " << itsFFrameFrequencies.size());

    for (unsigned int ch = 0; ch < itsTopoFrequencies.size(); ++ch) {

        ASKAPLOG_DEBUG_STR(logger,"Topocentric Channel " << ch << ":" << itsTopoFrequencies[ch]);
        ASKAPLOG_DEBUG_STR(logger,"Converted Channel " << ch << ":" << itsFFrameFrequencies[ch]);
        unsigned int allocation_index = floor(ch / nchanpercore);
        /// We allocate the frequencies based upon the topocentric range.
        /// We do this becuase it is easier for the user to understand.
        /// Plus - all beams will have the same allocation. Which will produce cubes/images
        /// that will easily merge.

        /// Beware the syntactic confusion here - we are allocating a frequency that is from
        /// the Topocentric list. But will match a channel based upon the barycentric frequency

        /// need to trim if itsChannels has been set
        ASKAPLOG_DEBUG_STR(logger,"Allocating frequency "<< itsTopoFrequencies[ch].getValue() \
        << " to worker " << allocation_index+1);

        itsAllocatedFrequencies[allocation_index].push_back(itsTopoFrequencies[ch].getValue());
    }


    // Now for each allocated workunit we need to fill in the rest of the workunit
    // we now have a workUnit for each channel in the allocation - but not
    // for each Epoch.

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
                int lc = 0;

                lc = match(set,thisAllocation[frequency]);
                if (lc >= 0) {
                    // there is a channel of this frequency in the measurement set


                    cp::ContinuumWorkUnit wu;

                    wu.set_payloadType(cp::ContinuumWorkUnit::WORK);
                    wu.set_channelFrequency(thisAllocation[frequency]);
                    wu.set_beam(myBeam);

                    if (itsTopoFrequencies.size() > 1)
                        wu.set_channelWidth(fabs(itsTopoFrequencies[1].getValue() - itsTopoFrequencies[0].getValue()));
                    else
                        wu.set_channelWidth(fabs(chanWidth[0][0]));

                    wu.set_localChannel(lc);
                    wu.set_globalChannel(globalChannel);
                    wu.set_dataset(ms[set]);
                    itsAllocatedWork[work].push_back(wu);
                    itsWorkUnitCount++;
                    ASKAPLOG_DEBUG_STR(logger,"MATCH Allocating barycentric freq " << thisAllocation[frequency] \
                    << " with local channel number " << lc << " ( " << chanFreq[set][lc] << " ) of width " << wu.get_channelWidth()  \
                    << " in set: " << ms[set] <<  " to rank " << work+1 << " this rank has " \
                    << itsAllocatedWork[work].size() << " of a total count " << itsWorkUnitCount \
                    << " the global channel is " << globalChannel);

                    allocated = true;
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


    unsigned int nWorkersPerWriter = floor(itsAllocatedWork.size() / nwriters);
    int mywriter = 0;
    for (int wrk = 0; wrk < itsAllocatedWork.size(); wrk++) {
        if (nwriters>1) {
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
                ASKAPCHECK(mywriter < itsAllocatedWork.size(),"Ran out of eligible writers");
            }
        }

        for (int unit = 0; unit < itsAllocatedWork[wrk].size(); unit++) {
                    itsAllocatedWork[wrk][unit].set_writer(mywriter+1); // plus 1 for rank
                    ASKAPLOG_DEBUG_STR(logger,"Set rank " << wrk+1 << " writer to be rank " << mywriter+1);
        }

    }
    for (int grp = 1; grp < itsComms.nGroups(); grp++) {
        for (int wrk = 0; wrk < nWorkersPerGroup; wrk++) {
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
    cp::ContinuumWorkUnit rtn;
    if (itsAllocatedWork[id].empty() == true) {
        ASKAPLOG_WARN_STR(logger, "Stack is empty for " << id+1);
        rtn.set_payloadType(cp::ContinuumWorkUnit::DONE);
        return rtn;
    }
    else {
        rtn = itsAllocatedWork[id].back();
        itsAllocatedWork[id].pop_back();
        itsWorkUnitCount--;
    }
    if (itsAllocatedWork[id].empty() == true) {
        // this is the last unitParset
        ASKAPLOG_WARN_STR(logger, "Final job for " << id+1);
        if (rtn.get_payloadType() != cp::ContinuumWorkUnit::NA){
            rtn.set_payloadType(cp::ContinuumWorkUnit::LAST);
        }
        else {
            ASKAPLOG_WARN_STR(logger, "Final job is bad for " << id+1);
            rtn.set_payloadType(cp::ContinuumWorkUnit::DONE);
        }
    }
    return rtn;
}

int AdviseDI::match(int ms_number, casa::MVFrequency testFreq) {
    /// Which channel does the frequency correspond to.
    /// IF the barycentr flag has been set then this will match
    /// the barycentred channel to it.
    vector<double>::iterator it_current = chanFreq[ms_number].begin();
    vector<double>::iterator it_end = chanFreq[ms_number].end()-1;
    double testVal = testFreq.getValue();

    if (in_range(*it_current,*it_end,testVal)) {
        int ch = 0;
        it_current=chanFreq[ms_number].begin();
        for (ch=0 ; ch < chanFreq[ms_number].size(); ++ch) {
            ASKAPLOG_DEBUG_STR(logger, "looking for " << testVal << \
            " in test frequency channel " << *it_current << \
                " width " << chanWidth[ms_number][ch]);
            double one_edge = (*it_current) - chanWidth[ms_number][ch]/2.;
            double other_edge = (*it_current) + chanWidth[ms_number][ch]/2.;

            if (in_range(one_edge,other_edge,testVal)) {

                return ch;
            }
            it_current++;

        }
    }

    return -1;


}
void AdviseDI::addMissingParameters() {
    this->addMissingParameters(this->itsParset);
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
      std::ostringstream pstr;
      // Only J2000 is implemented at the moment.
      pstr<<"["<<printLon(itsTangent[n])<<", "<<printLat(itsTangent[n])<<", J2000]";
      ASKAPLOG_INFO_STR(logger, "  updating parameter " << param << ": " << pstr.str().c_str());
      parset.replace(param, pstr.str().c_str());

      for (size_t img = 0; img < imageNames.size(); ++img) {
        param ="Images."+imageNames[img]+".direction";

        std::ostringstream pstr;
        // Only J2000 is implemented at the moment.
        pstr<<"["<<printLon(itsTangent[n])<<", "<<printLat(itsTangent[n])<<", J2000]";
        ASKAPLOG_INFO_STR(logger, "  updating parameter " << param << ": " << pstr.str().c_str());
        parset.replace(param, pstr.str().c_str());

      }
    }
  }
}

void AdviseDI::addMissingParameters(LOFAR::ParameterSet& parset)
{

    ASKAPLOG_INFO_STR(logger,"Adding missing params ");

    if (isPrepared == true) {
        ASKAPLOG_INFO_STR(logger,"Prepared therefore can add frequency label for the output image");
        std::vector<casa::MFrequency>::iterator begin_it;
        std::vector<casa::MFrequency>::iterator end_it;
        if (barycentre) {
            begin_it = itsFFrameFrequencies.begin();
            end_it = itsFFrameFrequencies.end()-1;
        }
        else {
            begin_it = itsTopoFrequencies.begin();
            end_it = itsTopoFrequencies.end()-1;

        }
        this->minFrequency = (*begin_it).getValue();
        this->maxFrequency = (*end_it).getValue();
        ASKAPLOG_INFO_STR(logger,"Min:Max frequency -- " << this->minFrequency << ":" << this->maxFrequency);


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
       std::ostringstream pstr;
       // Only J2000 is implemented at the moment.
       pstr<<"["<<printLon(itsTangent[0])<<", "<<printLat(itsTangent[0])<<", J2000]";
       ASKAPLOG_INFO_STR(logger, "  Advising on parameter (getting from the tangent point on the first measurement set) " << param << ": " << pstr.str().c_str());
       itsParset.add(param, pstr.str().c_str());
   }
   param = "Images.restFrequency";

   if ( !parset.isDefined(param) ) {
       std::ostringstream pstr;
       // Only J2000 is implemented at the moment.
       pstr<<"HI";
       ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << ": " << pstr.str().c_str());
       parset.add(param, pstr.str().c_str());
   }

   for (size_t img = 0; img < imageNames.size(); ++img) {

     param = "Images."+imageNames[img]+".cellsize";
     if ( !parset.isDefined(param) ) {
       cellsizeNeeded = true;
     }
     else {
       param = "Images.cellsize";
       if (!parset.isDefined(param) ) {
         const vector<string> cellSizeVector = parset.getStringVector("Images.cellsize");
         std::ostringstream pstr;
         pstr<<"["<< cellSizeVector[0].c_str() <<"arcsec,"<<cellSizeVector[1].c_str() <<"arcsec]";
         ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
         parset.add(param, pstr.str().c_str());
       }
     }
     param = "Images."+imageNames[img]+".shape";
     if ( !parset.isDefined(param) ) shapeNeeded = true;

     param = "Images."+imageNames[img]+".frequency";
     if ( !parset.isDefined(param) && isPrepared == true) {

       const string key="Images."+imageNames[img]+".frequency";
       char tmp[64];
       // changing this to match adviseParallel
       const double aveFreq = 0.5*(minFrequency+maxFrequency);
       sprintf(tmp,"[%f,%f]",aveFreq,aveFreq);
       string val = string(tmp);
       ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << val);
       parset.add(key,val);

     }
     param ="Images."+imageNames[img]+".direction";
     if ( !parset.isDefined(param) ) {

       if (parset.isDefined("Images.direction") ) {
         const std::vector<std::string> direction = parset.getStringVector("Images.direction");

         const double ra = SynthesisParamsHelper::convertQuantity(direction[0],"rad");
         const double dec = SynthesisParamsHelper::convertQuantity(direction[1],"rad");
         const casa::MVDirection itsDirection = casa::MVDirection(ra,dec);

         std::ostringstream pstr;
         // Only J2000 is implemented at the moment.
         pstr<<"["<<printLon(itsDirection)<<", "<<printLat(itsDirection)<<", J2000]";
         ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << " (obtained from Images.direction) : " << pstr.str().c_str());
         itsParset.add(param, pstr.str().c_str());
       }
       else {

         std::ostringstream pstr;
         // Only J2000 is implemented at the moment.
         pstr<<"["<<printLon(itsTangent[0])<<", "<<printLat(itsTangent[0])<<", J2000]";
         ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << "(obtained from tangent point of first MS): " << pstr.str().c_str());
         parset.add(param, pstr.str().c_str());
       }
     }
     param = "Images."+imageNames[img]+".nterms"; // if nterms is set, store it for later
     if (parset.isDefined(param)) {
       if ((nTerms>1) && (nTerms!=itsParset.getInt(param))) {
         ASKAPLOG_WARN_STR(logger, "  Imaging with different nterms may not work");
       }
       nTerms = itsParset.getInt(param);
     }

     if ( !parset.isDefined("Images."+imageNames[img]+".nchan") ) {

     }
   }

   if (nTerms > 1) { // check required MFS parameters
       param = "visweights"; // set to "MFS" if unset and nTerms > 1
       if (!parset.isDefined(param)) {
           std::ostringstream pstr;
           pstr<<"MFS";
           ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<" (obtained by default - we know no other weighting) : " << pstr.str().c_str());
           parset.add(param, pstr.str().c_str());
       }

       param = "visweights.MFS.reffreq"; // set to average frequency if unset and nTerms > 1
       if ((parset.getString("visweights")=="MFS")) {
           if (!itsParset.isDefined(param)) {
               char tmp[64];
               const double aveFreq = 0.5*(minFrequency+maxFrequency);
               sprintf(tmp,"%f",aveFreq);
               string val = string(tmp);
               ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<" (using average frequency):  " << val);
               parset.add(param,val);
           }

       }
   }

   // test for general missing parameters:
   if ( cellsizeNeeded && !parset.isDefined("nUVWMachines") ) {

   } else if ( cellsizeNeeded && !parset.isDefined("Images.cellsize") ) {

   } else if ( shapeNeeded && !parset.isDefined("Images.shape") ) {

   }
   ASKAPLOG_INFO_STR(logger,"Done adding missing params ");

}
// Utility function to get dataset names from parset.
std::vector<std::string> AdviseDI::getDatasets()
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

            LOFAR::ostringstream ss;
            ss << "dataset" << idx + 1;
            key = ss.str();
            ++idx;
        }
    }

    return ms;
}

/// the adviseDI should be smart enough to tell the difference between a straight list and 
/// actual frequencies .... or is that too different.
std::vector<double> AdviseDI::getFrequencies() {
    std::vector<double> f(3,0);
    std::vector<string> fstr(3,"0");

    if (!itsParset.isDefined("Frequencies")) {
    ASKAPLOG_WARN_STR(logger,
        "Frequencies keyword is not defined");
        
    }
    else {
        fstr = itsParset.getStringVector("Frequencies",true);
        f[0] = atof(fstr[0].c_str());
        f[1] = atof(fstr[1].c_str());
        f[2] = atof(fstr[2].c_str());
         // just a start and stop - assume no averaging
        if (f[2] == 0) {
           ASKAPLOG_WARN_STR(logger,
           "Channel width not specified this setting will be ignored"); 
        }
        
    }
    return f;
}

std::vector<int> AdviseDI::getChannels() {

    // channels should now behave more like the historical version
    // <nchan> <start> ...
    //

    std::vector<int> c(3,0);
    std::vector<string> cstr(3,"0");

    if (!itsParset.isDefined("Channels")) {
    ASKAPLOG_WARN_STR(logger,
        "Channels keyword is not defined");
        c[2] = -1;
    }
    else {
        cstr = itsParset.getStringVector("Channels",true);
        c[0] = atof(cstr[0].c_str());
        string wild = "%w";
        if (cstr[1].compare(wild) == 0) {
            c[1] = 0;
            c[2] = -1;
        }
        else {
            c[1] = atoi(cstr[1].c_str());
            c[2] = atoi(cstr[2].c_str());
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
        int last_channel = itsAllocatedWork[worker][0].get_globalChannel();
        for (int alloc=0; alloc < itsAllocatedWork[worker].size() ; alloc++) {
            int current_channel = itsAllocatedWork[worker][alloc].get_globalChannel();
            // need to test whether this is a distinct channel or a different epoch for
            // the same epoch
            if (current_channel != last_channel || alloc == 0) {
                if (itsAllocatedWork[worker][alloc].get_payloadType() != cp::ContinuumWorkUnit::NA) {
                    itsCubeComms.addWriter(itsAllocatedWork[worker][alloc].get_writer());

                    itsCubeComms.addChannelToWriter(itsAllocatedWork[worker][alloc].get_writer(),worker+1);
                    itsCubeComms.addChannelToWorker(worker+1);
                }
                last_channel = current_channel;
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
std::vector<int> AdviseDI::getBeams()
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
double AdviseDI::getBaseFrequencyAllocation(int workerNumber) {
    return itsAllocatedFrequencies[workerNumber][0];
}
} // namespace synthesis

} // namespace askap
