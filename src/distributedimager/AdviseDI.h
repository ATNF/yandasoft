/// @file
///
/// Support for parallel statistics accululation to advise on imaging parameters
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

#ifndef ASKAP_IMAGER_ADVISE_DI_H
#define ASKAP_IMAGER_ADVISE_DI_H

#include <messages/ContinuumWorkUnit.h>
#include <distributedimager/CubeComms.h>

#include <Common/ParameterSet.h>
#include <parallel/MEParallelApp.h>
#include <measurementequation/VisMetaDataStats.h>
#include <casacore/casa/Quanta/MVDirection.h>
#include <parallel/AdviseParallel.h>

#include <boost/shared_ptr.hpp>
#include <string>

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


namespace askap {

    namespace synthesis {

        /// @brief parallel helper for the advise utility
        /// @details This class does the core operation to run statistics estimators on every measurement set
        /// and aggregating the result. Most non-trivial actions happen in the parallel mode.
        /// @note It may be a bit untidy to derive this class from MEParallelApp just to reuse a bunch of existing code,
        /// but some subtle features like frequency conversion setup may come handy in the future. The goal is that it should
        /// work with only the single parameter present in the parset which describes the measurement set(s).
        /// @ingroup parallel
        class AdviseDI : public AdviseParallel
        {
        public:
            /// @brief Constructor from ParameterSet
            /// @details The parset is used to construct the internal state. We could
            /// also support construction from a python dictionary (for example).
            /// The command line inputs are needed solely for MPI - currently no
            /// application specific information is passed on the command line.
            /// @param comms communication object
            /// @param parset ParameterSet for inputs

            AdviseDI(askap::cp::CubeComms& comms, LOFAR::ParameterSet& parset);

            /// @brief Add the missing parameters
            /// @details Add whatever details we require for both master and
            /// worker implementations

            void addMissingParameters(LOFAR::ParameterSet& parset);

            void addMissingParameters();


            void updateDirectionFromWorkUnit(LOFAR::ParameterSet& parset, askap::cp::ContinuumWorkUnit& wu);

            /// @brief Access to the parset
            /// @details Returns the current parset after advice

            LOFAR::ParameterSet getParset() { return itsParset; };

            /// @brief the datasets
            std::vector<std::string> getDatasets();


            /// @brief get the channels
            std::vector<int> getChannels();

            /// @brief get the frequencies
            std::vector<double> getFrequencies();
            
            casa::MVDirection getTangent(int ms=0) {return itsTangent[ms];};

            casa::MVEpoch getEpoch(int ms=0) {return itsEpoch[ms]; };

            casa::MPosition getPosition(int ms=0) {return itsPosition[ms]; };

            vector<casa::MFrequency> getBaryFrequencies() {return itsFFrameFrequencies;};

            vector<casa::MFrequency> getTopoFrequencies() {return itsTopoFrequencies;};

            cp::ContinuumWorkUnit getAllocation(int id);

            double getBaseFrequencyAllocation(int workerNumber);

            void updateComms();

            int getWorkUnitCount() { return itsWorkUnitCount;};

            void prepare();

            bool isPrepared;

            bool barycentre;
            /// obtain frequency reference frame
            inline casa::MFrequency::Ref getFreqRefFrame() const { return itsFreqRefFrame;}

        private:

            int itsWorkUnitCount;

            LOFAR::ParameterSet& itsParset;

            casa::uInt itsRef;

            vector<casa::MFrequency> itsFFrameFrequencies;

            vector<casa::MFrequency> itsTopoFrequencies;

            vector<casa::MFrequency> itsRequestedFrequencies;

            /// @brief reference frame for frequency
            /// @details We may want to simulate/image in different reference frames.
            /// This field contains the reference frame selected in the parset.
            casa::MFrequency::Ref itsFreqRefFrame;  
            
            casa::MFrequency::Types itsFreqType;

            double minFrequency;

            double maxFrequency;



            std::vector<casa::MVDirection> itsTangent;

            std::vector<casa::Vector<casa::MDirection> > itsDirVec;

            std::vector<casa::MVEpoch> itsEpoch;

            std::vector<casa::MPosition> itsPosition;

            std::vector< std::vector<double> > chanFreq;
            std::vector< std::vector<double> > chanWidth;
            std::vector< std::vector<double> > effectiveBW;
            std::vector< std::vector<double> > resolution;
            std::vector< std::vector<double> > centre;

            std::vector< std::vector<double> > itsAllocatedFrequencies;
            std::vector< std::vector<cp::ContinuumWorkUnit> > itsAllocatedWork;

            int match(int ms_number,  casa::MVFrequency freq);

            std::vector<int> getBeams();


        };

    } // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_ADVISE_PARALLEL_H
