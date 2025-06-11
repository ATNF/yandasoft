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

#include <askap/messages/ContinuumWorkUnit.h>
#include <askap/distributedimager/CubeComms.h>

#include <Common/ParameterSet.h>
#include <askap/parallel/MEParallelApp.h>
#include <askap/measurementequation/VisMetaDataStats.h>
#include <casacore/casa/Quanta/MVDirection.h>
#include <askap/parallel/AdviseParallel.h>

#include <boost/shared_ptr.hpp>
#include <string>

#include <casacore/casa/BasicSL.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>
#include <casacore/ms/MSOper/MSReader.h>
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
            /// @param[in] extra - if true check for additional parameters besides frequencies and directions
            void addMissingParameters(bool extra=false);

            /// @brief Use workunit to set Image direction parameters in parset
            void updateDirectionFromWorkUnit(const askap::cp::ContinuumWorkUnit& wu);

            /// @brief the datasets
            std::vector<std::string> getDatasets() const;

            /// @brief get the channels
            std::vector<int> getChannels() const;

            /// @brief get the frequencies
            std::vector<double> getFrequencies() const;

            /// @brief get the epoch for each MS
            inline const casacore::MVEpoch& getEpoch(int ms=0) const {return itsEpoch[ms]; };

            inline const vector<casacore::MFrequency>& getTopoFrequencies() const {return itsInputFrequencies;};

            cp::ContinuumWorkUnit getAllocation(int id);

            double getBaseFrequencyAllocation(int workerNumber) const;

            void updateComms();

            inline int getWorkUnitCount() const { return itsWorkUnitCount;};

            void prepare();

            /// obtain frequency reference frame
            inline const casacore::MFrequency::Ref& getFreqRefFrame() const { return itsFreqRefFrame;}

        protected:

            /// @brief obtain metadata stats for one dataset
            /// @details Hopefully this is a temporary method until this class (or the whole imager) can be
            /// redesigned, at least to avoid accessing measurement set from first principles. Currently, this
            /// method adds to uglyness of the code - it uses the functionality of the base class (i.e. the originally
            /// written estimator) to get the basic info for a single measurement set without relying on the correct
            /// parallel distribution. The returned statistics is valid until the next call to this method or until
            /// the full parallel advise is done.
            /// @param[in] ms measurement set to work with
            /// @return const reference to the object populated with resulted metadata statistics
            const VisMetaDataStats& computeVisMetaDataStats(const std::string &ms);

        private:

            bool itsPrepared;

            bool itsDetailedLog;

            int itsWorkUnitCount;

            LOFAR::ParameterSet& itsParset;

            casacore::uInt itsRef;

            vector<casacore::MFrequency> itsInputFrequencies;

            vector<casacore::MFrequency> itsRequestedFrequencies;

            /// @brief reference frame for frequency
            /// @details We may want to simulate/image in different reference frames.
            /// This field contains the reference frame selected in the parset.
            casacore::MFrequency::Ref itsFreqRefFrame;

            casacore::MFrequency::Types itsFreqType;

            double itsMinFrequency;

            double itsMaxFrequency;

            std::vector<casacore::MVDirection> itsTangent;

            std::vector<casacore::MDirection> itsDirMeas;

            std::vector<casacore::MVEpoch> itsEpoch;

            std::vector<casacore::MPosition> itsPosition;

            std::vector< std::vector<double> > itsChanFreq;
            std::vector< std::vector<double> > itsChanWidth;

            std::vector< std::vector<double> > itsAllocatedFrequencies;
            std::vector< std::vector<cp::ContinuumWorkUnit> > itsAllocatedWork;

            vector<int> matchAll(int,  casacore::MVFrequency, casacore::MVFrequency) const;

            std::vector<int> getBeams() const;


        };

    } // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_ADVISE_DI_H
