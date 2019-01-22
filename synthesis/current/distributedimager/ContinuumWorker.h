/// @file ContinuumWorker.h
///
/// @copyright (c) 2009 CSIRO
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>
/// @author Stephen Ord <stephen.ord@csiro.au>

#ifndef ASKAP_CP_SIMAGER_CONTINUUMWORKER_H
#define ASKAP_CP_SIMAGER_CONTINUUMWORKER_H

// System includes
#include <string>

// ASKAPsoft includes
#include "boost/shared_ptr.hpp"
#include <Common/ParameterSet.h>
#include <fitting/INormalEquations.h>
#include <fitting/Params.h>
#include <dataaccess/TableDataSource.h>
#include <gridding/IVisGridder.h>

// Local includes
#include "distributedimager/AdviseDI.h"
#include "distributedimager/MSSplitter.h"
#include "distributedimager/CalcCore.h"
#include "messages/ContinuumWorkUnit.h"
#include "distributedimager/CubeBuilder.h"
#include "distributedimager/CubeComms.h"
namespace askap {
namespace cp {

class ContinuumWorker

    {
    public:
        ContinuumWorker(LOFAR::ParameterSet& parset,
                           CubeComms& comms);
        ~ContinuumWorker();

        void run(void);




    private:

        // My Advisor
        boost::shared_ptr<synthesis::AdviseDI> itsAdvisor;
         // The work units
        vector<ContinuumWorkUnit> workUnits;
        // cached files
        vector<string> cached_files;

        // Whether preconditioning has been requested
        bool itsDoingPreconditioning;

        // Whether the gridder is a Mosaicking one
        bool itsGridderCanMosaick;

        // Cache a workunit to a different location
        void cacheWorkUnit(ContinuumWorkUnit& wu, LOFAR::ParameterSet& unitParset);
        // Process a workunit
        void processWorkUnit(ContinuumWorkUnit& wu);
        // Delete a workunit from the cache
        void deleteWorkUnitFromCache(ContinuumWorkUnit& wu, LOFAR::ParameterSet& unitParset);
        // clear the current cached files
        void clearWorkUnitCache();

        // Vector of the stored parsets of the work allocations
        vector<LOFAR::ParameterSet> itsParsets;


        //For all workunits .... process

        void processChannels();

        // For a given workunit, just process a single snapshot - the channel is specified
        // in the parset ...
        void processSnapshot(LOFAR::ParameterSet& parset);


        // Setup the image specified in itsParset and add it to the Params instance.
        void setupImage(const askap::scimath::Params::ShPtr& params,
                    double channelFrequency, bool shapeOveride = false);

        void buildSpectralCube();

        // Root Parameter set good for information common to all workUnits
        LOFAR::ParameterSet& itsParset;

        // Communications class
        CubeComms& itsComms;

        // No support for assignment
        ContinuumWorker& operator=(const ContinuumWorker& rhs);

        // No support for copy constructor
        ContinuumWorker(const ContinuumWorker& src);

        // ID of the master process
        static const int itsMaster = 0;

        // List of measurement sets to work on
        vector<std::string> datasets;

        // the basechannel number assigned to this worker
        unsigned int baseChannel;

        // the baseFrequency associated with this channel
        double baseFrequency;
        // the baseFrequency associated with the cube if being built
        double baseCubeFrequency;
        // the global channel associated with this part of the cube
        int baseCubeGlobalChannel;
        // the number of channels in this cube (if writer)
        int nchanCube;

        boost::shared_ptr<CubeBuilder> itsImageCube;
        boost::shared_ptr<CubeBuilder> itsPSFCube;
        boost::shared_ptr<CubeBuilder> itsResidualCube;
        boost::shared_ptr<CubeBuilder> itsWeightsCube;
        boost::shared_ptr<CubeBuilder> itsPSFimageCube;
        boost::shared_ptr<CubeBuilder> itsRestoredCube;

        void handleImageParams(askap::scimath::Params::ShPtr params, unsigned int chan);

        void copyModel(askap::scimath::Params::ShPtr SourceParams, askap::scimath::Params::ShPtr SinkParams);

        void recordBeam(const askap::scimath::Axes &axes, const unsigned int globalChannel);
        void storeBeam(const unsigned int cubeChannel);

        std::map<unsigned int, casa::Vector<casa::Quantum<double> > > itsBeamList;

        unsigned int itsBeamReferenceChannel;
        void logBeamInfo();


};

};
};

#endif
