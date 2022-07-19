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
#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/gridding/IVisGridder.h>
#include <askap/askap/StatReporter.h>

// Local includes
#include "askap/distributedimager/AdviseDI.h"
#include "askap/distributedimager/MSSplitter.h"
#include "askap/distributedimager/CalcCore.h"
#include "askap/messages/ContinuumWorkUnit.h"
#include "askap/distributedimager/CubeBuilder.h"
#include "askap/distributedimager/CubeComms.h"
#include <askap/utils/StatsAndMask.h>

namespace askap {
namespace cp {

class ContinuumWorker

    {
    public:
        ContinuumWorker(LOFAR::ParameterSet& parset,
                           CubeComms& comms, StatReporter& stats);
        ~ContinuumWorker();

        void run(void);

        void writeCubeStatistics();

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
        void cacheWorkUnit(ContinuumWorkUnit& wu);
        // Process a workunit
        void preProcessWorkUnit(ContinuumWorkUnit& wu);
        // Compress all continuous channel allocations into individual workunits
        void compressWorkUnits();

        // Delete a workunit from the cache
        void deleteWorkUnitFromCache(ContinuumWorkUnit& wu);
        // clear the current cached files
        void clearWorkUnitCache();

        //For all workunits .... process

        void processChannels();

        // For a given workunit, just process a single snapshot - the channel is specified
        // in the parset ...
        void processSnapshot();


        // Setup the image specified in parset and add it to the Params instance.
        void setupImage(const askap::scimath::Params::ShPtr& params,
                    double channelFrequency, bool shapeOveride = false);

        void buildSpectralCube();

        // Root Parameter set good for information common to all workUnits
        LOFAR::ParameterSet& itsParset;

        // Communications class
        CubeComms& itsComms;

        // statistics
        StatReporter& itsStats;

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

        boost::shared_ptr<CubeBuilder<casacore::Float> > itsImageCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPSFCube;
        boost::shared_ptr<CubeBuilder<casacore::Complex> > itsPCFGridCube;
        boost::shared_ptr<CubeBuilder<casacore::Complex> > itsPSFGridCube;
        boost::shared_ptr<CubeBuilder<casacore::Complex> > itsVisGridCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPCFGridCubeReal;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPSFGridCubeReal;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsVisGridCubeReal;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPCFGridCubeImag;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPSFGridCubeImag;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsVisGridCubeImag;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsResidualCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsWeightsCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPSFimageCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsRestoredCube;
        boost::shared_ptr<askap::utils::StatsAndMask> itsRestoredStatsAndMask;
        boost::shared_ptr<askap::utils::StatsAndMask> itsResidualStatsAndMask;
        std::string itsWeightsName;
        std::string itsGridType;

        // calculate the statistic of the per plane image
        void calculateImageStats(boost::shared_ptr<askap::utils::StatsAndMask> statsAndMask,
                                 boost::shared_ptr<CubeBuilder<casacore::Float> > imgCube,
                                 int channel, const casacore::Array<float>& arr);

        void handleImageParams(askap::scimath::Params::ShPtr params, unsigned int chan);

        void copyModel(askap::scimath::Params::ShPtr SourceParams, askap::scimath::Params::ShPtr SinkParams);

        void initialiseBeamLog(const unsigned int numChannels);
        void recordBeam(const askap::scimath::Axes &axes, const unsigned int globalChannel);
        //void storeBeam(const unsigned int cubeChannel);

        std::map<unsigned int, casacore::Vector<casacore::Quantum<double> > > itsBeamList;
        unsigned int itsBeamReferenceChannel;
        void logBeamInfo();

        void initialiseWeightsLog(const unsigned int numChannels);
        void recordWeight(float wt, const unsigned int globalChannel);
        std::map<unsigned int, float> itsWeightsList;
        void logWeightsInfo();

        /// @brief Do we want a restored image?
        bool itsRestore;

        /// @brief Do we want a residual image
        bool itsWriteResidual;

        /// @brief write 'raw', unnormalised, natural weight psf
        bool itsWritePsfRaw;

        /// @brief write normalised, preconditioned psf
        bool itsWritePsfImage;

        /// @brief write weights image
        bool itsWriteWtImage;

        /// @brief write a weights log
        bool itsWriteWtLog;

        /// @brief write out the (clean) model image
        bool itsWriteModelImage;

        /// @brief write out the gridded data, pcf and psf
        bool itsWriteGrids;

        /// @brief write out the grids with UV coordinate grid
        bool itsGridCoordUV;

        /// @brief write out the FFT of the grids
        bool itsGridFFT;

        /// @brief the number of rank that can write to the cube
        int itsNumWriters;

};

};
};

#endif
