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

// boost includes
#include <boost/noncopyable.hpp>
#include "boost/shared_ptr.hpp"

// ASKAPsoft includes
#include <Common/ParameterSet.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/askap/StatReporter.h>

// Local includes
#include "askap/distributedimager/AdviseDI.h"
#include "askap/messages/ContinuumWorkUnit.h"
#include "askap/distributedimager/WorkUnitContainer.h"
#include "askap/distributedimager/CubeBuilder.h"
#include "askap/distributedimager/CubeComms.h"
#include "askap/distributedimager/DataSourceManager.h"
#include "askap/distributedimager/CalcCore.h"
#include <askap/utils/StatsAndMask.h>
#include <askap/gridding/IUVWeightCalculator.h>

namespace askap {
namespace cp {

class ContinuumWorker : public boost::noncopyable

    {
    public:
        ContinuumWorker(LOFAR::ParameterSet& parset,
                           CubeComms& comms, StatReporter& stats);

        void run(void);

        void writeCubeStatistics();

    private:
        /// @brief check if sample density grid needs to be built
        /// @details For now, use the shared pointer carrying uv weight calculator as a flag that we need to
        /// build grid of weights (this is what uv weight calculator works with). The appropriate factory method
        /// is called in the constructor, so technically the output of this method is predetermined for all 
        /// instances of this class.
        /// @return true if sample density grid needs to be built
        bool inline isSampleDensityGridNeeded() const { return static_cast<bool>(itsUVWeightCalculator); }

        /// @brief configure allocation in channels
        /// @details This method sets up channel allocation for cube writing (in local solver mode)
        /// or combines channels in the global solver mode if configured in the parset.
        void configureChannelAllocation();

        /// @brief configure reference channel used for the restoring beam
        /// @details This method populates itsBeamReferenceChannel based on the parset and
        /// the number of channels allocated to the cube handled by this rank
        void configureReferenceChannel();

        /// @brief initialise cube writing if in local solver mode
        /// @details This method encapsulates the code which handles cube writing in the local solver mode
        /// (i.e. when it is done from the worker). Safe to call in continuum mode too (as itsComms.isWriter
        /// would return false in this case)
        void initialiseCubeWritingIfNecessary();

        /// @brief check stopping thresholds in the model 
        /// @details This method is used at the end of minor cycle deconvolution to check whether to continue iterations.
        /// @param[in] model shared pointer to the scimath::Params object with the model
        /// @return true if stopping is required
        /// @note We do similar checks in both the master and in workers. So this method can be moved somewhere else to be shared.
        bool checkStoppingThresholds(const boost::shared_ptr<scimath::Params> &model) const;


        /// @brief add current image as a model
        /// @details This method adds fullres (if present) or ordinary image as model.slice in the given params object.
        /// It is expected that model.slice will be absent and either fullres or Nyquist resolution image should be
        /// present.
        /// @param[in] params shared pointer to the params object to work with (should be non-empty)
        /// @note I (MV) think there could be untidy design here - we probably make an extra copy which could be avoided
        static void addImageAsModel(const boost::shared_ptr<scimath::Params> &params);

        /// @brief perform one write job for a remote client
        /// @details This method is expected to be used for cube writing ranks only. It receives a single
        /// write job and performs it.
        void performSingleWriteJob(); 

        /// @brief cleanup outstanding write jobs
        /// @details This method is used for cube writing ranks (does nothing for non-writers), it loops over all outstanding
        /// work units and performs write operation assigned to this rank.
        /// @param[in] targetOutstanding desired number of outstanding write jobs at the end of execution
        ///                              (to spread writing across the iteration), default is all jobs
        /// @param[in] minOutstanding    minimal number of outstanding jobs to remain (default - none)
        /// @note (MV) I didn't fully understand the logic behind targetOutstanding and minOutstanding (one should be
        /// sufficient), the same behaviour as we had prior to refactoring has been implemented.
        void performOutstandingWriteJobs(int targetOutstanding = 0, int minOutstanding = -1);

        /// @brief perform write job allocated to this rank
        /// @details unlike performOutstandingWriteJobs or performSingleWriteJob this method deals with the write
        /// job handled entirely by this rank (i.e. its own write job) and, hence, provided explicitly rather than
        /// received from another rank.
        /// @param[in] globalChannel global channel (i.e. channel in the whole cube) to write
        /// @param[in] params shared pointer to the model with required info (should not be empty)
        void performOwnWriteJob(unsigned int globalChannel, const boost::shared_ptr<scimath::Params> &params);

        /// @brief send blank image to writer
        /// @details This method is expected to be used when calculation of a spectral plane is failed for some reason,
        /// but some other rank is responsible for writing it. Essentially it sends a blank image with parameters 
        /// (like frequency and channel) filled from the work unit. 
        /// @param[in] wu work unit to take the information from
        /// @note (MV:) This doesn't seem like a good design, but the behaviour is left the same as it was prior to
        /// the refactoring.
        void sendBlankImageToWriter(const cp::ContinuumWorkUnit &wu) const;

        /// @brief figure out if preconditioning is to be done
        /// @details This method encapsulates checks of the parset indicating that preconditioning is going to be done.
        /// It is necessary to configure writing of additional data products, although preconditioning itself 
        /// is enabled and done by the appropriate solver class.
        /// @param[in] parset parset to use (this method is expected to be used in the constructor, so it is handy not
        /// to rely on the itsParset data field).
        /// @return true if preconditioning is to be done, false otherwise
        static bool doingPreconditioning(const LOFAR::ParameterSet &parset);

        /// @brief helper method to obtain the number of writers for the cube
        /// @details This method obtains the number of writers from the parset and adjusts it if necessary.
        /// It is intended to be used in the constructor to fill itsNumWriters data field and requires 
        /// itsGridType and itsParset to be valid.
        /// @return number of writer ranks for grid export
        int configureNumberOfWriters();

        /// @brief helper method to process a single work unit
        /// @details This method encapsulates a part of the old processChannel calculating and merging NE for a single work 
        /// unit. The resulting NE is either added to the root imager passed as the parameter or a new CalcCore object is
        /// created and returned if the passed shared pointer is empty.
        /// @param[inout] rootImagerPtr shared pointer to CalcCore object to update or create (if empty shared pointer is passed)
        /// @param[in] wu work unit to process
        /// @param[in] lastcycle if this parameter is true and itsWriteGrids is true as well, the grids are extracted into root imager 
        /// for writing later on. We only do this in the last major cycle, hence the name. In the central solver case this option has
        /// no effect.
        void processOneWorkUnit(boost::shared_ptr<CalcCore> &rootImagerPtr, const cp::ContinuumWorkUnit &wu, bool lastcycle) const;

        /// @brief helper method to accumulate uv-weights for a single work unit
        /// @details It is alalogous to processOneWorkUnit, but is used in the section responsible for traditional weighting.
        /// The resulting sample density grid (wrapped into an NE-like object) is either added to the root imager passed as a parameter or
        /// a new CalcCore object is created and returned if the passed shared pointer is null.
        /// @param[inout] rootImagerPtr shared pointer to CalcCore object to update or create (if empty shared pointer is passed)
        /// @param[in] wu work unit to process
        void accumulateUVWeightsForOneWorkUnit(boost::shared_ptr<CalcCore> &rootImagerPtr, const cp::ContinuumWorkUnit &wu) const;


        /// @brief helper method to create and configure work and (optionally) root imagers
        /// @details This method encapsulates the part of single work unit processing where the work and root imagers are created. 
        /// Using two imager objects is a bit of the technical debt - ideally, one has to merge normal equations or models directly.
        /// But this is deeply in the design of this application and left as is for now. Normally, all gridding of data is taken place
        /// in the 'work imager' and the results are merged into 'root imager' when ready. If the root imager is not defined, the work imager
        /// becomes one for the subsequent data merge. However, the logic of creating these imagers depends on the mode (e.g. itsUpdateDir, itsLocalSolver).
        /// This method encapsulates all the logic, so it can be repeated easily for both normal gridding and sample grid calculation for traditional weighting.
        /// @note in the case of central solver (and no joint deconvolution), this method expects to receive model from the master rank if root imager is not
        /// defined (i.e. when the work imager created inside this method would become a new root imager for future accumulation)
        /// @param[in] wu work unit to work with (parameters like frequency, channel and the dataset may be used). Note, access to data currently happens
        /// in the joint imaging mode where the full image is created through dummy iteration hack.
        /// @param[inout] rootImagerPtr shared pointer to CalcCore object to use as the root imager (if empty, it is created in the joint imaging mode)
        /// @return shared pointer to CalcCore object to be used as a work imager for the given work unit
        boost::shared_ptr<CalcCore> createImagers(const cp::ContinuumWorkUnit &wu, boost::shared_ptr<CalcCore> &rootImagerPtr) const;

        /// @brief helper method to perform minor cycle activities
        /// @details This method encapsulates running the solver at the conclusion of each major cycle and
        /// associated data transfers, if necessary. It can be viewed at the place where the minor cycle is
        /// performed in the case of the local solver (or the interface part, if the master perofms it as
        /// it happens in the continuum mode).
        /// @param[in] rootImagerPtr shared pointer to the CalcCore object containing the result of the current 
        ///                          major cycle for the given worker.
        /// @param[in] haveMoreMajorCycles flag that more major cycles are to be done (subject to thresholds). In the
        ///                          central solver mode we expect to receive the new model in this flag is true or
        ///                          perform new solution ourselves in the case of the local solver.
        /// @return true if major cycles have to be terminated due to thresholds being reached
        /// @note empty rootImagerPtr causes an exception. This may happen if no data are processed.
        bool runMinorCycleSolver(const boost::shared_ptr<CalcCore> &rootImagerPtr, bool haveMoreMajorCycles) const;

        // My Advisor
        boost::shared_ptr<synthesis::AdviseDI> itsAdvisor;

        /// @brief the work units
        WorkUnitContainer itsWorkUnits;

        /// @brief whether preconditioning has been requested
        /// @details It is populated in the constructor based on the parset using doingPreconditioning method
        const bool itsDoingPreconditioning;

        // Process a workunit
        void preProcessWorkUnit(ContinuumWorkUnit& wu);

        // Compress all continuous channel allocations into individual workunits
        void compressWorkUnits();

        //For all workunits .... process
        void processChannels();
        // temporary keep the old version
        void processChannelsOld();

        // Setup the image specified in parset and add it to the Params instance.
        void setupImage(const askap::scimath::Params::ShPtr& params,
                    double channelFrequency, bool shapeOveride = false) const;

        /// @brief Load an image plane into a parameter
        /// @details Load the specified channel from the image specified in parset and add it to the Params instance.
        /// @params[in] params shared pointer to params object with image parameters
        /// @params[in] channel, the (global) channel number for the image plane to load
        void loadImage(const askap::scimath::Params::ShPtr& params, int channel) const;

        /// @brief Load an image model into a parameter
        /// @details Evaluate the MFS model at the specified frequency and add it to the Params instance.
        /// @params[in] params shared pointer to params object with image parameters
        /// @params[in] freq, the frequency to evaluate the model at
        /// @params[in] channel, the (global) channel number for the image plane to fill with the model
        void loadImageFromMFSModel(const askap::scimath::Params::ShPtr& params, double freq, int channel) const;

        /// @brief Load a starting model if required, or setup empty model
        /// @details Initialise the starting model image from a model cube, MFS model, or empty image
        /// @params[in] channel, the (global) channel number for the image plane to fill with the model
        /// @params[in] frequency, the frequency to evaluate the MFS model at, if supplied
        /// @return bool, true if a starting model was loaded, false if zero model was set up
        bool loadStartingModel(const askap::scimath::Params::ShPtr& params, uInt channel, double frequency) const;

        /// @brief Mask image with NaNs when corresponding weight is below cutoff
        /// @details When using jont deconvolution we want to mask the output mosaic below
        /// a weights image cutoff level. The cutoff is specified as a fraction of the maximum weight
        /// using the solver.Clean.tolerance parameter to match how the Clean solvers handle this.
        /// @params[inout] arr, the image plane to be masked
        /// @params[in] wts, the weights plane
        void maskOutput(casacore::Array<float>& arr, const casacore::Array<float>& wts);

        // Root Parameter set good for information common to all workUnits
        LOFAR::ParameterSet& itsParset;

        // Communications class
        CubeComms& itsComms;

        // statistics
        StatReporter& itsStats;

        /// @brief shared pointer to data source manager reponsible for dealing with measurement sets
        boost::shared_ptr<DataSourceManager> itsDSM;

        // ID of the master process
        static const int itsMaster = 0;

        // the basechannel number assigned to this worker
        unsigned int itsBaseChannel;

        // the baseFrequency associated with this channel
        double itsBaseFrequency;

        // the baseFrequency associated with the cube if being built
        double itsBaseCubeFrequency;

        // the global channel associated with this part of the cube
        int itsBaseCubeGlobalChannel;

        // the number of channels in this cube (if writer)
        int itsNChanCube;

        /// @brief true if solver is run locally (spectral line mode), 
        /// false for central solver (continuum). 
        const bool itsLocalSolver;

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

        void handleImageParams(askap::scimath::Params::ShPtr params, unsigned int chan);

        void copyModel(askap::scimath::Params::ShPtr SourceParams, askap::scimath::Params::ShPtr SinkParams) const;

        void initialiseBeamLog(const unsigned int numChannels);
        void recordBeam(const askap::scimath::Axes &axes, const unsigned int globalChannel);

        std::map<unsigned int, casacore::Vector<casacore::Quantum<double> > > itsBeamList;
        unsigned int itsBeamReferenceChannel;
        void logBeamInfo() const;

        void initialiseWeightsLog(const unsigned int numChannels);
        void recordWeight(float wt, const unsigned int globalChannel);
        std::map<unsigned int, float> itsWeightsList;
        void logWeightsInfo() const;

        /// @brief read starting model cube
        const bool itsReadStartingModelCube;

        /// @brief Do we want a restored image?
        const bool itsRestore;

        /// @brief Do we want a residual image
        const bool itsWriteResidual;

        /// @brief write 'raw', unnormalised, natural weight psf
        const bool itsWritePsfRaw;

        /// @brief write normalised, preconditioned psf
        const bool itsWritePsfImage;

        /// @brief write weights image
        const bool itsWriteWtImage;

        /// @brief write a weights log
        const bool itsWriteWtLog;

        /// @brief write out the (clean) model image
        const bool itsWriteModelImage;

        /// @brief write out the gridded data, pcf and psf
        const bool itsWriteGrids;

        /// @brief grid image type
        const std::string itsGridType;

        /// @brief write out the grids with UV coordinate grid
        const bool itsGridCoordUV;

        /// @brief write out the FFT of the grids
        const bool itsGridFFT;

        /// @brief the number of rank that can write to the cube
        const int itsNumWriters;

        /// @brief updatedirection option (switching on joint gridding)
        const bool itsUpdateDir;

        /// @brief do we have an MFS starting model for spectral imaging
        const bool itsMFSStartingModel;

        /// @brief masking level (as fraction of peak weight) for mosaic output
        const float itsMaskLevel;

        /// @brief do we mask mosaic output using the weight image?
        const bool itsMaskOutput;

        /// @brief shared pointer to the uv-weight calculator object
        /// @details it can also be used as a flag that the sample density grid is needed (and that traditional weighting is done)
        /// If defined, the sample density grid needs to be constructed and traditional weighting
        /// needs to be done. Note, this serves the same purpose as the similarly named data member in ImagerParallel class. However, 
        /// it is handy to get this object available in this class before the imager is constructed (because we can have many 
        /// imagers and it makes the code messy). Ideally, we need to clean up some technical debt and probably avoid having 
        /// unnecessary responsibilities assign to the imager class.
        const boost::shared_ptr<IUVWeightCalculator> itsUVWeightCalculator;
};

};
};

#endif
