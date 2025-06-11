/// @file StokesVFlagger.h
///
/// @copyright (c) 2012-2014 CSIRO
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

#ifndef ASKAP_SYNTHESIS_STOKESVFLAGGER_H
#define ASKAP_SYNTHESIS_STOKESVFLAGGER_H

// ASKAPsoft includes
#include "Common/ParameterSet.h"

// Local package includes
#include "askap/flagging/IFlagger.h"
#include "askap/flagging/FlaggingStats.h"

namespace askap {
namespace synthesis {

/// @brief Performs flagging based on Stokes-V thresholding. For each row
/// the mean and standard deviation for all Stokes-V correlations (i.e. all
/// channels within a given row). Then, where the stokes-V correlation exceeds
/// the average plus (stddev * threshold) all correlations for that channel in
/// that row will be flagged.
///
/// The one parameter that is read from the parset passed to the constructor is
/// "threshold". To flag at the five-sigma point specify a valud of "5.0".
class StokesVFlagger : public IFlagger {

    public:

        /// @brief Constructs zero or more instances of the StokesVFlagger.
        /// The flagger is responsible for reading the "parset" and constructing
        /// zero or more instances of itself, depending on the configuration.
        static std::vector<std::shared_ptr<IFlagger> > build(
                const LOFAR::ParameterSet& parset);

        /// @brief Constructor
        StokesVFlagger(float threshold, bool robustStatistics,
                       bool integrateSpectra, float spectraThreshold,
                       bool integrateTimes, float timesThreshold,
                       float highLimit);

        /// @see IFlagger::processRows()
        virtual void processRows(const accessors::IDataSharedIter& di,
            const casacore::Vector<bool>& rowFlag,
            const casacore::uInt pass, const bool dryRun) override;

        /// @see IFlagger::stats()
        virtual FlaggingStats stats(void) const override;

        /// @see IFlagger::stats()
        virtual casacore::Bool processingRequired(const casacore::uInt pass) const override;

    private:

        // Flagging statistics
        FlaggingStats itsStats;

        // Flagging threshold (in standard deviations)
        float itsThreshold;

        // Use the median and interquartile range to estimate the mean and stddev
        bool itsRobustStatistics;

        // Generate averaged spectra and search these for peaks to flag
        bool itsIntegrateSpectra;
        // Flagging threshold
        casacore::Float itsSpectraThreshold;

        // Generate averaged time series and search these for peaks to flag
        bool itsIntegrateTimes;
        // Flagging threshold
        casacore::Float itsTimesThreshold;

        // When integrating, used to limit flag generation to a single call to
        // "processRow"
        bool itsAverageFlagsAreReady;

        // absolute limit flagging
        float itsHighLimit;

        // Calculate the median, the interquartile range, the min and the max
        // of a simple array without masking
        casacore::Vector<casacore::Float> getRobustStats(casacore::Vector<casacore::Float> amplitudes);

        // Calculate the median, the interquartile range, the min and the max
        // of a masked array
        casacore::Vector<casacore::Float>getRobustStats(casacore::MaskedArray<casacore::Float> maskedAmplitudes);

        // Generate a key for a given row and polarisation
        rowKey getRowKey(const accessors::IDataSharedIter& di, const casacore::uInt row) const;

        // Maps of accumulation std::vectors for averaging spectra and generating flags
        std::map<rowKey, casacore::Vector<casacore::Double> > itsAveSpectra;
        std::map<rowKey, casacore::Vector<casacore::Bool> > itsMaskSpectra;
        std::map<rowKey, casacore::Vector<casacore::Int> > itsCountSpectra;

        // Maps of accumulation std::vectors for averaging time series and generating flags
        std::map<rowKey, casacore::Vector<casacore::Float> > itsAveTimes;
        std::map<rowKey, casacore::Vector<casacore::Bool> > itsMaskTimes;
        std::map<rowKey, casacore::Int> itsCountTimes;

        // Functions to handle accumulation std::vectors and indices
        void updateTimeVectors(const rowKey &key, const casacore::uInt pass);
        void initSpectrumVectors(const rowKey &key, const casacore::IPosition &shape);

        // Set flags based on integrated quantities
        void setFlagsFromIntegrations(void);

};

}
}

#endif
