/// @file AmplitudeFlagger.h
///
/// @copyright (c) 2013,2014 CSIRO
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

#ifndef ASKAP_SYNTHESIS_AMPLITUDEFLAGGER_H
#define ASKAP_SYNTHESIS_AMPLITUDEFLAGGER_H

// ASKAPsoft includes
#include "Common/ParameterSet.h"
#include "casacore/measures/Measures/Stokes.h"

// Local package includes
#include "askap/flagging/IFlagger.h"
#include "askap/flagging/FlaggingStats.h"

namespace askap {
namespace synthesis {

/// @brief Applies flagging based on amplitude thresholding.
class AmplitudeFlagger : public IFlagger {

    public:

        /// @brief Constructs zero or more instances of the AmplitudeFlagger.
        /// The flagger is responsible for reading the "parset" and constructing
        /// zero or more instances of itself, depending on the configuration.
        ///
        /// @throw AskapError   If an upper or lower threshold is not specified
        ///                     in the parset.
        static std::vector<std::shared_ptr<IFlagger> > build(
                const LOFAR::ParameterSet& parset);

        /// @brief Constructor
        /// @throw AskapError   If an upper or lower threshold is not specified
        ///                     in the parset.
        AmplitudeFlagger(const LOFAR::ParameterSet& parset);

        /// @see IFlagger::processRows()
        virtual void processRows(const accessors::IDataSharedIter& di,
            const casacore::Vector<casacore::Bool>& rowFlag,
            const casacore::uInt pass, const bool dryRun) override;

        /// @see IFlagger::stats()
        virtual FlaggingStats stats(void) const override;

        /// @see IFlagger::processingRequired()
        virtual casacore::Bool processingRequired(const casacore::uInt pass) const override;

    private:
        // load and log relevant parset parameters
        void loadParset(const LOFAR::ParameterSet& parset);
        void logParsetSummary(const LOFAR::ParameterSet& parset) const;

        // Flagging statistics
        FlaggingStats itsStats;

        // True if an upper amplitude limit has been set, otherwise false
        bool itsHasHighLimit;
        // True if lower amplitude limit has been set, otherwise false
        bool itsHasLowLimit;
        // The upper amplitude limit
        casacore::Float itsHighLimit;
        // The lower amplitude limit
        casacore::Float itsLowLimit;

        // Automatically set either of these limits that are unset
        bool itsAutoThresholds;
        // sigma multiplier used to set cutoffs
        casacore::Float itsThresholdFactor;

        // Generate averaged spectra and search these for peaks to flag
        bool itsIntegrateSpectra;
        // sigma multiplier used to set cutoffs
        casacore::Float itsSpectraFactor;

        // Generate averaged time series and search these for peaks to flag
        bool itsIntegrateTimes;
        // sigma multiplier used to set cutoffs
        casacore::Float itsTimesFactor;

        // When integrating, do not separate spectra based on baseline, etc.
        bool itsAveAll;
        // When integrating, do separate spectra for different polarisations
        bool itsAveAllButPol;
        // When integrating, do separate spectra for different beams
        bool itsAveAllButBeam;

        // When integrating, used to limit flag generation to a single call to
        // "processRow"
        bool itsAverageFlagsAreReady;

        // The set of correlation products for which these flagging rules should
        // be applied. An empty list means apply to all correlation products.
        std::set<casacore::Stokes::StokesTypes> itsStokes;

        // Maps of accumulation std::vectors for averaging spectra and generating flags
        std::map<rowKey, casacore::Vector<casacore::Double> > itsAveSpectra;
        std::map<rowKey, casacore::Vector<casacore::Bool> > itsMaskSpectra;
        std::map<rowKey, casacore::Vector<casacore::Int> > itsCountSpectra;

        // Maps of accumulation std::vectors for averaging time series and generating flags
        std::map<rowKey, casacore::Vector<casacore::Float> > itsAveTimes;
        std::map<rowKey, casacore::Vector<casacore::Bool> > itsMaskTimes;
        std::map<rowKey, casacore::Int> itsCountTimes;

        // Generate a key for a given row and polarisation
        rowKey getRowKey(const accessors::IDataSharedIter& di,
            const casacore::uInt row, const casacore::uInt corr) const;

        // Functions to handle accumulation std::vectors and indices
        void updateTimeVectors(const rowKey &key, const casacore::uInt pass);
        void initSpectrumVectors(const rowKey &key, const casacore::IPosition &shape);

        // Calculate the median, the interquartile range, the min and the max
        // of a masked array
        casacore::Vector<casacore::Float>
            getRobustStats(casacore::MaskedArray<casacore::Float> maskedAmplitudes);

        // Set flags based on integrated quantities
        void setFlagsFromIntegrations(void);

};

}
}

#endif
