/// @file StokesVFlagger.cc
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

// Include own header file first
#include "StokesVFlagger.h"

// Include package level header file
#include "askap/askap_synthesis.h"

// ASKAPsoft includes
#include "askap/askap/AskapLogging.h"
#include "askap/dataaccess/IFlagDataAccessor.h"
#include "askap/dataaccess/TableConstDataIterator.h"
#include "askap/scimath/utils/PolConverter.h"

// Local package includes
#include "askap/flagging/FlaggingStats.h"

ASKAP_LOGGER(logger, ".StokesVFlagger");

using namespace std;
using namespace askap;
using namespace casacore;
using namespace askap::synthesis;
using namespace askap::accessors;

vector<shared_ptr<IFlagger> > StokesVFlagger::build(
        const LOFAR::ParameterSet& parset)
{
    vector< shared_ptr<IFlagger> > flaggers;
    const string key = "stokesv_flagger.enable";
    if (parset.isDefined(key) && parset.getBool(key)) {
        const LOFAR::ParameterSet subset = parset.makeSubset("stokesv_flagger.");

        const float threshold = subset.getFloat("threshold", 5.0);
        const bool robustStatistics = subset.getBool("useRobustStatistics", false);
        const bool integrateSpectra = subset.getBool("integrateSpectra", false);
        const float spectraThreshold = subset.getFloat("integrateSpectra.threshold", 5.0);
        const bool integrateTimes = subset.getBool("integrateTimes", false);
        const float timesThreshold = subset.getFloat("integrateTimes.threshold", 5.0);
        const float highLimit = subset.getFloat("high", 0.0);

        ASKAPLOG_INFO_STR(logger, "Parameter Summary:");
        if (highLimit > 0) {
            ASKAPLOG_INFO_STR(logger, "Searching for outliers with absolute value > "<<highLimit<<" Jy");
        }
        ASKAPLOG_INFO_STR(logger, "Searching for outliers with a "<<threshold<<"-sigma cutoff");
        if (robustStatistics) {
            ASKAPLOG_INFO_STR(logger, "Using robust statistics");
        }
        if (integrateSpectra) {
            ASKAPLOG_INFO_STR(logger,
                "Searching for outliers in integrated spectra with a "
                <<spectraThreshold<<"-sigma cutoff");
        }
        if (integrateTimes) {
            ASKAPLOG_INFO_STR(logger,
                "Searching for outliers in integrated time series with a "
                <<timesThreshold<<"-sigma cutoff");
        }

        flaggers.push_back(shared_ptr<IFlagger>
            (new StokesVFlagger(threshold,robustStatistics,integrateSpectra,
                                spectraThreshold,integrateTimes,timesThreshold,
                                highLimit)));
    }
    return flaggers;
}

StokesVFlagger::StokesVFlagger(float threshold, bool robustStatistics,
                                bool integrateSpectra, float spectraThreshold,
                                bool integrateTimes, float timesThreshold,
                                float highLimit)
    : itsStats("StokesVFlagger"),
      itsThreshold(threshold), itsRobustStatistics(robustStatistics),
      itsIntegrateSpectra(integrateSpectra), itsSpectraThreshold(spectraThreshold),
      itsIntegrateTimes(integrateTimes), itsTimesThreshold(timesThreshold),
      itsAverageFlagsAreReady(true), itsHighLimit(highLimit)
{
    ASKAPCHECK(itsThreshold > 0.0, "Threshold must be greater than zero");
}

FlaggingStats StokesVFlagger::stats(void) const
{
    return itsStats;
}

casacore::Bool StokesVFlagger::processingRequired(const casacore::uInt pass) const
{
    if (itsIntegrateSpectra || itsIntegrateTimes) {
        return (pass<2);
    } else {
        return (pass<1);
    }
}

void StokesVFlagger::processRows(const IDataSharedIter& di,
    const casacore::Vector<bool>& rowFlag,
    const casacore::uInt pass, const bool dryRun)
{
    // Convert data to Stokes V (imag(data(2,i))-imag(data(3,i)))
    const casacore::Vector<casacore::Stokes::StokesTypes> target(1, Stokes::V);
    scimath::PolConverter polConv(di->stokes(),target);
    const Cube<casacore::Complex>& data = di->visibility();
    Cube<casacore::Bool> flags = di->flag();
    const casacore::uInt nRow = di->nRow();
    const casacore::uInt nPol = di->nPol();
    const casacore::uInt nChan = di->nChannel();
    // does it help to keep this around?
    casacore::Matrix<casacore::Complex> vdata(nChan, nRow, 0);
    casacore::Vector<casacore::Complex> in(nPol);
    casacore::Vector<casacore::Complex> out(1);

    // Only need to write out the flag matrix if it was updated
    bool wasUpdated = false;

    for (casacore::uInt row = 0; row < nRow;  row++) {

        //do the conversion to Stokes V using PolConverter
        bool allFlagged = true;
        std::vector<casacore::Float> tmpamps;
        for (size_t chan = 0; chan < nChan; ++chan) {
            bool anyFlagged = false;
            for (casacore::uInt pol=0; pol < nPol; pol++) {
                if (flags(pol, chan, row)) anyFlagged = true;
            }
            if (!anyFlagged) {
                if (pass == 0) {
                    in(0) = data(0,chan,row);
                    in(1) = data(1,chan,row);
                    in(2) = data(2,chan,row);
                    in(3) = data(3,chan,row);
                    polConv.convert(out,in);
                    vdata(chan, row) = out(0);
                }
                allFlagged = false;
            }
        }

        // Now go through and flag based on absolute levels
        // also build a vector with the remaining amplitudes
        if (pass == 0) {
            for (size_t chan = 0; chan < nChan; ++chan) {
                const casacore::Float amp = abs(vdata(chan,row));
                // Apply threshold based flagging
                if (itsHighLimit > 0 && amp > itsHighLimit) {
                    for (casacore::uInt pol = 0; pol < flags.nrow(); ++pol) {
                        if (flags(pol, chan, row)) {
                            itsStats.visAlreadyFlagged++;
                            continue;
                        }
                        flags(pol, chan, row) = true;
                        wasUpdated = true;
                        itsStats.visFlagged++;
                    }
                } else {
                    tmpamps.push_back(amp);
                }
            }
            if (tmpamps.size()==0) {
                allFlagged = true;
            }
        }


        // normalise averages and search them for peaks to flag
        if ( !itsAverageFlagsAreReady && (pass==1) ) {
            ASKAPLOG_INFO_STR(logger, "Finalising averages at the start of pass "
                <<pass+1);
            setFlagsFromIntegrations();
        }

        // return a key that indicates which integration this row is in
        const rowKey key = getRowKey(di, row);
        // update a counter for this row and the storage vectors
        // do it before any processing that is dependent on "pass"
        if ( itsIntegrateTimes ) {
            updateTimeVectors(key, pass);
        }

        // if this is the first instance of this key, initialise storage vectors
        if ( itsIntegrateSpectra && (pass==0) &&
               (itsAveSpectra.find(key) == itsAveSpectra.end()) ) {
            initSpectrumVectors(key, IPosition(1,nChan));
        }

        // If all visibilities are flagged, nothing to do
        if (allFlagged) continue;

        bool wasUpdatedRow = false;

        if ( pass==0 ) {

            // Convert to a casacore::Vector so we can use ArrayMath functions
            // to determine the mean and stddev
            casacore::Vector<casacore::Float> amps(tmpamps);

            // Flag all correlations where the Stokes V product
            // is greater than the threshold
            casacore::Float sigma, avg;
            if (itsRobustStatistics) {
                const casacore::Vector<casacore::Float> statsVector = getRobustStats(amps);
                avg = statsVector[0];
                sigma = statsVector[1];
                // if min and max are bounded, they all are.
                // so skip if there is no other reason to loop over frequencies
                if ((statsVector[2] >= (avg - (sigma * itsThreshold))) &&
                    (statsVector[3] <= (avg + (sigma * itsThreshold))) &&
                    !itsIntegrateSpectra && !itsIntegrateTimes) {
                    continue;
                }
            }
            else {
                avg = mean(amps);
                sigma = stddev(amps);
            }

            // If stokes-v can't be formed due to lack of the necessary input products
            // then vdata will contain all zeros. In this case, no flagging can be done.
            const casacore::Float epsilon = std::numeric_limits<casacore::Float>::epsilon();
            if (near(sigma, 0.0, epsilon) && near(avg, 0.0, epsilon)) {
                continue;
            }

            // Apply threshold based flagging and accumulate any averages
            // only need these if itsIntegrateTimes
            casacore::Double aveTime = 0.0;
            casacore::uInt countTime = 0;
            for (size_t chan = 0; chan < nChan; ++chan) {
                const casacore::Float amp = abs(vdata(chan, row));
                // Apply threshold based flagging
                if (amp > (avg + (sigma * itsThreshold))) {
                    for (casacore::uInt pol = 0; pol < nPol; ++pol) {
                        if (flags(pol, chan, row)) {
                            // no double counting
                            if (itsHighLimit == 0 || amp < itsHighLimit) {
                                itsStats.visAlreadyFlagged++;
                            }
                            continue;
                        }
                        flags(row, chan, pol) = true;
                        wasUpdatedRow = true;
                        itsStats.visFlagged++;
                    }
                }
                // Accumulate any averages
                else if ( itsIntegrateSpectra || itsIntegrateTimes ) {
                    if ( itsIntegrateSpectra ) {
                        if (itsAverageFlagsAreReady) ASKAPLOG_INFO_STR(logger,"integrate spectra ");
                        // do spectra integration
                        itsAveSpectra[key][chan] += amp;
                        itsCountSpectra[key][chan]++;
                        itsAverageFlagsAreReady = casacore::False;
                    }
                    if ( itsIntegrateTimes ) {
                        // do time-series integration
                        aveTime += amp;
                        countTime++;
                    }
                }
            }
            if ( itsIntegrateTimes ) {
                // normalise this average
                if ( countTime>0 ) {
                    itsAveTimes[key][itsCountTimes[key]] =
                        aveTime/casacore::Double(countTime);
                    itsMaskTimes[key][itsCountTimes[key]] = casacore::True;
                    itsAverageFlagsAreReady = casacore::False;
                } else {
                    itsMaskTimes[key][itsCountTimes[key]] = casacore::False;
                }
            }
        }
        else if ( (pass==1) &&  ( itsIntegrateSpectra || itsIntegrateTimes ) ) {
            // only flag unflagged data, so that new flags can be counted.
            // "flags" is true for flags, "mask*" are false for flags
            bool rowFlagged = false;
            if ( itsIntegrateTimes ) {
                // apply itsMaskTimes flags. Could just use flagRow,
                // but not sure that all applications support flagRow
                if ( !itsMaskTimes[key][itsCountTimes[key]] ) {
                    rowFlagged = true;
                    itsStats.rowsFlagged++;
                    for (size_t chan = 0; chan < nChan; ++chan) {
                        for (casacore::uInt pol = 0; pol < nPol; ++pol) {
                            if (!flags(pol, chan, row)) {
                                flags(pol, chan, row) = true;
                                wasUpdatedRow = true;
                                itsStats.visFlagged++;
                            }
                        }
                    }
                }
            }
            // apply itsIntegrateSpectra flags
            if ( itsIntegrateSpectra && !rowFlagged ) {
                for (size_t chan = 0; chan < nChan; ++chan) {
                    if ( !itsMaskSpectra[key][chan] ) {
                        for (casacore::uInt pol = 0; pol < nPol; ++pol) {
                            if (!flags(pol, chan, row)) {
                                flags(pol, chan, row) = true;
                                wasUpdatedRow = true;
                                itsStats.visFlagged++;
                            }
                        }
                    }
                }
            }
        }
        if (wasUpdatedRow) wasUpdated = true;
    }
    if (wasUpdated && !dryRun) {
        IFlagDataAccessor &fda=dynamic_cast<IFlagDataAccessor&>(*di);
        fda.rwFlag() = flags;
    }
}


// return the median, the interquartile range, and the min/max of a masked array
casacore::Vector<casacore::Float>StokesVFlagger::getRobustStats(
    casacore::MaskedArray<casacore::Float> maskedAmplitudes)
{
    // extract all of the unflagged amplitudes
    casacore::Vector<casacore::Float> amplitudes = maskedAmplitudes.getCompressedArray();

    // return with zeros if all of the data are flagged
    if (amplitudes.nelements() == 0) {
        const casacore::Vector<casacore::Float> statsVector(4,0.0);
        return(statsVector);
    }
    return getRobustStats(amplitudes);
}

// return the median, the interquartile range, and the min/max of a non-masked array
casacore::Vector<casacore::Float>StokesVFlagger::getRobustStats(
    casacore::Vector<casacore::Float> amplitudes)
{
    casacore::Float minVal, maxVal;
    casacore::minMax(minVal,maxVal,amplitudes);

    // Now find median and IQR
    // Use the fact that nth_element does a partial sort:
    // all elements before the selected element will be smaller
    const casacore::uInt n = amplitudes.nelements();
    std::vector<casacore::Float> vamp = amplitudes.tovector();
    const casacore::uInt Q1 = n / 4;
    const casacore::uInt Q2 = n / 2;
    const casacore::uInt Q3 = 3 * n /4;
    std::nth_element(vamp.begin(),        vamp.begin() + Q2, vamp.end());
    std::nth_element(vamp.begin(),        vamp.begin() + Q1, vamp.begin() + Q2);
    std::nth_element(vamp.begin() + Q2+1, vamp.begin() + Q3, vamp.end());

    casacore::Vector<casacore::Float> statsVector(4);
    statsVector[0] = vamp[Q2]; // median
    statsVector[1] = (vamp[Q3] - vamp[Q1]) / 1.34896; // sigma from IQR
    statsVector[2] = minVal; // min
    statsVector[3] = maxVal; // max
    return(statsVector);
}

// Generate a key for a given row and polarisation
rowKey StokesVFlagger::getRowKey(const IDataSharedIter& di,
    const casacore::uInt row) const
{
    auto tdi = di.dynamicCast<TableConstDataIterator>();
    const casacore::Int field = tdi->currentFieldID();
    const casacore::Int feed1 = di->feed1()(row);
    const casacore::Int ant1  = di->antenna1()(row);
    const casacore::Int ant2  = di->antenna2()(row);
#ifdef TUPLE_INDEX
    // casacore::Int feed2 = di->feed2()(row);
    // looking for outliers in a single polarisation, so set the corr key to zero
    return std::make_tuple(field, feed1, /*feed2,*/ ant1, ant2, 0);
#else
    const casacore::uLong nant = tdi->subtableInfo().getAntenna().getNumberOfAntennas();
    const casacore::uLong nfeed = 36; // too hard to find number of feeds through iterator/accessor
    return feed1+nfeed*(ant1+nant*(ant2+nant*field));
#endif
}

void StokesVFlagger::updateTimeVectors(const rowKey &key, const casacore::uInt pass)
{
    if (itsCountTimes.find(key) == itsCountTimes.end()) {
        itsCountTimes[key] = 0; // init counter for this key
    }
    else {
        itsCountTimes[key]++;
    }
    if ( pass==0 ) {
        const int index = itsCountTimes[key];
        // allocate 10 at a time to reduce reallocations with copy
        if (itsAveTimes[key].size() <= index) {
            itsAveTimes[key].resize(index + 10, casacore::True);
            itsMaskTimes[key].resize(index + 10,casacore::True);
            for (int i = index; i < index + 10; i++) {
                itsAveTimes[key][i] = 0;
                itsMaskTimes[key][i] = casacore::False;
            }
        }
     }
}

void StokesVFlagger::initSpectrumVectors(const rowKey &key, const casacore::IPosition &shape)
{
    itsAveSpectra[key].resize(shape);
    itsAveSpectra[key].set(0.0);
    itsCountSpectra[key].resize(shape);
    itsCountSpectra[key].set(0);
    itsMaskSpectra[key].resize(shape);
    itsMaskSpectra[key].set(casacore::False);
}

// Set flags based on integrated quantities
void StokesVFlagger::setFlagsFromIntegrations(void)
{
    if ( itsIntegrateSpectra ) {

        for (auto el : itsAveSpectra) {

            // get the spectra
            casacore::Vector<casacore::Float> aveSpectrum(el.second.shape());
            casacore::Vector<casacore::Int>& countSpectrum = itsCountSpectra[el.first];
            casacore::Vector<casacore::Bool>& maskSpectrum = itsMaskSpectra[el.first];
            //std::vector<casacore::Float> tmpamps; // use instead of MaskedArray?

            for (size_t chan = 0; chan < aveSpectrum.size(); ++chan) {
                if (countSpectrum[chan]>0) {
                    aveSpectrum[chan] = el.second[chan] /
                                        casacore::Double(countSpectrum[chan]);
                    //tmpamps.push_back(aveSpectrum[chan]);
                    countSpectrum[chan] = 1;
                    maskSpectrum[chan] = casacore::True;
                }
                else {
                    maskSpectrum[chan] = casacore::False;
                }
            }

            // generate the flagging stats. could fill the unflagged spectrum
            // directly in the preceding loop, but the full vector is needed below
            casacore::MaskedArray<casacore::Float>
                maskedAmplitudes(aveSpectrum, maskSpectrum);
            const casacore::Vector<casacore::Float> statsVector =
                getRobustStats(maskedAmplitudes);
            const casacore::Float median = statsVector[0];
            const casacore::Float sigma_IQR = statsVector[1];

            // check min and max relative to thresholds.
            // do not loop over data again if all unflagged channels are good
            if ((statsVector[2] < median-itsSpectraThreshold*sigma_IQR) ||
                (statsVector[3] > median+itsSpectraThreshold*sigma_IQR)) {

                for (size_t chan = 0; chan < aveSpectrum.size(); ++chan) {
                    if (maskSpectrum[chan]==casacore::False) continue;
                    if ((aveSpectrum[chan]<median-itsSpectraThreshold*sigma_IQR) ||
                        (aveSpectrum[chan]>median+itsSpectraThreshold*sigma_IQR)) {
                        maskSpectrum[chan]=casacore::False;
                    }
                }

            }

        }

    }

    if ( itsIntegrateTimes ) {

        for (auto el : itsAveTimes) {
            // reset the counter for this key
            itsCountTimes[el.first] = -1;

            // get the spectra
            casacore::Vector<casacore::Float>& aveTime = el.second;
            casacore::Vector<casacore::Bool>& maskTime = itsMaskTimes[el.first];

            // generate the flagging stats
            casacore::MaskedArray<casacore::Float> maskedAmplitudes(aveTime, maskTime);
            const casacore::Vector<casacore::Float>
                statsVector = getRobustStats(maskedAmplitudes);
            const casacore::Float median = statsVector[0];
            const casacore::Float sigma_IQR = statsVector[1];
            // check min and max relative to thresholds.
            // do not loop over data again if all unflagged times are good
            if ((statsVector[2] < median-itsTimesThreshold*sigma_IQR) ||
                (statsVector[3] > median+itsTimesThreshold*sigma_IQR)) {

                for (size_t t = 0; t < aveTime.size(); ++t) {
                    if (maskTime[t] == casacore::False) continue;
                    if ((aveTime[t] < median-itsTimesThreshold*sigma_IQR) ||
                        (aveTime[t] > median+itsTimesThreshold*sigma_IQR)) {
                        maskTime[t] = casacore::False;
                    }
                }

            }

        }

    }

    itsAverageFlagsAreReady = casacore::True;

}
