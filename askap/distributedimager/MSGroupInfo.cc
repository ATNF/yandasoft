/// @file MSGroupInfo.cc
///
/// @copyright (c) 2013 CSIRO
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
#include <askap/distributedimager/MSGroupInfo.h>

// Include package level header file
#include <askap/askap_synthesis.h>

// System includes
#include <string>
#include <vector>
#include <numeric>
#include <limits>
#include <iomanip>

// ASKAPsoft includes
#include <askap/AskapError.h>
#include <askap/AskapLogging.h>
#include <casacore/casa/Quanta/Unit.h>
#include <casacore/measures/Measures/MFrequency.h>
#include <askap/dataaccess/IConstDataSource.h>
#include <askap/dataaccess/TableConstDataSource.h>
#include <askap/dataaccess/IConstDataIterator.h>
#include <askap/dataaccess/IDataConverter.h>
#include <askap/dataaccess/IDataSelector.h>
#include <askap/dataaccess/IDataIterator.h>
#include <askap/dataaccess/SharedIter.h>

ASKAP_LOGGER(logger, ".MSGroupInfo");

using namespace askap::cp;
using namespace casa;
using namespace std;

MSGroupInfo::MSGroupInfo() : itsTotalNumChannels(0)
{
}

MSGroupInfo::MSGroupInfo(const std::vector<std::string>& ms)
{
    // NOTE: This function makes the assumption that each iteration will have
    // the same number of channels. This may not be true, but reading through the
    // entire dataset to validate this assumption is going to be too slow.

    // Pre-conditions
    ASKAPCHECK(!ms.empty(), "No measurement sets specified");

    // Frequency for each channel
    vector<casacore::Quantity> freqinfo;

    // Iterate over measurement sets
    const casacore::Unit frequnit("Hz");
    for (size_t i = 0; i < ms.size(); ++i) {
        // Open dataset and get access to the first row
        askap::accessors::TableConstDataSource ds(ms[i]);
        askap::accessors::IDataSelectorPtr sel = ds.createSelector();
        askap::accessors::IDataConverterPtr conv = ds.createConverter();
        conv->setFrequencyFrame(casacore::MFrequency::Ref(casacore::MFrequency::TOPO), frequnit);
        conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));
        const askap::accessors::IConstDataSharedIter it = ds.createConstIterator(sel, conv);

        // Extract # of channels
        itsNumChannels.push_back(it->nChannel());

        // Extract frequency info
        for (size_t i = 0; i < it->nChannel(); ++i) {
            freqinfo.push_back(casacore::Quantity(it->frequency()(i), frequnit));
        }
    }

    // Calculate aggregate number of channels
    itsTotalNumChannels = std::accumulate(itsNumChannels.begin(), itsNumChannels.end(), 0);
    ASKAPCHECK(itsTotalNumChannels > 1,
               "Only one channel, need at least two to calculate frequency increment");

    // Calcuate first freq and freq increment
    itsFirstFreq = freqinfo[0];
    const Quantity secondfreq = freqinfo[1];
    itsFreqInc = secondfreq - itsFirstFreq;

    // Ensure the frequency increment is valid for all frequencies
    Quantity fi = itsFirstFreq;
    for (size_t i = 0; i < freqinfo.size(); ++i) {
        const double expected = fi.getValue(frequnit);
        const double actual = freqinfo[i].getValue(frequnit);
        const double tolerance = 1; // i.e. one Hz
        if (abs(expected - actual) > tolerance) {
            ASKAPLOG_ERROR_STR(logger,
                               "Expected: " << setprecision(16) << expected <<
                               ", Actual: " << actual);
            ASKAPTHROW(AskapError, "Non-constant frequency increment is not supported");
        }
        fi += itsFreqInc;
    }
}

MSGroupInfo::~MSGroupInfo()
{
}

casacore::uInt MSGroupInfo::getNumChannels(const casacore::uInt n) const
{
    return itsNumChannels.at(n);
}

casacore::uInt MSGroupInfo::getTotalNumChannels() const
{
    return itsTotalNumChannels;
}

casacore::Quantity MSGroupInfo::getFirstFreq() const
{
    return itsFirstFreq;
}

casacore::Quantity MSGroupInfo::getFreqInc() const
{
    return itsFreqInc;
}
