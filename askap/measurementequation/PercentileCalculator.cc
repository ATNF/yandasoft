/// @file
///
/// @brief A percentile value estimator
/// @details Calculate a quick estimate the value at a given percentile
/// in a single pass, without storing the data
///
/// @copyright (c) 2021 CSIRO
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
/// @author Mark Wieringa <mark.wieringa@csiro.au>
///

#include <askap/measurementequation/PercentileCalculator.h>
#include <askap/askap/AskapError.h>

namespace askap {

namespace synthesis {


    /// Good value for step is same scale as std dev of data
    /// Percentile value should be between 0 and 1
    PercentileCalculator::PercentileCalculator(double percentile, double step):
        itsPerc(percentile),itsStep(step),itsCount(0),itsValue(0)
    {}

    /// @brief (re)initialise
    void PercentileCalculator::init()
    {
       itsCount = 0;
       itsValue = 0;
    }

    /// @brief add a datapoint to the stats
    void PercentileCalculator::add(double data)
    {
        if (itsCount++ == 0) {
            itsValue = data;
            return;
        }
        if (itsValue > data) {
            itsValue -= itsStep * (1.0 - itsPerc);
        } else if (itsValue < data) {
            itsValue += itsStep * itsPerc;
        }
        if (abs(data-itsValue) < itsStep) {
            itsStep /= 2;
        }
    }
    /// @brief merge in another dataset using weighted mean
    void PercentileCalculator::merge(const PercentileCalculator & other)
    {
        ASKAPCHECK(itsPerc == other.itsPerc,"Inconsistent percentile values in merge");
        itsValue = itsValue * itsCount + other.value() * other.count();
        itsCount += other.count();
        if (itsCount>0) {
            itsValue /= itsCount;
        }
    }
    /// @brief serialise
    void PercentileCalculator::writeToBlob(LOFAR::BlobOStream& os) const
    {
        os << itsPerc << itsStep << (LOFAR::TYPES::uint64)itsCount << itsValue;
    }
    /// @brief deserialise
    void PercentileCalculator::readFromBlob(LOFAR::BlobIStream& is)
    {
        LOFAR::TYPES::uint64 count;
        is >> itsPerc >> itsStep >> count >> itsValue;
        itsCount = count;
    }

} // namespace synthesis

} // namespace askap
