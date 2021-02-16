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

#ifndef SYNTHESIS_PERCENTILE_CALCULATOR_H
#define SYNTHESIS_PERCENTILE_CALCULATOR_H

#include <askap/scimath/fitting/ISerializable.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>

#include <utility>

namespace askap {

namespace synthesis {
/// @brief W percentile calculator
/// @brief helper class to calculate running percentile (approximate, no sort or storage)
class PercentileCalculator
{
public:
    /// Good value for step is same scale as std dev of data
    /// Percentile value should be between 0 and 1
    explicit PercentileCalculator(double percentile=0.95, double step=1000.0);

    /// @brief (re)initialise
    void init();

    /// @brief return percentile value of the data
    double value() const {
        return itsValue;
    }

    /// @brief return count of values
    size_t count() const {
        return itsCount;
    }

    /// @brief add a datapoint to the stats
    void add(double data);

    /// @brief merge in another dataset using weighted mean
    void merge(const PercentileCalculator & other);

    /// @brief serialise
    void writeToBlob(LOFAR::BlobOStream& os) const;

    /// @brief deserialise
    void readFromBlob(LOFAR::BlobIStream& is);

private:
    double itsPerc;
    double itsStep;
    size_t itsCount;
    double itsValue;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_PERCENTILE_CALCULATOR_H
