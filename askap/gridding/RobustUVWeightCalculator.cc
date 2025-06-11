/// @file
/// @brief Class calculating robust weights from the grid of weights
/// @details This is an implementation of UVWeight calculator interface making robust weights.
/// The math follows the original code by Daniel Mitchell, who presumably checked it against casa.
///
/// @copyright (c) 2023 CSIRO
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

// own includes
#include <askap/gridding/RobustUVWeightCalculator.h>

#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.robustuvweightcalculator");
#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>


namespace askap {

namespace synthesis {

/// @brief perform processing for the given weight (single grid slice along the 3rd axis)
/// @details For performance reasons, slices along the 3rd axis are taken inside finalise method
/// of the builder (this can be changed if we ever had any effect where frequency dependence matter).
/// At this stage, we can guarantee that supplied matrix has contiguous storage.
/// @param[in] wt weight to work with (it is modified in situ).
/// @note The shape is supposed to stay intact.
void RobustUVWeightCalculator::process(casacore::Matrix<float> &wt) const
{
    ASKAPTRACE("RobustUVWeightCalculator::process");
    ASKAPLOG_DEBUG_STR(logger, "Before robustness sum of grid = " << sum(wt) << " robustness = "<<itsRobustness );

    // MV: do processing in double precision as in the original Mitch's hack. It may be worth checking whether we can
    // get away with single precision float here (it seems possible, perhaps with specialised sum and sum of squares methods)
    casa::Matrix<double> dblWt(wt.shape());
    casa::convertArray<double,float>(dblWt, wt);
    const double aveWgtSum = sum(dblWt*dblWt) / sum(dblWt);
    ASKAPLOG_DEBUG_STR(logger, "Average weight sum estimate: " << aveWgtSum);

    const float noisePower = (1.0/aveWgtSum)*25.0*std::pow(10., -2.0*itsRobustness);
    ASKAPLOG_DEBUG_STR(logger, "noisePower estimate: " << noisePower);

    for (casacore::uInt iv=0; iv < wt.ncolumn(); ++iv) {
        for (casacore::uInt iu=0; iu < wt.nrow(); ++iu) {
            // note, the weight is reciprocal to that in the original hack - this is because we apply the weight directly to visibilities
            // by multiplication, not via noise method (sigma estimate) which enters as inverse variance weight resulting in a division.
            // The user can use "Reciprocal" weight calculator in the chain if for some reason 1 / wt is desired instead.
            wt(iu, iv) = 1.f / (noisePower * wt(iu, iv) + 1.f);
        }
    }

}

} // namespace synthesis

} // namespace askap

