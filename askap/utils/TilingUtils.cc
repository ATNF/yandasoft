/// @file TilingUtils.cc
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
/// @author Mark Wieringa <mark.wieringa@csiro.au>
///

#include <askap/utils/TilingUtils.h>

#include <askap/askap/AskapLogging.h>
#include <askap/dataaccess/TableDataSelector.h>
#include <askap/dataaccess/DataAccessError.h>
#include <string>
ASKAP_LOGGER(logger, ".utils");


namespace askap {
namespace utils {


/// @brief set selection to distribute work over ranks by tiles
/// @param[in] IDataSelectorPtr sel : shared pointer to data selector
/// @param[in] string column : name of column to inspect (usually DATA or FLAG)
/// @param[in] uint nRanks: number of MPI ranks
/// @param[in] uint rank : MPI rank of current process
void distributeByTile(accessors::IDataSelectorPtr& sel, const std::string& column, uint nRanks, uint rank)
{
    ASKAPCHECK(rank < nRanks, "rank value has to be less than number of ranks");
    // Tiles are Table specific, so cast to TableDataSelector
    boost::shared_ptr<accessors::TableDataSelector> tsel = boost::dynamic_pointer_cast<accessors::TableDataSelector>(sel);
    if (!tsel) {
        ASKAPTHROW(accessors::DataAccessError,"Tile selection can only be used for data stored in Tables");
    }
    casacore::IPosition tileShape, hyperCubeShape;
    const uint nTiles = tsel->getTiling(column,tileShape,hyperCubeShape);
    if (nTiles == 0) {
        ASKAPLOG_WARN_STR(logger,"Column is not tiled or has multiple tiling schemes, cannot distribute by tile");
        return;
    }
    const uint div = nTiles / nRanks;
    const uint rem = nTiles % nRanks;
    // Simple round-robin: the first `rem` ranks receive an extra item
    const uint firstTile = rank * div + (rank < rem ? rank : rem);
    const uint numTiles = div + (rank < rem);

    if (nRanks > nTiles) {
        ASKAPLOG_WARN_STR(logger, "More ranks than tiles - some workers will be idle");
    }
    tsel->chooseDataTiles(numTiles, firstTile);
}

}
}
