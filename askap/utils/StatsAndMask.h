/// @file StatsAndMask.h
/// @brief 
/// @details 
///
///
/// @copyright (c) 2016 CSIRO
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
/// @author Minh Vuong <minh.vvuong@csiro.au>
///

#ifndef ASKAP_ACCESSORS_STATS_AND_MASK_H
#define ASKAP_ACCESSORS_STATS_AND_MASK_H

#include <askap/imageaccess/IImageAccess.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/images/Regions/ImageRegion.h>
#include <casacore/images/Regions/RegionHandler.h>
#include <casacore/casa/Arrays/Array.h>

#include <string>
#include <vector>
#include <map>
#include <string>

#include <boost//weak_ptr.hpp>
#include <boost//shared_ptr.hpp>

#include <askap/askapparallel/AskapParallel.h>

namespace askap {
namespace utils {
using Channel = uint;

struct Stats {
    Channel channel;    
    float rms;
    float std;
    float mean;
    float onepc;
    float median;
    float madfm;
    float maxval;
    float minval;
};
class StatsAndMask {
    public:
        using Channel = uint;

        StatsAndMask(askapparallel::AskapParallel &comms, const std::string& cubeName, boost::shared_ptr<askap::accessors::IImageAccess<>> imageCube);
        void setScaleFactor(float scaleFactor);
        void setUnits(const std::string& unit); 
        void calculate(const std::string& name, Channel channel,const casacore::IPosition& blc, const casacore::IPosition& trc);
        const std::map<Channel,Stats>& statsPerChannelMap() const
        {
            return itsStatsPerChannelMap;
        }
        const std::map<Channel,bool>& masksPerChannelMap() const
        {
            return itsMasksPerChannelMap;
        }
        void writeStatsToImageTable(const std::string& name);
        // @TODO. This method is just a placeholder for now.
        void receiveStats();
        void sendStats(int destRank = 0);
    private:
        boost::weak_ptr<askap::accessors::IImageAccess<>> itsImageCube;
        std::string itsImageName;
        std::string itsUnit;
        float itsScaleFactor;
        askapparallel::AskapParallel& itsComms;
        std::map<Channel,Stats> itsStatsPerChannelMap;
        std::map<Channel,bool>  itsMasksPerChannelMap;
};
}
}
#endif
