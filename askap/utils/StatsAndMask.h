/// @file StatsAndMask.h
///
/// @copyright (c) 2022 CSIRO
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

#ifndef ASKAP_YANDASOFT_STATS_AND_MASK_H
#define ASKAP_YANDASOFT_STATS_AND_MASK_H

#include <askap/imageaccess/IImageAccess.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/images/Regions/ImageRegion.h>
#include <casacore/images/Regions/RegionHandler.h>
#include <casacore/casa/Arrays/Array.h>

#include <string>
#include <set>
#include <vector>
#include <map>
#include <string>

#include <boost//weak_ptr.hpp>
#include <boost//shared_ptr.hpp>

#include <askap/askapparallel/AskapParallel.h>

namespace askap {
namespace utils {
using Channel = unsigned int;

/// @brief a structure to store the statistics
struct Stats {
    Channel channel;
    double freq;
    float rms;
    float std;
    float mean;
    float onepc;
    float median;
    float madfm;
    float maxval;
    float minval;
};

/// @details A class that calculates the image cube's statistics and performs the masking of the
///          image's channels based on some user defined thresholds. It does similar tasks as the
///          the current cubestatsHelpers.py in the askap-pipeline package and the maskBadChannels
///          in this package.
class StatsAndMask {
    public:
        using Channel = uint;

        /// @brief - constructor
        /// @param[in] comms - MPI comms
        /// @param[in] cubeName - name of image cube
        /// @param[in] imageCube - a boost shared pointer to the image access instance
        StatsAndMask(askapparallel::AskapParallel &comms, const std::string& cubeName, boost::shared_ptr<askap::accessors::IImageAccess<>> imageCube);
        StatsAndMask(askapparallel::AskapParallel &comms)
        : itsComms(comms)
        {}

        /// @brief - disallowed default, copy and move constructors
        StatsAndMask() = delete;
        StatsAndMask(const StatsAndMask&) = delete;
        StatsAndMask(StatsAndMask&&) = delete;
        StatsAndMask& operator=(const StatsAndMask&) = delete;
        StatsAndMask& operator=(StatsAndMask&&) = delete;
        ~StatsAndMask() {}

        /// @brief - set the scaling factor
        /// @param[in] scaleFactor - factor to be set
        void setScaleFactor(float scaleFactor);
        /// @brief - set the unit
        /// @param[in] unit - unit to be set
        void setUnits(const std::string& unit);

        /// @brief calculates the per plane statistics of the image cube
        /// @param[in] blc - bottom left corner of the image plane
        /// @param[in] trc - top right corner of the image plane
        void calculate(const casacore::IPosition& blc, const casacore::IPosition& trc);
        /// @brief calculates the per plane statistics
        /// @param[in] channel - channel of the image where the statistics are to be calculated
        /// @param[in] arr - the channel image where the statistics are calculated
        void calculate(Channel channel, const casacore::Array<float>& arr);

        /// @brief returns the image cube's statistics
        /// @return A map of the per plane statistics of the image cube.
        ///         The key of the map is the image channel and the value
        ///         contains the iimage statistics of the channel.
        const std::map<Channel,Stats>& statsPerChannelMap() const
        {
            return itsStatsPerChannelMap;
        }

        /// @brief this method masks the channels in the image cube if they dont meet the user
        ///        defined thresholds. It is based very much on the code in the maskBadChannels.cc
        /// @param[in] image - name of the image
        /// @param[in] threshold - user defined threshold (significanceLevel) in the parset.
        /// @param[in] badMADFM - user defined  madfm threshold (madfmThreshold) in the parset.
        /// @param[in] maskBlank - a boolean defined in the parset to indicate (with other flags)
        ///                        if the channel is masked.
        /// @param[in] useSignificance - a boolean defined in the parset
        /// @param[in] useNoise - a boolean defined in the parset. If set and the channel has considered
        ///                       to have bad noise then the channel is masked
        /// @param[in] editStats - a boolean defined in the parset. If set, the statistics are written to
        ///                        to a file (see outputStats)
        /// @param[in] editImage - a boolean defined in the parset. If set and the channel is
        ///                        considered "bad" then the image plane for that channel is masked.
        /// @param[in] outputStats - the file name that contains the image cube's  statistics
        void maskBadChannels(const std::string& image, float threshold, float badMADFM, bool maskBlank,
                             bool useSignificance, bool useNoise, bool editStats, bool editImage,
                             const std::string& outputStats, int master = 0);

        /// @brief This method writes the statistics to the image cube
        /// @param[in] name - name of the image
        void writeStatsToImageTable(const std::string& name);

        /// @brief This method writes the statistics to a file
        /// @param[in] catalogue - name of the file to be written to
        void writeStatsToFile(const std::string& catalogue);

        /// @brief This method receives the statistics from the sendStats() method
        /// @details This method receives the statistics from the sendStats() method
        ///          and adds them to the itsStatsPerChannelMap. The idea is that the
        ///          master rank/process calls this method to collect all the statistics
        ///          from the workers.
        void receiveStats(const std::set<unsigned int>& excludedRanks = std::set<unsigned int> {0});

        /// @brief This method sends the statistics of the channels that it collects back
        ///        to the master (i.e destRank)
        /// @details This method sends the statistics of the channels that it collects back
        ///        to the master (i.e destRank). The idea is that the worker rank(s) calls
        ///        this method to send all the statistics it collects back to the master rank.
        /// @param[in] destRank - process (master rank) that receives the statistics.
        void sendStats(int destRank = 0);
        void print() const;
        void set(const std::string& cubeName, boost::shared_ptr<askap::accessors::IImageAccess<>> imageCube)
        {
            itsImageCube = imageCube;
            itsImageName = cubeName;
        }
        /// @brief This static function  writes the statistics to the image cube
        /// @detail This static function writes the statistics to the image cube.
        ///         In the imager.cc code, it is called by the ImagerParallel object
        ///         after the image cube is created/saved i.e ( after SynthesisParamsHelper::saveImageParameter)
        ///         is invoked. The caller can stop the statistics collection/calculation by
        ///         setting the Cimager.calcstats = false
        /// @param[in] name - name of the image
        /// @param[in] comms - MPI comms
        /// @param[in] iacc - image access object
        /// @param[in] imgName - name of image cube
        /// @param[in] parset - configuration parameters object
        static void writeStatsToImageTable(askapparallel::AskapParallel &comms,
                                           accessors::IImageAccess<float>& iacc,
                                           const std::string& imgName,
                                           const LOFAR::ParameterSet& parset);
    private:
        void updateParams();

        /// @brief calculates the per plane statistics
        /// @param[in] channel - chanel of the image where the statistics are to be calculated
        /// @param[in] arr - the channel image where the statistics are calculated
        /// @return Stats - the statistics of the channel image
        Stats calculateImpl(Channel channel, const casacore::Array<float>& arr);

        /// A weak pointer to the image access object
        boost::weak_ptr<askap::accessors::IImageAccess<>> itsImageCube;
        /// Name of the image cube
        std::string itsImageName;
        std::string itsUnit;
        /// The reference frequency of the image cube
        double itsRefFrequency;
        /// The increment frequency of the image cube
        double itsFreqIncrement;
        float itsScaleFactor;
        /// MPI comms
        askapparallel::AskapParallel& itsComms;
        /// A map contains the per plane/chanel image statistics
        std::map<Channel,Stats> itsStatsPerChannelMap;
        //std::map<Channel,bool>  itsMasksPerChannelMap;
};
}
}
#endif
