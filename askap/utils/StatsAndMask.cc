/// @file StatsAndMask.cc
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

#include <askap/utils/StatsAndMask.h>

#include <askap_accessors.h>

#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapUtil.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/coordinates/Coordinates.h>
#include <boost/shared_array.hpp>

#include <fstream>
#include <iomanip>

ASKAP_LOGGER(logger, ".statsAndMask");

using namespace askap;
using namespace askap::utils;
using namespace askap::accessors;

/// @brief - constructor
/// @param[in] comms - MPI comms
/// @param[in] cubeName - name of image cube
/// @param[in] imageCube - a boost shared pointer to the image access instance
StatsAndMask::StatsAndMask(askapparallel::AskapParallel &comms, const std::string& cubeName,
                           boost::shared_ptr<askap::accessors::IImageAccess<>> imageCube)
    : itsComms(comms), itsImageName(cubeName), itsImageCube(imageCube)
{
#ifndef HAVE_MPI
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
#endif
}

/// @brief - set the scaling factor
/// @param[in] scaleFactor - factor to be set
void StatsAndMask::setScaleFactor(float scaleFactor)
{
    itsScaleFactor = scaleFactor;
}

/// @brief - set the unit
/// @param[in] unit - unit to be set
void StatsAndMask::setUnits(const std::string& unit)
{
    itsUnit = unit;
}

/// @brief calculates the per plane statistics of the image cube
/// @param[in] blc - bottom left corner of the image plane
/// @param[in] trc - top right corner of the image plane
void StatsAndMask::calculate(const casacore::IPosition& blc, const casacore::IPosition& trc)
{
#ifdef HAVE_MPI
    if ( boost::shared_ptr<IImageAccess<>> imageCube = itsImageCube.lock() ) {
        casacore::Array<float> imgPerPlane = imageCube->read(itsImageName,blc,trc);
        // remove all the NaN from the imgPerPlane because it messes it the calculation
        // of the statistics
        casacore::Vector<float> nonMaskArray(imgPerPlane.size());
        unsigned long index = 0;
        for(size_t i=0; i < imgPerPlane.size(); i++) {
            if ( ! casacore::isNaN(imgPerPlane.data()[i]) ) {
                nonMaskArray[index] = imgPerPlane.data()[i];
                index += 1;
            }
        }
        nonMaskArray.resize(index,true);
        // The client code that uses this method normally wants to collect the 
        // statistics of one particular channel. However, the code here is generic
        // enough that it is able to collect statistics for a range of channels.
        auto startChan = trc(3);
        auto endChan = blc(3);
        for(auto channel = startChan; channel <= endChan; channel++) {
            Stats stats;
            if ( nonMaskArray.size() > 0 ) {
                stats = calculateImpl(channel,nonMaskArray);
            } else {
                // nonMaskArray.size() == 0 if the image plane is masked or contains NaN pixels
                // in this case, dont use nonMaskArray
                stats = calculateImpl(channel,imgPerPlane);
            }
            itsStatsPerChannelMap.insert(std::make_pair(channel,stats));
        }
    }
#else
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
#endif
}

/// @brief calculates the per plane statistics
/// @param[in] channel - chanel of the image where the statistics are to be calculated
/// @param[in] arr - the channel image where the statistics are calculated
void StatsAndMask::calculate(Channel channel, const casacore::Array<float>& arr)
{
#ifdef HAVE_MPI
    // remove all the NaN from the input arr because it messes it the calculation
    // of the statistics
    casacore::Vector<float> nonMaskArray(arr.size());
    unsigned long index = 0;
    for(size_t i=0; i < arr.size(); i++) {
        if ( ! casacore::isNaN(arr.data()[i]) ) {
            nonMaskArray[index] = arr.data()[i];
            index += 1;
        }
    }
    nonMaskArray.resize(index,true);

    Stats stats;
    if ( nonMaskArray.size() > 0 ) {
      stats = calculateImpl(channel,nonMaskArray);
    } else {
        // nonMaskArray.size() == 0 if the image plane is masked or contains NaN pixels
        // in this case, dont use nonMaskArray
        stats = calculateImpl(channel,arr);
    }

    // check if the channel is in the map
    auto search = itsStatsPerChannelMap.find(channel);
    if ( search != itsStatsPerChannelMap.end() ) {
        // erase it
        itsStatsPerChannelMap.erase(search);
    }
    itsStatsPerChannelMap.insert(std::make_pair(channel,stats));
#else
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
#endif
}

Stats StatsAndMask::calculateImpl(Channel channel, const casacore::Array<float>& imgPerPlane)
{
#ifdef HAVE_MPI
    updateParams();

    Stats stats;

    if ( boost::shared_ptr<IImageAccess<>> imageCube = itsImageCube.lock() ) {
        stats.rms = itsScaleFactor * casacore::rms(imgPerPlane);
        stats.std = itsScaleFactor * casacore::stddev(imgPerPlane);
        stats.mean = itsScaleFactor * casacore::mean(imgPerPlane);
        stats.median = itsScaleFactor * casacore::median(imgPerPlane);
        stats.madfm = itsScaleFactor * casacore::madfm(imgPerPlane);
        stats.maxval = itsScaleFactor * casacore::max(imgPerPlane);
        stats.minval = itsScaleFactor * casacore::min(imgPerPlane);
        stats.onepc =  itsScaleFactor * casacore::fractile(imgPerPlane,0.01);
        stats.channel = channel;
        stats.freq = (itsRefFrequency + (channel*itsFreqIncrement))/1.0e6;


        if ( std::isnan(stats.rms) )
            stats.rms = 0.0;

        if (  std::isnan(stats.std) )
            stats.std = 0.0;

        if ( std::isnan(stats.mean) )
            stats.mean = 0.0;

        if ( std::isnan(stats.onepc) )
            stats.onepc = 0.0;

        if ( std::isnan(stats.median) )
            stats.median = 0.0;

        if ( std::isnan(stats.madfm) )
            stats.madfm = 0.0;

        if ( std::isnan(stats.maxval) )
            stats.maxval = 0.0;

        if ( std::isnan(stats.minval) )
            stats.minval = 0.0;

        ASKAPLOG_INFO_STR(logger,"channel: " << stats.channel << ", rms: " << stats.rms << ", std: " << stats.std
                            << ", mean: " << stats.mean << ", median: " << stats.median << ", madfm: " << stats.madfm
                            << ", maxval: " << stats.maxval << ", stats.minval: " << stats.minval
                            << ", onepc: " << stats.onepc << ", #points: "<<imgPerPlane.size());
    }
    return stats;
#else
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
    Stats stats;
    return stats;
#endif
}

/// @brief returns the image cube's statistics
/// @return A map of the per plane statistics of the image cube.
///         The key of the map is the image channel and the value
///         contains the iimage statistics of the channel.
void StatsAndMask::receiveStats(const std::set<unsigned int>& excludedRanks)
{
#ifdef HAVE_MPI
    MPI_Status status;
    // this method is called by the master i.e rank = 0
    // number of workers not including the master
    int numWorkers = itsComms.nProcs() - 1;
    int ranks = itsComms.nProcs();
    for ( int rank = 0; rank < ranks; rank++ ) {
        if ( excludedRanks.find(rank) != excludedRanks.end() ) {
            ASKAPLOG_INFO_STR(logger,"Skipping rank: " << rank);
            continue;
        }
        int source = rank;
        // first read the message size
        unsigned char* msgSizebuffer[sizeof(unsigned long)];
        itsComms.receive(msgSizebuffer,sizeof(unsigned long),source);
        unsigned long msgSize;
        std::memcpy(reinterpret_cast<void *>(&msgSize),msgSizebuffer,sizeof(unsigned long));
        // check that msgSize is a multiple of sizeof(Stats)
        ASKAPCHECK((msgSize % sizeof(Stats)) == 0,
                    "StatsAndMask::receiveStats: msgSize is a multiple of sizeof(Stats)");
        // ASKAPLOG_INFO_STR(logger,"node: " << itsComms.nodeName() << ", rank: " << itsComms.rank() << " - number of stats objects received: " << msgSize/sizeof(Stats));
        // now read the stats but nothing to read if msgSize == 0
        if ( msgSize > 0 ) {
            boost::shared_array<unsigned char> buffer(new unsigned char[msgSize]);
            itsComms.receive(buffer.get(),msgSize,source);
            // add the stats to the map
            unsigned long numberOfStatsObjsReceived = (msgSize/sizeof(Stats));
            unsigned char* ptr = buffer.get();
            for (unsigned long i = 0; i < numberOfStatsObjsReceived; i++) {
                Stats s;
                //std::memcpy(reinterpret_cast<void *>(&s),buffer.get(),sizeof(Stats));
                std::memcpy(reinterpret_cast<void *>(&s),ptr,sizeof(Stats));
                ASKAPLOG_INFO_STR(logger,"node: " << itsComms.nodeName() << ", rank: " << itsComms.rank()
                                << " - Received stats from channel: " << s.channel);

                itsStatsPerChannelMap.insert(std::make_pair(s.channel,s));
                ptr = ptr + sizeof(Stats);
            }
        }
    }
#else
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
#endif
}

/// @brief This method sends the statistics of the channels that it collects back
///        to the master (i.e destRank)
/// @details This method sends the statistics of the channels that it collects back
///        to the master (i.e destRank). The idea is that the worker rank(s) calls
///        this method to send all the statistics it collects back to the master rank.
/// @param[in] destRank - process (master rank) that receives the statistics.
void StatsAndMask::sendStats(int destRank)
{
#ifdef HAVE_MPI
    auto mapSize = itsStatsPerChannelMap.size();
    auto msgSize = static_cast<unsigned long> (mapSize * sizeof(Stats));

    // send the message size
    unsigned char* msgSizebuffer[sizeof(unsigned long)];
    std::memcpy(msgSizebuffer,reinterpret_cast<void *>(&msgSize),sizeof(unsigned long));
    itsComms.send(msgSizebuffer,sizeof(unsigned long),destRank);

    ASKAPLOG_INFO_STR(logger,"rank: " << itsComms.rank() << " sends  "
                        << msgSize << " bytes to rank " << destRank);
    // if the map is empty then got nothing to send
    if ( msgSize > 0 ) {
        // put the stats in the buffer
        boost::shared_array<unsigned char> buffer(new unsigned char[msgSize]);
        unsigned char* ptr = buffer.get();
        for (auto& kvp : itsStatsPerChannelMap) {
            auto& stats = kvp.second; // stat's type is Stats
            std::memcpy(ptr,reinterpret_cast<void *> (&stats),sizeof(Stats));
            ptr = ptr + sizeof(Stats);
        }
        // now send the stats back to the destRank
        itsComms.send(buffer.get(),msgSize,destRank);
    }
#else
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
#endif
}

void StatsAndMask::print() const
{
    for (auto& kvp : itsStatsPerChannelMap) {
        auto& stats = kvp.second; // stat's type is Stats
        ASKAPLOG_INFO_STR(logger,"=====================================");
        ASKAPLOG_INFO_STR(logger, "Rank: " << itsComms.rank()
                                << "chan: " << stats.channel
                                << ", rms: " << stats.rms
                                << ", std: " << stats.std
                                << ", mean: " << stats.mean
                                << ", onepc: " << stats.onepc
                                << ", median: " << stats.median
                                << ", madfm: " << stats.madfm
                                << ", maxval: " << stats.maxval
                                << ", minval: " << stats.minval);
    }
}

/// @brief This method writes the statistics to the image cube
/// @param[in] catalogue - name of the file to be written to
void StatsAndMask::writeStatsToImageTable(const std::string& name)
{
#ifdef HAVE_MPI
    ASKAPLOG_INFO_STR(logger,"writeStatsToImageTable - name: " << name);
    casacore::Record statsRecord;
    casacore::Record statsTable;
    statsRecord.define("IMGSTATS","Image Statistics");

    const auto rows = itsStatsPerChannelMap.size();
    casacore::IPosition shape(1);
    shape(0) = rows;
    casacore::Vector<casacore::uInt> chanVect(rows);
    casacore::Vector<casacore::Double> freqVect(rows);
    casacore::Vector<casacore::Float> rmsVect(rows);
    casacore::Vector<casacore::Float> stdVect(rows);
    casacore::Vector<casacore::Float> meanVect(rows);
    casacore::Vector<casacore::Float> onepcVect(rows);
    casacore::Vector<casacore::Float> medianVect(rows);
    casacore::Vector<casacore::Float> madfmVect(rows);
    casacore::Vector<casacore::Float> maxvalVect(rows);
    casacore::Vector<casacore::Float> minvalVect(rows);
    auto index = static_cast<unsigned long> (0);
    ASKAPLOG_DEBUG_STR(logger,"writeStatsToImageTable - rows: " << rows);
    for ( const auto& kvp : itsStatsPerChannelMap ) {
        const auto& statsPerChan = kvp.second;
        chanVect[index] = statsPerChan.channel;
        freqVect[index] = statsPerChan.freq;

        rmsVect[index] = statsPerChan.rms;
        stdVect[index] = statsPerChan.std;
        meanVect[index] = statsPerChan.mean;
        onepcVect[index] = statsPerChan.onepc;
        medianVect[index] = statsPerChan.median;
        madfmVect[index] = statsPerChan.madfm;
        maxvalVect[index] = statsPerChan.maxval;
        minvalVect[index] = statsPerChan.minval;

        index += 1;
    }
    statsTable.define("#Channel",chanVect);
    statsTable.define("Frequency",freqVect);
    statsTable.define("Mean",meanVect);
    statsTable.define("Std",stdVect);
    statsTable.define("Median",medianVect);
    statsTable.define("MADFM",madfmVect);
    statsTable.define("1%ile",onepcVect);
    statsTable.define("Min",minvalVect);
    statsTable.define("Max",maxvalVect);
    statsTable.define("rms",rmsVect);

    const auto columns = static_cast<unsigned int> (10); // number of columns
    casacore::Vector<casacore::String> units(columns);
    units[0] = " ";
    units[1] = "MHz";
    units[2] = "mJy/beam";
    units[3] = "mJy/beam";
    units[4] = "mJy/beam";
    units[5] = "mJy/beam";
    units[6] = "mJy/beam";
    units[7] = "mJy/beam";
    units[8] = "mJy/beam";
    units[9] = "mJy/beam";
    statsTable.define("Units",units);

    statsRecord.defineRecord("ImgStats",statsTable);

    if ( boost::shared_ptr<IImageAccess<>> imageCube = itsImageCube.lock() ) {
        imageCube->setInfo(name,statsRecord);
    }
    ASKAPLOG_DEBUG_STR(logger,"writeStatsToImageTable - exit: " << rows);
#else
    ASKAPLOG_WARN_STR(logger,"No image stats written to table as StatsAndMask class requires MPI library which is not installed");
#endif
}

/// @brief This method writes the statistics to the image cube
/// @param[in] catalogue - name of the file to be written to
void StatsAndMask::writeStatsToFile(const std::string& catalogue)
{
#ifdef HAVE_MPI
    std::ofstream ofile;

    ofile.setf(std::ios::left);
    ofile.open(catalogue,std::ios::out|std::ios::trunc);
    if ( ofile ) {
        std::ios_base::fmtflags oflags = ofile.flags();
        ofile.setf(std::ios::right);
        ofile << std::setw(10) << "#channel";
        ofile << std::setw(15) << "frequency";
        ofile << std::setw(10) << "Mean";
        ofile << std::setw(10) << "Std";
        ofile << std::setw(10) << "Median";
        ofile << std::setw(10) << "MADFM";
        ofile << std::setw(10) << "1%ile";
        ofile << std::setw(10) << "Min";
        ofile << std::setw(10)  << "Max";
        ofile << std::endl;
        ofile << std::setw(10) << " ";
        ofile << std::setw(15) << "MHz";
        ofile << std::setw(10) << "mJy/beam";
        ofile << std::setw(10) << "mJy/beam";
        ofile << std::setw(10) << "mJy/beam";
        ofile << std::setw(10) << "mJy/beam";
        ofile << std::setw(10) << "mJy/beam";
        ofile << std::setw(10) << "mJy/beam";
        ofile << std::setw(10) << "mJy/beam";
        ofile << std::endl;

        ofile.setf(std::ios::right);
        for ( const auto& kvp : itsStatsPerChannelMap ) {
            const auto& statsPerChan = kvp.second;
            ofile << std::setw(10) << statsPerChan.channel;
            ofile << std::setw(15) << std::fixed <<  std::setprecision(6) << statsPerChan.freq;
            ofile << std::setw(10) << std::fixed << std::setprecision(3) << statsPerChan.mean;
            ofile << std::setw(10) << std::fixed << std::setprecision(3) << statsPerChan.std;
            ofile << std::setw(10) << std::fixed << std::setprecision(3) << statsPerChan.median;
            ofile << std::setw(10) << std::fixed << std::setprecision(3) << statsPerChan.madfm;
            ofile << std::setw(10) << std::fixed << std::setprecision(3) << statsPerChan.onepc;
            ofile << std::setw(10) << std::fixed << std::setprecision(3) << statsPerChan.minval;
            ofile << std::setw(10) << std::fixed << std::setprecision(3) << statsPerChan.maxval;
            ofile << std::endl;
        }

        ofile.close();
    }
#else
    ASKAPLOG_WARN_STR(logger,"No image stats written to file as StatsAndMask class requires MPI library which is not installed");
#endif
}

/// @brief this method masks the channels in the image cube if they dont meet the user
///        defined thresholds. It is based very much on the code in the maskBadChannels.cc
void StatsAndMask::maskBadChannels(const std::string& image, float threshold, float badMADFM,
                                   bool maskBlank, bool useSignificance, bool useNoise,
                                   bool editStats, bool editImage,
                                   const std::string& outputStats, int master)
{
#ifdef HAVE_MPI
    // only do the masking on the master rank
    // this code is more or less copied from the maskBadChannels.cc in yandasoft package
    if ( itsComms.rank() == master ) {
        casacore::Vector<double> ratio(itsStatsPerChannelMap.size());
        unsigned long index = 0;
        for ( const auto& kvp : itsStatsPerChannelMap ) {
            const auto& statsPerChan = kvp.second;
            if ( statsPerChan.madfm > 0.0 ) {
                ratio[index] = statsPerChan.onepc/statsPerChan.madfm;
            } else {
                ratio[index] = 0.0;
            }
            index += 1;
        }
        double med = casa::median(ratio);
        double mad = casa::median(abs(ratio-med));
        casa::Vector<double> significance=(ratio-med)/mad;
        ASKAPLOG_INFO_STR(logger, "Ratio median = " << med << " and madfm = " << mad);
        ASKAPLOG_INFO_STR(logger, "Acceptable channels have ratio values between "
                                  << med - threshold * mad << " and "
                                  << med + threshold * mad );

        if ( boost::shared_ptr<IImageAccess<>> iacc = itsImageCube.lock() ) {
            casacore::IPosition shape = iacc->shape(image);
            ASKAPLOG_INFO_STR(logger, "Shape of input image = " << shape);
            casa::CoordinateSystem coo = iacc->coordSys(image);
            int specAxis = coo.spectralAxisNumber();
            ASKAPLOG_INFO_STR(logger, "Spectral axis = " << specAxis);
            casa::IPosition chanShape=shape;
            ASKAPLOG_INFO_STR(logger, "Shape of a single channel = " << chanShape);
            chanShape[specAxis]=1;
            casa::Vector<float> datavec(chanShape.product(),0.);
            for(size_t i=0;i<datavec.size();i++) {
                casa::setNaN(datavec[i]);
            }
            casa::Array<float> data(chanShape,datavec.data());
            unsigned int numBad = 0;
            unsigned long i = 0;
            for ( auto& kvp : itsStatsPerChannelMap ) {
                auto& statsPerChan = kvp.second;
                ASKAPLOG_INFO_STR(logger, "Channel " << i << " has values "
                                      << statsPerChan.onepc << " and "
                                      << statsPerChan.madfm << " for ratio of "
                                      << ratio[i] << " and significance of "
                                      << significance[i]);
                bool isBlank = !(statsPerChan.madfm > 0.);
                bool isHighSig = (abs(significance[i]) > threshold);
                bool isBadNoise = (statsPerChan.madfm > badMADFM);

                bool maskChannel = false;
                maskChannel = maskChannel || (maskBlank && isBlank);
                maskChannel = maskChannel || (useSignificance && isHighSig);
                maskChannel = maskChannel || (useNoise && isBadNoise);

                if ( maskChannel ) {
                    if (isBlank) {
                        ASKAPLOG_INFO_STR(logger, "Blank Channel #" << statsPerChan.channel << ": std = "<<statsPerChan.std);
                    } else {
                        ASKAPLOG_INFO_STR(logger, "Bad Channel #" << statsPerChan.channel
                                              << ": MADFM = "<< statsPerChan.madfm
                                              << " 1%ile = " << statsPerChan.onepc
                                              << " 1%ile/MADFM = " << ratio[i]);
                    }

                    numBad++;
                    if (editImage) {
                        casa::IPosition loc(shape.size(), 0);
                        loc[specAxis] = statsPerChan.channel;
                        ASKAPLOG_INFO_STR(logger, "Writing data of shape " << data.shape() << " to position " << loc);
                        iacc->write(image,data,loc);
                        std::stringstream ss;
                        if ( statsPerChan.madfm > 0.) {
                            ss << "Masked blank channel " << statsPerChan.channel;
                        } else {
                            ss << "Masked channel " << statsPerChan.channel << " : madfm=" << statsPerChan.madfm
                               <<", 1%ile/madfm=" << ratio[i];
                        }
                    }
                    if (editStats) {
                        // mean[i] = std[i] = median[i] = madfm[i] = onepc[i] = minval[i] = maxval[i] = 0.0;
                        statsPerChan.mean = 0.0;
                        statsPerChan.std = 0.0;
                        statsPerChan.median = 0.0;
                        statsPerChan.madfm = 0.0;
                        statsPerChan.onepc = 0.0;
                        statsPerChan.minval = 0.0;
                        statsPerChan.maxval = 0.0;
                    }
                }
                i += 1;
            }
            ASKAPLOG_INFO_STR(logger, "Detected " << numBad << " bad channels to be masked");
            if (editStats) {
                ASKAPLOG_INFO_STR(logger, "Writing new stats file to " << outputStats);
                writeStatsToFile(outputStats);
            }
        }
    }
#else
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
#endif
}

void StatsAndMask::updateParams()
{
#ifdef HAVE_MPI
    if ( boost::shared_ptr<IImageAccess<>> imageCube = itsImageCube.lock() ) {
        // NOTE: dont put these statements in the constructor because in some cases, the imager which
        // uses this class for statistics collection passes in the IImageAccess to the constructor without calling the create()
        // method. If this method is called then the create() has been invoked and hence the IImageAccess is valid.
        itsUnit = imageCube->getUnits(itsImageName);
        itsScaleFactor = 1.0;

        if ( itsUnit == "Jy/beam" ) {
            itsScaleFactor = 1000.0;
        }
        ASKAPLOG_INFO_STR(logger,"unit: " << itsUnit << ", scale: " << itsScaleFactor);


        const casacore::CoordinateSystem coord = imageCube->coordSys(itsImageName);
        const casacore::SpectralCoordinate& spectralCoord = coord.spectralCoordinate();
        const casacore::Vector<casacore::Double> refFreq = spectralCoord.referenceValue();
        const casacore::Vector<casacore::Double> increment = spectralCoord.increment();
        ASKAPCHECK(refFreq.capacity() == 1, "Reference freq vector is NOT 1");
        ASKAPCHECK(increment.capacity() == 1, "Freq increment vector is NOT 1");
        itsRefFrequency = refFreq[0];
        itsFreqIncrement = increment[0];
        if ( itsComms.rank() == 0 ) {
            ASKAPLOG_INFO_STR(logger,"refFreq size: " << refFreq.capacity() << ", value: " << refFreq[0]);
            ASKAPLOG_INFO_STR(logger,"increment size: " << increment.capacity() << ", value: " << increment[0]);
        }
        // END NOTE
    }
#else
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
#endif
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
void StatsAndMask::writeStatsToImageTable(askapparallel::AskapParallel &comms,
                                          accessors::IImageAccess<float>& iacc,
                                          const std::string& imgName,
                                          const LOFAR::ParameterSet& parset)
{
#ifdef HAVE_MPI
    ASKAPLOG_INFO_STR(logger,"writeStatsToImageTable using static method. imgName: " << imgName);
    //const bool singleoutputfile = parset.getBool("singleoutputfile", false);
    const std::string statsFile = parset.getString("outputStats","");
    const bool calcstats = parset.getBool("calcstats", false);
    if ( !calcstats ) {
        ASKAPLOG_INFO_STR(logger,"calcstats param is false - StatsAndMask is not used");
        return;
    }

    if ( (imgName.find("taylor.0") != std::string::npos) ||
         (imgName.find("taylor") == std::string::npos) ) {
        // Stats is only done for taylor.0 image if nterms > 1
        // the iacc is created on the stack so we dont want iaccPtr to delete it when
        // it goes out of scope and hence contruct it with NullDeleter
        boost::shared_ptr<accessors::IImageAccess<float>> iaccPtr(&iacc,askap::utility::NullDeleter{});
        casa::IPosition cubeShape = iaccPtr->shape(imgName);
        casa::CoordinateSystem coo = iaccPtr->coordSys(imgName);
        const int specAxis = coo.spectralAxisNumber();
        const auto numberOfChansInImage = cubeShape[specAxis];
        ASKAPLOG_INFO_STR(logger, "numberOfChansInImage: " << numberOfChansInImage);
        StatsAndMask stats(comms,imgName,iaccPtr);
        for (unsigned int chan = 0; chan < numberOfChansInImage; chan++) {
            casacore::IPosition blc(4,0,0,0,chan);
            casacore::IPosition trc = cubeShape - 1;
            trc(3) = chan;
            casacore::Array<float> imagePerPlane = iaccPtr->read(imgName,blc,trc);
            stats.calculate(chan,imagePerPlane);
        }
        stats.writeStatsToImageTable(imgName);
        if ( statsFile != "" ) {
           stats.writeStatsToFile(imgName+"_"+statsFile);
        }
    }
#else
    ASKAPLOG_WARN_STR(logger,"No image stats calculated as StatsAndMask class requires MPI library which is not installed");
#endif
}
