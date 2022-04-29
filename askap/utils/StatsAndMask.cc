#include <askap/utils/StatsAndMask.h>

#include <askap_accessors.h>

#include <askap/askap/AskapLogging.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <boost/shared_array.hpp>

ASKAP_LOGGER(logger, ".statsAndMask");

using namespace askap;
using namespace askap::utils;
using namespace askap::accessors;

void StatsAndMask::setScaleFactor(float scaleFactor)
{
}

void StatsAndMask::setUnits(const std::string& unit)
{
}

void StatsAndMask::calculate(const std::string& name, Channel channel,const casacore::IPosition& blc, const casacore::IPosition& trc)
{
    if ( boost::shared_ptr<IImageAccess<>> imageCube = itsImageCube.lock() ) {
        casacore::Array<float> imgPerPlane = imageCube->read(name,blc,trc);
        
        Stats stats;
        stats.rms = casacore::rms(imgPerPlane);
        stats.std = casacore::stddev(imgPerPlane);
        stats.mean = casacore::mean(imgPerPlane);
        stats.median = casacore::median(imgPerPlane);
        stats.madfm = casacore::madfm(imgPerPlane);
        stats.maxval = casacore::max(imgPerPlane);
        stats.minval = casacore::min(imgPerPlane);
        stats.onepc = casacore::fractile(imgPerPlane,0.01);
        stats.channel = channel;

        itsStatsPerChannelMap.insert(std::make_pair(channel,stats));
    }
}

void StatsAndMask::receiveStats()
{
    MPI_Status status;
    // this method is called by the master i.e rank = 0
    // number of workers not including the master
    int numWorkers = itsComms.nProcs() - 1;
    for ( int workerRank = 1; workerRank <= numWorkers; workerRank++ ) {
        int source = workerRank;
        // first read the message size    
        unsigned char* msgSizebuffer[sizeof(unsigned long)];
        itsComms.receive(msgSizebuffer,sizeof(unsigned long),source);
        unsigned long msgSize;
        std::memcpy(reinterpret_cast<void *>(&msgSize),msgSizebuffer,sizeof(unsigned long));
        // check that msgSize is a multiple of sizeof(Stats)
        ASKAPCHECK((msgSize % sizeof(Stats)) == 0, 
                    "StatsAndMask::receiveStats: msgSize is a multiple of sizeof(Stats)");
        // now read the stats
        boost::shared_array<unsigned char> buffer(new unsigned char[msgSize]);
        itsComms.receive(buffer.get(),msgSize,source);
        
        // add the stats to the map
        unsigned long numberOfStatsObjsReceived = (msgSize % sizeof(Stats));
        unsigned char* ptr = buffer.get();
        for (unsigned long i = 0; i < numberOfStatsObjsReceived; i++) {
            Stats s;
            std::memcpy(reinterpret_cast<void *>(&s),buffer.get(),sizeof(Stats));
            itsStatsPerChannelMap.insert(std::make_pair(s.channel,s));
        }
    }

}

void StatsAndMask::sendStats(int destRank)
{
    auto mapSize = itsStatsPerChannelMap.size();
    auto msgSize = static_cast<unsigned long> (mapSize * sizeof(Stats)); 
    
    // send the message size
    unsigned char* msgSizebuffer[sizeof(unsigned long)];
    std::memcpy(msgSizebuffer,reinterpret_cast<void *>(&msgSize),sizeof(unsigned long));
    itsComms.send(msgSizebuffer,sizeof(unsigned long),destRank);

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
