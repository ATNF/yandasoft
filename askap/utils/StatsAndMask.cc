#include <askap/utils/StatsAndMask.h>

#include <askap_accessors.h>

#include <askap/askap/AskapLogging.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Containers/Record.h>
#include <boost/shared_array.hpp>

ASKAP_LOGGER(logger, ".statsAndMask");

using namespace askap;
using namespace askap::utils;
using namespace askap::accessors;

StatsAndMask::StatsAndMask(askapparallel::AskapParallel &comms, const std::string& cubeName, 
                           boost::shared_ptr<askap::accessors::IImageAccess<>> imageCube)
    : itsComms(comms), itsImageName(cubeName), itsImageCube(imageCube),
    itsUnit(imageCube->getUnits(cubeName)), itsScaleFactor(1.0)
{
    if ( itsUnit == "Jy/beam" ) {
        itsScaleFactor = 1000.0;
    }
    ASKAPLOG_INFO_STR(logger,"unit: " << itsUnit << ", scale: " << itsScaleFactor);
}

void StatsAndMask::setScaleFactor(float scaleFactor)
{
    itsScaleFactor = scaleFactor;
}

void StatsAndMask::setUnits(const std::string& unit)
{
    itsUnit = unit;
}

void StatsAndMask::calculate(const std::string& name, Channel channel,const casacore::IPosition& blc, const casacore::IPosition& trc)
{
    if ( boost::shared_ptr<IImageAccess<>> imageCube = itsImageCube.lock() ) {
        casacore::Array<float> imgPerPlane = imageCube->read(name,blc,trc);
        
        if ( itsComms.rank() == 0 ) {
                
        }
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
        ASKAPLOG_INFO_STR(logger,"node: " << itsComms.nodeName() << ", rank: " << itsComms.rank() << " - number of stats objects received: " << msgSize/sizeof(Stats));
        // now read the stats
        boost::shared_array<unsigned char> buffer(new unsigned char[msgSize]);
        itsComms.receive(buffer.get(),msgSize,source);
        // add the stats to the map
        unsigned long numberOfStatsObjsReceived = (msgSize/sizeof(Stats));
        ASKAPLOG_INFO_STR(logger,"node: " << itsComms.nodeName() << ", rank: " << itsComms.rank() 
                                << " - msgSize: " << msgSize << "; sizeof(Stats): " << sizeof(Stats)
                          << "; numberOfStatsObjsReceived: " << numberOfStatsObjsReceived);
        unsigned char* ptr = buffer.get();
        for (unsigned long i = 0; i < numberOfStatsObjsReceived; i++) {
            Stats s;
            std::memcpy(reinterpret_cast<void *>(&s),buffer.get(),sizeof(Stats));
            ASKAPLOG_INFO_STR(logger,"node: " << itsComms.nodeName() << ", rank: " << itsComms.rank() 
                                << " - Received stats from channel: " << s.channel);
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

void StatsAndMask::writeStatsToImageTable(const std::string& name)
{
    ASKAPLOG_INFO_STR(logger,"writeStatsToImageTable - name: " << name);
    casacore::Record statsRecord;
    casacore::Record statsTable;
    statsRecord.define("IMGSTATS","Image Statistics");

    const auto rows = itsStatsPerChannelMap.size();
    casacore::IPosition shape(1);
    shape(0) = rows;
    casacore::Vector<casacore::uInt> chanVect(rows);
    casacore::Vector<casacore::Float> rmsVect(rows);
    casacore::Vector<casacore::Float> stdVect(rows);
    casacore::Vector<casacore::Float> meanVect(rows);
    casacore::Vector<casacore::Float> onepcVect(rows);
    casacore::Vector<casacore::Float> medianVect(rows);
    casacore::Vector<casacore::Float> madfmVect(rows);
    casacore::Vector<casacore::Float> maxvalVect(rows);
    casacore::Vector<casacore::Float> minvalVect(rows);
    auto index = static_cast<unsigned long> (0);
    ASKAPLOG_INFO_STR(logger,"writeStatsToImageTable - rows: " << rows);
    for ( const auto& kvp : itsStatsPerChannelMap ) {
        const auto& statsPerChan = kvp.second;
        chanVect[index] = statsPerChan.channel;
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
    statsTable.define("channel",chanVect);
    statsTable.define("rms",rmsVect);
    statsTable.define("std",stdVect);
    statsTable.define("mean",meanVect);
    statsTable.define("onepc",onepcVect);
    statsTable.define("median",medianVect);
    statsTable.define("madfm",madfmVect);
    statsTable.define("max",maxvalVect);
    statsTable.define("min",minvalVect);

    const auto columns = static_cast<unsigned int> (9); // number of columns
    casacore::Vector<casacore::String> units(columns);
    units[0] = "uint";
    units[1] = "float";
    units[2] = "float";
    units[3] = "float";
    units[4] = "float";
    units[5] = "float";
    units[6] = "float";
    units[7] = "float";
    units[8] = "float";
    statsTable.define("Units",units); 

    statsRecord.defineRecord("ImgStats",statsTable);

    if ( boost::shared_ptr<IImageAccess<>> imageCube = itsImageCube.lock() ) {
        imageCube->setInfo(name,statsRecord);
    }
    
}

