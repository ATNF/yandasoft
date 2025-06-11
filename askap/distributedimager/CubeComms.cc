/// @file CubeComms.cc
///
/// Class to provide extra MPI communicator functionality to manage the writing
/// of distributed spectral cubes.
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
/// @author Stephen Ord <Stephen.Ord@csiro.au>
///

#include <limits>
#include <map>
#include <list>
/// out own header first


#include "askap/distributedimager/CubeComms.h"
#include "askap/messages/IMessage.h"

///ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askap/AskapError.h>
#include "Blob/BlobIStream.h"
#include "Blob/BlobIBufVector.h"
#include "Blob/BlobOStream.h"
#include "Blob/BlobOBufVector.h"


#include "casacore/casa/OS/Timer.h"



ASKAP_LOGGER(logger, ".CubeComms");

namespace askap {
namespace cp {

CubeComms::CubeComms(int argc, const char** argv) : AskapParallel(argc, const_cast<const char **>(argv))
{

    ASKAPLOG_DEBUG_STR(logger, "Constructor");
    writerCount = 1; // always at least one
    singleSink = false;

}
CubeComms::~CubeComms()
{
    ASKAPLOG_DEBUG_STR(logger, "Destructor");
}

size_t CubeComms::buildCommIndex()
{

    const int nWorkers = nProcs() - 1;
    std::vector<int> ranks(nWorkers, -1);

    for (size_t wrk = 0; wrk < nWorkers; ++wrk) {
        ranks[wrk] = 1 + wrk;
    }
    itsComrades = createComm(ranks);
    //ASKAPLOG_DEBUG_STR(logger, "Interworker communicator index is " << itsComrades);
    return itsComrades;
}
size_t CubeComms::buildWriterIndex(size_t comm)
{

    std::map<int, int>::iterator it = writerMap.begin();
    std::vector<int> ranks(writerMap.size());
    std::vector<int> comm_ranks(writerMap.size());
    int wrt = 0;
    while (it != writerMap.end()) {
        ranks[wrt] = it->first;
        it++;
        wrt++;
    }
    if (comm != 0){
    comm_ranks = translateRanks(ranks, comm);
    } else {
        comm_ranks = ranks;
    }
    itsWriters = createComm(comm_ranks, comm);
    //ASKAPLOG_DEBUG_STR(logger, "Interwriter communicator index is " << itsComrades);
    return itsWriters;
}
#ifdef HAVE_MPI
std::vector<int> CubeComms::translateRanks(std::vector<int> ranks, size_t comm)
{
    int n = ranks.size();
    std::vector<int> ranksOut(n);

    MPI_Group world_group;
    MPI_Group comm_group;

    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Comm_group(this->getComm(comm), &comm_group);

    MPI_Group_translate_ranks(world_group, n, ranks.data(), comm_group, ranksOut.data());
    return ranksOut;
}
#else
std::vector<int> CubeComms::translateRanks(std::vector<int> ranks, size_t comm)
{
    ASKAPTHROW(AskapError, "CubeComms::translateRanks() cannot be used - configured without MPI");
}
#endif

void CubeComms::initWriters(int nwriters, int nchanpercore)
{

    unsigned int nWorkers = nProcs() - 1;
    unsigned int nWorkersPerGroup = nWorkers / nGroups();
    unsigned int nWorkersPerWriter = floor(nWorkers / nwriters);

    if (nWorkersPerWriter < 1) {
        nWorkersPerWriter = 1;
    }

    int mywriter = 1;
    for (int wrk = 0; wrk < nWorkersPerGroup; wrk = wrk + nWorkersPerWriter) {
        if (nwriters > 1) {
            mywriter = floor(wrk / nWorkersPerWriter) * nWorkersPerWriter + 1;
        }

        std::pair<std::map<int, int>::iterator, bool> ret;
        ret = writerMap.insert(std::pair<int, int> (mywriter, writerCount));
        if (ret.second == false) {
            //ASKAPLOG_DEBUG_STR(logger, "element " << mywriter << " already existed");
        } else {
            writerCount++;
        }

    }

}
int CubeComms::isWriter()
{
    ASKAPLOG_DEBUG_STR(logger, "Providing writer status");
    /// see if my rank is in the writers list
    std::map<int, int>::iterator it = writerMap.begin();
    for (it = writerMap.begin(); it != writerMap.end(); ++it) {
        if (itsRank == it->first) {
            return it->second;
        }
    }
    return 0;
}
void CubeComms::addWriter(unsigned int writerRank)
{

    std::pair<std::map<int, int>::iterator, bool> ret;
    ret = writerMap.insert(std::pair<int, int> (writerRank, 0));

    if (ret.second == false) {
        //ASKAPLOG_DEBUG_STR(logger, "element " << writerRank << " already existed");
    } else {
        writerCount++;
    }


    // ret2 = clientMap.insert(std::pair<int,std::map<int,int> > (writerRank,clientlist) );
    // if (ret2.second==false) {
    //     ASKAPLOG_WARN_STR(logger, "element " << writerRank << " already existed");
    // }

}

void CubeComms::showWorkerMap() {
  std::map<int, int>::iterator it;
  it = workerMap.begin();
  while (it != workerMap.end()) {
    ASKAPLOG_INFO_STR(logger,"workerMap " << it->first << ":" << it->second);
    it++;
  }
}

void CubeComms::addWorker(unsigned int workerRank)
{

    std::pair<std::map<int, int>::iterator, bool> ret;
    ret = workerMap.insert(std::pair<int, int> (workerRank, 0));

    if (ret.second == false) {
        // ASKAPLOG_WARN_STR(logger, "worker " << workerRank << " already existed");
        //ASKAPLOG_DEBUG_STR(logger, "worker " << workerRank << " already existed");
    } else {
        //ASKAPLOG_DEBUG_STR(logger, "added worker rank " << " to workerMap");
    }

}
bool CubeComms::isCubeCreator()
{
    std::list<int>::iterator it;
    it = cubeCreators.begin();
    while (it != cubeCreators.end()) {
        if (*it == rank())
            return true;
        it++;
    }
    return false;
}

void CubeComms::addCubeCreator(int creator_rank)
{
    cubeCreators.push_back(creator_rank);
}

std::list<int> CubeComms::getCubeCreators()
{
    return cubeCreators;
}

void CubeComms::setSingleSink()
{
    std::map<int, int>::iterator it;
    it = writerMap.begin();
    addCubeCreator(it->first);
    this->singleSink = true;
}
void CubeComms::setMultiSink()
{
    std::map<int, int>::iterator it;
    it = writerMap.begin();
    while (it != writerMap.end()) {
        addCubeCreator(it->first);
        it++;
    }
    this->singleSink = false;
}
void CubeComms::addChannelToWorker(unsigned int workerRank)
{

    int rtn = addChannelToMap(workerMap, workerRank);
    if (rtn < 1) {
        ASKAPLOG_WARN_STR(logger, "Adding channel to non-existent worker");

    } else {
        //ASKAPLOG_DEBUG_STR(logger, "Added channel to worker rank " << workerRank);
    }

}
void CubeComms::removeChannelFromWorker(unsigned int workerRank)
{

    int rtn = removeChannelFromMap(workerMap, workerRank);

    if (rtn < 0) {
        ASKAPLOG_WARN_STR(logger, "Removing channel from non-existent writer");
    } else {
        ASKAPLOG_DEBUG_STR(logger, "Removed channel from worker rank" << workerRank);
    }

}
void CubeComms::addChannelToWriter(unsigned int writerRank, unsigned int workerRank)
{

    int rtn = addChannelToMap(writerMap, writerRank);
    if (rtn < 1) {
        ASKAPLOG_WARN_STR(logger, "Adding channel to non-existent writer");

    } else {
        //ASKAPLOG_DEBUG_STR(logger, "Added channel to writer rank " << writerRank);
    }

    //ASKAPLOG_DEBUG_STR(logger, "Adding " << writerRank << " to clientMap");

    this->clientMap[writerRank];
    //ASKAPLOG_DEBUG_STR(logger, "Pushing back " << workerRank);
    this->clientMap[writerRank].push_back(workerRank);
    //ASKAPLOG_DEBUG_STR(logger, "clientMap " << clientMap);
    this->clientMap[writerRank].sort();
    //ASKAPLOG_DEBUG_STR(logger, "Sorted clientMap " << clientMap);
    this->clientMap[writerRank].unique();
    //ASKAPLOG_DEBUG_STR(logger, "Uniquified clientMap " << clientMap);



}
void CubeComms::removeChannelFromWriter(unsigned int writerRank)
{

    int rtn = removeChannelFromMap(writerMap, writerRank);

    if (rtn < 0) {
        ASKAPLOG_WARN_STR(logger, "Removing channel from non-existent writer");
    } else {
        ASKAPLOG_DEBUG_STR(logger, "Removed channel from writer rank" << writerRank);
    }

}
int CubeComms::removeChannelFromMap(std::map<int, int>& theMap, unsigned int theRank)
{

    std::map<int, int>::iterator it;
    it = theMap.find(theRank);
    if (it != theMap.end()) {
        int oldcount = it->second;
        int newcount = oldcount - 1;
        theMap.erase(it);
        std::pair<std::map<int, int>::iterator, bool> ret;
        ret = theMap.insert(std::pair<int, int> (theRank, newcount));
        return newcount;
    } else {
        return -1;

    }

}
int CubeComms::addChannelToMap(std::map<int, int>& theMap, unsigned int theRank)
{

    std::map<int, int>::iterator it;
    it = theMap.find(theRank);
    if (it != theMap.end()) {
        int oldcount = it->second;
        int newcount = oldcount + 1;
        theMap.erase(it);
        std::pair<std::map<int, int>::iterator, bool> ret;
        ret = theMap.insert(std::pair<int, int> (theRank, newcount));
        //ASKAPLOG_DEBUG_STR(logger, "added a channel to rank " << theRank  \
        //                   << " old count was " << oldcount << " current count is " << newcount);
        return newcount;

    } else {
        return -1;

    }

}
/// Dont like putting this here - would rather make the MPIComms version public
void* CubeComms::addByteOffset(const void *ptr, size_t offset) const
{
    char *cptr = static_cast<char*>(const_cast<void*>(ptr));
    cptr += offset;

    return cptr;
}
int CubeComms::getOutstanding()
{
    return writerMap[rank()];
}
std::list<int> CubeComms::getClients()
{
    // return list of clients with at least one piece of work outstanding
    std::list<int> clients_with_work;
    std::list<int>::iterator it;
    it = clientMap[rank()].begin();

    while (it != clientMap[rank()].end()) {
        if (workerMap[*it] > 0) {
            if (*it != rank())
                clients_with_work.push_back(*it);
        }
        it++;
    }
    return clients_with_work;
}

} // end namespace cp
} // end namespace askap
