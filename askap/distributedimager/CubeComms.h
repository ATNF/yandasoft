/// @file CubeComms.h
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

#ifndef ASKAP_CP_IMAGER_CUBECOMMS_H
#define ASKAP_CP_IMAGER_CUBECOMMS_H

#include <map>
#include <list>
#ifdef HAVE_MPI
#include "mpi.h"
#endif
///ASKAP includes ...
#include <askap/askapparallel/AskapParallel.h>
#include "askap/messages/IMessage.h"

namespace askap {
namespace cp {

class CubeComms: public askapparallel::AskapParallel {
    public:
        /// @brief Constructor
        /// @details The command line inputs are needed solely for MPI - currently no
        /// application specific information is passed on the command line.
        /// @param argc Number of command line inputs
        /// @param argv Command line inputs
        CubeComms(int argc, const char** argv);

        /// @brief isWriter
        /// @details This bool can be tested to find out whether the current
        /// rank is a writer
        /// returns writerCount 1->nwriters
        int isWriter();
        /// @brief adds a writer to the list by rank
        /// @details This adds a rank to a vector of ranks
        /// each is the rank of a writer
        void addWriter(unsigned int writer_rank);
        /// @brief adds a worker to the list by rank
        /// @details This adds a rank to a vector of ranks
        /// each is the rank of a writer
        void addWorker(unsigned int writer_rank);

        /// @brief isCubeCreator
        /// @details This bool can be tested to find out whether the current
        /// rank creates an output cube
        ///
        bool isCubeCreator();
        /// @brief adds a creator to the list by rank
        /// @details This adds a rank to a vector of ranks
        /// each is the rank of a creator of a cube
        void addCubeCreator(int creator_rank);
        /// @brief Return the list of creator ranks
        std::list<int> getCubeCreators();

        /// @brief initialises the writer list
        /// @details By evenly dividing the writing across
        /// the ranks

        void initWriters(int nwriters, int nchanpercore);
        /// @brief increments a counter (one for each rank)
        /// @details Takes the rank of the writer
        void addChannelToWriter(unsigned int writer_rank, unsigned int worker);
        void addClientToWriter(unsigned int writer_rank, unsigned int client_rank);

        void removeChannelFromWriter(unsigned int writer_rank);

        void addChannelToWorker(unsigned int worker_rank);
        void removeChannelFromWorker(unsigned int worker_rank);

        int getOutstanding();
        std::list<int> getClients();

        int anyWork();

        /// @brief its communicator for its fellow workers
        size_t buildCommIndex();
        /// @brief its communicator for its fellow writers
        size_t buildWriterIndex(size_t comm = 0);

        /// @brief translate vector of ranks from MPI_COMM_WORLD to specified comm
        std::vector<int> translateRanks(std::vector<int> ranks, size_t comm);
      
        /// @brief sets the cubecreator to be the first writer
        void setSingleSink();
        bool isSingleSink()
        {
            return singleSink;
        };
        /// @brief sets the cubecreators to be the all writers
        void setMultiSink();

        size_t theWorkers() {return itsComrades;};
        size_t theBusyComrades() {return itsBusyComrades;};
        size_t theWriters() {return itsWriters;};

        void showWorkerMap();

        ~CubeComms();

    private:
        // Add a byte offset to the  specified pointer, returning the result
        void* addByteOffset(const void *ptr, size_t offset) const;

        // MAP for each writer and the number of work allocations (channels)
        // it has to write
        std::map<int, int> writerMap;
        // MAP for each client with the number of good channels it intends to
        // write

        std::map<int, int > workerMap;

        std::map<int, std::list<int> > clientMap;

        // List of ranks that are creating a cube
        std::list<int> cubeCreators;

        int addChannelToMap(std::map<int, int>& theMap, unsigned int theRank);
        int removeChannelFromMap(std::map<int, int>& theMap, unsigned int theRank);

        int writerCount;
        size_t itsBusyComrades;
        size_t itsComrades;
        size_t itsWriters;
        bool singleSink;
};
}
}
#endif
