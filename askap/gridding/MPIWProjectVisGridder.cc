/// @file WProjectVisGridder.cc
///
/// @copyright (c) 2007,2016 CSIRO
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

// Package level header file
#include <askap/askap_synthesis.h>

// System includes
#include <cmath>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askap/StatReporter.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <askap/scimath/fft/FFTWrapper.h>
#include <askap/profile/AskapProfiler.h>

// Local package includes
#include <askap/gridding/MPIWProjectVisGridder.h>
#include <askap/gridding/SupportSearcher.h>

ASKAP_LOGGER(logger, ".gridding.mpiwprojectvisgridder");

namespace askap {
namespace synthesis {

// Initialise of class static variables
MPI_Aint MPIWProjectVisGridder::itsWindowSize;
int      MPIWProjectVisGridder::itsWindowDisp;
MPI_Win  MPIWProjectVisGridder::itsWindowTable = MPI_WIN_NULL;
MPI_Comm MPIWProjectVisGridder::itsNodeComms;
MPI_Comm MPIWProjectVisGridder::itsNonRankZeroComms;
MPI_Group MPIWProjectVisGridder::itsWorldGroup = MPI_GROUP_NULL;
MPI_Group MPIWProjectVisGridder::itsGridderGroup = MPI_GROUP_NULL;
int MPIWProjectVisGridder::itsNodeSize;
int MPIWProjectVisGridder::itsNodeRank;
int MPIWProjectVisGridder::itsWorldRank;
imtypeComplex* MPIWProjectVisGridder::itsMpiSharedMemory = nullptr;
bool MPIWProjectVisGridder::itsMpiMemSetup = false;
unsigned int MPIWProjectVisGridder::ObjCount = 0;
std::mutex MPIWProjectVisGridder::ObjCountMutex;

MPIWProjectVisGridder::MPIWProjectVisGridder(const double wmax,
                                       const int nwplanes,
                                       const double cutoff,
                                       const int overSample,
                                       const int maxSupport,
                                       const int limitSupport,
                                       const std::string& name,
                                       const float alpha,
                                       const bool shareCF,
                                       const bool mpipresetup) :
        WProjectVisGridder(wmax, nwplanes, cutoff,overSample,maxSupport,limitSupport,name,alpha,shareCF),
        itsMpiMemPreSetup(mpipresetup)
{
    ASKAPCHECK(overSample > 0, "Oversampling must be greater than 0");
    ASKAPCHECK(maxSupport > 0, "Maximum support must be greater than 0")
    
    std::lock_guard<std::mutex> lk(ObjCountMutex);
    ObjCount += 1;
}

MPIWProjectVisGridder::~MPIWProjectVisGridder()
{
    std::lock_guard<std::mutex> lk(ObjCountMutex);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ASKAPLOG_DEBUG_STR(logger,"~MPIWProjectVisGridder ObjCount: " << ObjCount << ", itsNodeRank: " << itsNodeRank << ", rank: " << rank);
    ObjCount -= 1;

    if ( ObjCount == 0 ) {
        if ( itsWorldGroup != MPI_GROUP_NULL )
            MPI_Group_free(&itsWorldGroup);

        if ( itsGridderGroup != MPI_GROUP_NULL )
            MPI_Group_free(&itsGridderGroup);

        if ( itsNodeComms != MPI_COMM_NULL )
            MPI_Comm_free(&itsNodeComms);

        if ( itsNonRankZeroComms != MPI_COMM_NULL )
            MPI_Comm_free(&itsNonRankZeroComms);

        if ( itsWindowTable != MPI_WIN_NULL )
            MPI_Win_free(&itsWindowTable);

        itsMpiMemSetup = false;
    }

}

/// @brief copy constructor
/// @details It is required to decouple internal arrays in the input
/// object and the copy
/// @param[in] other input object
MPIWProjectVisGridder::MPIWProjectVisGridder(const MPIWProjectVisGridder &other) :
        IVisGridder(other), WProjectVisGridder(other),
        itsMpiMemPreSetup(other.itsMpiMemPreSetup)
{
	std::lock_guard<std::mutex> lk(ObjCountMutex);
    ASKAPLOG_DEBUG_STR(logger, "copy constructor");
	ObjCount += 1;
}


/// Clone a copy of this Gridder
IVisGridder::ShPtr MPIWProjectVisGridder::clone()
{
    ASKAPLOG_DEBUG_STR(logger, "clone()");
    return IVisGridder::ShPtr(new MPIWProjectVisGridder(*this));
}

/// Initialize the convolution function into the cube. If necessary this
/// could be optimized by using symmetries.
void MPIWProjectVisGridder::initConvolutionFunction(const accessors::IConstDataAccessor&)
{
    ASKAPTRACE("MPIWProjectVisGridder::initConvolutionFunction");
    /// We have to calculate the lookup function converting from
    /// row and channel to plane of the w-dependent convolution
    /// function

    if (itsSupport > 0) {
        return;
    }

    itsSupport = 0;

    if (isOffsetSupportAllowed()) {
        // this command is executed only once when itsSupport is not set.
        initConvFuncOffsets(nWPlanes());
    }

    if (isPCFGridder()) {
        doPCFGridder();
        return;
    }

    ASKAPLOG_DEBUG_STR(logger, "theirCFCache size: " << theirCFCache.size() << ", itsShareCF: " << itsShareCF);
    if (itsShareCF && theirCFCache.size()>0) {
        // we already have what we need
        itsSupport = 1;
        size_t size = deepRefCopyOfSTDVector(WProjectVisGridder::theirCFCache,itsConvFunc)/1024/1024;
        ASKAPLOG_INFO_STR(logger, "Using cached convolution functions ("<<size<<" MB)");
        if (isOffsetSupportAllowed()) {
            for (size_t i=0; i<WProjectVisGridder::theirConvFuncOffsets.size(); i++) {
                setConvFuncOffset(i,WProjectVisGridder::theirConvFuncOffsets[i].first,WProjectVisGridder::theirConvFuncOffsets[i].second);
            }
        }
        return;
    }


    // only rank 0 of each node does the calculation and then copy itsConvFunc to the MPI
    // shared memory then other ranks within the node fill their itsConvFunc vector from the
    // shared memory.
    if ( itsNodeRank == 0 ) {
        generate();
    } // End of rank 0 in a given node
    // ranks > 0 of a given node wait here
    //MPI_Barrier(itsNodeComms);
    MPI_Barrier(itsNonRankZeroComms);

    // Copy the offsets from rank 0 to other ranks on a per node basis
    if (isOffsetSupportAllowed()) {
        for (int nw=0; nw<nWPlanes(); nw++) {
            int offsetPerPlane[3] = {-1,-1,-1};
            if ( itsNodeRank == 0 ) {
                std::pair<int,int> offset = getConvFuncOffset(nw);
                offsetPerPlane[0] = nw;
                offsetPerPlane[1] = offset.first;
                offsetPerPlane[2] = offset.second;
            }
            MPI_Bcast(offsetPerPlane,3,MPI_INT,0,itsNodeComms);
            if ( itsNodeRank != 0 ) {
                setConvFuncOffset(offsetPerPlane[0],offsetPerPlane[1],offsetPerPlane[2]);
            }
        }
    }

    //MPI_Barrier(itsNodeComms);
    MPI_Barrier(itsNonRankZeroComms);
    // Save the CF to the cache
    if (itsShareCF) {
        if ( itsNodeRank == 0 ) {
            ASKAPLOG_DEBUG_STR(logger, "Copy to shared memory etc ...");
            ASKAPLOG_DEBUG_STR(logger, "number of elements in itsConvFunc: " << itsConvFunc.size());
        }

        size_t total = 0; // in bytes
        if ( itsNodeRank == 0 ) {
            // workout the size of itsConvFunc
            for (auto it = itsConvFunc.begin();
                it != itsConvFunc.end(); ++it) {
                total += it->nelements() * sizeof(imtypeComplex);
            }
        } 
        // create MPI shared memory
	    MPI_Barrier(itsNodeComms);
        // only rank 0 has a total value > 0 and other ranks have a value of 0 but there is
        // nothing wrong with this
        if ( ! itsMpiMemPreSetup )
        {
            std::lock_guard<std::mutex> lk(ObjCountMutex);
            setupMpiMemory(total);
        } else {
            ASKAPLOG_INFO_STR(logger, "itsNodeRank: " << itsNodeRank << ".  memory has been presetup.");
        }

        // copy itsConvFunc to shared memory
        // itsConvFuncMatSize keeps an array of pairs whose values are number of rows and columns
        // of the matrixes of the itsConvFunc vector. This variable is required by ranks > 0 because
        // their itsConvFunc is empty up until now.
        std::vector<std::pair<int,int> > itsConvFuncMatSize;
        copyToSharedMemory(itsConvFuncMatSize);
	    MPI_Barrier(itsNodeComms);

    	//ASKAPLOG_DEBUG_STR(logger, "copy shared memory back to itsConvFunc");
        copyFromSharedMemory(itsConvFuncMatSize);

        //MPI_Barrier(itsNodeComms);
        // this is the original code from WProjectVisGridder
        total = deepRefCopyOfSTDVector(itsConvFunc,WProjectVisGridder::theirCFCache);
        if (isOffsetSupportAllowed()) {
            theirConvFuncOffsets.resize(nWPlanes());
            for (int nw=0; nw<nWPlanes(); nw++) {
                theirConvFuncOffsets[nw]=getConvFuncOffset(nw);
            }
        }
    }
    //ASKAPLOG_INFO_STR(logger, "itsWorldRank: " << itsWorldRank << ", itsNodeRank: " << itsNodeRank << ".  Waiting for all ranks to get here");
    //MPI_Barrier(itsNonRankZeroComms);
    //ASKAPLOG_INFO_STR(logger, "itsWorldRank: " << itsWorldRank << ", itsNodeRank: " << itsNodeRank << ".  All ranks arrived here");
}

/// @brief static method to create gridder
/// @details Each gridder should have a static factory method, which is
/// able to create a particular type of the gridder and initialise it with
/// the parameters taken form the given parset. It is assumed that the
/// method receives a subset of parameters where the gridder name is already
/// taken out.
/// @param[in] parset input parset file
/// @return a shared pointer to the gridder instance
IVisGridder::ShPtr MPIWProjectVisGridder::createGridder(const LOFAR::ParameterSet& parset)
{

    const double wmax = parset.getDouble("wmax", 35000.0);
    const int nwplanes = parset.getInt32("nwplanes", 65);
    const double cutoff = parset.getDouble("cutoff", 1e-3);
    const int oversample = parset.getInt32("oversample", 8);
    const int maxSupport = parset.getInt32("maxsupport", 256);
    const int limitSupport = parset.getInt32("limitsupport", 0);
    const string tablename = parset.getString("tablename", "");
    const float alpha=parset.getFloat("alpha", 1.);
    const bool useDouble = parset.getBool("usedouble",false);

    const bool mpipresetup = parset.getBool("mpipresetup", false);

    ASKAPLOG_INFO_STR(logger, "---> Gridding using W projection with " << nwplanes << " w-planes");
    ASKAPLOG_INFO_STR(logger, "Gridding using maxsupport: " << maxSupport );
    ASKAPLOG_INFO_STR(logger, "Using " << (useDouble ? "double":"single")<<
                      " precision to calculate convolution functions");
    boost::shared_ptr<MPIWProjectVisGridder> gridder(new MPIWProjectVisGridder(wmax, nwplanes,
            cutoff, oversample, maxSupport, limitSupport, tablename, alpha, useDouble,mpipresetup));
    gridder->configureGridder(parset);
    gridder->configureWSampling(parset);

    return gridder;
}

/// @brief additional operations to configure gridder
/// @details This method is supposed to be called from createGridder and could be
/// used in derived classes to avoid too much duplication of the code. For this
/// particular class it configures variable/offset support and cutoff behavior.
/// @param[in] parset input parset file
void MPIWProjectVisGridder::configureGridder(const LOFAR::ParameterSet& parset)
{
    std::lock_guard<std::mutex> lk(ObjCountMutex);
    ASKAPLOG_INFO_STR(logger, "configureGridder");
    const bool planeDependentSupport = parset.getBool("variablesupport", false);

    WProjectVisGridder::planeDependentSupport(planeDependentSupport);

    const bool offsetSupport = parset.getBool("offsetsupport", false);

    ASKAPCHECK((!offsetSupport && !planeDependentSupport) || planeDependentSupport,
               "offsetsupport option of the gridder should only be used together with variablesupport option");
    WProjectVisGridder::offsetSupport(offsetSupport);

    const bool absCutoff = parset.getBool("cutoff.absolute", false);

    if (absCutoff) {
        ASKAPLOG_INFO_STR(logger, "Cutoff value of " << itsCutoff <<
            " will be treated as an absolute threshold during CF generation");
    } else {
        ASKAPLOG_INFO_STR(logger, "Cutoff value of " << itsCutoff <<
            " will be treated as a threshold relative to the peak during CF generation");
        ASKAPCHECK(itsCutoff > 0.0, "Cutoff must be positive");
        ASKAPCHECK(itsCutoff < 1.0, "Cutoff must be less than 1.0");
    }

    setAbsCutoffFlag(absCutoff);

    itsShareCF = parset.getBool("sharecf",false);

    //std::lock_guard<std::mutex> lk(ObjCountMutex);

    MPI_Comm_rank(MPI_COMM_WORLD, &itsWorldRank);

    // Each rank should only run the code below once.
    if ( ObjCount > 1 ) return;

    ASKAPLOG_INFO_STR(logger, "Setup MPI subgroup communicator for MPI Shared Memory");

    // The code below setup a MPI communicator for each node where the 
    // rank 0 (in the COMM_WORLD) of the first node is not included.
    int r;
    r = MPI_Comm_group(MPI_COMM_WORLD, &itsWorldGroup);
    ASKAPASSERT(r == MPI_SUCCESS);
    const int exclude_ranks[1] = {0};
    r = MPI_Group_excl(itsWorldGroup, 1, exclude_ranks, &itsGridderGroup);
    ASKAPCHECK(r == MPI_SUCCESS,"rank: " << itsWorldRank << " - MPI_Group_excl() failed");
    r = MPI_Comm_create_group(MPI_COMM_WORLD, itsGridderGroup, 0, &itsNonRankZeroComms);
    ASKAPCHECK(r == MPI_SUCCESS,"rank: " << itsWorldRank << " - MPI_Comm_create_group() failed");
    r = MPI_Comm_split_type(itsNonRankZeroComms, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &itsNodeComms);
	ASKAPCHECK(r == MPI_SUCCESS,"rank: " << itsWorldRank << " - MPI_Comm_split_type() failed");
    // number of ranks within the node
    MPI_Comm_size(itsNodeComms, &itsNodeSize);
    // rank within a node
    MPI_Comm_rank(itsNodeComms, &itsNodeRank);
    ASKAPCHECK(itsNodeComms != MPI_COMM_NULL, "rank: " << itsWorldRank << ", itsNodeRank: " << itsNodeRank << " has itsNodeComms = MPI_COMM_NULL");
    MPI_Barrier(itsNodeComms);

    ASKAPLOG_INFO_STR(logger,"rank: " << itsWorldRank << ", itsNodeSize: " << itsNodeSize << ", itsNodeRank: " << itsNodeRank);

    if ( itsMpiMemPreSetup ) {
        std::string mpiMem = parset.getString("mpipresetmemory","");
        ASKAPCHECK(mpiMem != "", "mpipresetmemory option is not set");
        unsigned long memoryInBytes = std::stoul(mpiMem);;
        ASKAPLOG_INFO_STR(logger,"Presetup the shared memory. " << memoryInBytes);
        setupMpiMemory(memoryInBytes); 
    }
}



/// @brief assignment operator
/// @details Defined as private, so it can't be called (to enforce usage of the
/// copy constructor
/// @param[in] other input object
/// @return reference to itself
MPIWProjectVisGridder& MPIWProjectVisGridder::operator=(const MPIWProjectVisGridder &)
{
    ASKAPTHROW(AskapError, "This method is not supposed to be called!");
    return *this;
}

void  MPIWProjectVisGridder::setupMpiMemory(size_t bufferSize /* in bytes */)
{
    //std::lock_guard<std::mutex> lk(ObjCountMutex);
    if ( itsMpiMemSetup  ) {
	    ASKAPLOG_INFO_STR(logger,"itsNodeRank: " << itsNodeRank << " - mpi shared memory already setup. ObjCount: " << ObjCount);
        return;
    }

    itsMpiMemSetup = true;


    //size_t memSizeInBytes = 0;
    MPI_Aint memSizeInBytes = 0;
    if ( itsNodeRank == 0 ) {
        memSizeInBytes = bufferSize;
        // memSizeInBytes = 41032802304;
    }
    char estring[MPI_MAX_ERROR_STRING];
    int elen = 0;
    ASKAPLOG_INFO_STR(logger,"itsNodeRank: " << itsNodeRank << ", bufferSize: " << memSizeInBytes);
    ASKAPCHECK(itsNodeComms != MPI_COMM_NULL,"itsNodeRank: " << itsNodeRank << " - itsNodeComms is Null");
    int r = MPI_Win_allocate_shared(memSizeInBytes,sizeof(imtypeComplex),
                                MPI_INFO_NULL, itsNodeComms, &itsMpiSharedMemory,
                                &itsWindowTable);
    if ( r != MPI_SUCCESS ) {
        MPI_Error_string(r, estring, &elen);
        ASKAPLOG_INFO_STR(logger,"MPI_Win_allocate_shared - " << estring);    
    }
    ASKAPCHECK(r == MPI_SUCCESS, "itsNodeRank: " << itsNodeRank << " - MPI_Win_allocate_shared() failed.");
    // For itsNodeRanks != 0, get their itsMpiSharedMemory pointer variable to point the
    // start of the MPI shared memory allocated by itsNodeRank = 0
    //MPI_Barrier(itsNodeComms);
    if ( itsNodeRank != 0 ) {
        int r = MPI_Win_shared_query(itsWindowTable, 0, &itsWindowSize, &itsWindowDisp, &itsMpiSharedMemory);
        ASKAPCHECK(r == MPI_SUCCESS, "MPI_Win_shared_query failed.");
    }
        
    MPI_Barrier(itsNodeComms);
	//unsigned long* val;
	//int flag = 1;
	//r = MPI_Win_get_attr(itsWindowTable,MPI_WIN_SIZE,&val,&flag);
}

void MPIWProjectVisGridder::copyToSharedMemory(std::vector<std::pair<int,int>>& itsConvFuncMatSize)
{
    itsConvFuncMatSize.clear();

    imtypeComplex* shareMemPtr = itsMpiSharedMemory; // a contiguous chunk of shared memory
    // itsConvFuncMatSize keeps an array of pairs whose values are number of rows and columns
    // of the matrixes of the itsConvFunc vector. This variable is required by ranks > 0 because
    // their itsConvFunc is empty up until now.
    // only itsNodeRank 0 does the copy 
    unsigned int numberOfElements = 0;

	if ( itsNodeRank == 0 ) {
        ASKAPLOG_DEBUG_STR(logger, "copy itsConvFunc to shared memory");
        numberOfElements = itsConvFunc.size();

        int matrixSize[2];
        for ( auto it = itsConvFunc.begin();
        	it != itsConvFunc.end(); ++it) {
           	std::copy(it->data(),it->data() + it->nelements(),shareMemPtr);
           	shareMemPtr += it->nelements();
       	}
	}
    // Send the nrows and ncolumns of each matrix in the itsConvFunc vector of rank 0
    // to other ranks
    MPI_Bcast(&numberOfElements,1,MPI_UNSIGNED_LONG,0,itsNodeComms);
    ASKAPLOG_DEBUG_STR(logger, "itsNodeRank: " << itsNodeRank << ", numberOfElements: " << numberOfElements);
    for ( unsigned int elem = 0; elem < numberOfElements; elem++) {
        unsigned long matrixSize[2];
        if ( itsNodeRank == 0 ) {
            matrixSize[0] = itsConvFunc[elem].nrow();
            matrixSize[1] = itsConvFunc[elem].ncolumn();
        }
        MPI_Bcast(matrixSize,2,MPI_UNSIGNED_LONG,0,itsNodeComms);
        itsConvFuncMatSize.push_back(std::make_pair(matrixSize[0],matrixSize[1]));
    }
}

void MPIWProjectVisGridder::copyFromSharedMemory(const std::vector<std::pair<int,int>>& itsConvFuncMatSize)
{
    ASKAPLOG_DEBUG_STR(logger, "itsNodeRank: " << itsNodeRank << " - copy shared memory back to itsConvFunc");
    unsigned int numOfElems = itsConvFuncMatSize.size();
    imtypeComplex* shareMemPtr = itsMpiSharedMemory;
    for (unsigned int elem = 0; elem < numOfElems; elem++) {
        casacore::IPosition pos(2);
        pos(0) = itsConvFuncMatSize[elem].first;
        pos(1) = itsConvFuncMatSize[elem].second;
        
        casacore::Matrix<imtypeComplex> m(pos,shareMemPtr,casacore::SHARE);
        itsConvFunc[elem].reference(m);
        shareMemPtr += m.nelements();
    }
}

} // namespace askap
} // namespace synthesis
