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
                                       const int cfRank,
                                       const bool mpipresetup) :
        WProjectVisGridder(wmax, nwplanes, cutoff,overSample,maxSupport,limitSupport,name,alpha,shareCF),
        itsMpiMemPreSetup(mpipresetup), itsCFRank(cfRank)
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
        itsMpiMemPreSetup(other.itsMpiMemPreSetup),
        itsCFRank(other.itsCFRank)
{
	std::lock_guard<std::mutex> lk(ObjCountMutex);
	ObjCount += 1;
}


/// Clone a copy of this Gridder
IVisGridder::ShPtr MPIWProjectVisGridder::clone()
{
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

    // generating the convolution plane for the cache can be time consumming so the idea
    // here is to let each rank within a node to generate a portion the convoulution cache
    // (ie each rank does a portion of nWPlanes()*itsOverSample*itsOverSample)
    ASKAPCHECK(itsCFRank > 0,"CF Rank (i.e itsCFRank) is <= 0");
    // itsCFRank => number of ranks that participate in the CF calculation
    if ( itsNodeRank < itsCFRank ) {
        int numberOfPlanePerNodeRank;
        int startPlane;
        int endPlane;
        int rem = nWPlanes() % itsCFRank;
        if ( rem == 0 ) {
            numberOfPlanePerNodeRank = nWPlanes()/itsCFRank;
            startPlane = itsNodeRank * numberOfPlanePerNodeRank;
            endPlane = numberOfPlanePerNodeRank*(itsNodeRank + 1);
        } else {
            numberOfPlanePerNodeRank = (nWPlanes() - rem)/itsCFRank;
            if ( itsNodeRank == itsCFRank - 1 ) {
                startPlane = numberOfPlanePerNodeRank*(itsCFRank - 1);
                endPlane = nWPlanes();
            } else {
                startPlane = itsNodeRank * numberOfPlanePerNodeRank;
                endPlane = (itsNodeRank * numberOfPlanePerNodeRank) + numberOfPlanePerNodeRank;
            }
        }

        ASKAPLOG_INFO_STR(logger,"itsNodeRank: " << itsNodeRank 
                            << ", itsCFRank: " << itsCFRank
                            << ", startPlane: " << startPlane
                            << ", endPlane: " << endPlane);

        // Only ranks smaller than itsCFRank participate in the CF calculation
        generate(startPlane,endPlane);
    }
    // ranks > itsCFRank of a given node wait here
    MPI_Barrier(itsNodeComms);
    //MPI_Barrier(itsNonRankZeroComms);

    // calculate the total size of the matrix data in itsConvFunc vector. the itsConvFunc vector
    // of each rank within the node only contains a partial of the data.
    size_t total = 0;
    for ( auto it = itsConvFunc.begin();
        it != itsConvFunc.end(); ++it) {
        total += it->nelements() * sizeof(imtypeComplex);
    }
    // now collect the total number of bytes from all the ranks to create the shared memory
    unsigned long totalFromAllRanks = 0;
    MPI_Reduce(&total,&totalFromAllRanks,1,MPI_UNSIGNED_LONG,MPI_SUM,0,itsNodeComms);
    if ( ! itsMpiMemPreSetup ) {
        std::lock_guard<std::mutex> lk(ObjCountMutex);
        setupMpiMemory(totalFromAllRanks);
    } else {
        ASKAPLOG_INFO_STR(logger, "itsNodeRank: " << itsNodeRank << ".  memory has been presetup.");
    }
    // Copy the offsets from rank 0 to other ranks on a per node basis
    if (isOffsetSupportAllowed()) {
        copyConvFuncOffset();
    }

    MPI_Barrier(itsNodeComms);
    //MPI_Barrier(itsNonRankZeroComms);
    // Save the CF to the cache
    if (itsShareCF) {
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
    const int cfRank = parset.getInt32("cfrank",1);
    const bool variablesupport  = parset.getBool("variablesupport", false);

    // if we split CF calculation among ranks, variablesupport must be true
    // otherwise we get this error :
    // Askap error in: Peak of PSF(0) is at [1, 5]: not at centre pixel: [1024,1024] (thrown in yandasoft/askap/deconvolution/DeconvolverBase.tcc:99)
    if ( cfRank > 1 ) {
        ASKAPCHECK(variablesupport, "Can only do CF calculation among ranks if variablesupport = true");
    }
    

    ASKAPLOG_INFO_STR(logger, "Gridding using maxsupport: " << maxSupport );
    ASKAPLOG_INFO_STR(logger, "Using " << (useDouble ? "double":"single")<<
                      " precision to calculate convolution functions");
    boost::shared_ptr<MPIWProjectVisGridder> gridder(new MPIWProjectVisGridder(wmax, nwplanes,
            cutoff, oversample, maxSupport, limitSupport, tablename, alpha, useDouble,cfRank,mpipresetup));
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
    //ASKAPLOG_INFO_STR(logger, "configureGridder");
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


    MPI_Aint memSizeInBytes = 0;
    if ( itsNodeRank == 0 ) {
        memSizeInBytes = bufferSize;
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
    if ( itsNodeRank != 0 ) {
        int r = MPI_Win_shared_query(itsWindowTable, 0, &itsWindowSize, &itsWindowDisp, &itsMpiSharedMemory);
        ASKAPCHECK(r == MPI_SUCCESS, "MPI_Win_shared_query failed.");
    }
        
    MPI_Barrier(itsNodeComms);
}

void MPIWProjectVisGridder::copyToSharedMemory(std::vector<std::pair<int,int>>& itsConvFuncMatSize)
{
    // itsConvFuncMatSize keeps an array of pairs whose values are number of rows and columns
    // of the matrixes of the itsConvFunc vector.
    unsigned int numberOfElements = itsConvFunc.size();;
    itsConvFuncMatSize.resize(numberOfElements);

    ASKAPLOG_INFO_STR(logger, "copyToSharedMemory itsNodeRank: " << itsNodeRank 
                            << ", number of planes: " << numberOfElements);
    // the itsConvFunc vector of each rank only has a portion of the data
    for (unsigned int iw = 0; iw < numberOfElements; iw++) {
        ASKAPLOG_DEBUG_STR(logger, "copyToSharedMemory itsNodeRank: " << itsNodeRank << 
                            " - " << iw << "/" << itsConvFuncMatSize.size());
        unsigned long matrixSize[3] = {0,0,0};
        // rank that has data (i.e casacore::Matrix nelements() != 0, sends
        // its plane index, nrow and ncol to other ranks
        if (itsConvFunc[iw].nelements() != 0) {
            matrixSize[0] = iw;
            matrixSize[1] = itsConvFunc[iw].nrow();
            matrixSize[2] = itsConvFunc[iw].ncolumn();
            itsConvFuncMatSize[matrixSize[0]] = std::make_pair(matrixSize[1],matrixSize[2]);
            for (int rank=0; rank < itsNodeSize; rank++) {
                if ( rank != itsNodeRank ) { // dont send to ourselves
                    MPI_Send(matrixSize,3,MPI_UNSIGNED_LONG,rank,0,itsNodeComms);
                }
            }
        } else {
            // Since we dont know the sender, we use MPI_ANY_SOURCE
            MPI_Recv(matrixSize,3,MPI_UNSIGNED_LONG,MPI_ANY_SOURCE,0,itsNodeComms,MPI_STATUS_IGNORE);
            if ( matrixSize[1] != 0 && matrixSize[2] != 0 ) {
                itsConvFuncMatSize[matrixSize[0]] = std::make_pair(matrixSize[1],matrixSize[2]);
            }
        }
    }
    //MPI_Barrier(itsNodeComms);
    // if we get here, itsConvFuncMatSize of each rank knows the nrow and ncol and contains the full
    // size of the itsConvFunc vector

    // now each rank copies its CF to the shared memory
    imtypeComplex* shareMemPtr = itsMpiSharedMemory; // a contiguous chunk of shared memory
    for (unsigned int iw = 0; iw < itsConvFuncMatSize.size(); iw++) {
        if (itsConvFunc[iw].nelements() != 0) {
            std::copy(itsConvFunc[iw].data(),itsConvFunc[iw].data() + itsConvFunc[iw].nelements(),shareMemPtr);
            // std::memcpy(shareMemPtr,itsConvFunc[iw].data(),sizeof(imtypeComplex) * itsConvFunc[iw].nelements());
        }
        shareMemPtr += itsConvFuncMatSize[iw].first *  itsConvFuncMatSize[iw].second;;
    }
    ASKAPLOG_INFO_STR(logger, "copyToSharedMemory itsNodeRank: " << itsNodeRank << " - DONE");
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

void MPIWProjectVisGridder::copyConvFuncOffset()
{
    for (int nw=0; nw<nWPlanes(); nw++) {
        std::pair<int,int> offset = getConvFuncOffset(nw);
        if ( offset.first != 0 || offset.second != 0 ) {
            int offsetPerPlane[3] = {-1,-1,-1};
            offsetPerPlane[0] = nw;
            offsetPerPlane[1] = offset.first;
            offsetPerPlane[2] = offset.second;
            for (unsigned int rank = 0; rank < itsNodeSize; rank++) {
                MPI_Bcast(offsetPerPlane,3,MPI_INT,rank,itsNodeComms);
                setConvFuncOffset(offsetPerPlane[0],offsetPerPlane[1],offsetPerPlane[2]);
            }
        }
    }
}

} // namespace askap
} // namespace synthesis
