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
#include <askap/profile/AskapProfiler.h>

// Local package includes
#include <askap/gridding/MPIWProjectVisGridder.h>
#include <askap/gridding/SupportSearcher.h>

ASKAP_LOGGER(logger, ".gridding.mpiwprojectvisgridder");

namespace askap {
namespace synthesis {
#ifdef HAVE_MPI
// Initialise of class static variables
MPI_Aint MPIWProjectVisGridder::itsWindowSize;
int      MPIWProjectVisGridder::itsWindowDisp;
MPI_Win  MPIWProjectVisGridder::itsWindowTable = MPI_WIN_NULL;
MPI_Comm MPIWProjectVisGridder::itsNodeComms;
MPI_Comm MPIWProjectVisGridder::itsNonRankZeroComms;
MPI_Group MPIWProjectVisGridder::itsWorldGroup = MPI_GROUP_NULL;
MPI_Group MPIWProjectVisGridder::itsGridderGroup = MPI_GROUP_NULL;
#endif
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
                                       const bool masterDoesWork,
                                       const bool mpipresetup) :
        WProjectVisGridder(wmax, nwplanes, cutoff,overSample,maxSupport,limitSupport,name,alpha,shareCF),
        itsMpiMemPreSetup(mpipresetup), itsCFRank(cfRank), itsSerial(false), itsMasterDoesWork(masterDoesWork)
{
#ifdef HAVE_MPI
    ASKAPCHECK(overSample > 0, "Oversampling must be greater than 0");
    ASKAPCHECK(maxSupport > 0, "Maximum support must be greater than 0")

    std::lock_guard<std::mutex> lk(ObjCountMutex);
    ObjCount += 1;
#else
    ASKAPTHROW(AskapError, "Cant use MPIWProject gridder without MPI library");
#endif
}

MPIWProjectVisGridder::~MPIWProjectVisGridder()
{
#ifdef HAVE_MPI
    std::lock_guard<std::mutex> lk(ObjCountMutex);
    ObjCount -= 1;

    ASKAPLOG_DEBUG_STR(logger,"MPIWProjectVisGridder::destructor. itsSerial: " << itsSerial << ", ObjCount: " << ObjCount);
    if ( ObjCount == 0 && !itsSerial ) {
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
#endif
}

/// @brief copy constructor
/// @details It is required to decouple internal arrays in the input
/// object and the copy
/// @param[in] other input object
MPIWProjectVisGridder::MPIWProjectVisGridder(const MPIWProjectVisGridder &other) :
        IVisGridder(other), WProjectVisGridder(other),
        itsMpiMemPreSetup(other.itsMpiMemPreSetup),
        itsCFRank(other.itsCFRank),
        itsSerial(other.itsSerial),
        itsMasterDoesWork(other.itsMasterDoesWork)
{
// the ifdef here is just for completeness and is not needed because the MPIWProjectVisGridder
// object cant be created without having MPI in the first place
#ifdef HAVE_MPI
	std::lock_guard<std::mutex> lk(ObjCountMutex);
	ObjCount += 1;
#else
    ASKAPTHROW(AskapError, "Cant use MPIWProject gridder without MPI library");
#endif
}


/// Clone a copy of this Gridder
IVisGridder::ShPtr MPIWProjectVisGridder::clone()
{
#ifdef HAVE_MPI
    return IVisGridder::ShPtr(new MPIWProjectVisGridder(*this));
#else
    ASKAPTHROW(AskapError, "Cant use MPIWProject gridder without MPI library");
    return IVisGridder::ShPtr(new MPIWProjectVisGridder(*this));
#endif
}

/// Initialize the convolution function into the cube. If necessary this
/// could be optimized by using symmetries.
void MPIWProjectVisGridder::initConvolutionFunction(const accessors::IConstDataAccessor& acc)
{
#ifdef HAVE_MPI
    ASKAPTRACE("MPIWProjectVisGridder::initConvolutionFunction");

    if ( itsSerial ) {
        ASKAPLOG_DEBUG_STR(logger,"MPI WPRoject gridder runs in serial mode. Delegate the call to WProjectVisGridder::initConvolutionFunction");
        WProjectVisGridder::initConvolutionFunction(acc);
    } else {
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
        ASKAPCHECK(itsCFRank <= itsNodeSize,"CF Rank (i.e itsCFRank) is more than the number of ranks per node");
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
        ASKAPLOG_DEBUG_STR(logger,"itsNodeRank: " << itsNodeRank << " waits at mpi barrier");
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
#else
    ASKAPTHROW(AskapError, "Cant use MPIWProject gridder without MPI library");
#endif
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

    const bool masterDoesWork = parset.getBool("masterDoesWork",false);

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
            cutoff, oversample, maxSupport, limitSupport, tablename, alpha, useDouble,cfRank,masterDoesWork,mpipresetup));
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
#ifdef HAVE_MPI
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

    // use the MPI_Initialized() function to test if the gridder is running in parallel or serial
    // assuming that all multi ranks MPI task must first call MPI_init() and so if the MPI_Initialized()
    // function returns false (0) then the task is not an MPI task
    itsSerial = false;
    int mpiInitalised = 0;
    MPI_Initialized(&mpiInitalised);
    if ( mpiInitalised == 1 ) {
        int numRanks = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        if ( numRanks == 1 ) {
            // MPI has only one rank so it must be running in serial
            ASKAPLOG_INFO_STR(logger, "MPI WProject gridder runs in serial");
            itsSerial = true;
        }
    }

    if ( !itsSerial ) {
        ASKAPLOG_INFO_STR(logger, "MPI WProject gridder runs in parallel");
        // Each rank should only run the code below once.
        if ( ObjCount > 1 ) return;

        ASKAPLOG_INFO_STR(logger, "Setup MPI subgroup communicator for MPI Shared Memory");

        // The code below setup a MPI communicator for each node where the
        // rank 0 (in the COMM_WORLD) of the first node is not included.
        int r;
        r = MPI_Comm_group(MPI_COMM_WORLD, &itsWorldGroup);
        ASKAPASSERT(r == MPI_SUCCESS);
        ASKAPLOG_DEBUG_STR(logger,"itsMasterDoesWork: " << itsMasterDoesWork);
        if ( itsMasterDoesWork )  {
            // We have a master that wants to get its hand dirty and do some work. Good on you
            MPI_Comm_group(MPI_COMM_WORLD,&itsGridderGroup);
        } else {
            // We have a lazy master that likes to sit around and farm out the work and get the
            // reward and credit of its workers
            const int exclude_ranks[1] = {0};
            r = MPI_Group_excl(itsWorldGroup, 1, exclude_ranks, &itsGridderGroup);
        }

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
    } else {
        // The gridder is running in serial so it should behave as if it is a WProject gridder object
        ASKAPLOG_INFO_STR(logger,"MPI gridder is running in serial so it should behave as if it is a WProject gridder");
    }
#else
    ASKAPTHROW(AskapError, "Cant use MPI WProject gridder without MPI library");
#endif
}



/// @brief assignment operator
/// @details Defined as private, so it can't be called (to enforce usage of the
/// copy constructor
/// @param[in] other input object
/// @return reference to itself
//MPIWProjectVisGridder& MPIWProjectVisGridder::operator=(const MPIWProjectVisGridder &)
//{
//    ASKAPTHROW(AskapError, "This method is not supposed to be called!");
//    return *this;
//}

void  MPIWProjectVisGridder::setupMpiMemory(size_t bufferSize /* in bytes */)
{
#ifdef HAVE_MPI
    //std::lock_guard<std::mutex> lk(ObjCountMutex);
    if ( itsMpiMemSetup  ) {
	    ASKAPLOG_INFO_STR(logger,"itsNodeRank: " << itsNodeRank << " - mpi shared memory already setup. ObjCount: " << ObjCount);
        return;
    }

    ASKAPLOG_INFO_STR(logger,"itsNodeRank: " << itsNodeRank << " - Setup mpi shared memory. ObjCount: " << ObjCount);

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
#else
    ASKAPTHROW(AskapError, "Cant use MPIWProject gridder without MPI library" );
#endif
}

void MPIWProjectVisGridder::copyToSharedMemory(std::vector<std::pair<int,int>>& itsConvFuncMatSize)
{
#ifdef HAVE_MPI
    // itsConvFuncMatSize keeps an array of pairs whose values are number of rows and columns
    // of the matrixes of the itsConvFunc vector.
    unsigned int numberOfElements = itsConvFunc.size();;
    itsConvFuncMatSize.resize(numberOfElements);

    ASKAPLOG_INFO_STR(logger, "copyToSharedMemory itsNodeRank: " << itsNodeRank
                            << ", number of planes: " << numberOfElements);
    // the itsConvFunc vector of each rank only has a portion of the data
    // matrixSizeVect variable stores the planes and shapes of the matrices.
    std::vector<std::tuple<unsigned long, unsigned long, unsigned long>> matrixSizeVect;
    for (unsigned int iw = 0; iw < numberOfElements; iw++) {
        if (itsConvFunc[iw].nelements() != 0) {
            matrixSizeVect.push_back(std::make_tuple(iw,itsConvFunc[iw].nrow(),itsConvFunc[iw].ncolumn()));
        }
    }
    unsigned long howmany = matrixSizeVect.size();
    ASKAPLOG_DEBUG_STR(logger, "copyToSharedMemory itsNodeRank: " << itsNodeRank <<
                        " - size of matrixSizeVect = " << howmany);

    // for each rank, broadcast its matrixSizeVect to other ranks
    // First, broadcast the size of the matrixSizeVect and if the size
    // is not 0, then broadcast the payload (i.e the content of matrixSizeVect)
    // as one large message.
    for (int rank=0; rank < itsNodeSize; rank++) {
        unsigned long len = 0;
        unsigned long* buffer = nullptr;
        if ( rank == itsNodeRank ) {
            len = howmany;
        }

        MPI_Bcast(&len,1,MPI_UNSIGNED_LONG,rank,itsNodeComms);
        if ( len != 0 ) {
            unsigned long* buffer = new unsigned long[len * 3];
            if ( rank == itsNodeRank ) {
                unsigned int count = 0;
                for (auto i = 0; i < len; i++) {
                    buffer[count] = std::get<0>(matrixSizeVect[i]);
                    buffer[count+1] = std::get<1>(matrixSizeVect[i]);
                    buffer[count+2] = std::get<2>(matrixSizeVect[i]);
                    count += 3;
                }
            }
            MPI_Bcast(buffer,len*3,MPI_UNSIGNED_LONG,rank,itsNodeComms);
            for (unsigned long i = 0; i < len*3; i += 3) {
                itsConvFuncMatSize[buffer[i]] = std::make_pair(buffer[i+1],buffer[i+2]);
            }

            delete []buffer;
        }
    }
    // if we get here, the itsConvFuncMatSize variable contains the shapes of all the matrices
    // from 0 to (nwplane*oversample^2 - 1)
    ASKAPLOG_DEBUG_STR(logger, "copyToSharedMemory itsNodeRank: " << itsNodeRank <<
                        " writes/saves itsConvFunc data to shared memory");
    //MPI_Barrier(itsNodeComms);

    // now each rank copies its CF to the shared memory
    imtypeComplex* shareMemPtr = itsMpiSharedMemory; // a contiguous chunk of shared memory
    for (unsigned int iw = 0; iw < itsConvFuncMatSize.size(); iw++) {
        if (itsConvFunc[iw].nelements() != 0) {
            std::copy(itsConvFunc[iw].data(),itsConvFunc[iw].data() + itsConvFunc[iw].nelements(),shareMemPtr);
        }
        shareMemPtr += itsConvFuncMatSize[iw].first *  itsConvFuncMatSize[iw].second;;
    }
    ASKAPLOG_INFO_STR(logger, "copyToSharedMemory itsNodeRank: " << itsNodeRank << " - DONE");
#else
    ASKAPTHROW(AskapError, "Cant use MPIWProject gridder without MPI library" );
#endif
}

void MPIWProjectVisGridder::copyFromSharedMemory(const std::vector<std::pair<int,int>>& itsConvFuncMatSize)
{
    ASKAPLOG_DEBUG_STR(logger, "itsNodeRank: " << itsNodeRank <<
                                " - copy shared memory back to itsConvFunc. number of CFs = " << itsConvFuncMatSize.size());
    unsigned int numOfElems = itsConvFuncMatSize.size();
    imtypeComplex* shareMemPtr = itsMpiSharedMemory;
    for (unsigned int elem = 0; elem < numOfElems; elem++) {
        casacore::IPosition pos(2);
        pos(0) = itsConvFuncMatSize[elem].first;
        pos(1) = itsConvFuncMatSize[elem].second;
        ASKAPCHECK(itsConvFuncMatSize[elem].first != 0 || itsConvFuncMatSize[elem].second != 0, "shape of CF is zero");

        casacore::Matrix<imtypeComplex> m(pos,shareMemPtr,casacore::SHARE);
        itsConvFunc[elem].reference(m);
        shareMemPtr += m.nelements();
    }
}

void MPIWProjectVisGridder::copyConvFuncOffset()
{
#ifdef HAVE_MPI
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
#else
    ASKAPTHROW(AskapError, "Cant use MPIWProject gridder without MPI library");
#endif
}

} // namespace askap
} // namespace synthesis
