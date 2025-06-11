/// @file WProjectVisGridder.h
///
/// WProjectVisGridder: W projection gridding

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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///
#ifndef ASKAP_SYNTHESIS_MPIWPROJECTVISGRIDDER_H_
#define ASKAP_SYNTHESIS_MPIWPROJECTVISGRIDDER_H_

// ASKAPsoft includes
#include <askap/gridding/WProjectVisGridder.h>
#include <askap/askapparallel/AskapParallel.h>

// Local package includes
#include <askap/dataaccess/IConstDataAccessor.h>

namespace askap
{
    namespace synthesis
    {
        /// @brief Visibility gridder using W projection by utilising MPI
        ///        shared memory to reduce memory footprint.
        /// @details The visibilities are gridded using a convolution
        /// function that implements a Fresnel transform. This corrects
        /// for the w term in the full synthesis measurement equation.
        ///
        /// The convolution function is calculated straightforwardly
        /// by constructing an image of the complex w-dependent
        /// phasor and Fourier transforming. The calculation is
        /// done using a coarse but large grid in image space so
        /// that it is sub-sampled in uv space.
        ///
        /// The scaling is slow in data points, fast in w planes.
        ///
        /// NOTE: At the moment, this class only works for the imager
        ///       but not cimager because it is setup in a way that rank
        ///       0 of the first node does not have the gridder.
        ///
        /// @ingroup gridding
        class MPIWProjectVisGridder : public WProjectVisGridder
        {
            public:
                /// @brief Construct a gridder for W projection
                /// @param wmax Maximum baseline (wavelengths)
                /// @param nwplanes Number of w planes
                /// @param cutoff Cutoff in determining support e.g. 10^-3 of the peak
                /// @param overSample Oversampling (currently limited to <=1)
                /// @param maxSupport Maximum support to be allowed
                /// @param limitSupport Upper limit of support
                /// @param name Name of table to save convolution function into
                /// @param alpha spheroidal function alpha value
                /// @param useDouble Use Double precision for the convolution functions
                MPIWProjectVisGridder(/* const askap::askapparallel::AskapParallel& comms, */
                        const double wmax, const int nwplanes,
                        const double cutoff, const int overSample,
                        const int maxSupport, const int limitSupport,
                        const std::string& name=std::string(""),
                        const float alpha=1.,
                        const bool shareCF=false,
                        const int  cfRank = 1,
                        const bool masterDoesWork = false,
                        const bool mpipresetup = false);

                /// @brief dont allow default constructor
                MPIWProjectVisGridder() = delete;
                /// @brief assignment operator
                /// @details Defined as private, so it can't be called (to enforce usage of the
                /// copy constructor
                /// @param[in] other input object
                /// @return reference to itself
                MPIWProjectVisGridder& operator=(const MPIWProjectVisGridder &other) = delete;

                /// @brief destructor
                virtual ~MPIWProjectVisGridder();

                /// @brief copy constructor
                /// @details It is required to decouple internal arrays in the input
                /// object and the copy.
                /// @param[in] other input object
                MPIWProjectVisGridder(const MPIWProjectVisGridder &other);

                /// Clone a copy of this Gridder
                IVisGridder::ShPtr clone() override;

                /// @brief static method to get the name of the gridder
                /// @details We specify parameters per gridder type in the parset file.
                /// This method returns the gridder name which should be used to extract
                /// a subset of parameters for createGridder method.
                static inline std::string gridderName() { return "MPIWProject";}

                /// @brief static method to create gridder
                /// @details Each gridder should have a static factory method, which is
                /// able to create a particular type of the gridder and initialise it with
                /// the parameters taken from the given parset. It is assumed that the
                /// method receives a subset of parameters where the gridder name is already
                /// taken out.
                /// @param[in] parset input parset file
                /// @return a shared pointer to the gridder instance
                static IVisGridder::ShPtr createGridder(const LOFAR::ParameterSet& parset);

            protected:
                /// @brief additional operations to configure gridder
                /// @details This method is supposed to be called from createGridder and could be
                /// used in derived classes to avoid too much duplication of the code. For this
                /// particular class it configures variable/offset support and cutoff behavior.
                /// @param[in] parset input parset file
                void configureGridder(const LOFAR::ParameterSet& parset) override;

                /// Initialize convolution function
                /// @param[in] acc const data accessor to work with
                void initConvolutionFunction(const accessors::IConstDataAccessor& acc) override;

            private:

                void setupMpiMemory(size_t bufferSize /* in bytes */);
                void copyToSharedMemory(std::vector<std::pair<int,int> >& itsConvFuncMatSize);
                void copyFromSharedMemory(const std::vector<std::pair<int,int> >& itsConvFuncMatSize);
                void copyConvFuncOffset();

#ifdef HAVE_MPI
                /// @brief MPI variables used to setup MPI shared memory
		        /// @details - MPI_Win_get_attr(itsWindowTable,MPI_WIN_SIZE,&val,&flag)
		        ///            should return a val which is the size of the shared memory segment.
                static MPI_Aint itsWindowSize;
                static int      itsWindowDisp;
                static MPI_Win  itsWindowTable;
		        /// @brief - get the group associated with the COMM_WORLD communicator
		        static MPI_Group itsWorldGroup;
		        /// @brief this communicator contains all the ranks except rank 0 of COMM_WORLD
                ///        if itsMasterDoesWork flag is false otherwise it is the same as COMM_WORLD
                static MPI_Comm itsNonRankZeroComms;
		        /// @details - get the group associated with the itsNonRankZeroComms communicator.
		        ///            This group does not contain rank 0
		        static MPI_Group itsGridderGroup;
		        /// @brief - communicator for a given node
		        static MPI_Comm itsNodeComms;
#endif

                /// @brief number of ranks within a node
                static int     itsNodeSize;
                /// @brief current rank in the itsNodeComms
                static int      itsNodeRank;
                static int      itsWorldRank;
                /// @brief a pointer to the MPI shared memory
                static imtypeComplex* itsMpiSharedMemory;

		        /// @details - These are used to synchronise and keep track of how many
		        ///            gridder objects are instantiated. The ObjCount member is
		        ///            employed to determined when the process should cleanup the
		        ///            shared memory segment. This is done when ObjCount is 0.
                static std::mutex ObjCountMutex;
                static unsigned int ObjCount;
                static bool     itsMpiMemSetup;
                /// mpi shared memory on startup
                bool itsMpiMemPreSetup;

                /// number of ranks to do the CF calculation
                int itsCFRank;

                bool itsSerial;
                bool itsMasterDoesWork;
        };
    }
}
#endif
