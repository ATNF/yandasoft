/// @file 
/// @brief iterator implementing parallel write
/// @details This is an implementation of data iterator
/// (see Base/accessors) which runs in a particular worker to
/// allow parallel writing of visibilities. Read operation is not
/// supported for simplicity. The server has to be executed at the
/// master side at the same time. It gathers the data (and distributes jobs
/// between workers). The decision was made to have this class in synthesis/parallel
/// rather than Base/accessors because it uses master-working specific code and
/// is not a general purpose class. The master (server iterator) is implemented as 
/// a static method of this class, so the communication protocol is encapsulated here.
///
/// @copyright (c) 2007 CSIRO
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>
///

#ifndef ASKAP_SYNTHESIS_PARALLEL_WRITE_ITERATOR_H
#define ASKAP_SYNTHESIS_PARALLEL_WRITE_ITERATOR_H

#include <dataaccess/IDataIterator.h>
#include <askap/parallel/ParallelAccessor.h>
#include <dataaccess/SharedIter.h>
#include <askapparallel/AskapParallel.h>

namespace askap {

namespace synthesis {

/// @brief client iterator implementing parallel write
/// @details This is an implementation of data iterator
/// (see Base/accessors) which runs in a particular worker to
/// allow parallel writing of visibilities. Read operation is not
/// supported for simplicity. The server has to be executed at the
/// master side at the same time. It gathers the data (and distributes jobs
/// between workers). The decision was made to have this class in synthesis/parallel
/// rather than Base/accessors because it uses master-working specific code and
/// is not a general purpose class.
/// @ingroup parallel
class ParallelWriteIterator : public accessors::IDataIterator {
public:
   /// @brief optional extensions to operation intended to be or'ed
   /// @details NONE is the write-only operation as envisaged in the original design,
   /// READ enables reading (and distributing) visivilities before update, otherwise 
   /// the visiibility cube is initialised with zeros for all workers, SYNCFLAGANDNOISE 
   /// enables synchronisation of flag and noise cubes back to the master. Finally,
   /// ALL enables all such extensions
   /// @note The mode of operations is controlled by the master and is distributed to
   /// all workers via MPI.
   enum OpExtension {
       NONE = 0,
       READ = 1,
       SYNCFLAG = 2,
       SYNCNOISE = 4,
       SYNCFLAGANDNOISE = 6,
       ALL = 7
   };

   /// @brief constructor
   /// @details 
   /// @param comms communication object
   /// @param[in] cacheSize uvw-machine cache size
   /// @param[in] tolerance pointing direction tolerance in radians, exceeding
   /// which leads to initialisation of a new UVW machine and recompute of the rotated uvws/delays  
   explicit ParallelWriteIterator(askap::askapparallel::AskapParallel& comms, size_t cacheSize = 1, double tolerance = 1e-6);
    
    
   // Return the data accessor (current chunk) in various ways	

   /// @brief reference to data accessor (current chunk)
   /// @return a reference to the current chunk
   /// @note constness of the return type is changed to allow read/write
   /// operations.
   virtual accessors::IDataAccessor& operator*() const;
	
   /// Switch the output of operator* and operator-> to one of 
   /// the buffers. This is meant to be done to provide the same 
   /// interface for a buffer access as exists for the original 
   /// visibilities (e.g. it->visibility() to get the cube).
   /// It can be used for an easy substitution of the original 
   /// visibilities to ones stored in a buffer, when the iterator is
   /// passed as a parameter to mathematical algorithms. 
   /// 
   /// The operator* and operator-> will refer to the chosen buffer
   /// until a new buffer is selected or the chooseOriginal() method
   /// is executed to revert operators to their default meaning
   /// (to refer to the primary visibility data).
   ///
   /// @param[in] bufferID  the name of the buffer to choose
   ///
   virtual void chooseBuffer(const std::string &bufferID);

   /// Switch the output of operator* and operator-> to the original
   /// state (present after the iterator is just constructed) 
   /// where they point to the primary visibility data. This method
   /// is indended to cancel the results of chooseBuffer(casacore::uInt)
   ///
   virtual void chooseOriginal();

   /// return any associated buffer for read/write access. The 
   /// buffer is identified by its bufferID. The method 
   /// ignores a chooseBuffer/chooseOriginal setting.
   /// 
   /// @param[in] bufferID the name of the buffer requested
   /// @return a reference to writable data accessor to the
   ///         buffer requested
   ///
   /// Because IDataAccessor has both const and non-const visibility()
   /// methods defined separately, it is possible to detect when a
   /// write operation took place and implement a delayed writing
   virtual accessors::IDataAccessor& buffer(const std::string &bufferID) const;

   /// Restart the iteration from the beginning
   virtual void init();

   /// Checks whether there are more data available.
   /// @return True if there are more data available
   casacore::Bool hasMore() const throw(); 

   /// advance the iterator one step further
   /// @return True if there are more data (so constructions like
   ///         while(it.next()) {} are possible)
   casacore::Bool next();

   /// @brief server method
   /// @details It iterates through the given iterator, serves metadata
   /// to client iterators and combines visibilities in a single cube.
   /// @param comms communication object
   /// @param iter shared iterator to use
   /// @param mode operating mode, a flag or flags allowing extensions to be enabled
   /// Supported modes are:  NONE is the write-only operation as envisaged in the original design,
   /// READ enables reading (and distributing) visivilities before update, otherwise 
   /// the visiibility cube is initialised with zeros for all workers, SYNCFLAGANDNOISE 
   /// enables synchronisation of flag and noise cubes back to the master. Finally,
   /// ALL enables all such extensions
   /// @note The mode of operations is controlled by the master (i.e. this method) and is distributed to
   /// all workers via MPI.
   static void masterIteration(askap::askapparallel::AskapParallel& comms, const accessors::IDataSharedIter &iter, OpExtension mode = NONE);


   /// @brief obtain the channel offset
   /// @details The method is specific to this particular adapter. It is sometimes handly to get the channel offset of
   /// the chunk dealt with by this worker to avoid re-implementing this functionality in the user code as it would require
   /// additional communication.
   /// @return channel offset of the first channel dealt with by this worker w.r.t. the first channel of the whole spectrum
   casacore::uInt inline chanOffset() const { return itsChanOffset; }
	
protected:
    
   /// @brief obtain metadata for the next iteration
   /// @details This is a core method of the class. It receives the
   /// status message from the master and reads the metadata if not at
   /// the last iteration. If not at the first iteration, it also syncronises
   /// the visibility cube with the master before advancing to the next iteration.
   void advance();

private:
   /// @brief communicator
   askap::askapparallel::AskapParallel& itsComms;
	
   /// @brief true, if at least one chunk of data have been received
   /// @details We use this flag to assert that no attempts to restart
   /// iteration have been performed.
   bool itsNotAtOrigin;

   /// @brief buffer for the accessor
   mutable ParallelAccessor itsAccessor;
	
   /// @brief true if current accessor contains valid data
   bool itsAccessorValid;

   /// @brief current channel offset w.r.t the full axis
   /// @details Workers may need to know which part of the full spectrum they are working with. Although ideally, we should've used
   /// the accessor itself (e.g. the frequency axis) to figure this out.
   casacore::uInt itsChanOffset;

   /// @brief true if flags need to be synchronised from the current accessor
   /// @note This field is only checked/makes sense if accessor is valid and iterator is not at origin
   bool itsFlagNeedSync;
};

/// @brief or operator for operational mode flags
/// @param[in] a first operand
/// @param[in] b second operand
/// @return or'ed result
inline ParallelWriteIterator::OpExtension operator|(const ParallelWriteIterator::OpExtension a, const ParallelWriteIterator::OpExtension b) 
     {return static_cast<ParallelWriteIterator::OpExtension>(static_cast<int>(a) | static_cast<int>(b));}

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_PARALLEL_WRITE_ITERATOR_H


