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

#include <askap/parallel/ParallelWriteIterator.h>
#include <askap/parallel/ParallelIteratorStatus.h>
#include <askap/dataaccess/IFlagDataAccessor.h>
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/RangePartition.h>

#include <Blob/BlobString.h>
#include <Blob/BlobIBufString.h>
#include <Blob/BlobOBufString.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobArray.h>


ASKAP_LOGGER(logger, ".parallel");

namespace askap {

namespace synthesis {

/// @brief constructor
/// @details 
/// @param comms communication object
/// @param[in] cacheSize uvw-machine cache size
/// @param[in] tolerance pointing direction tolerance in radians, exceeding
/// which leads to initialisation of a new UVW machine and recompute of the rotated uvws/delays  
ParallelWriteIterator::ParallelWriteIterator(askap::askapparallel::AskapParallel& comms, size_t cacheSize, double tolerance) : 
   itsComms(comms), itsNotAtOrigin(false), itsAccessor(cacheSize, tolerance), itsAccessorValid(false), itsChanOffset(0u),
   itsFlagNeedSync(false)
{
  ASKAPCHECK(itsComms.isWorker() && itsComms.isParallel(), 
      "ParallelWriteIterator class is supposed to be used only in workers in the parallel mode");
  advance();
}    
    
// Return the data accessor (current chunk) in various ways	

/// @brief reference to data accessor (current chunk)
/// @return a reference to the current chunk
/// @note constness of the return type is changed to allow read/write
/// operations.
accessors::IDataAccessor& ParallelWriteIterator::operator*() const
{
  ASKAPCHECK(itsAccessorValid, "An attempt to obtain accessor following the end of iteration");
  return itsAccessor;
}
		
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
void ParallelWriteIterator::chooseBuffer(const std::string &bufferID)
{
  ASKAPTHROW(AskapError, "An attempt to choose the buffer "<<bufferID<<
             ". Operation is not supported by the parallel iterator");
}

/// Switch the output of operator* and operator-> to the original
/// state (present after the iterator is just constructed) 
/// where they point to the primary visibility data. This method
/// is indended to cancel the results of chooseBuffer(casacore::uInt)
///
void ParallelWriteIterator::chooseOriginal() {}

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
accessors::IDataAccessor& ParallelWriteIterator::buffer(const std::string &bufferID) const
{
  ASKAPTHROW(AskapError, "An attempt to access the buffer "<<bufferID<<
             ". Operation is not supported by the parallel iterator");  
}
	
/// Restart the iteration from the beginning
void ParallelWriteIterator::init()
{
  ASKAPCHECK(!itsNotAtOrigin, "Restart of the iteration is not supported by the parallel iterator");
}
	
/// Checks whether there are more data available.
/// @return True if there are more data available
casacore::Bool ParallelWriteIterator::hasMore() const throw()
{
  return itsAccessorValid;
}
	
/// advance the iterator one step further
/// @return True if there are more data (so constructions like
///         while(it.next()) {} are possible)
casacore::Bool ParallelWriteIterator::next()
{
  itsNotAtOrigin = true;
  advance();
  return hasMore();
}

/// @brief obtain metadata for the next iteration
/// @details This is a core method of the class. It receives the
/// status message from the master and reads the metadata if not at
/// the last iteration. If not at the first iteration, it also syncronises
/// the visibility cube with the master before advancing to the next iteration.
void ParallelWriteIterator::advance()
{
  ASKAPDEBUGASSERT(itsComms.isWorker());
  if (itsNotAtOrigin) {
      // sync the result
      //ASKAPLOG_INFO_STR(logger, "About to send visibilities from rank "<<itsComms.rank());
      ASKAPDEBUGASSERT(itsAccessor.itsVisibility.shape() == itsAccessor.itsFlag.shape()); 
      LOFAR::BlobString bs;
      bs.resize(0);
      LOFAR::BlobOBufString bob(bs);
      LOFAR::BlobOStream out(bob);
      out.putStart("AccessorVisibilities",itsFlagNeedSync ? 2 : 1);
      if (itsFlagNeedSync) {
          out<<itsAccessor.itsFlag;
      }
      out<<itsAccessor.itsVisibility;
      out.putEnd();
      itsComms.sendBlob(bs, 0);      
  }
  // get status
  // update itsAccessorValid from status
  ParallelIteratorStatus status;
  {
    LOFAR::BlobString bs;
    bs.resize(0);
    itsComms.broadcastBlob(bs,0);
    LOFAR::BlobIBufString bib(bs);
    LOFAR::BlobIStream in(bib);
    in>>status;
    itsAccessorValid = status.itsHasMore;
    //ASKAPLOG_INFO_STR(logger, "Received status "<<itsAccessorValid<<" (rank "<<itsComms.rank()<<")");
  }
  const bool readVis = (status.itsMode & READ) == READ;
  itsFlagNeedSync = (status.itsMode & SYNCFLAG) == SYNCFLAG;

  // channel distribution is static, but the number of channels in the measurement set can change.
  // although the latter case is not supported, we have to cater for this possbility at some minimal level 
  // to avoid nasty bugs.
  // This functionality is only needed for high-level user, the code below is agnostic on which data it works with
  ASKAPDEBUGASSERT(itsComms.nProcs() > 1);
  ASKAPDEBUGASSERT(itsComms.rank() > 0);
  itsChanOffset = utility::RangePartition(status.itsTotalNChan, static_cast<unsigned int>(itsComms.nProcs()) - 1u).first(itsComms.rank() - 1);
        
  if (itsAccessorValid) {
      // receive common metadata
      {
        LOFAR::BlobString bs;
        bs.resize(0);
        itsComms.broadcastBlob(bs,0);
        LOFAR::BlobIBufString bib(bs);
        LOFAR::BlobIStream in(bib);
        casacore::Matrix<casacore::Double> uvwBuf;
        casacore::Matrix<casacore::Vector<casacore::Double> > dirBuf;
        casacore::Vector<casacore::Int> stokesBuf;
        const int version = in.getStart("AccessorMetadata");
        ASKAPCHECK(version == 1, "Version mismatch for AccessorMetadata stream, you have version="<<version);
        in >> itsAccessor.itsAntenna1 >> itsAccessor.itsAntenna2 >> itsAccessor.itsFeed1 >> itsAccessor.itsFeed2 >> 
              itsAccessor.itsFeed1PA >> itsAccessor.itsFeed2PA >> dirBuf >> uvwBuf >> itsAccessor.itsTime >> stokesBuf;        
        in.getEnd();        
        ASKAPASSERT(dirBuf.nrow() == status.itsNRow);
        ASKAPASSERT(uvwBuf.nrow() == status.itsNRow);
        ASKAPASSERT(dirBuf.ncolumn() == 4);
        ASKAPASSERT(uvwBuf.ncolumn() == 3);
        itsAccessor.itsUVW.resize(uvwBuf.nrow());
        itsAccessor.itsPointingDir1.resize(dirBuf.nrow());
        itsAccessor.itsPointingDir2.resize(dirBuf.nrow());
        itsAccessor.itsDishPointing1.resize(dirBuf.nrow());
        itsAccessor.itsDishPointing2.resize(dirBuf.nrow());
        for (casacore::uInt row = 0; row<dirBuf.nrow(); ++row) {
             itsAccessor.itsPointingDir1[row] = casacore::MVDirection(dirBuf(row,0));                  
             itsAccessor.itsPointingDir2[row] = casacore::MVDirection(dirBuf(row,1));
             itsAccessor.itsDishPointing1[row] = casacore::MVDirection(dirBuf(row,2));                  
             itsAccessor.itsDishPointing2[row] = casacore::MVDirection(dirBuf(row,3));
             for (casacore::uInt col = 0; col<3; ++col) {
                  itsAccessor.itsUVW[row](col) = uvwBuf(row,col);
             }
        }
        itsAccessor.itsStokes.resize(stokesBuf.nelements());
        ASKAPASSERT(stokesBuf.nelements() == status.itsNPol);
        for (casacore::uInt pol = 0; pol<status.itsNPol; ++pol) {
             itsAccessor.itsStokes[pol] = casacore::Stokes::StokesTypes(stokesBuf[pol]);
        }
        // consistency check
        ASKAPASSERT(itsAccessor.itsAntenna1.nelements() == status.itsNRow);
        ASKAPASSERT(itsAccessor.itsAntenna2.nelements() == status.itsNRow);
        ASKAPASSERT(itsAccessor.itsFeed1.nelements() == status.itsNRow);
        ASKAPASSERT(itsAccessor.itsFeed2.nelements() == status.itsNRow);
        ASKAPASSERT(itsAccessor.itsFeed1PA.nelements() == status.itsNRow);
        ASKAPASSERT(itsAccessor.itsFeed2PA.nelements() == status.itsNRow);        
      }
      // receive unique metadata, fill itsAccessor
      {
        //ASKAPLOG_INFO_STR(logger, "About to receive rank-specific metadata in rank "<<itsComms.rank());        
        LOFAR::BlobString bs;
        bs.resize(0);
        itsComms.receiveBlob(bs, 0);
        LOFAR::BlobIBufString bib(bs);
        LOFAR::BlobIStream in(bib);
        const int version = in.getStart("AccessorVariableMetadata");
        ASKAPCHECK(version == (readVis ? 2 : 1), "Version mismatch receiving rank-specific metadata");
        in>>itsAccessor.itsFlag>>itsAccessor.itsNoise>>itsAccessor.itsFrequency;
        if (readVis) {
            in>>itsAccessor.itsVisibility;
        }
        in.getEnd();
        if (readVis) {
            ASKAPASSERT(itsAccessor.itsVisibility.nrow() == itsAccessor.itsFlag.nrow());
            ASKAPASSERT(itsAccessor.itsVisibility.ncolumn() == itsAccessor.itsFlag.ncolumn());
            ASKAPASSERT(itsAccessor.itsVisibility.nplane() == itsAccessor.itsFlag.nplane());
        } else {
            itsAccessor.itsVisibility.resize(itsAccessor.itsFlag.nrow(), itsAccessor.itsFlag.ncolumn(), itsAccessor.itsFlag.nplane());
            itsAccessor.itsVisibility.set(0.);    
        }
        // consistency checks
        ASKAPASSERT(itsAccessor.nRow() == itsAccessor.itsVisibility.nplane());
        ASKAPASSERT(itsAccessor.nChannel() == itsAccessor.itsVisibility.ncolumn());
        ASKAPASSERT(itsAccessor.nPol() == itsAccessor.itsVisibility.nrow());
        ASKAPASSERT(itsAccessor.nChannel() == itsAccessor.itsFrequency.nelements());            
      }
  }
}


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
void ParallelWriteIterator::masterIteration(askap::askapparallel::AskapParallel& comms, const accessors::IDataSharedIter &iter, ParallelWriteIterator::OpExtension mode)
{
  ASKAPCHECK(mode < SYNCNOISE, "Requested extension mode is not supported at the moment");
  ASKAPDEBUGASSERT(comms.isMaster());
  ASKAPDEBUGASSERT(comms.nProcs() > 1);
  accessors::IDataSharedIter it(iter);
  bool contFlag = true;
  ParallelIteratorStatus status;
  do {
    status.itsHasMore = it.hasMore();
    if (status.itsHasMore) {
       ASKAPCHECK(it->nChannel() >= static_cast<casacore::uInt>(comms.nProcs() - 1), 
                  "Idle workers are not currently supported. Number of spectral channels("<<it->nChannel()<<
                  ") should not be less than the number of workers ("<<(comms.nProcs() - 1)<<")");
       
       status.itsNRow = it->nRow();
       status.itsNPol = it->nPol();
       status.itsMode = mode;
       if (status.itsTotalNChan > 0) {
           ASKAPCHECK(status.itsTotalNChan == it->nChannel(), "The current code does not support change to the spectral axis throughout the dataset");
       }
       status.itsTotalNChan = it->nChannel();
    } else {
      contFlag = false;
    }
    {
      LOFAR::BlobString bs;
      bs.resize(0);
      LOFAR::BlobOBufString bob(bs);
      LOFAR::BlobOStream out(bob);
      out << status;
      //ASKAPLOG_INFO_STR(logger, "About to send status "<<status.itsHasMore<<" (rank "<<comms.rank()<<")");      
      comms.broadcastBlob(bs,0);
    }
    
    if (contFlag) {
        // broadcast common metadata
        {
          LOFAR::BlobString bs;
          bs.resize(0);
          LOFAR::BlobOBufString bob(bs);
          LOFAR::BlobOStream out(bob);
          out.putStart("AccessorMetadata", 1);
          casacore::Matrix<casacore::Double> uvwBuf(it->nRow(),3);
          casacore::Matrix<casacore::Vector<casacore::Double> > dirBuf(it->nRow(),4);
          for (casacore::uInt row = 0; row<it->nRow(); ++row) {
               for (casacore::uInt col = 0; col<3; ++col) {
                    uvwBuf(row,col) = it->uvw()[row](col);
               }
               dirBuf(row,0) = it->pointingDir1()[row].get();
               dirBuf(row,1) = it->pointingDir2()[row].get();
               dirBuf(row,2) = it->dishPointing1()[row].get();
               dirBuf(row,3) = it->dishPointing2()[row].get();               
          }
          casacore::Vector<casacore::Int> stokesBuf(it->stokes().nelements());
          for (casacore::uInt pol = 0; pol<stokesBuf.nelements(); ++pol) {
               stokesBuf[pol] = casacore::Int(it->stokes()[pol]);
          }
          out << it->antenna1() << it->antenna2() << it->feed1() << it->feed2() << it->feed1PA() <<
                 it->feed2PA() << dirBuf << uvwBuf << it->time() << stokesBuf;
          out.putEnd();   
          comms.broadcastBlob(bs,0);
        }
        // point-to-point transfer of data which differ
        for (int worker = 0; worker < comms.nProcs() - 1; ++worker) {
             //ASKAPLOG_INFO_STR(logger, "About to send rank-specific metadata to rank "<<worker + 1);
             // start and stop of the slice
             casacore::IPosition start(3,0);
             ASKAPDEBUGASSERT((it->nRow()!=0) && (it->nChannel()!=0) && (it->nPol()));

             utility::RangePartition rp(it->nChannel(), static_cast<unsigned int>(comms.nProcs()) - 1u);

             //casacore::IPosition end(3,int(it->nRow()) - 1, static_cast<int>(rp.last(worker)), int(it->nPol()) - 1);
             casacore::IPosition end(3,int(it->nPol()) - 1, static_cast<int>(rp.last(worker)), int(it->nRow()) - 1);
             start(1) = static_cast<int>(rp.first(worker));
             ASKAPASSERT(end(1) < int(it->nChannel()));
             ASKAPDEBUGASSERT(start(1)<=end(1));
             const casacore::IPosition vecStart(1, start(1));
             const casacore::IPosition vecEnd(1, end(1));
             // send slices of flags, noise and frequency. Visibility can be assumed to be zero or read as well.
             {
               const bool readVis = ((mode & READ) == READ);
               LOFAR::BlobString bs;
               bs.resize(0);
               LOFAR::BlobOBufString bob(bs);
               LOFAR::BlobOStream out(bob);
               out.putStart("AccessorVariableMetadata", readVis ? 2 : 1);
               casacore::Cube<casacore::Bool> flagBuf(it->flag());
               casacore::Cube<casacore::Complex> noiseBuf(it->noise());
               casacore::Vector<casacore::Double> freqBuf(it->frequency());
               out<<flagBuf(start,end)<<noiseBuf(start,end)<<freqBuf(vecStart,vecEnd);
               if (readVis) {
                   casacore::Cube<casacore::Complex> visBuf(it->visibility());
                   out<<visBuf(start,end);
               }
               out.putEnd();
               comms.sendBlob(bs, worker + 1);
             }
        }
        
        // receive the result and store it in rwVisibility
        for (int worker = 0; worker < comms.nProcs() - 1; ++worker) {
             //ASKAPLOG_INFO_STR(logger, "About to receive visibilities from rank "<<worker + 1);
             // start and stop of the slice
             casacore::IPosition start(3,0);
             ASKAPDEBUGASSERT((it->nRow()!=0) && (it->nChannel()!=0) && (it->nPol()));

             utility::RangePartition rp(it->nChannel(), static_cast<unsigned int>(comms.nProcs()) - 1u);

             //casacore::IPosition end(3,int(it->nRow()) - 1, static_cast<int>(rp.last(worker)), int(it->nPol()) - 1);
             casacore::IPosition end(3,int(it->nPol()) - 1, static_cast<int>(rp.last(worker)), int(it->nRow()) - 1);
             start(1) = static_cast<int>(rp.first(worker));
             ASKAPASSERT(end(1) < int(it->nChannel()));

             ASKAPDEBUGASSERT(start(1)<=end(1));
             // receive a slice of visibility
             {
               const bool recvFlags = (status.itsMode & SYNCFLAG) == SYNCFLAG;
               LOFAR::BlobString bs;
               bs.resize(0);
               comms.receiveBlob(bs, worker + 1);               
               LOFAR::BlobIBufString bib(bs);
               LOFAR::BlobIStream in(bib);
               const int version = in.getStart("AccessorVisibilities");
               ASKAPCHECK(version == recvFlags ? 2 : 1, "Version mismatch in serialising of visibilities");
               if (recvFlags) {
                   casacore::Cube<casacore::Bool> flagBuf;
                   in>>flagBuf;
                   const boost::shared_ptr<accessors::IFlagDataAccessor> fda = boost::dynamic_pointer_cast<accessors::IFlagDataAccessor>(boost::shared_ptr<accessors::IDataAccessor>(it.operator->(), utility::NullDeleter()));
                   ASKAPCHECK(fda, "Flag write operation is requested, but supplied data accessor doesn't support writing flags");
                   casacore::Cube<casacore::Bool> flagSlice = fda->rwFlag()(start,end);
                   ASKAPCHECK(flagSlice.shape() == flagBuf.shape(), "Shape mismatch of the flag cube, received has shape="<<
                        flagBuf.shape()<<" expected shape="<<flagSlice.shape());
                   flagSlice = flagBuf;
               }
               casacore::Cube<casacore::Complex> visBuf;
               in>>visBuf;
               in.getEnd();
               casacore::Cube<casacore::Complex> visSlice = it->rwVisibility()(start,end);
               ASKAPCHECK(visSlice.shape() == visBuf.shape(), "Shape mismatch of the visibility cube, received has shape="<<
                        visBuf.shape()<<" expected shape="<<visSlice.shape());
               visSlice = visBuf;
             }
        }
        
        it.next();
    }
  } while (contFlag);  
}


} // namespace synthesis

} // namespace askap
