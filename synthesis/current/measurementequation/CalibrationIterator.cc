/// @file
/// 
/// @brief An iterator adapter which applies some calibration
/// @details The measurement equations used for imaging and calibration are
/// written in a different way: calibration equations work with individual
/// accessors, while imaging equations work with iterator as a whole. I hope
/// this will change in the future, but at this stage adapters are required.
/// One of such adapters is FakeSingleStepIterator (in the dataaccess subpackage).
/// It is an iterator over a single accessor and can be used to plug an imaging
/// equation into the calibration framework (as the source of uncorrupted 
/// visibilities). This class is doing a reverse. It allows to use calibration
/// equation as a source of visibilities for imaging equation. As a result,
/// the measurement equation formed this way would deal with calibrated data.
/// @note I hope this adapter is temporary, and a better way of handling 
/// composite equations will be adopted in the future.
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

#include <measurementequation/CalibrationIterator.h>
#include <dataaccess/OnDemandNoiseAndFlagDA.h>
#include <askap/AskapError.h>

#include <measurementequation/ICalibrationApplicator.h>


using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

/// @brief construct iterator
/// @details The input iterator is remembered and switched to the 
/// original visibilities (can be switched to a buffer later, but 
/// via this adapter interface). Note, direct manipulation with this 
/// input iterator after it is assigned to CalibrationIterator can
/// lead to an unpredictable result.
/// @param[in] iter input iterator
/// @param[in] calME calibration measurement equation
/// @note If calME belongs to a type initialized with an iterator, it doesn't
/// matter which iterator it has been initialized with. This class always
/// uses accessor-based methods.
CalibrationIterator::CalibrationIterator(const IDataSharedIter &iter, 
              const boost::shared_ptr<ICalibrationApplicator> &calME) :
    itsWrappedIterator(iter), itsCalibrationME(calME), itsBufferFlag(false)
{
  ASKAPASSERT(itsWrappedIterator);
  itsWrappedIterator.chooseOriginal();
  ASKAPASSERT(itsCalibrationME);
}
	
/// Return the data accessor (current chunk) in various ways
/// operator* delivers a reference to data accessor (current chunk)
///
/// @return a reference to the current chunk
///
/// constness of the return type is changed to allow read/write
/// operations.
///
IDataAccessor& CalibrationIterator::operator*() const
{
  ASKAPDEBUGASSERT(itsWrappedIterator);
  if (itsBufferFlag) {
      // buffer is active, just return what wrapped iterator points to
      return *itsWrappedIterator;
  } 
  if (!itsDataAccessor) {
      // need to fill the local buffer with the calibrated visibilities
      
      // note, accessor is remembered by reference in OnDemandNoiseAndFlagDA!
      // But it's OK because we invalidate the cache each time before the iterator
      // changes.
      itsDataAccessor.reset(new OnDemandNoiseAndFlagDA(*itsWrappedIterator));

      ASKAPDEBUGASSERT(itsDataAccessor);
      // copy perfect data
      itsDataAccessor->rwVisibility() = itsWrappedIterator->visibility().copy();
      
      // correct data
      ASKAPDEBUGASSERT(itsCalibrationME);              
      itsCalibrationME->correct(*itsDataAccessor);
  }
  return *itsDataAccessor;
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
void CalibrationIterator::chooseBuffer(const std::string &bufferID)
{
  itsBufferFlag = true; // a buffer is active
  itsDataAccessor.reset(); // invalidate cache of calibrated visibilities
  ASKAPDEBUGASSERT(itsWrappedIterator)
  itsWrappedIterator.chooseBuffer(bufferID);
}

/// Switch the output of operator* and operator-> to the original
/// state (present after the iterator is just constructed) 
/// where they point to the primary visibility data. This method
/// is indended to cancel the results of chooseBuffer(casa::uInt)
///
void CalibrationIterator::chooseOriginal()
{
  itsBufferFlag = false; // original is active
  itsDataAccessor.reset(); // invalidate cache of calibrated visibilities
  ASKAPDEBUGASSERT(itsWrappedIterator)
  itsWrappedIterator.chooseOriginal();
}

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
IDataAccessor& CalibrationIterator::buffer(const std::string &bufferID) const
{
  ASKAPDEBUGASSERT(itsWrappedIterator)
  return itsWrappedIterator.buffer(bufferID);
}

/// Restart the iteration from the beginning
void CalibrationIterator::init()
{
  ASKAPDEBUGASSERT(itsWrappedIterator)
  itsDataAccessor.reset(); // invalidate cache of calibrated visibilities
  itsWrappedIterator.init();
}
	
/// Checks whether there are more data available.
/// @return True if there are more data available
casa::Bool CalibrationIterator::hasMore() const throw()
{
  #ifdef ASKAP_DEBUG
  if (!itsWrappedIterator) {
      return false;
  }
  #endif 
  return itsWrappedIterator.hasMore();
}
	
/// advance the iterator one step further 
/// @return True if there are more data (so constructions like 
///         while(it.next()) {} are possible)
casa::Bool CalibrationIterator::next()
{
  itsDataAccessor.reset(); // invalidate cache of calibrated visibilities
  ASKAPDEBUGASSERT(itsWrappedIterator)
  return itsWrappedIterator.next();
}


