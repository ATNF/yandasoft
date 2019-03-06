/// @file
/// 
/// @brief This class is used inside the measurement equation object implementing
/// pre-averaging (or pre-summing to be exact) algorithm for calibration. Strictly
/// speaking it is not an adapter and it doesn't behave as an accessor. However, it
/// mimics the accessor interface, so we can reuse the existing code to a greater
/// extent. In addition, we can extend the code to a more complicated types of calibration 
/// later (i.e. with equations using more metadata). Current implementation is derived 
/// from DataAccessorAdapter just to speed up the development. None of the functionality
/// of this base class is used (except throwing exceptions if methods which are not
/// intended to be used are called). The plan is to always keep the DataAccessorAdapter
/// in the detached state.
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

#ifndef SYNTHESIS_PREAVGCALBUFFER_H
#define SYNTHESIS_PREAVGCALBUFFER_H

#include <dataaccess/DataAccessorAdapter.h>
#include <dataaccess/IConstDataAccessor.h>
#include <measurementequation/IMeasurementEquation.h>
#include <fitting/PolXProducts.h>

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief This class is used inside the measurement equation object implementing
/// pre-averaging (or pre-summing to be exact) algorithm for calibration. Strictly
/// speaking it is not an adapter and it doesn't behave as an accessor. However, it
/// mimics the accessor interface, so we can reuse the existing code to a greater
/// extent. In addition, we can extend the code to a more complicated types of calibration 
/// later (i.e. with equations using more metadata). Current implementation is derived 
/// from DataAccessorAdapter just to speed up the development. None of the functionality
/// of this base class is used (except throwing exceptions if methods which are not
/// intended to be used are called). The plan is to always keep the DataAccessorAdapter
/// in the detached state.
/// @note At the moment all frequency channels are summed up together. Later we may want
/// to implement a partial averaging in frequency.
/// @ingroup measurementequation
class PreAvgCalBuffer : public accessors::DataAccessorAdapter {
public:
   /// @brief default constructor
   /// @details preaveraging is initialised based on the first encountered accessor
   PreAvgCalBuffer();
   
   /// @brief constructor with explicit averaging parameters
   /// @details This version of the constructor explicitly defines the number of 
   /// antennas and beams to initialise the buffer appropriately.
   /// @param[in] nAnt number of antennas, indices are expected to run from 0 to nAnt-1
   /// @param[in] nBeam number of beams, indices are expected to run from 0 to nBeam-1
   /// @param[in] nChan number of channels to buffer, 1 (default) is a special case
   /// assuming that measurement equation is frequency-independent
   PreAvgCalBuffer(casa::uInt nAnt, casa::uInt nBeam, casa::uInt nChan = 1);
   
   /// @brief constructor with explicit averaging parameters
   /// @details This version of the constructor explicitly defines the number of 
   /// antennas to initialise the buffer appropriately. Unlike the version with 
   /// explicitly given number of beams with nBeam set to 1, this constructor configures
   /// the buffer to ignore the beam index (i.e. assuming the measurement equation is beam-independent)
   /// @param[in] nAnt number of antennas, indices are expected to run from 0 to nAnt-1
   /// @note frequency-independent case is implied
   PreAvgCalBuffer(casa::uInt nAnt); 
   
   /// @brief configure beam-independent accumulation
   /// @param[in] flag if true, accumulation is beam-independent
   void beamIndependent(bool flag);
   
   /// @brief initialise accumulation via an accessor
   /// @details This method resets the buffers and sets the shape using the given accessor
   /// as a template.
   /// @param[in] acc template accessor
   /// @param[in] fdp frequency dependency flag, if true a separate buffer is created for every
   /// presented spectral channel. Otherwise (default), all channels contribute to the same buffer
   /// (which is appropriate for frequency-independent effects)
   void initialise(const IConstDataAccessor &acc, const bool fdb = false);
   
   /// @brief initialise accumulation explicitly
   /// @details This method resets the buffers and sets the shape to accommodate the given
   /// number of antennas and beams (i.e. the buffer size is nBeams*nAnt*(nAnt-1)/2)
   /// @param[in] nAnt number of antennas, indices are expected to run from 0 to nAnt-1
   /// @param[in] nBeam number of beams, indices are expected to run from 0 to nBeam-1
   /// @param[in] nChan number of channels to buffer, 1 (default) is a special case
   /// assuming that measurement equation is frequency-independent
   void initialise(casa::uInt nAnt, casa::uInt nBeam, casa::uInt nChan = 1);
   
   // implemented accessor methods
   
   /// The number of rows in this chunk
   /// @return the number of rows in this chunk
   virtual casa::uInt nRow() const throw();
  	
   /// The number of spectral channels (equal for all rows)
   /// @return the number of spectral channels
   virtual casa::uInt nChannel() const throw();

   /// The number of polarization products (equal for all rows)
   /// @return the number of polarization products (can be 1,2 or 4)
   virtual casa::uInt nPol() const throw();

   /// First antenna IDs for all rows
   /// @return a vector with IDs of the first antenna corresponding
   /// to each visibility (one for each row)
   virtual const casa::Vector<casa::uInt>& antenna1() const;

   /// Second antenna IDs for all rows
   /// @return a vector with IDs of the second antenna corresponding
   /// to each visibility (one for each row)
   virtual const casa::Vector<casa::uInt>& antenna2() const;
  
   /// First feed IDs for all rows
   /// @return a vector with IDs of the first feed corresponding
   /// to each visibility (one for each row)
   virtual const casa::Vector<casa::uInt>& feed1() const;

   /// Second feed IDs for all rows
   /// @return a vector with IDs of the second feed corresponding
   /// to each visibility (one for each row)
   virtual const casa::Vector<casa::uInt>& feed2() const;

   /// Cube of flags corresponding to the output of visibility() 
   /// @return a reference to nRow x nChannel x nPol cube with flag 
   ///         information. If True, the corresponding element is flagged bad.
   virtual const casa::Cube<casa::Bool>& flag() const;

   /// @brief polarisation type for each product
   /// @return a reference to vector containing polarisation types for
   /// each product in the visibility cube (nPol() elements).
   /// @note All rows of the accessor have the same structure of the visibility
   /// cube, i.e. polarisation types returned by this method are valid for all rows.
   virtual const casa::Vector<casa::Stokes::StokesTypes>& stokes() const;

   // access to accumulated statistics
   
   /// @brief obtain weighted sum of polarisation products
   /// @details We accumulate cross-products of polarisation
   /// components of model and measured visibilities 
   /// (and model visibilities alone) using the helper class of
   /// type PolXProducts. This method allows to get a reference
   /// to the buffer class (to be used to construct normal equations)
   /// managed by this instance of PreAvgCalBuffer.
   /// @return const reference to the helper class
   inline const scimath::PolXProducts& polXProducts() const { return itsPolXProducts; }
   
   
   // the actual summing in of an accessor with data
   
   /// @brief process one accessor
   /// @details This method processes the given accessor and updates the internal 
   /// buffers. The measurement equation is used to calculate model visibilities 
   /// corresponding to measured visibilities.
   /// @param[in] acc input accessor with measured data
   /// @param[in] me shared pointer to the measurement equation
   /// @param[in] fdp frequency dependency flag (see initialise). It is used if initialisation from accessor
   /// is required. Otherwise, it is just checked for consistency (i.e. more than one channel is defined, if it is true)
   /// @note only predict method of the measurement equation is used.
   void accumulate(const IConstDataAccessor &acc, const boost::shared_ptr<IMeasurementEquation const> &me, const bool fdp = false);

   // access stats
   
   /// @brief number of visibilities ignored due to type
   /// @details This includes auto-correlations and cross-beam cross-correlations.
   /// @return number of visibilities ignored due to visibility type
   inline casa::uInt ignoredDueToType() const { return itsVisTypeIgnored;}
   
   /// @brief number of visibilities ignored due to lack of match
   /// @details This covers the visibilities which fell outside the range of
   /// indices covered by this buffer
   /// @return number of ignored visibilities due to lack of match
   inline casa::uInt ignoredNoMatch() const { return itsNoMatchIgnored;}
   
   /// @brief number of visibilities ignored due to flags
   /// @details This includes visibilities flagged in the input accessor
   /// and visibilities which correspond to polarisation products not managed 
   /// by the buffer.
   /// @return number of visibilities ignored due to flags
   inline casa::uInt ignoredDueToFlags() const { return itsFlagIgnored;}
    
protected:
   /// @brief helper method to find a match row in the buffer
   /// @details It goes over antenna and beam indices and finds a buffer row which 
   /// corresponds to the given indices.
   /// @param[in] ant1 index of the first antenna
   /// @param[in] ant2 index of the second antenna
   /// @param[in] beam beam index
   /// @return row number in the buffer corresponding to the given (ant1,ant2,beam) or -1 if 
   /// there is no match
   int findMatch(casa::uInt ant1, casa::uInt ant2, casa::uInt beam); 
      
private:
   /// @brief indices of the first antenna for all rows
   casa::Vector<casa::uInt> itsAntenna1;   
   
   /// @brief indices of the second antenna for all rows
   casa::Vector<casa::uInt> itsAntenna2;   
   
   /// @brief indices of the beam for all rows
   /// @note beam cross-products are not supported here
   casa::Vector<casa::uInt> itsBeam;
   
   /// @brief flags for all rows, channels and polarisations
   casa::Cube<casa::Bool> itsFlag;
   
   /// @brief types of polarisation products
   casa::Vector<casa::Stokes::StokesTypes> itsStokes;
   
   /// @brief buffer for accumulated cross-products
   /// @details This helper object handles two buffers for
   /// cross-products of conj(Vmodel)*Vmodel and conj(Vmodel)*Vmeasured.
   scimath::PolXProducts itsPolXProducts;
   
   // statistics on the number of samples excluded due to various criteria
   
   /// @brief number of ignored samples due to visibility type
   /// @details (autocorrelation or cross-beam visibility)
   casa::uInt itsVisTypeIgnored;
   
   /// @brief number of ignored samples due to lack of match
   /// @details We increment this number everytime a visibility sample
   /// is ignored because it doesn't match any row in the buffer
   casa::uInt itsNoMatchIgnored;
   
   /// @brief number of ignored samples due to flags in the input accessor
   /// @details This number also accounts for the ignored samples due to 
   /// polarisation index being too large (beyond the range of the buffer)
   casa::uInt itsFlagIgnored;
   
   /// @brief if true, beam index is ignored
   bool itsBeamIndependent;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SYNTHESIS_PREAVGCALBUFFER_H

