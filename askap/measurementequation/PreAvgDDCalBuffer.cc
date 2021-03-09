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

#include <askap/measurementequation/PreAvgDDCalBuffer.h>
#include <askap/measurementequation/ComponentEquation.h>
#include <askap/measurementequation/ImagingEquationAdapter.h>
#include <askap/askap/AskapError.h>
#include <askap/dataaccess/DDCalBufferDataAccessor.h>
#include <askap/scimath/utils/PolConverter.h>

#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.preavgcalbuffer");

using namespace askap;
using namespace askap::synthesis;

/// @brief default constructor
/// @details preaveraging is initialised based on the first encountered accessor
PreAvgDDCalBuffer::PreAvgDDCalBuffer() : itsPolXProducts(0), // set nPol = 0 for now as a proper initialisation is pending
    itsVisTypeIgnored(0), itsNoMatchIgnored(0), itsFlagIgnored(0), itsBeamIndependent(false), itsNDir(1) {}
   
/// @brief constructor with explicit averaging parameters
/// @details This version of the constructor explicitly defines the number of 
/// antennas and calibrator directions to initialise the buffer appropriately.
/// @param[in] nAnt number of antennas, indices are expected to run from 0 to nAnt-1
/// @param[in] nDir number of beams, indices are expected to run from 0 to nDir-1
/// @param[in] nChan number of channels to buffer, 1 (default) is a special case
/// assuming that measurement equation is frequency-independent
PreAvgDDCalBuffer::PreAvgDDCalBuffer(casacore::uInt nAnt, casacore::uInt nDir, casacore::uInt nChan) :
      itsAntenna1(nDir*nAnt*(nAnt-1)/2), itsAntenna2(nDir*nAnt*(nAnt-1)/2), itsBeam(nDir*nAnt*(nAnt-1)/2),
      itsFlag(nDir*nAnt*(nAnt-1)/2,casacore::Int(nChan),4), itsStokes(4), itsNDir(1),
      itsPolXProducts(4,casacore::IPosition(2,int(nDir*nAnt*(nAnt-1)/2),casacore::Int(nChan))),
      itsVisTypeIgnored(0), itsNoMatchIgnored(0), itsFlagIgnored(0), itsBeamIndependent(false)
{
  initialise(nAnt,nDir,nChan);
}  

/// @brief configure beam-independent accumulation
/// @param[in] flag if true, accumulation is beam-independent
void PreAvgDDCalBuffer::beamIndependent(bool flag)
{
  itsBeamIndependent = flag;
}  
   
/// @brief initialise accumulation via an accessor
/// @details This method resets the buffers and sets the shape using the given accessor
/// as a template.
/// @param[in] acc template accessor
/// @param[in] fdp frequency dependency flag, if true a separate buffer is created for every
/// presented spectral channel. Otherwise (default), all channels contribute to the same buffer
/// (which is appropriate for frequency-independent effects)
void PreAvgDDCalBuffer::initialise(const IConstDataAccessor &acc, const bool fdp)
{
  // resize buffers
  const casacore::uInt numberOfRows = acc.nRow();
  const casacore::uInt numberOfPol = acc.nPol();
  const casacore::uInt numberOfChan = fdp ? acc.nChannel() : 1;
  if (itsFlag.shape() != casacore::IPosition(3,numberOfRows, numberOfChan, numberOfPol)) {
      // resizing buffers
      itsAntenna1.resize(numberOfRows);
      itsAntenna2.resize(numberOfRows);
      itsBeam.resize(numberOfRows);
      itsFlag.resize(numberOfRows, numberOfChan, numberOfPol);
      itsPolXProducts.resize(numberOfPol,
          casacore::IPosition(2,casacore::Int(numberOfRows),casacore::uInt(numberOfChan)),false);
      itsStokes.resize(numberOfPol);      
  }
  // initialise buffers
  itsAntenna1 = acc.antenna1();
  itsAntenna2 = acc.antenna2();
  if (itsBeamIndependent) {
      itsBeam.set(0);
  } else {
      itsBeam = acc.feed1();
  }
  const casacore::uInt unusedBeamId = casacore::max(itsBeam)*10;
  const casacore::Vector<casacore::uInt>& feed2 = acc.feed2();
  for (casacore::uInt row=0; row<numberOfRows; ++row) {
       if ((itsBeam[row] != feed2[row]) || (itsBeamIndependent && (itsBeam[row] != 0))) {
           // so it is kept flagged (not very tidy way of doing the check,
           // we can introduce a separate vector to track this condition)
           itsBeam[row] = unusedBeamId;
       }
  }
  itsStokes = acc.stokes();
  // all elements are flagged until at least something is averaged in
  itsFlag.set(true); 
  itsPolXProducts.reset();
  // initialise stats
  itsVisTypeIgnored = 0;
  itsNoMatchIgnored = 0;
  itsFlagIgnored = 0;
}
   
/// @brief initialise accumulation explicitly
/// @details This method resets the buffers and sets the shape to accommodate the given
/// number of antennas and beams (i.e. the buffer size is nDirs*nAnt*(nAnt-1)/2)
/// @param[in] nAnt number of antennas, indices are expected to run from 0 to nAnt-1
/// @param[in] nDir number of beams, indices are expected to run from 0 to nDir-1
/// @param[in] nChan number of channels to buffer, 1 (default) is a special case
/// assuming that measurement equation is frequency-independent
void PreAvgDDCalBuffer::initialise(casacore::uInt nAnt, casacore::uInt nDir, casacore::uInt nChan)
{
  ASKAPDEBUGASSERT(nChan > 0);
  ASKAPASSERT(nAnt > 0);
  const casacore::uInt numberOfRows = nDir*nAnt*(nAnt-1)/2;
  itsNDir = nDir;
  if (itsFlag.shape() != casacore::IPosition(3,int(numberOfRows),int(nChan),4)) {
     // resizing buffers
     itsAntenna1.resize(numberOfRows);
     itsAntenna2.resize(numberOfRows);
     itsBeam.resize(numberOfRows);
     // npol=4
     itsFlag.resize(numberOfRows,casacore::Int(nChan),4);
     itsPolXProducts.resize(4,casacore::Int(nDir),
         casacore::IPosition(2,casacore::Int(numberOfRows), casacore::Int(nChan)), false);
     itsStokes.resize(4);
  }
  // initialising buffers
  itsFlag.set(true); // everything is bad, unless at least one sample is summed into the buffer
  itsPolXProducts.reset();
  
  for (casacore::uInt beam=0,row=0; beam<nDir; ++beam) {
       for (casacore::uInt ant1=0; ant1<nAnt; ++ant1) {
            for (casacore::uInt ant2 = ant1 + 1; ant2<nAnt; ++ant2,++row) {
                 ASKAPDEBUGASSERT(row<numberOfRows);
                 itsAntenna1[row] = ant1;
                 itsAntenna2[row] = ant2;
                 itsBeam[row] = beam;
            }
       }
  }
  
  // we don't track polarisation at this stage leaving this up to the user of this class
  // just fill the vector with Linear Stokes
  for (casacore::uInt pol = 0; pol<itsStokes.nelements(); ++pol) {
       itsStokes[pol] = scimath::PolConverter::stokesFromIndex(pol, casacore::Stokes::XX);
  }
  
  // initialise stats
  itsVisTypeIgnored = 0;
  itsNoMatchIgnored = 0;
  itsFlagIgnored = 0;
}
   
// implemented accessor methods
   
/// The number of rows in this chunk
/// @return the number of rows in this chunk
casacore::uInt PreAvgDDCalBuffer::nRow() const throw()
{
  return itsBeam.nelements();
}
  	
/// The number of spectral channels (equal for all rows)
/// @return the number of spectral channels
casacore::uInt PreAvgDDCalBuffer::nChannel() const throw()
{
  return itsFlag.ncolumn();
}

/// The number of polarization products (equal for all rows)
/// @return the number of polarization products (can be 1,2 or 4)
casacore::uInt PreAvgDDCalBuffer::nPol() const throw()
{
  return itsFlag.nplane();
}

/// First antenna IDs for all rows
/// @return a vector with IDs of the first antenna corresponding
/// to each visibility (one for each row)
const casacore::Vector<casacore::uInt>& PreAvgDDCalBuffer::antenna1() const
{
  return itsAntenna1;
}

/// Second antenna IDs for all rows
/// @return a vector with IDs of the second antenna corresponding
/// to each visibility (one for each row)
const casacore::Vector<casacore::uInt>& PreAvgDDCalBuffer::antenna2() const
{
  return itsAntenna2;
}
  
/// First feed IDs for all rows
/// @return a vector with IDs of the first feed corresponding
/// to each visibility (one for each row)
const casacore::Vector<casacore::uInt>& PreAvgDDCalBuffer::feed1() const
{
  return itsBeam;
}

/// Second feed IDs for all rows
/// @return a vector with IDs of the second feed corresponding
/// to each visibility (one for each row)
const casacore::Vector<casacore::uInt>& PreAvgDDCalBuffer::feed2() const
{
  return itsBeam;
}

/// Cube of flags corresponding to the output of visibility() 
/// @return a reference to nRow x nChannel x nPol cube with flag 
///         information. If True, the corresponding element is flagged bad.
const casacore::Cube<casacore::Bool>& PreAvgDDCalBuffer::flag() const
{
  return itsFlag;
}

/// @brief polarisation type for each product
/// @return a reference to vector containing polarisation types for
/// each product in the visibility cube (nPol() elements).
/// @note All rows of the accessor have the same structure of the visibility
/// cube, i.e. polarisation types returned by this method are valid for all rows.
const casacore::Vector<casacore::Stokes::StokesTypes>& PreAvgDDCalBuffer::stokes() const
{
  return itsStokes;
}

/// @brief helper method to find a match row in the buffer
/// @details It goes over antenna and beam indices and finds a buffer row which 
/// corresponds to the given indices.
/// @param[in] ant1 index of the first antenna
/// @param[in] ant2 index of the second antenna
/// @param[in] beam beam index
/// @return row number in the buffer corresponding to the given (ant1,ant2,beam) or -1 if 
/// there is no match
int PreAvgDDCalBuffer::findMatch(casacore::uInt ant1, casacore::uInt ant2, casacore::uInt beam)
{
  ASKAPDEBUGASSERT(itsAntenna1.nelements() == itsAntenna2.nelements());
  ASKAPDEBUGASSERT(itsAntenna1.nelements() == itsBeam.nelements());
  // we can probably implement a more clever search algorithm here because the 
  // metadata are almost always ordered
  for (casacore::uInt row=0; row<itsAntenna1.nelements(); ++row) {
       if ((itsAntenna1[row] == ant1) && (itsAntenna2[row] == ant2) && (itsBeam[row] == (itsBeamIndependent ? 0 : beam))) {
           return int(row);
       }
  } 
  return -1;
}

/// @brief process one accessor
/// @details This method processes the given accessor and updates the internal 
/// buffers. The measurement equation is used to calculate model visibilities 
/// corresponding to measured visibilities.
/// @param[in] acc input accessor with measured data
/// @param[in] me shared pointer to the measurement equation
/// @param[in] fdp frequency dependency flag (see initialise). It is used if initialisation from accessor
/// is required. Otherwise, it is just checked for consistency (i.e. more than one channel is defined, if it is true)
/// @note only predict method of the measurement equation is used.
void PreAvgDDCalBuffer::accumulate(const IConstDataAccessor &acc,
         const boost::shared_ptr<IMeasurementEquation const> &me, const bool fdp)
{

  if (acc.nRow() == 0) {
      // nothing to process
      return;
  }
  ASKAPCHECK(me, "Uninitialised shared pointer to the measurement equation has been encountered");
  if (itsFlag.nrow() == 0) {
      // initialise using the given accessor as a template
      initialise(acc,fdp);
  } else {
     if (fdp) {
         ASKAPCHECK(nChannel() == acc.nChannel(), 
             "Number of channels in the accessor passed to PreAvgDDCalBuffer::accumulate "
             "doesn't match the number of frequency buffers"); 
     } 
  }
  ASKAPDEBUGASSERT(itsPolXProducts.nPol() > 0);

  // set up the model vis DataAccessor
  accessors::DDCalBufferDataAccessor modelAcc(acc);

  // DDCALTAG
  // let the DataAccessor and MeasurementEquation know how many directions are being calibrated
  if (itsNDir > 1) {
      ASKAPLOG_DEBUG_STR(logger, "DDCalBufferDataAccessor::setNDir("<<itsNDir<<")");
      modelAcc.setNDir(itsNDir);
  
      // Currently only set up for the ComponentEquation and ImagingEquationAdapter classes
      if (typeid(*me) == typeid(ComponentEquation)) {
          try {
              ASKAPLOG_DEBUG_STR(logger, "Casting MeasEq to CompEq to call setNDir("<<itsNDir<<")");
              const boost::shared_ptr<ComponentEquation const>
                  &ce = boost::dynamic_pointer_cast<ComponentEquation const>(me);
              ce->setNDir(itsNDir);
          }
          catch (const std::bad_cast&) {}
      } else if (typeid(*me) == typeid(ImagingEquationAdapter)) {
          try {
              ASKAPLOG_DEBUG_STR(logger, "Casting MeasEq to CompEq to call setNDir("<<itsNDir<<")");
              const boost::shared_ptr<ImagingEquationAdapter const>
                  &ie = boost::dynamic_pointer_cast<ImagingEquationAdapter const>(me);
              ie->setNDir(itsNDir);
          }
          catch (const std::bad_cast&) {}
      }
  }
  
  me->predict(modelAcc);
  const casacore::Cube<casacore::Complex> &modelVis = modelAcc.visibility();
  const casacore::Cube<casacore::Complex> &measuredVis = acc.visibility();
  const casacore::Cube<casacore::Complex> &measuredNoise = acc.noise();
  const casacore::Cube<casacore::Bool> &measuredFlag = acc.flag();
  ASKAPDEBUGASSERT(measuredFlag.nrow() == acc.nRow());
  ASKAPDEBUGASSERT(measuredFlag.ncolumn() == acc.nChannel());
  ASKAPDEBUGASSERT(measuredFlag.nplane() == acc.nPol());
  ASKAPDEBUGASSERT(measuredVis.shape() == measuredNoise.shape());
  ASKAPDEBUGASSERT(measuredVis.shape() == measuredFlag.shape());
  const casacore::uInt bufferNPol = nPol();
  ASKAPDEBUGASSERT(bufferNPol == itsPolXProducts.nPol());
  ASKAPASSERT(modelVis.nrow() == itsNDir*acc.nRow());
  ASKAPASSERT(modelVis.ncolumn() == acc.nChannel());
  ASKAPASSERT(modelVis.nplane() == acc.nPol());
  
  // references to metadata
  const casacore::Vector<casacore::uInt> &beam1 = acc.feed1();
  const casacore::Vector<casacore::uInt> &beam2 = acc.feed2();
  const casacore::Vector<casacore::uInt> &antenna1 = acc.antenna1();
  const casacore::Vector<casacore::uInt> &antenna2 = acc.antenna2(); 

  // bufLenPerDir will often equal acc.nRow(), but not always (e.g. when are auto-correlations are present)
  const int bufLenPerDir = itsPolXProducts.getModelMeasProductShape()(0) / itsPolXProducts.nDir();
  // do a quick sanity check
  ASKAPASSERT(bufLenPerDir*itsPolXProducts.nDir()*(itsPolXProducts.nDir()+1)/2 == 
              itsPolXProducts.getModelProductShape()(0));

  ASKAPCHECK(fdp || (nChannel() == 1), 
     "Only single spectral channel is supported by the pre-averaging calibration buffer in the frequency-independent mode");
  for (casacore::uInt row = 0; row<acc.nRow(); ++row) {
       // could simply force beam1 and beam2 to zero...
       ASKAPCHECK(beam1[row] == 0, "DD calibration is only set up for single-beam data");
       if ((beam1[row] != beam2[row]) || (antenna1[row] == antenna2[row])) {
           // cross-beam correlations and auto-correlations are not supported
           itsVisTypeIgnored += acc.nChannel() * acc.nPol();
           continue;
       }
       // search which row of the buffer corresponds to the same metadata
       const int matchRow = findMatch(antenna1[row],antenna2[row],beam1[row]);
       if (matchRow<0) {
           // there is no match, skip this sample
           itsNoMatchIgnored += acc.nChannel() * acc.nPol();
           continue;
       }
       const casacore::uInt bufRow = casacore::uInt(matchRow);
       ASKAPDEBUGASSERT(bufRow < itsFlag.nrow());

       // the code below works updates itsPolxProds for a row/channel at a time
       for (casacore::uInt chan = 0; chan<acc.nChannel(); ++chan) {
            // in the frequency-independent mode do averaging of all frequency channels together
            const casacore::uInt bufChan = fdp ? chan : 0;

            // if any polarisations are flagged, ignore this visibility
            casacore::Bool JonesFlag = casacore::False;
            for (casacore::uInt pol = 0; pol<acc.nPol(); ++pol) {
                if (measuredFlag(row,chan,pol)) {
                    JonesFlag = casacore::True;
                    break;
                }
            }
            if (JonesFlag) {
                itsFlagIgnored += acc.nPol();
                continue;
            }

            for (casacore::uInt pol = 0; pol<acc.nPol(); ++pol) {
                 if (pol < bufferNPol) {
                     const float visNoise =
                         casacore::square(casacore::real(measuredNoise(row,chan,pol)));
                     const float weight = (visNoise > 0.) ? 1./visNoise : 0.;
                     for (casacore::uInt pol2 = 0; pol2<acc.nPol(); ++pol2) {
                          // different polarisations can have different weight?
                          // ignoring for now

                          // DDCALTAG -- loop over dir may be more efficient on the outside
                          for (casacore::uInt dir = 0; dir<itsNDir; ++dir) {
                              itsPolXProducts.addModelMeasProduct(bufRow+dir*bufLenPerDir, bufChan, pol, pol2,
                                  weight * std::conj(modelVis(row+dir*acc.nRow(),chan,pol)) *
                                                     measuredVis(row,chan,pol2));
                          }
    
                          if (pol2<=pol) {
                              for (casacore::uInt dir = 0; dir<itsNDir; ++dir) {
                                  itsPolXProducts.addModelProduct(bufRow+dir*bufLenPerDir, bufChan, pol, pol2,
                                      weight * std::conj(modelVis(row+dir*acc.nRow(),chan,pol)) *
                                                         modelVis(row+dir*acc.nRow(),chan,pol2));
                              }
                              // generate any model cross products and put after the main buffers
                              casacore::uInt rowOffset = itsNDir*bufLenPerDir;
                              for (casacore::uInt dir = 0; dir<itsNDir-1; ++dir) {
                                   for (casacore::uInt dir2 = dir+1; dir2<itsNDir; ++dir2) {
                                       itsPolXProducts.addModelProduct(bufRow+rowOffset, bufChan, pol, pol2,
                                           weight * std::conj(modelVis(row+dir*acc.nRow(),chan,pol)) *
                                                              modelVis(row+dir2*acc.nRow(),chan,pol2));
                                       rowOffset += bufLenPerDir;
                                   }
                              }
                          }

                     }
                     //std::cout<<"accumulated ("<<bufRow<<","<<pol<<"): "<<model<<
                     //    " "<<measuredVis(row,chan,pol)<<std::endl;
                     // unflag this row because it now has some data
                     itsFlag(bufRow,bufChan,pol) = false;
                 } else {
                     ++itsFlagIgnored;
                 }
            }
       }
  }
}


