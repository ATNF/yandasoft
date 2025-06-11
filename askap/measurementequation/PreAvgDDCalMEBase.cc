/// @file
/// 
/// @brief Base class for generic measurement equation for calibration with pre-averaging.
/// @details This is a base class for a template designed to represent any 
/// possible measurement equation we expect to encounter in calibration. 
/// It is similar to CalibrationMEBase, but implements pre-averaging (or pre-summing to be
/// exact) using PreAvgDDCalBuffer, so that only one iteration over the data is required.
/// Because of this, the method to calculate normal equations without parameters is the one
/// which is supposed to be used.
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

#include <algorithm>
#include <cstdlib>

#include <askap/measurementequation/PreAvgDDCalMEBase.h>
#include <askap/dataaccess/IDataIterator.h>
#include <askap/askap/AskapError.h>
#include <askap/scimath/fitting/ComplexDiffMatrix.h>
#include <askap/scimath/fitting/ComplexDiff.h>
#include <askap/scimath/fitting/DesignMatrix.h>
#include <askap/scimath/fitting/PolXProducts.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/CalParamNameHelper.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
#include <askap/measurementequation/IonosphericTerm.h>
ASKAP_LOGGER(logger, ".measurementequation.preavgddcalmebase");

using namespace askap;
using namespace askap::synthesis;

/// @brief constructor setting up only parameters
/// @param[in] ip Parameters
PreAvgDDCalMEBase::PreAvgDDCalMEBase(const askap::scimath::Params& ip) :
    scimath::GenericEquation(ip), itsNoDataProcessedFlag(true),
    itsIsIonospheric(false), itsMinTime(0.), itsMaxTime(0.)
     {}

/// @brief Standard constructor using the parameters and the
/// data iterator.
/// @param[in] ip Parameters
/// @param[in] idi data iterator
/// @param[in] ime measurement equation describing perfect visibilities
/// @note This version does iteration over the dataset and all accumulation.
PreAvgDDCalMEBase::PreAvgDDCalMEBase(const askap::scimath::Params& ip,
        const accessors::IDataSharedIter& idi, 
        const boost::shared_ptr<IMeasurementEquation const> &ime) :
        scimath::GenericEquation(ip), itsNoDataProcessedFlag(true), 
        itsIsIonospheric(false), itsMinTime(0.), itsMaxTime(0.)
{
  accumulate(idi,ime);
}

/// @brief configure beam-independent case
/// @details This method configures the underlying buffer to apply 
/// optimisations if the measurement equation is known to be beam-independent.
/// By default, the general case is assumed (but it should work for the 
/// beam-independent case as well)
/// @param[in] flag if true, the ME is assumed to be beam-independent
void PreAvgDDCalMEBase::beamIndependent(const bool flag)
{
  itsBuffer.beamIndependent(flag);
}

          
/// @brief accumulate one accessor
/// @details This method processes one accessor and accumulates the data.
/// It is essentially a proxy for the accumulate method of the buffer.
/// @param[in] acc data accessor
/// @param[in] me measurement equation describing perfect visibilities
void PreAvgDDCalMEBase::accumulate(const accessors::IConstDataAccessor &acc,  
          const boost::shared_ptr<IMeasurementEquation const> &me)
{
  itsBuffer.accumulate(acc,me,isFrequencyDependent());
  accumulateStats(acc);
}
          
/// @brief accumulate all data
/// @details This method iterates over the whole dataset and accumulates all
/// the data.
/// @param[in] idi data iterator
/// @param[in] ime measurement equation describing perfect visibilities
void PreAvgDDCalMEBase::accumulate(const accessors::IDataSharedIter& idi, 
        const boost::shared_ptr<IMeasurementEquation const> &ime)
{
  const bool fdp = isFrequencyDependent();
  accessors::IDataSharedIter iter(idi);
  for (; iter.hasMore(); iter.next()) {
       itsBuffer.accumulate(*iter,ime,fdp);
       accumulateStats(*iter);
  }
}        
                    
/// @brief Predict model visibilities for one accessor (chunk).
/// @details This class cannot be used for prediction 
/// (use CalibrationMEBase instead). Therefore this method just 
/// throws an exception.
void PreAvgDDCalMEBase::predict() const
{
  ASKAPTHROW(AskapError, "PreAvgDDCalMEBase::predict() is not supposed to be called");
}

/// @brief a helper method to manage dataset-related statistics   
/// @details It manages statistics data fields and processes one data accessor.
/// @param[in] acc input data accessor
void PreAvgDDCalMEBase::accumulateStats(const accessors::IConstDataAccessor &acc)
{
  if (itsNoDataProcessedFlag) {
      itsNoDataProcessedFlag = false;
      itsMinTime = acc.time();
      itsMaxTime = itsMinTime;
  } else {
      const double time = acc.time();
      if (time < itsMinTime) {
          itsMinTime = time;
      }
      if (time > itsMaxTime) {
          itsMaxTime = time;
      }
  }
}

/// @brief a helper method to update metadata associated with the normal equations
/// @details This method manipulates metadata stored in the normal equations indexed
/// by the given keyword. If itsNoDataProcessedFlag is true, the given item is removed,
/// otherwise it is updated with the given value.
/// @param[in] ne normal equations to work with
/// @param[in] keyword keyword of the metadata of interest
/// @param[in] val new value
void PreAvgDDCalMEBase::updateMetadata(scimath::GenericNormalEquations &ne, const std::string &keyword, 
                      const double val) const
{
  scimath::Params& metadata = ne.metadata();  
  if (itsNoDataProcessedFlag) {
      if (metadata.has(keyword)) {
          metadata.remove(keyword);
      }
  } else {
      if (metadata.has(keyword)) {
          metadata.update(keyword,val);
      } else {
          metadata.add(keyword,val);
      }
  }
}

void PreAvgDDCalMEBase::initIndexedNormalMatrixAndParameterIndex(scimath::GenericNormalEquations &ne) const
{
    ASKAPLOG_INFO_STR(logger, "Initializing indexed normal matrix and parameter index.");

    const casacore::uInt chanOffset = static_cast<casacore::uInt>(rwParameters()->has("chan_offset") ? rwParameters()->scalarValue("chan_offset") : 0);

    std::set<std::string> baseParamNames;
    for (const auto &name : rwParameters()->freeNames()) {
        std::string baseParamName = scimath::CalParamNameHelper::extractBaseParamName(name);
        baseParamNames.insert(baseParamName);
    }
    size_t nChannelsLocal = itsBuffer.nChannel();
    size_t nBaseParameters = baseParamNames.size();

    // Allocate and initialize the indexed normal matrix.
    ne.initIndexedNormalMatrix(nBaseParameters, nChannelsLocal, chanOffset);

    // Allocate and initialize the indexed data vector.
    ne.initIndexedDataVector(nBaseParameters, nChannelsLocal);

    // Building parameter index map.
    for (const auto& baseName: baseParamNames) {
        for (casa::uInt chan = 0; chan < itsBuffer.nChannel(); ++chan) {
            size_t trueChanNumber = chan + chanOffset;
            std::string parName = scimath::CalParamNameHelper::addChannelInfo(baseName, trueChanNumber);
            ne.addParameterNameToIndexMap(parName);
        }
    }
    ASKAPCHECK(nBaseParameters == ne.getNumberBaseParameters(), "Wrong number of base parameters!");
}

/// @brief calculate normal equations in the general form
/// @details This method calculates normal equations for the
/// given set of parameters. It is assumed that some data have already
/// been accumulated.
/// @param[in] ne normal equations to update
void PreAvgDDCalMEBase::calcGenericEquations(scimath::GenericNormalEquations &ne) const
{
    ASKAPLOG_INFO_STR(logger, "Calculating Generic Equations.");

    const scimath::PolXProducts &polXProducts = itsBuffer.polXProducts();
    ASKAPDEBUGASSERT(itsBuffer.nChannel() > 0);

    if (ne.useIndexedNormalMatrix()) {
        initIndexedNormalMatrixAndParameterIndex(ne);
    }

    ASKAPLOG_INFO_STR(logger, "Building the normal matrix.");

    const casacore::uInt nDir = polXProducts.nDir();
    const casacore::uInt rowsPerDir = itsBuffer.nRow() / nDir;

    for (casacore::uInt row = 0; row < itsBuffer.nRow(); ++row) {
        scimath::ComplexDiffMatrix cdm = buildComplexDiffMatrix(itsBuffer, row);
        ASKAPDEBUGASSERT(cdm.nRow() == itsBuffer.nPol());
        ASKAPDEBUGASSERT(cdm.nColumn() == itsBuffer.nPol() * itsBuffer.nChannel());

        const casacore::uInt dirRow = row % rowsPerDir;
        const casacore::uInt dir = row / rowsPerDir;

        for (casa::uInt chan = 0; chan < itsBuffer.nChannel(); ++chan) {

            // take a slice, this takes care of indices along the first two axes (row and channel)
            const scimath::PolXProducts pxpSlice = polXProducts.roSlice(row, chan);
            if (isFrequencyDependent()) {
                // cdm is a block matrix
                const size_t columnOffset = chan * itsBuffer.nPol();
                const size_t neChan = isIonospheric() ? 1 : chan;
                ne.add(cdm, pxpSlice, columnOffset, neChan);
            } else {
                // cdm is a normal matrix
                ne.add(cdm, pxpSlice);
            }
        }
        for (casa::uInt dir2 = 0; dir2 < nDir; ++dir2) {

            // in the model product buffer, store the nDir self products then the nDir*(nDir-1)/2 cross products
            casa::uInt modelBlock, doConj;
            if (dir == dir2) {
                continue;
            }
            else if (dir < dir2) {
                modelBlock = (nDir-1)*(dir+1) - (dir+1)*dir/2 + dir2;
                doConj = 0;
            }
            else {
                modelBlock = (nDir-1)*(dir2+1) - (dir2+1)*dir2/2 + dir;
                doConj = 1;
            }
            const casa::uInt blockRow = modelBlock * rowsPerDir + dirRow;
            const casa::uInt row2 = dir2 * rowsPerDir + dirRow;

            scimath::ComplexDiffMatrix cdm2 = buildComplexDiffMatrix(itsBuffer, row2);

            for (casa::uInt chan = 0; chan < itsBuffer.nChannel(); ++chan) {
                const scimath::PolXProducts pxpSlice2 = polXProducts.roModelSlice(blockRow, chan);
                if (isFrequencyDependent()) {
                    // cdm is a block matrix
                    const size_t columnOffset = chan * itsBuffer.nPol();
                    const size_t neChan = isIonospheric() ? 1 : chan;
                    ne.DDupdate(cdm, cdm2, pxpSlice2, doConj, columnOffset, neChan);
                } else {
                    // cdm is a normal matrix
                    ne.DDupdate(cdm, cdm2, pxpSlice2, doConj, 0, 0);
                }

            }
        }
    }
    updateMetadata(ne, "min_time", itsMinTime);
    updateMetadata(ne, "max_time", itsMaxTime);
}
  
/// @brief initialise accumulation
/// @details Resets the buffer and configure it to the given number of
/// antennas, beams and channels
/// @param[in] nAnt number of antennas
/// @param[in] nDir number of beams
/// @param[in] nChan number of channels
void PreAvgDDCalMEBase::initialise(casacore::uInt nAnt, casacore::uInt nDir, casacore::uInt nChan)
{
  itsBuffer.initialise(nAnt,nDir,nChan);
  itsNoDataProcessedFlag = true;
  itsMinTime = 0.;
  itsMaxTime = 0.;
}

/// @brief destructor 
/// @details This method just prints statistics on the number of
/// visibilities not accumulated due to various reasons
PreAvgDDCalMEBase::~PreAvgDDCalMEBase()
{
  ASKAPLOG_DEBUG_STR(logger, "PreAvgDDCalMEBase statistics on ignored visibilities");
  ASKAPLOG_DEBUG_STR(logger, "   ignored due to type (e.g. autocorrelations): "<<itsBuffer.ignoredDueToType());
  ASKAPLOG_DEBUG_STR(logger, "   no match found for baseline/beam: "<<itsBuffer.ignoredNoMatch());
  ASKAPLOG_DEBUG_STR(logger, "   ignored because of flags: "<<itsBuffer.ignoredDueToFlags());
  if (!itsNoDataProcessedFlag) {
      ASKAPLOG_DEBUG_STR(logger, "Last solution calculated using the data from time range ("<<itsMinTime<<", "<<itsMaxTime<<")");
  }
}

/// @brief check that some data were accumulated for the given antenna, beam and channel
////@details check flags in the buffer that at least one element is unflagged for the 
/// given antenna, beam and channel. This allows us to fix parameters for which we don't have
/// data.
/// @param[in] ant antenna index to query
/// @param[in] beam beam index to query
/// @param[in] pol polarisation index to check
/// @param[in] chan channel number to query or 0 for channel-independent case
/// @return true if there are some data accumulated for the given antenna, beam and channel
bool PreAvgDDCalMEBase::hasDataAccumulated(casa::uInt ant, casa::uInt beam, casa::uInt pol, casa::uInt chan)
{
   ASKAPCHECK(chan < itsBuffer.nChannel(), "Requested channel id = "<<chan<<" is outside the buffer size");
   ASKAPCHECK(pol < itsBuffer.nPol(), "Requested polarisation id = "<<pol<<" is outside the buffer size");
   const casa::Cube<casa::Bool> &flag = itsBuffer.flag();
   const casa::Vector<casa::uInt> &antenna1 = itsBuffer.antenna1();
   const casa::Vector<casa::uInt> &antenna2 = itsBuffer.antenna2();
   const casa::Vector<casa::uInt> &beam1= itsBuffer.feed1();
   ASKAPDEBUGASSERT(chan < flag.ncolumn());
   ASKAPDEBUGASSERT(pol < flag.nrow());
   ASKAPDEBUGASSERT(antenna1.nelements() == antenna2.nelements());
   ASKAPDEBUGASSERT(antenna1.nelements() == beam1.nelements());
   ASKAPDEBUGASSERT(flag.nplane() == beam1.nelements());
   for (casa::uInt row = 0; row < flag.nplane(); ++row)  {
        if (beam1[row] == beam  && (antenna1[row] == ant || antenna2[row] == ant)) {
            if (!flag(pol, chan, row)) {
                return true;
            }
        }
   }
   return false;
}

/// @brief polarisation type for each product
/// @return a reference to vector containing polarisation types for
/// each product the pre-averaging buffer is setup to handle
const casa::Vector<casa::Stokes::StokesTypes>& PreAvgDDCalMEBase::stokes() const
{
   return itsBuffer.stokes();
}

  
