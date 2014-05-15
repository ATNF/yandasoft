/// @file
/// 
/// @brief Base class for generic measurement equation for calibration with pre-averaging.
/// @details This is a base class for a template designed to represent any 
/// possible measurement equation we expect to encounter in calibration. 
/// It is similar to CalibrationMEBase, but implements pre-averaging (or pre-summing to be
/// exact) using PreAvgCalBuffer, so that only one iteration over the data is required.
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

#ifndef PRE_AVG_CAL_ME_BASE_H
#define PRE_AVG_CAL_ME_BASE_H

#include <measurementequation/PreAvgCalBuffer.h>
#include <dataaccess/IDataAccessor.h>
#include <dataaccess/SharedIter.h>
#include <fitting/GenericEquation.h>
#include <fitting/Params.h>
#include <fitting/GenericNormalEquations.h>
#include <fitting/ComplexDiffMatrix.h>
#include <fitting/ComplexDiff.h>

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief Base class for generic measurement equation for calibration with pre-averaging.
/// @details This is a base class for a template designed to represent any 
/// possible measurement equation we expect to encounter in calibration. 
/// It is similar to CalibrationMEBase, but implements pre-averaging (or pre-summing to be
/// exact) using PreAvgCalBuffer, so that only one iteration over the data is required.
/// Because of this, the method to calculate normal equations without parameters is the one
/// which is supposed to be used.
/// @ingroup measurementequation
class PreAvgCalMEBase : public scimath::GenericEquation 
{
public:
  /// @brief constructor setting up only parameters
  /// @param[in] ip Parameters
  PreAvgCalMEBase(const askap::scimath::Params& ip);

  /// @brief Standard constructor using the parameters and the
  /// data iterator.
  /// @param[in] ip Parameters
  /// @param[in] idi data iterator
  /// @param[in] ime measurement equation describing perfect visibilities
  /// @note This version does iteration over the dataset and all accumulation.
  PreAvgCalMEBase(const askap::scimath::Params& ip,
          const accessors::IDataSharedIter& idi, 
          const boost::shared_ptr<IMeasurementEquation const> &ime);
          
  /// @brief accumulate one accessor
  /// @details This method processes one accessor and accumulates the data.
  /// It is essentially a proxy for the accumulate method of the buffer.
  /// @param[in] acc data accessor
  /// @param[in] me measurement equation describing perfect visibilities
  void accumulate(const accessors::IConstDataAccessor &acc,  
          const boost::shared_ptr<IMeasurementEquation const> &me);
          
  /// @brief accumulate all data
  /// @details This method iterates over the whole dataset and accumulates all
  /// the data.
  /// @param[in] idi data iterator
  /// @param[in] ime measurement equation describing perfect visibilities
  void accumulate(const accessors::IDataSharedIter& idi, 
          const boost::shared_ptr<IMeasurementEquation const> &ime);
                    
  /// @brief Predict model visibilities for one accessor (chunk).
  /// @details This class cannot be used for prediction 
  /// (use CalibrationMEBase instead). Therefore this method just 
  /// throws an exception.
  virtual void predict() const;

  /// @brief calculate normal equations in the general form 
  /// @details This method calculates normal equations for the
  /// given set of parameters. It is assumed that some data have already 
  /// been accumulated.
  /// @param[in] ne normal equations to update
  virtual void calcGenericEquations(scimath::GenericNormalEquations &ne) const;
  
  /// @brief initialise accumulation
  /// @details Resets the buffer and configure it to the given number of
  /// antennas and beams.
  /// @param[in] nAnt number of antennas
  /// @param[in] nBeam number of beams
  /// @param[in] nChan number of channels, 1 channel is a special case of
  /// frequency-independent buffering
  void initialise(casa::uInt nAnt, casa::uInt nBeam, casa::uInt nChan = 1);

  /// @brief destructor 
  /// @details This method just prints statistics on the number of
  /// visibilities not accumulated due to various reasons
  ~PreAvgCalMEBase();
  
  /// @brief configure beam-independent case
  /// @details This method configures the underlying buffer to apply 
  /// optimisations if the measurement equation is known to be beam-independent.
  /// By default, the general case is assumed (but it should work for the 
  /// beam-independent case as well)
  /// @param[in] flag if true, the ME is assumed to be beam-independent
  void beamIndependent(const bool flag);
  
protected:  
  /// @brief a helper method to form a ComplexDiffMatrix for a given row
  /// @details This is the only method which depends on the template type.
  /// Therefore in this class it is just declared pure virtual. This method
  /// is used on the most outer level of the measurement equation chain. Therefore,
  /// making it virtual doesn't cause problems with the compile time building of
  /// the measurement equation.
  /// @param[in] acc input data accessor (to define metadata for a given row)
  /// @param[in] row the row number to work with
  /// @return ComplexDiffMatrix encapsulating information about measurement 
  ///         equation corresponding to the given row
  virtual scimath::ComplexDiffMatrix buildComplexDiffMatrix(const accessors::IConstDataAccessor &acc,
                    casa::uInt row) const = 0;
  
  /// @brief check whether the measurement equation is frequency-dependent
  /// @details For frequency-dependent effects the buildComplexDiffMatrix method returns block matrix with 
  /// one block corresponding to every channel (i.e. the size is nPol x nPol*nChannel 
  /// instead of simply nPol x nPol). We need this flag to unroll the matrix multiplication
  /// correctly.
  /// @return true, if the effect is frequency-dependent
  virtual bool isFrequencyDependent() const = 0;
  
  /// @brief a helper method to manage dataset-related statistics   
  /// @details It manages statistics data fields and processes one data accessor.
  /// @param[in] acc input data accessor
  void accumulateStats(const accessors::IConstDataAccessor &acc);
  
  /// @brief a helper method to update metadata associated with the normal equations
  /// @details This method manipulates metadata stored in the normal equations indexed
  /// by the given keyword. If itsNoDataProcessedFlag is true, the given item is removed,
  /// otherwise it is updated with the given value.
  /// @param[in] ne normal equations to work with
  /// @param[in] keyword keyword of the metadata of interest
  /// @param[in] val new value
  void updateMetadata(scimath::GenericNormalEquations &ne, const std::string &keyword, 
                      const double val) const;
private:
  /// @brief buffer with partial sums
  PreAvgCalBuffer itsBuffer;    
  
  // some statistics useful to tag the resulting calibration solution
  // (workers manipulate data, so we need to add these statistics to the normal equations to pass
  // them to the master)
  
  /// @brief true if no data have been accumulated
  bool itsNoDataProcessedFlag;
  
  /// @brief minimal time encountered in the data
  /// @details (the units are the same as returned by time() method of the data accessor).
  /// This data field is undefined if itsNoDataProcessedFlag is true.
  double itsMinTime; 

  /// @brief maximal time encountered in the data
  /// @details (the units are the same as returned by time() method of the data accessor)
  /// This data field is undefined if itsNoDataProcessedFlag is true.
  double itsMaxTime;   
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef PRE_AVG_CAL_ME_BASE_H

