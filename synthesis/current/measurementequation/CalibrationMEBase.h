/// @file
/// 
/// @brief Base class for generic measurement equation for calibration.
/// @details This is a base class for a template designed to represent any 
/// possible measurement equation we expect to encounter in calibration. 
/// It is a result of evolution of the former GainCalibrationEquation, which 
/// will probably be completely substituted by this template in the future. 
/// See CalibrationME template for more details. This class contains all
/// functionality, which doesn't depend on the template parameter.
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

#ifndef CALIBRATION_ME_BASE_H
#define CALIBRATION_ME_BASE_H


// own includes
#include <fitting/GenericEquation.h>
#include <fitting/Params.h>
#include <fitting/GenericNormalEquations.h>
#include <measurementequation/ICalibrationApplicator.h>

#include <measurementequation/GenericMultiChunkEquation.h>
#include <dataaccess/IDataAccessor.h>
#include <fitting/ComplexDiffMatrix.h>
#include <fitting/ComplexDiff.h>
#include <fitting/Equation.h>

// boost includes
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {


/// @brief Base class for generic measurement equation for calibration.
/// @details This is a base class for a template designed to represent any 
/// possible measurement equation we expect to encounter in calibration. 
/// It is a result of evolution of the former GainCalibrationEquation, which 
/// will probably be completely substituted by this template in the future. 
/// See CalibrationME template for more details. This class contains all
/// functionality, which doesn't depend on the template parameter.
/// @ingroup measurementequation
class CalibrationMEBase : public GenericMultiChunkEquation,
                          virtual public ICalibrationApplicator
{
public:

  /// @brief Standard constructor using the parameters and the
  /// data iterator.
  /// @param[in] ip Parameters
  /// @param[in] idi data iterator
  /// @param[in] ime measurement equation describing perfect visibilities
  /// @note In the future, measurement equations will work with accessors
  /// only, and, therefore, the dependency on iterator will be removed
  CalibrationMEBase(const askap::scimath::Params& ip,
          const accessors::IDataSharedIter& idi, 
          const boost::shared_ptr<IMeasurementEquation const> &ime);
  
  /// @brief Predict model visibilities for one accessor (chunk).
  /// @details This version of the predict method works with
  /// a single chunk of data only. It seems that all measurement
  /// equations should work with accessors rather than iterators
  /// (i.e. the iteration over chunks should be moved to the higher
  /// level, outside this class). In the future, I expect that
  /// predict() without parameters will be deprecated.
  /// @param[in] chunk a read-write accessor to work with
  virtual void predict(accessors::IDataAccessor &chunk) const;

  /// @brief correct model visibilities for one accessor (chunk).
  /// @details This method corrects the data in the given accessor
  /// (accessed via rwVisibility) for the calibration errors 
  /// represented by this measurement equation (i.e. an inversion of
  /// the matrix has been performed). 
  /// @param[in] chunk a read-write accessor to work with
  /// @note Need to think what to do in the inversion is unsuccessful
  /// e.g. amend flagging information? This is not yet implemented as
  /// existing accessors would throw an exception if flagging info is 
  /// changed.
  virtual void correct(accessors::IDataAccessor &chunk) const;

  /// @brief Calculate the normal equation for one accessor (chunk).
  /// @details This version of the method works on a single chunk of
  /// data only (one iteration).It seems that all measurement
  /// equations should work with accessors rather than iterators
  /// (i.e. the iteration over chunks should be moved to the higher
  /// level, outside this class). In the future, I expect that
  /// the variant of the method without parameters will be deprecated.
  /// @param[in] chunk a read-write accessor to work with
  /// @param[in] ne Normal equations
  virtual void calcGenericEquations(const accessors::IConstDataAccessor &chunk,
                         askap::scimath::GenericNormalEquations& ne) const;
  
    
  using GenericMultiChunkEquation::predict;
  using GenericMultiChunkEquation::calcEquations;
  using GenericMultiChunkEquation::calcGenericEquations;
 
  /// @brief determines whether to scale the noise estimate
  /// @details This is one of the configuration methods, it controlls
  /// whether the noise estimate is scaled aggording to applied calibration
  /// factors or not.
  /// @param[in] scale if true, the noise will be scaled
  virtual void scaleNoise(bool scale);
  
  /// @brief determines whether to allow data flagging
  /// @details This is one of the configuration methods, it controlls
  /// the behavior of the correct method in the case when the matrix inversion
  /// fails. If data flagging is allowed, corresponding visibilities are flagged
  /// otherwise an exception is thrown.
  /// @param[in] flag if true, the flagging is allowed
  virtual void allowFlag(bool flag);
  
  /// @brief determines whether beam=0 calibration is used for all beams or not
  /// @details It is handy to be able to apply the same solution for all beams. 
  /// With this flag set, beam=0 solution will be used for all beams.
  /// @param[in] flag if true, beam=0 calibration is applied to all beams
  virtual void beamIndependent(bool flag);

  /// @brief normalise gain update amplitudes before application
  /// @details It may be useful to only update phases during calibration, for
  /// instance during self-calibration with a limited number of antennas.
  /// @param[in] flag if true, only phases will be updated when applying
  /// calibration solutions
  virtual void phaseOnly(bool flag);
       
 
protected:  

  /// @brief form correction matrix
  /// @details This is a helper method which forms correction matrix out of
  /// ComplexDiffMatrix (to encapsulate the code used both in frequency-dependent
  /// and frequency-independent cases).
  /// @param[in] cdm a square ComplexDiffMatrix describing the effect
  /// @return correction matrix
  static casa::Matrix<casa::Complex> getCorrectionMatrix(const scimath::ComplexDiffMatrix &cdm);

  /// @brief a helper method to form a ComplexDiffMatrix for a given row
  /// @details This is the only method which depends on the template type.
  /// Therefore in this class it is just declared pure virtual. This method
  /// is used on the most outer level of the measurement equation chain. Therefore,
  /// making it virtual doesn't cause problems with the compile time building of
  /// the measurement equation.
  /// @param[in] acc input data accessor with the perfect visibilities
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
      
private:
  /// @brief measurement equation giving perfect visibilities
  boost::shared_ptr<IMeasurementEquation const> itsPerfectVisME;  
};


} // namespace synthesis

} // namespace askap

#endif // #ifndef CALIBRATION_ME_BASE_H
