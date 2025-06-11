/// @file
///
/// @brief Interface for a measurement equation which can apply calibration
/// @details This interface defines a single 'correct' method which is
/// supposed to correct a chunk of visibilities for calibration errors. Initially
/// we had this functionality in the CalibrationMEBase, but it is worth while to
/// have a separate interface, so we could do a generic correction based on the
/// parameters supplied by the Calibration Solution accessor.
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

#ifndef I_CALIBRATION_APPLICATOR_H
#define I_CALIBRATION_APPLICATOR_H

// own includes
#include <askap/dataaccess/IDataAccessor.h>

// boost includes
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief Interface for a measurement equation which can apply calibration
/// @details This interface defines a single 'correct' method which is
/// supposed to correct a chunk of visibilities for calibration errors. Initially
/// we had this functionality in the CalibrationMEBase, but it is worth while to
/// have a separate interface, so we could do a generic correction based on the
/// parameters supplied by the Calibration Solution accessor.
/// @ingroup measurementequation
struct ICalibrationApplicator {
  /// @brief virtual destructor to keep the compiler happy
  virtual ~ICalibrationApplicator() {}

  /// @brief correct visibilities for one accessor (chunk).
  /// @details This method corrects the data in the given accessor
  /// (accessed via rwVisibility) for the calibration errors
  /// represented by this measurement equation (i.e. an inversion of
  /// the matrix has been performed).
  /// @param[in] chunk a read-write accessor to work with
  virtual void correct(accessors::IDataAccessor &chunk) const = 0;

  /// @brief predict visibilities for one accessor (chunk).
  /// @details This method corrupts the data in the given accessor
  /// (accessed via rwVisibility) with the calibration errors
  /// represented by this measurement equation (i.e. an inversion of
  /// the matrix has not been performed).
  /// @param[in] chunk a read-write accessor to work with
  virtual void predict(accessors::IDataAccessor &chunk) const = 0;

  /// @brief determines whether to scale the noise estimate
  /// @details This is one of the configuration methods, it controlls
  /// whether the noise estimate is scaled aggording to applied calibration
  /// factors or not.
  /// @param[in] scale if true, the noise will be scaled
  virtual void scaleNoise(bool scale) = 0;

  /// @brief determines whether to allow data flagging
  /// @details This is one of the configuration methods, it controlls
  /// the behavior of the correct method in the case when the matrix inversion
  /// fails. If data flagging is allowed, corresponding visibilities are flagged
  /// otherwise an exception is thrown.
  /// @param[in] flag if true, the flagging is allowed
  virtual void allowFlag(bool flag) = 0;

  /// @brief determines whether beam=0 calibration is used for all beams or not
  /// @details It is handy to be able to apply the same solution for all beams. With
  /// this flag set, beam=0 solution will be used for all beams.
  /// @param[in] flag if true, beam=0 calibration is applied to all beams
  virtual void beamIndependent(bool flag) = 0;

  /// @brief determines whether channel dependent calibration is used or not
  /// @details It is handy to be able to apply the same solution for all channels. With
  /// this flag set, channel=0 solution will be used for all channels.
  /// @param[in] flag if true, channel=0 calibration is applied to all channels
  virtual void channelIndependent(bool flag) {};

  /// @brief determines whether leakage calibration is used or not
  /// @details It is handy to know if leakages are present. With
  /// this flag set, the parallel hand solutions will be used.
  /// @param[in] flag if true, leakage free calibration is applied
  virtual void leakageFree(bool flag) {};

  /// @brief determines whether time interpolation is used
  /// @details Interpolating solutions in time improves dynamic range
  /// @param[in] flag, if true interpolate between solutions
  virtual void interpolateTime(bool flag) {};

};

} // namespace synthesis

} // namespace askap

#endif // #ifndef I_CALIBRATION_APPLICATOR_H
