/// @file
///
/// @brief measurement equation to apply calibration.
/// @details This is a special type of the measurement equation (i.e. it is not
/// even derived from the scimath::Equation class because it is not solvable). It
/// corrects a chunk of visibilities for calibration, leakages and bandpasses
/// obtained via the solution access interface. Unlike CalibrationMEBase and
/// PreAvgCalMEBase this class has the full measurement equation built in
/// (essentially implemented by the solution access class returning a complete
/// jones matrix for each antenna/beam combination).
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

#ifndef CALIBRATION_APPLICATOR_ME_H
#define CALIBRATION_APPLICATOR_ME_H

// own includes
#include <askap/measurementequation/ICalibrationApplicator.h>
#include <askap/calibaccess/ICalSolutionConstSource.h>
#include <askap/calibaccess/ICalSolutionConstAccessor.h>
#include <askap/measurementequation/CalibrationSolutionHandler.h>
#include <askap/dataaccess/IDataAccessor.h>

// boost includes
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief measurement equation to apply calibration.
/// @details This is a special type of the measurement equation (i.e. it is not
/// even derived from the scimath::Equation class because it is not solvable). It
/// corrects a chunk of visibilities for calibration, leakages and bandpasses
/// obtained via the solution access interface. Unlike CalibrationMEBase and
/// PreAvgCalMEBase this class has the full measurement equation built in
/// (essentially implemented by the solution access class returning a complete
/// jones matrix for each antenna/beam combination). This class handles time-dependence
/// properly provided the solution source interface supports it as well.
/// @ingroup measurementequation
class CalibrationApplicatorME : virtual public ICalibrationApplicator,
                                protected CalibrationSolutionHandler {
public:

  /// @brief constructor
  /// @details It initialises ME for a given solution source.
  /// @param[in] src calibration solution source to work with
  CalibrationApplicatorME(const boost::shared_ptr<accessors::ICalSolutionConstSource> &src);

  /// @brief correct model visibilities for one accessor (chunk).
  /// @details This method corrects the data in the given accessor
  /// (accessed via rwVisibility) for the calibration errors
  /// represented by this measurement equation (i.e. an inversion of
  /// the matrix has been performed).
  /// @param[in] chunk a read-write accessor to work with
  virtual void correct(accessors::IDataAccessor &chunk) const;

  /// @brief predict visibilities for one accessor (chunk).
  /// @details This method corrupts the data in the given accessor
  /// (accessed via rwVisibility) with the calibration errors
  /// represented by this measurement equation (i.e. an inversion of
  /// the matrix has not been performed).
  /// @param[in] chunk a read-write accessor to work with
  virtual void predict(accessors::IDataAccessor &chunk) const;

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
  /// @details It is handy to be able to apply the same solution for all beams. With
  /// this flag set, beam=0 solution will be used for all beams.
  /// @param[in] flag if true, beam=0 calibration is applied to all beams
  virtual void beamIndependent(bool flag);

  /// @brief determines whether channel dependent calibration is used or not
  /// @details It is handy to be able to apply the same solution for all channels. With
  /// this flag set, channel=0 solution will be used for all channels.
  /// @param[in] flag if true, channel=0 calibration is applied to all channels
  virtual void channelIndependent(bool flag);

  /// @brief determines whether time interpolation is used
  /// @details Interpolating solutions in time improves dynamic range
  /// @param[in] flag, if true interpolate between solutions
  virtual void interpolateTime(bool flag);

private:
  /// @brief correct or corrupt visibilities for one accessor
  /// @details This method corrects or corrupts the data in the given accessor
  /// (accessed via rwVisibility) for the calibration errors
  /// represented by this measurement equation (inversion of
  /// the matrix is done if correct=true).
  /// @param[in] chunk a read-write accessor to work with
  /// @param[in] correct correct data if true, corrupt if false
  void generic(accessors::IDataAccessor &chunk, bool correct) const;

  /// @brief correct or corrupt visibilities for one accessor
  /// @details This method corrects or corrupts the data in the given accessor
  /// (accessed via rwVisibility) for the calibration errors
  /// represented by this measurement equation (inversion of
  /// the matrix is done if correct=true).
  /// This is an optimized version for data with exactly 4 polarizations
  /// @param[in] chunk a read-write accessor to work with
  /// @param[in] correct correct data if true, corrupt if false
  void generic4(accessors::IDataAccessor &chunk, bool correct) const;

  /// @brief true, if correct method is to scale the noise estimate
  bool itsScaleNoise;

  /// @brief true, if the data corresponding to invalid solution are flagged
  bool itsFlagAllowed;

  /// @brief true, if beam index should be ignored and beam=0 corrections applied to all beams
  bool itsBeamIndependent;
  /// @brief true, if channel index can be ignored and channel=0 corrections applied to all channels
  bool itsChannelIndependent;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef CALIBRATION_APPLICATOR_ME_H
