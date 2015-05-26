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

#include <measurementequation/CalibrationApplicatorME.h>
#include <casa/Arrays/MatrixMath.h>
#include <scimath/Mathematics/MatrixMathLA.h>
#include <scimath/Mathematics/RigidVector.h>
#include <askap/AskapError.h>
#include <utils/PolConverter.h>
#include <dataaccess/IFlagAndNoiseDataAccessor.h>


#include <askap/AskapUtil.h>
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation");

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief constructor 
/// @details It initialises ME for a given solution source.
/// @param[in] src calibration solution source to work with
CalibrationApplicatorME::CalibrationApplicatorME(const boost::shared_ptr<accessors::ICalSolutionConstSource> &src) :
     CalibrationSolutionHandler(src), itsScaleNoise(false), itsFlagAllowed(false), itsBeamIndependent(false) {}

/// @brief correct model visibilities for one accessor (chunk).
/// @details This method corrects the data in the given accessor
/// (accessed via rwVisibility) for the calibration errors 
/// represented by this measurement equation (i.e. an inversion of
/// the matrix has been performed). 
/// @param[in] chunk a read-write accessor to work with
void CalibrationApplicatorME::correct(accessors::IDataAccessor &chunk) const
{
  casa::Cube<casa::Complex> &rwVis = chunk.rwVisibility();
  ASKAPDEBUGASSERT(rwVis.nelements());
  updateAccessor(chunk.time());
  const casa::Vector<casa::uInt>& antenna1 = chunk.antenna1();
  const casa::Vector<casa::uInt>& antenna2 = chunk.antenna2();
  const casa::Vector<casa::uInt>& beam1 = chunk.feed1();
  const casa::Vector<casa::uInt>& beam2 = chunk.feed2();
  const casa::Vector<casa::Stokes::StokesTypes> stokes = chunk.stokes();   
  
  const casa::uInt nPol = chunk.nPol();
  ASKAPDEBUGASSERT(nPol <= 4);
  casa::Matrix<casa::Complex> mueller(nPol, nPol);
  casa::Matrix<casa::Complex> reciprocal(nPol, nPol);
  
  casa::RigidVector<casa::uInt, 4> indices(0u);
  for (casa::uInt pol = 0; pol<nPol; ++pol) {
       indices(pol) = scimath::PolConverter::getIndex(stokes[pol]);
  }
  
  boost::shared_ptr<accessors::IFlagAndNoiseDataAccessor> noiseAndFlagDA;
  // attempt to cast interface only if we need it
  if (itsScaleNoise || itsFlagAllowed) {
      boost::shared_ptr<accessors::IDataAccessor> chunkPtr(&chunk, utility::NullDeleter());
      ASKAPDEBUGASSERT(chunkPtr);
      noiseAndFlagDA = boost::dynamic_pointer_cast<accessors::IFlagAndNoiseDataAccessor>(chunkPtr);
  }
  
  // full 4x4 Mueller matrix
  casa::SquareMatrix<casa::Complex, 2> fullMueller(casa::SquareMatrix<casa::Complex, 2>::General);
  
  for (casa::uInt row = 0; row < chunk.nRow(); ++row) {
       casa::Matrix<casa::Complex> thisRow = rwVis.yzPlane(row);
       for (casa::uInt chan = 0; chan < chunk.nChannel(); ++chan) {
            casa::SquareMatrix<casa::Complex, 2> jones1 = calSolution().jones(antenna1[row], 
                             itsBeamIndependent ? 0 : beam1[row], chan);
            casa::SquareMatrix<casa::Complex, 2> jones2 = calSolution().jones(antenna2[row],
                             itsBeamIndependent ? 0 : beam2[row], chan);
            for (casa::uInt i = 0; i < nPol; ++i) {
                 for (casa::uInt j = 0; j < nPol; ++j) {
                      const casa::uInt index1 = indices(i);
                      const casa::uInt index2 = indices(j);
                      mueller(i,j) = jones1(index1 / 2, index2 / 2) * conj(jones2(index1 % 2, index2 % 2));
                 }
            }
            
            casa::Complex det = 0.;
            invert(reciprocal, det, mueller);

            casa::Vector<casa::Complex> thisChan = thisRow.row(chan);

            const float detThreshold = 1e-25;
            if (itsFlagAllowed) {
                if (casa::abs(det)<detThreshold) {
                    ASKAPCHECK(noiseAndFlagDA, "Accessor type passed to CalibrationApplicatorME does not support change of flags");
                    noiseAndFlagDA->rwFlag().yzPlane(row).row(chan).set(true);
                    thisChan.set(0.);
                    continue;
                }
            } else {
              ASKAPCHECK(casa::abs(det)>detThreshold, "Unable to apply calibration for (antenna1,beam1)=("<<antenna1[row]<<","<<beam1[row]<<") and (antenna2,beam2)=("<<antenna2[row]<<
                               ","<<beam2[row]<<"), time="<<chunk.time()/86400.-55000<<" determinate is too close to 0. D="<<casa::abs(det)<<" matrix="<<mueller
                       <<" jones1="<<jones1.matrix()<<" jones2="<<jones2.matrix()<<" dir="<<askap::printDirection(chunk.pointingDir1()[row]));           
            }           
            const casa::Vector<casa::Complex> origVis = thisChan.copy();
            ASKAPDEBUGASSERT(thisChan.nelements() == nPol);
            // matrix multiplication
            for (casa::uInt pol = 0; pol < nPol; ++pol) {
                 casa::Complex temp(0.,0.);
                 for (casa::uInt k = 0; k < nPol; ++k) {
                     temp += reciprocal(pol,k) * origVis[k];
                 }
                 thisChan[pol] = temp;
            }
            if (itsScaleNoise) {
                ASKAPCHECK(noiseAndFlagDA, "Accessor type passed to CalibrationApplicatorME does not support change of the noise estimate");
                casa::Vector<casa::Complex> thisChanNoise = noiseAndFlagDA->rwNoise().yzPlane(row).row(chan);
                const casa::Vector<casa::Complex> origNoise = thisChanNoise.copy();
                ASKAPDEBUGASSERT(thisChanNoise.nelements() == nPol);
                // propagating noise estimate through the matrix multiplication
                for (casa::uInt pol = 0; pol < nPol; ++pol) {
                     float tempRe = 0., tempIm = 0.;
                     for (casa::uInt k = 0; k < nPol; ++k) {
                         tempRe += casa::square(casa::real(reciprocal(pol,k)) * casa::real(origNoise[k])) + 
                                   casa::square(casa::imag(reciprocal(pol,k)) * casa::imag(origNoise[k]));
                         tempIm += casa::square(casa::real(reciprocal(pol,k)) * casa::imag(origNoise[k])) + 
                                   casa::square(casa::imag(reciprocal(pol,k)) * casa::real(origNoise[k]));                                   
                     }
                     thisChanNoise[pol] = casa::Complex(sqrt(tempRe), sqrt(tempIm));
                }
            }
       }       
  }
}

/// @brief determines whether to scale the noise estimate
/// @details This is one of the configuration methods, it controlls
/// whether the noise estimate is scaled aggording to applied calibration
/// factors or not.
/// @param[in] scale if true, the noise will be scaled
void CalibrationApplicatorME::scaleNoise(bool scale)
{
  itsScaleNoise = scale;
  if (itsScaleNoise) {
     ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will scale noise estimate when applying calibration solution");
  } else {
     ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will not change the noise estimate");
  }
}
  
/// @brief determines whether to allow data flagging
/// @details This is one of the configuration methods, it controlls
/// the behavior of the correct method in the case when the matrix inversion
/// fails. If data flagging is allowed, corresponding visibilities are flagged
/// otherwise an exception is thrown.
/// @param[in] flag if true, the flagging is allowed
void CalibrationApplicatorME::allowFlag(bool flag)
{
  itsFlagAllowed = flag;
  if (itsFlagAllowed) {
     ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will flag visibilites if inversion of the calibration solution fails");
  } else {
     ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will throw an exception if inversion of the calibration solution fails");
  }
}

/// @brief determines whether beam=0 calibration is used for all beams or not
/// @details It is handy to be able to apply the same solution for all beams. With 
/// this flag set, beam=0 solution will be used for all beams.
/// @param[in] flag if true, beam=0 calibration is applied to all beams
void CalibrationApplicatorME::beamIndependent(bool flag)
{
  itsBeamIndependent = flag;
  if (itsBeamIndependent) {
      ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will apply beam=0 calibration solutions to all beams encountered");
  } else {
      ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will apply beam-dependent calibration solutions");
  }
}


} // namespace synthesis

} // namespace askap

