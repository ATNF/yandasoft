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

#include <askap/measurementequation/CalibrationApplicatorME.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/scimath/Mathematics/MatrixMathLA.h>
#include <casacore/scimath/Mathematics/RigidVector.h>
#include <casacore/scimath/Mathematics/SquareMatrix.h>
#include <askap/askap/AskapError.h>
#include <askap/scimath/utils/PolConverter.h>
#include <askap/dataaccess/IFlagAndNoiseDataAccessor.h>


#include <askap/askap/AskapUtil.h>
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation");

#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief constructor
/// @details It initialises ME for a given solution source.
/// @param[in] src calibration solution source to work with
CalibrationApplicatorME::CalibrationApplicatorME(const boost::shared_ptr<accessors::ICalSolutionConstSource> &src) :
     CalibrationSolutionHandler(src), itsScaleNoise(false), itsFlagAllowed(false), itsBeamIndependent(false),
     itsChannelIndependent(false)
{}

void CalibrationApplicatorME::correct(accessors::IDataAccessor &chunk) const
{
    generic(chunk,true);
}

// Note that the CalibrationME class can also be used to do predict, but not (yet) DDCal type predict
// This predict method was created to allow the DDCal path to be implemented more easily
void CalibrationApplicatorME::predict(accessors::IDataAccessor &chunk) const
{
    generic(chunk,false);
}

void CalibrationApplicatorME::generic(accessors::IDataAccessor &chunk, bool correct) const
{
  const casacore::uInt nPol = chunk.nPol();
  casacore::RigidVector<casacore::uInt, 4> indices(0u);
  const casacore::Vector<casacore::Stokes::StokesTypes> stokes = chunk.stokes();
  for (casacore::uInt pol = 0; pol<nPol; ++pol) {
       indices(pol) = scimath::PolConverter::getIndex(stokes[pol]);
  }

 // Use the optimized version if we can: 4 pols in canonical order
  // MV: this seems like a hack/not the C++ way of doing it. Code duplication/technical dept
  if (nPol==4 && indices(0)==0 && indices(1)==1 && indices(2)==2 && indices(3)==3) {
      generic4(chunk, correct);
      return;
  }

  casa::Cube<casa::Complex> &rwVis = chunk.rwVisibility();
  ASKAPDEBUGASSERT(rwVis.nelements());
  bool solutionAccessorNeedsTime = true;
  const casa::Vector<casa::uInt>& antenna1 = chunk.antenna1();
  const casa::Vector<casa::uInt>& antenna2 = chunk.antenna2();
  const casa::Vector<casa::uInt>& beam1 = chunk.feed1();
  const casa::Vector<casa::uInt>& beam2 = chunk.feed2();

  ASKAPDEBUGASSERT(nPol <= 4);
  casacore::Matrix<casacore::Complex> mueller(nPol, nPol);
  casacore::Matrix<casacore::Complex> reciprocal(nPol, nPol);

  boost::shared_ptr<accessors::IFlagAndNoiseDataAccessor> noiseAndFlagDA;
  // attempt to cast interface only if we need it
  if (itsScaleNoise || itsFlagAllowed) {
      boost::shared_ptr<accessors::IDataAccessor> chunkPtr(&chunk, utility::NullDeleter());
      ASKAPDEBUGASSERT(chunkPtr);
      noiseAndFlagDA = boost::dynamic_pointer_cast<accessors::IFlagAndNoiseDataAccessor>(chunkPtr);
  }

  // full 4x4 Mueller matrix
  casacore::SquareMatrix<casacore::Complex, 2> fullMueller(casacore::SquareMatrix<casacore::Complex, 2>::General);

  // MV: we have to use rwFlag to avoid caching the wrong reference (it's a bit ugly)
  const casacore::Cube<casacore::Bool> &flag = noiseAndFlagDA ? noiseAndFlagDA->rwFlag() : chunk.flag();

  for (casacore::uInt row = 0; row < chunk.nRow(); ++row) {
       //casacore::Matrix<casacore::Complex> thisRow = rwVis.yzPlane(row);
       casacore::Matrix<casacore::Complex> thisRow = rwVis.xyPlane(row);
       //casacore::Matrix<casacore::Bool> thisRowFlag = flag.yzPlane(row);
       casacore::Matrix<casacore::Bool> thisRowFlag = flag.xyPlane(row);
       for (casacore::uInt chan = 0; chan < chunk.nChannel(); ++chan) {
            bool allFlagged = true;
            // we don't really support partial polarisation flagging, but to avoid nasty surprises it is better to flag such samples completely.
            bool needFlag = false;
            //ASKAPDEBUGASSERT(thisRowFlag.ncolumn() == nPol);
            ASKAPDEBUGASSERT(thisRowFlag.nrow() == nPol);
            for (casacore::uInt pol = 0; pol < nPol; ++pol) {
                 if (thisRowFlag(chan,pol)) {
                     needFlag = true;
                 } else {
                     allFlagged = false;
                 }
            }
            // don't bother with the rest of processing if the sample is flagged anyway in its entirety
            if (allFlagged) {
                continue;
            }
            // the whole accessor has a single timestamp. Putting the call here allows to avoid calling it for flagged data
            if (solutionAccessorNeedsTime) {
                solutionAccessorNeedsTime = false;
                updateAccessor(chunk.time());
            }
            const bool validSolution =  calSolution().jonesValid(antenna1[row],
                     itsBeamIndependent ? 0 : beam1[row], chan) &&
                     calSolution().jonesValid(antenna2[row], itsBeamIndependent ? 0 : beam2[row], chan);

            //ASKAPLOG_DEBUG_STR(logger, "row = "<<row<<" chan = "<<chan<<" ant1 = "<<antenna1[row]<<" ant2 = "<<antenna2[row]<<
            //            " beam = "<<beam1[row]<<" allFlagged: "<<allFlagged<<" needFlag: "<<needFlag<<" validSolution: "<<validSolution);
            casacore::Complex det = 0.;

            if (validSolution) {
                casacore::SquareMatrix<casacore::Complex, 2> jones1 = calSolution().jones(antenna1[row],
                                  itsBeamIndependent ? 0 : beam1[row], chan);
                casacore::SquareMatrix<casacore::Complex, 2> jones2 = calSolution().jones(antenna2[row],
                                  itsBeamIndependent ? 0 : beam2[row], chan);
                for (casacore::uInt i = 0; i < nPol; ++i) {
                     for (casacore::uInt j = 0; j < nPol; ++j) {
                          const casacore::uInt index1 = indices(i);
                          const casacore::uInt index2 = indices(j);
                          mueller(i,j) = jones1(index1 / 2, index2 / 2) * conj(jones2(index1 % 2, index2 % 2));
                     }
                }

                if (correct) {
                    invert(reciprocal, det, mueller);
                } else {
                    reciprocal = mueller;
                    det = 1.0;
                }
            }

            // mv: this is actually known to be slow (casacore slice from slice), need to change this to access given channel as matrix
            //casacore::Vector<casacore::Complex> thisChan = thisRow.row(chan);
            casacore::Vector<casacore::Complex> thisChan = thisRow.column(chan);

            const float detThreshold = 1e-25;
            if (itsFlagAllowed) {
                if (casacore::abs(det)<detThreshold || !validSolution || needFlag) {
                    ASKAPCHECK(noiseAndFlagDA, "Accessor type passed to CalibrationApplicatorME does not support change of flags");
                    noiseAndFlagDA->rwFlag().xyPlane(row).column(chan).set(true);
                    thisChan.set(0.);
                    continue;
                }
            } else {
              ASKAPCHECK(validSolution && !needFlag, "Encountered unflagged data and invalid solution, but flagging samples has not been allowed");
              ASKAPCHECK(casacore::abs(det)>detThreshold, "Unable to apply calibration for (antenna1,beam1)=("<<antenna1[row]<<","<<beam1[row]<<") and (antenna2,beam2)=("<<antenna2[row]<<
                               ","<<beam2[row]<<"), time="<<chunk.time()/86400.-55000<<" determinate is too close to 0. D="<<casacore::abs(det)<<" matrix="<<mueller
                       //<<" jones1="<<jones1.matrix()<<" jones2="<<jones2.matrix()
                       <<" dir="<<askap::printDirection(chunk.pointingDir1()[row]));
            }
            const casacore::Vector<casacore::Complex> origVis = thisChan.copy();
            ASKAPDEBUGASSERT(thisChan.nelements() == nPol);
            // matrix multiplication
            for (casacore::uInt pol = 0; pol < nPol; ++pol) {
                 casacore::Complex temp(0.,0.);
                 for (casacore::uInt k = 0; k < nPol; ++k) {
                     temp += reciprocal(pol,k) * origVis[k];
                 }
                 thisChan[pol] = temp;
            }
            if (itsScaleNoise) {
                ASKAPCHECK(noiseAndFlagDA, "Accessor type passed to CalibrationApplicatorME does not support change of the noise estimate");
                //casacore::Vector<casacore::Complex> thisChanNoise = noiseAndFlagDA->rwNoise().yzPlane(row).row(chan);
                casacore::Vector<casacore::Complex> thisChanNoise = noiseAndFlagDA->rwNoise().xyPlane(row).column(chan);
                const casacore::Vector<casacore::Complex> origNoise = thisChanNoise.copy();
                ASKAPDEBUGASSERT(thisChanNoise.nelements() == nPol);
                // propagating noise estimate through the matrix multiplication
                for (casacore::uInt pol = 0; pol < nPol; ++pol) {
                     float tempRe = 0., tempIm = 0.;
                     for (casacore::uInt k = 0; k < nPol; ++k) {
                         tempRe += casacore::square(casacore::real(reciprocal(pol,k)) * casacore::real(origNoise[k])) +
                                   casacore::square(casacore::imag(reciprocal(pol,k)) * casacore::imag(origNoise[k]));
                         tempIm += casacore::square(casacore::real(reciprocal(pol,k)) * casacore::imag(origNoise[k])) +
                                   casacore::square(casacore::imag(reciprocal(pol,k)) * casacore::real(origNoise[k]));
                     }
                     thisChanNoise[pol] = casacore::Complex(sqrt(tempRe), sqrt(tempIm));
                }
            }
       }
  }
}

void CalibrationApplicatorME::generic4(accessors::IDataAccessor &chunk, bool correct) const
{
  const casacore::uInt nPol = chunk.nPol();
  ASKAPDEBUGASSERT(nPol == 4);

  casacore::Cube<casacore::Complex>& rwVis = chunk.rwVisibility();
  casacore::Cube<casacore::Bool> rwFlag;
  ASKAPDEBUGASSERT(rwVis.nelements());
  bool solutionAccessorNeedsTime = true;
  const casa::Vector<casa::uInt>& antenna1 = chunk.antenna1();
  const casa::Vector<casa::uInt>& antenna2 = chunk.antenna2();
  const casa::Vector<casa::uInt>& beam1 = chunk.feed1();
  const casa::Vector<casa::uInt>& beam2 = chunk.feed2();

  boost::shared_ptr<accessors::IFlagAndNoiseDataAccessor> noiseAndFlagDA;
  // attempt to cast interface only if we need it
  if (itsScaleNoise || itsFlagAllowed) {
      boost::shared_ptr<accessors::IDataAccessor> chunkPtr(&chunk, utility::NullDeleter());
      ASKAPDEBUGASSERT(chunkPtr);
      noiseAndFlagDA = boost::dynamic_pointer_cast<accessors::IFlagAndNoiseDataAccessor>(chunkPtr);
      ASKAPCHECK(noiseAndFlagDA, "Accessor type passed to CalibrationApplicatorME does not support change of flags");
      rwFlag.reference(noiseAndFlagDA->rwFlag());
  }
  // MV: we can do something more clever, but accessing original read-only flags is not too bad
  // it will be the same as rwFlag above by reference
  const casacore::Cube<casacore::Bool> &flag = chunk.flag();

  casa::SquareMatrix<casa::Complex, 4> mueller;
  const float detThreshold = 1e-25;

  const casa::uInt nChan = chunk.nChannel();
  const casa::uInt nRow = chunk.nRow();
  casa::RigidVector<casa::Complex,4> vis;
  const casa::uInt nVisRow = rwVis.nplane();
  ASKAPDEBUGASSERT(nRow > 0);
  ASKAPCHECK(nVisRow % nRow == 0,"#visibility rows should be multiple of # rows in chunk");

  casa::uInt nDir = nVisRow / nRow;
  for (casa::uInt dir = 0; dir < nDir; ++dir) {
      for (casa::uInt row = 0; row < nRow; ++row) {
        const casa::uInt visRow =  row + dir * nRow;
        bool needMueller = true;
        bool validSolution = false;
        casa::Float det = 0.;
        for (casa::uInt chan = 0; chan < nChan; ++chan) {
            bool allFlagged = true;
            // we don't really support partial polarisation flagging, but to avoid nasty surprises it is better to flag such samples completely.
            bool needFlag = false;
            //ASKAPDEBUGASSERT(thisRowFlag.ncolumn() == nPol);
            for (casa::uInt pol = 0; pol < nPol; ++pol) {
                 if (flag(pol,chan,row)) {
                     needFlag = true;
                 } else {
                     allFlagged = false;
                 }
            }
            // don't bother with the rest of processing if the sample is flagged anyway in its entirety
            if (allFlagged) {
                continue;
            }

            if (!itsChannelIndependent || needMueller) {
                // we only need to fill the mueller matrix once in the non bandpass case
                needMueller = false;
                const int b1 = itsBeamIndependent ? 0 : (nDir > 1 ? dir : beam1[row]);
                const int b2 = itsBeamIndependent ? 0 : (nDir > 1 ? dir : beam2[row]);
                if (solutionAccessorNeedsTime) {
                    solutionAccessorNeedsTime = false;
                    updateAccessor(chunk.time());
                }
                std::pair<casa::SquareMatrix<casa::Complex, 2>, bool> jv1 =
                    calSolution().jonesAndValidity(antenna1[row], b1, chan);
                std::pair<casa::SquareMatrix<casa::Complex, 2>, bool> jv2 =
                    calSolution().jonesAndValidity(antenna2[row], b2, chan);

                validSolution = jv1.second && jv2.second;
                if (validSolution) {
                    if (getInterpolateTime()) {
                        // get second lot of calibration parameters
                        // a=after b=before 1=ant1 2=ant2
                        std::pair<casa::SquareMatrix<casa::Complex, 2>, bool> jv1a =
                            nextCalSolution().jonesAndValidity(antenna1[row], b1, chan);
                        std::pair<casa::SquareMatrix<casa::Complex, 2>, bool> jv2a =
                            nextCalSolution().jonesAndValidity(antenna2[row], b2, chan);
                        const bool valid = jv1a.second && jv2a.second;
                        if (valid) {
                            casa::SquareMatrix<casa::Complex, 2>& j1b = jv1.first;
                            casa::SquareMatrix<casa::Complex, 2>& j2b = jv2.first;
                            casa::SquareMatrix<casa::Complex, 2>& j1a = jv1a.first;
                            casa::SquareMatrix<casa::Complex, 2>& j2a = jv2a.first;
                            const double tb = calSolutionTime();
                            const double ta = nextCalSolutionTime();
                            if (ta > tb) {
                                const double t = chunk.time();
                                const float factor = (t - tb) / (ta - tb);
                                // simple interpolation scheme.. causes decorrelation
                                //j1b += (j1a - j1b) * factor;
                                //j2b += (j2a - j2b) * factor;
                                // Interpolate ampl and phase separately
                                // can only interpolate pure gain this way, not leakage
                                for (int i=0; i<=1; i++) {
                                    casa::Complex g1b = j1b(i,i);
                                    casa::Complex g2b = j2b(i,i);
                                    casa::Complex g1a = j1a(i,i);
                                    casa::Complex g2a = j2a(i,i);
                                    if (abs(g1a) > 0 && abs(g1b) > 0) {
                                        const casa::Complex g = g1a / g1b;
                                        const float mag = abs(g);
                                        j1b(i,i) = g1b * (1 + (mag-1)*factor) * pow(g/mag,factor);
                                    }
                                    if (abs(g2a) > 0 && abs(g2b) > 0) {
                                        const casa::Complex g = g2a / g2b;
                                        const float mag = abs(g);
                                        j2b(i,i) = g2b * (1 + (mag-1)*factor) * pow(g/mag,factor);
                                    }
                                }
                            }
                        }
                    }
                    const casa::SquareMatrix<casa::Complex, 2>& j1 = jv1.first;
                    const casa::SquareMatrix<casa::Complex, 2>& j2 = jv2.first;
                    const casa::Complex det1 = j1(0,0)*j1(1,1)-j1(0,1)*j1(1,0);
                    const casa::Complex det2 = j2(0,0)*j2(1,1)-j2(0,1)*j2(1,0);
                    det = casa::real(det1*conj(det1))*casa::real(det2*conj(det2));
                    // Inverse of 4x4 mueller is directProduct of 2x2 jones inverses
                    if (det > detThreshold) {
                        if (correct) {
                            directProduct(mueller, j1.inverse(),conj(j2).inverse());
                        } else {
                            directProduct(mueller, j1,conj(j2));
                        }
                    }
                }
            }
            if (validSolution) {
                for (casa::uInt pol = 0; pol < nPol; ++pol) {
                    //vis(pol) = rwVis(row,chan,pol);
                    vis(pol) = rwVis(pol,chan,visRow);
                }
            }

            if (itsFlagAllowed) {
                if (det <= detThreshold || !validSolution || needFlag) {
                    ASKAPCHECK(noiseAndFlagDA, "Accessor type passed to CalibrationApplicatorME does not support change of flags");
                    for (casa::uInt pol = 0; pol < nPol; ++pol) {
                        rwFlag(pol,chan,row)=true;
                        rwVis(pol,chan,visRow)=0.;
                    }
                    // go to next channel
                    continue;
                }
            } else {
              ASKAPCHECK(validSolution && !needFlag, "Encountered unflagged data and invalid solution, but flagging samples has not been allowed");
              ASKAPCHECK(det>detThreshold, "Unable to apply calibration for (antenna1,beam1)=("<<antenna1[row]<<","<<beam1[row]<<") and (antenna2,beam2)=("<<antenna2[row]<<
                               ","<<beam2[row]<<"), time="<<chunk.time()/86400.-55000<<" determinant is too close to 0. D="<<det<<" matrix="<<mueller.matrix()
                       //<<" jones1="<<jones1.matrix()<<" jones2="<<jones2.matrix()
                       <<" dir="<<askap::printDirection(chunk.pointingDir1()[row]));
            }

            // do the actual calibration - we have a valid invertable solution
            vis*=mueller;
            // write back to chunk
            for (casacore::uInt pol = 0; pol < nPol; ++pol) {
                rwVis(pol,chan,visRow) = vis(pol);
            }

            if (itsScaleNoise) {
                ASKAPCHECK(noiseAndFlagDA, "Accessor type passed to CalibrationApplicatorME does not support change of the noise estimate");
                //casacore::Vector<casacore::Complex> thisChanNoise = noiseAndFlagDA->rwNoise().yzPlane(row).row(chan);
                casacore::Vector<casacore::Complex> thisChanNoise = noiseAndFlagDA->rwNoise().xyPlane(row).column(chan);
                ASKAPDEBUGASSERT(thisChanNoise.nelements() == nPol);
                // propagating noise estimate through the matrix multiplication
                casa::RigidVector<casa::Complex,4> noise = thisChanNoise;
                // Somehow the indexing below calls the non const version and throws
                // an exception for a diagonal matrix. Avoid this with explicitly const version
                const casa::SquareMatrix<casa::Complex,4>& cmueller = mueller;

                for (casacore::uInt pol = 0; pol < nPol; ++pol) {
                    float tempRe = 0., tempIm = 0.;
                    for (casa::uInt k = 0; k < nPol; ++k) {
                        tempRe += casa::square(casa::real(cmueller(pol,k)) * casa::real(noise(k))) +
                                casa::square(casa::imag(cmueller(pol,k)) * casa::imag(noise(k)));
                        tempIm += casa::square(casa::real(cmueller(pol,k)) * casa::imag(noise(k))) +
                               casa::square(casa::imag(cmueller(pol,k)) * casa::real(noise(k)));
                    }
                    thisChanNoise(pol) = casacore::Complex(sqrt(tempRe),sqrt(tempIm));
                }
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
/// @brief determines whether channel dependent calibration is used or not
/// @details It is handy to be able to apply the same solution for all channels. With
/// this flag set, the same solution will be used for all channels.
/// @param[in] flag if true, the same calibration is applied to all channels
void CalibrationApplicatorME::channelIndependent(bool flag)
{
  itsChannelIndependent = flag;
  if (itsChannelIndependent) {
      ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will apply the same gain calibration solutions to all channels");
  } else {
      ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will apply bandpass calibration solutions");
  }
}

/// @brief determines whether time interpolation is used or not
/// @details when applying gains determined on a coarser time grid than the data
/// interpolating the gains in time improves dynamic range
/// @param[in] flag if true, interpolate the gains in time
void CalibrationApplicatorME::interpolateTime(bool flag)
{
  // make sure the parent CalibrationSolutionHandler knows about this too
  setInterpolateTime(flag);
  if (getInterpolateTime()) {
      ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will interpolate the gains in time");
  } else {
      ASKAPLOG_INFO_STR(logger, "CalibrationApplicatorME will not interpolate the gains in time");
  }
}

} // namespace synthesis

} // namespace askap
