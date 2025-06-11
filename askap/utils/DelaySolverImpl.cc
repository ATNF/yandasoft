/// @file
///
/// Actual algorithm implementing delay solver tool 
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
///

// a bit hacky way to get logs tagged with the cp-prefix
#define ASKAP_PACKAGE_NAME "cp"

#include <askap/utils/DelaySolverImpl.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapUtil.h>
ASKAP_LOGGER(logger, ".DelaySolverImpl");

#include <askap/askap/AskapError.h>
#include <casacore/casa/Arrays/MatrixMath.h>

// std
#include <set>
#include <fstream>
#include <iomanip>

using namespace askap;
using namespace askap::utils;

/// @brief constructor
/// @param[in] targetRes target spectral resolution in Hz, data are averaged to match the desired resolution 
/// note, integral number of channels are averaged.
/// @param[in] pol polarisation index to use
/// @param[in] ampCutoff if positive, amplitudes above ampCutoff will be flagged
/// @param[in] refAnt reference antenna index   
DelaySolverImpl::DelaySolverImpl(double targetRes, casa::Stokes::StokesTypes pol, float ampCutoff, casa::uInt refAnt) :
   itsTargetRes(targetRes), itsPol(pol), itsAmpCutoff(ampCutoff), itsRefAnt(refAnt), itsNAvg(0u), itsDelayEstimator(targetRes),
   itsChanToAverage(1u), itsVerbose(true), itsDataShapeWarningGiven(false), itsQualityThreshold(-1.)
{
  ASKAPCHECK(itsTargetRes > 0, "Target spectral resolution should be positive, you have "<<itsTargetRes<<" Hz");
} 
    
/// @brief set baselines to exclude
/// @details An empty vector configures the class to take all available baselines into account.
/// @param[in] bsln vector with pair of indices representing baselines to exclude from the solution
void DelaySolverImpl::excludeBaselines(const casa::Vector<std::pair<casa::uInt, casa::uInt> > &bsln) 
{
  itsExcludedBaselines.assign(bsln.copy());
}

/// @brief helper method to check that all channels/rows are flagged
/// @param[in] flags matrix with flags
/// @return true, if all channels and rows are flagged
bool DelaySolverImpl::checkAllFlagged(const casa::Matrix<bool> &flags)
{
   for (casa::uInt row=0; row<flags.nrow(); ++row) {
        for (casa::uInt chan=0; chan<flags.ncolumn(); ++chan) {
             if (!flags(row,chan)) {
                 return false;
             }
        }
   }
   return true;
}

/// @brief set verbose flag
/// @details Some important messages are set with high-severity (INFO and above), so they make it to the 
/// main system log given how we use this application now. This is the default behavior. However, setting false
/// to this flag allows to downgrade these messages to DEBUG severity and suppress these messages. We use this
/// in two-stage solution with initial estimation via lags. Otherwise, it is possible to have a duplication of
/// messages.
/// @param[in] flag true, to get messages published with high severity (default), false otherwise
void DelaySolverImpl::setVerboseFlag(bool flag)
{
   itsVerbose = flag;
}

/// @brief set optional quality threshold to fail the run if any unflagged baseline has a bad solution
/// @details LSQ fit is not very robust against data which have non-Gaussian errors (e.g. RFI), i.e.
/// problems at one baseline may affect the whole solution in a non-trivial way. Although there is no
/// good solution to this problem, in general, we have an option to abort the solution if any baseline
/// which is considered valid (i.e. is not flagged) returns solution of a bad quality. Negative value
/// for this threshold (default) essentially bypasses the check as the quality is from [0,1] range.
/// @param[in] threshold desired quality threshold 
/// @note The definition of quality is different for FFT and phase slope based algorithms.
void DelaySolverImpl::setQualityThreshold(double threshold)
{
   ASKAPCHECK(threshold <= 1., "Quality threshold greater than 1. will always cause a failure - such a useless setting of the parameter is unexpected");
   itsQualityThreshold = threshold;
}
     
/// @brief process one data accessor
/// @param[in] acc data accessor to process
void DelaySolverImpl::process(const accessors::IConstDataAccessor &acc)
{
   const casa::Vector<casa::Stokes::StokesTypes>& stokes = acc.stokes();
   casa::uInt pol2use = stokes.nelements();
   for (pol2use = 0; pol2use < stokes.nelements(); ++pol2use) {
        if (stokes[pol2use] == itsPol) {
            break;
        }
   }
   ASKAPCHECK(pol2use < stokes.nelements(), "Unable to find "<<casa::Stokes::name(itsPol)<<" polarisation product in the data");
   
   //const casa::Matrix<bool>& flags = acc.flag().xyPlane(pol2use);   
   const casa::Matrix<bool>& flags = acc.flag().yzPlane(pol2use);   
   if (checkAllFlagged(flags)) {
       return;
   }
   if (itsFreqAxis.nelements() == 0) {
       // this is the first time stamp, resize buffers, etc
       itsFreqAxis = acc.frequency();
       itsAnt1IDs = acc.antenna1();
       itsAnt2IDs = acc.antenna2();
       ASKAPDEBUGASSERT(itsAnt1IDs.nelements() == itsAnt2IDs.nelements());
       // determine the number of channels to average based on the target and actual resolution
       ASKAPCHECK(itsFreqAxis.nelements() > 1, "Need at least two spectral channels, you have "<<acc.nChannel());
       const double actualRes = (itsFreqAxis[itsFreqAxis.nelements() - 1] - itsFreqAxis[0]) / double(itsFreqAxis.nelements() - 1);
       ASKAPCHECK(abs(actualRes) > 0, "Unable to determine spectral resolution of the data");
       if (itsTargetRes > abs(actualRes)) {
           ASKAPDEBUGASSERT(itsTargetRes > 0);
           itsChanToAverage = casa::uInt(itsTargetRes / abs(actualRes));
       } 
       ASKAPLOG_DEBUG_STR(logger, "Averaging "<<itsChanToAverage<<" consecutive spectral channels");
       ASKAPDEBUGASSERT(itsChanToAverage > 0);
       itsDelayEstimator.setResolution(actualRes * itsChanToAverage);
       const casa::uInt targetNChan = acc.nChannel() / itsChanToAverage;
       ASKAPCHECK(targetNChan > 1, "Too few spectral channels remain after averaging: in="<<acc.nChannel()<<" out="<<targetNChan);
       itsSpcBuffer.resize(acc.nRow(), targetNChan);
       itsAvgCounts.resize(acc.nRow(), targetNChan);
       itsSpcBuffer.set(casa::Complex(0.,0.));
       itsAvgCounts.set(0u);
   } else {
       if (itsFreqAxis.nelements() != acc.nChannel()) {
           if (!itsDataShapeWarningGiven) {
               if (itsVerbose) {
                    ASKAPLOG_WARN_STR(logger, "The number of frequency channels has been changed, was "<<itsFreqAxis.nelements()<<" now "
                                      << acc.nChannel()<<", ignoring");
               } else {
                    ASKAPLOG_DEBUG_STR(logger, "The number of frequency channels has been changed, was "<<itsFreqAxis.nelements()<<" now "
                                      << acc.nChannel()<<", ignoring");
               }
           }
           itsDataShapeWarningGiven = true;
           return;
       }
       if (itsAnt1IDs.nelements() != acc.nRow()) {
           if (!itsDataShapeWarningGiven) {
               if (itsVerbose) {
                    ASKAPLOG_WARN_STR(logger, "The number of rows has been changed, was "<<itsAnt1IDs.nelements()<<" now "
                                      << acc.nRow()<<", ignoring");
               } else {
                    ASKAPLOG_DEBUG_STR(logger, "The number of rows has been changed, was "<<itsAnt1IDs.nelements()<<" now "
                                      << acc.nRow()<<", ignoring");
               }
               itsDataShapeWarningGiven = true;
           }
           return;
       }
       for (casa::uInt row = 0; row < acc.nRow(); ++row) {
            if (itsAnt1IDs[row] != acc.antenna1()[row]) {
                if (!itsDataShapeWarningGiven) {
                     if (itsVerbose) {
                         ASKAPLOG_WARN_STR(logger, "Antenna 1 index has been changed for row ="<<row<<", was "<<itsAnt1IDs[row]<<
                                           " now "<<acc.antenna1()[row]<<", ignoring");
                     } else {
                         ASKAPLOG_DEBUG_STR(logger, "Antenna 1 index has been changed for row ="<<row<<", was "<<itsAnt1IDs[row]<<
                                           " now "<<acc.antenna1()[row]<<", ignoring");
                     }
                     itsDataShapeWarningGiven = true;
                }
                return;                  
            }                     
            if (itsAnt2IDs[row] != acc.antenna2()[row]) {
                if (!itsDataShapeWarningGiven) {
                    if (itsVerbose) {
                        ASKAPLOG_WARN_STR(logger, "Antenna 2 index has been changed for row ="<<row<<", was "<<itsAnt2IDs[row]<<
                                          " now "<<acc.antenna2()[row]<<", ignoring");
                    } else {
                        ASKAPLOG_DEBUG_STR(logger, "Antenna 2 index has been changed for row ="<<row<<", was "<<itsAnt2IDs[row]<<
                                          " now "<<acc.antenna2()[row]<<", ignoring");
                    }
                    itsDataShapeWarningGiven = true;
                }
                return;                  
            }                     
       }
   }
   
   const casa::Matrix<casa::Complex> vis = acc.visibility().yzPlane(pol2use);
   for (casa::uInt row = 0; row < acc.nRow(); ++row) {
        const casa::Vector<casa::Complex> thisRowVis = vis.column(row);
        const casa::Vector<bool> thisRowFlags = flags.column(row);
        ASKAPDEBUGASSERT(itsSpcBuffer.nrow() == acc.nRow());
        ASKAPDEBUGASSERT(itsAvgCounts.nrow() == acc.nRow());        
        casa::Vector<casa::Complex> thisBufRowVis = itsSpcBuffer.row(row);
        casa::Vector<casa::uInt> thisRowCounts = itsAvgCounts.row(row);
        
        const double thisRowApproxDelay = delayApproximation(row);        
        
        ASKAPDEBUGASSERT(thisRowVis.nelements() == thisRowFlags.nelements());
        for (casa::uInt chan = 0, index = 0; chan < thisBufRowVis.nelements(); ++chan) {
             for (casa::uInt ch2 = 0; ch2 < itsChanToAverage; ++ch2,++index) {
                  ASKAPDEBUGASSERT(index < thisRowFlags.nelements());
                  
                  if (!thisRowFlags[index]) {
                      const casa::Complex curVis = thisRowVis[index];
                      if ((itsAmpCutoff < 0) || (abs(curVis) < itsAmpCutoff)) {
                          if (itsDelayApproximation.nelements() > 0) {
                              ASKAPDEBUGASSERT(index < itsFreqAxis.nelements());
                              // technically we don't need to subtract the start frequency, but it helps
                              // to avoid issues with subtracting two large numbers
                              const double phase = -casa::C::_2pi * (itsFreqAxis[index] - itsFreqAxis[0]) * thisRowApproxDelay;
                              const casa::Complex phasor(cos(phase),sin(phase));
                              thisBufRowVis[chan] += curVis * phasor;
                          } else {
                              thisBufRowVis[chan] += curVis;
                          }
                          ++thisRowCounts[chan];
                      }
                  }
             }   
        }
   }
   ++itsNAvg;
}

/// @brief helper method to obtain delay approximation for the given row
/// @details Zero is returned if itsDelayApproximation is empty. Otherwise,
/// delay for the given row is extracted (based on metadata contained in the buffers).
/// Note, delay units are the same as itsDelayApproximation (i.e. seconds).
/// @param[in] row row of interest
double DelaySolverImpl::delayApproximation(casa::uInt row) const
{
  if (itsDelayApproximation.nelements() == 0) {
      return 0.;
  }
  ASKAPDEBUGASSERT(row < itsAnt1IDs.nelements());
  ASKAPDEBUGASSERT(row < itsAnt2IDs.nelements());
  const casa::uInt ant1 = itsAnt1IDs[row];
  const casa::uInt ant2 = itsAnt2IDs[row];
  ASKAPCHECK(ant1 < itsDelayApproximation.nelements(), "Initial delay approximation not found for antenna "<<ant1);
  ASKAPCHECK(ant2 < itsDelayApproximation.nelements(), "Initial delay approximation not found for antenna "<<ant2);
  return itsDelayApproximation[ant1] - itsDelayApproximation[ant2];
}

/// @brief set target resolution
/// @details
/// @param[in] targetRes target spectral resolution in Hz, data are averaged to match the desired resolution 
/// note, integral number of channels are averaged.
void DelaySolverImpl::setTargetResolution(double targetRes)
{
  itsTargetRes = targetRes;
  ASKAPCHECK(itsTargetRes > 0, "Target spectral resolution should be positive, you have "<<itsTargetRes<<" Hz");
}

/// @brief set antenna names
/// @details If antenna names are set, these names will be used in the log messages related to the given antennas.
/// Otherwise (default), only antenna id will be shown.
/// @param[in] names vector with antenna names - index is id
/// @note It is allowed to have fewer named antennas than ids, if the array is too small antenna will be referred by id only.
void DelaySolverImpl::setAntennaNames(const std::vector<std::string> &names)
{
   itsAntennaNames = names;
}

/// @brief helper method to get string referring to antenna
/// @details Depending on whether the name is defined for the antenna with the 
/// given id, this method returns either "with id=..." or "name (id=..)". To be
/// used in log messages.
/// @param[in] id antenna id
/// @return string referring to the antenna with the given id
std::string DelaySolverImpl::antennaNameString(casa::uInt id) const
{
   if (id < itsAntennaNames.size()) {
       return itsAntennaNames[id] + " (id=" + utility::toString(id) + ")";
   } 
   return "with id=" + utility::toString(id);
}
    
/// @brief solve for antenna-based delays
/// @details This method estimates delays for all baselines and then solves for
/// antenna-based delays honouring baselines to be excluded.
/// @param[in] useFFT if true, FFT-based delay estimator is used. It is less accurate but
/// is more robust for large delays and less sensitive to flagged data.  
/// @return a vector with one delay per antenna (antennas are in the index-increasing order).
casa::Vector<double> DelaySolverImpl::solve(bool useFFT) const
{
  ASKAPCHECK(itsNAvg > 0, "No valid data found. At least one chunk of data have to be processed before delays can be estimated");
  ASKAPCHECK(itsFreqAxis.nelements() > 1, "Unable to estimate delays from monochromatic data");
  
  const casa::uInt nAnt = casa::max(casa::max(itsAnt1IDs), casa::max(itsAnt2IDs)) + 1;
  
  ASKAPLOG_DEBUG_STR(logger, "Using "<<itsNAvg<<" cycles to estimate delays for "<<nAnt<<" antennas; reference = "<<itsRefAnt);
  // build a set of baselines (rows) to exclude
  std::set<casa::uInt> rows2exclude;
  ASKAPDEBUGASSERT(itsAnt1IDs.nelements() == itsAnt2IDs.nelements());
  ASKAPDEBUGASSERT(itsAnt1IDs.nelements() == itsAvgCounts.nrow());
  for (casa::uInt row=0; row<itsAnt1IDs.nelements(); ++row) {
       // first, exclude explicitly listed baselines
       for (casa::uInt bsln = 0; bsln < itsExcludedBaselines.nelements(); ++bsln) {
            if ((itsExcludedBaselines[bsln].first == itsAnt1IDs[row]) && 
                (itsExcludedBaselines[bsln].second == itsAnt2IDs[row])) {
                 rows2exclude.insert(row);
            } 
       }
       // if this row is not excluded, do additional checks based on data
       if (rows2exclude.find(row) == rows2exclude.end()) {
           const casa::Vector<casa::uInt> thisRowCounts = itsAvgCounts.row(row);
           bool allFlagged = true;
           for (casa::uInt chan = 0; chan < thisRowCounts.nelements(); ++chan) {
                if (thisRowCounts[chan] > 0) {
                    allFlagged = false;
                    break;
                }
           }
           if (allFlagged) {
               rows2exclude.insert(row);
           }
       }
  }
  ASKAPCHECK(itsAnt1IDs.nelements() > rows2exclude.size(), "Looks like all data are flagged or excluded");
  ASKAPLOG_DEBUG_STR(logger, "Using "<<itsAnt1IDs.nelements() - rows2exclude.size()<<" rows(baselines) out of "<<
                             itsAnt1IDs.nelements()<<" available in the dataset");
  // build a list of excluded (flagged) antennas to ensure their
  // delays are set to zero
  std::set<casa::uInt> excludedAntennas;
  for (casa::uInt ant = 0; ant < nAnt; ++ant) {
       bool dataPresent = false;
       bool referencePresent = false;
       for (casa::uInt bsln=0; bsln<itsAnt1IDs.nelements(); ++bsln) {
            if (rows2exclude.find(bsln) == rows2exclude.end()) {
                if ((itsAnt1IDs[bsln] == ant) || (itsAnt2IDs[bsln] == ant)) {
                     dataPresent = true;
                     if ((itsAnt1IDs[bsln] == itsRefAnt) || (itsAnt2IDs[bsln] == itsRefAnt)) {
                         referencePresent = true;
                         break;
                     }
                }
            }
       }
       //ASKAPCHECK(dataPresent == referencePresent, "It looks like there are valid data for antenna "<<ant<<", but all baselines to reference="<<itsRefAnt<<" are flagged or missing.");
       if (dataPresent) {
           if (!referencePresent) {
               if (itsVerbose) {
                   ASKAPLOG_WARN_STR(logger, "Antenna "<<antennaNameString(ant)<<" has valid data, but not in baseline with the reference antenna "<<antennaNameString(itsRefAnt)<<", degeneracy possible");
               } else {
                   ASKAPLOG_DEBUG_STR(logger, "Antenna "<<antennaNameString(ant)<<" has valid data, but not in baseline with the reference antenna "<<antennaNameString(itsRefAnt)<<", degeneracy possible");
               }
           }
       } else {
           excludedAntennas.insert(ant);
       }
  }
  
  // build design equations
  ASKAPDEBUGASSERT(itsAnt1IDs.nelements() == itsSpcBuffer.nrow());
  casa::Vector<double> delays(itsSpcBuffer.nrow()+1,0.);
  casa::Vector<double> quality(itsSpcBuffer.nrow(),0.);
  casa::Matrix<double> dm(delays.nelements(),nAnt,0.);
  std::ofstream os("avgspectrum.dat");
  for (casa::uInt bsln = 0; bsln < itsSpcBuffer.nrow(); ++bsln) {
       if (rows2exclude.find(bsln) == rows2exclude.end()) {
           casa::Vector<casa::Complex> buf = itsSpcBuffer.row(bsln).copy();
           const casa::Vector<casa::uInt> thisRowCounts = itsAvgCounts.row(bsln);
           ASKAPDEBUGASSERT(buf.nelements() == thisRowCounts.nelements());           
           for (casa::uInt chan=0; chan < buf.nelements(); ++chan) {
                if (thisRowCounts[chan] > 0) {
                    buf[chan] /= float(thisRowCounts[chan]);
                } else {
                    if (chan > 0) {
                        // don't really have much better way to guess the phase of the flagged channel
                        // We could've interpolate, but unwrapping phase may be difficult, especially if
                        // more than one channel is flagged
                        buf[chan] = buf[chan - 1];
                    } else {
                        // can't really guess the first channel either
                        buf[chan] = 0;
                    }
                }
                os<<itsAnt1IDs[bsln]<<" "<<itsAnt2IDs[bsln]<<" "<<chan<<" "<<arg(buf[chan])/casa::C::pi*180.<<std::endl;
           }
           
           delays[bsln] = useFFT ? itsDelayEstimator.getDelayWithFFT(buf) : itsDelayEstimator.getDelay(buf);
           quality[bsln] = itsDelayEstimator.quality();
           ASKAPCHECK(quality[bsln] >= itsQualityThreshold, "Baseline between antenna "<<antennaNameString(itsAnt1IDs[bsln])<<
                      " and "<<antennaNameString(itsAnt2IDs[bsln])<<
                      " has bad solution quality ("<<quality[bsln]<<", threshold is "<<itsQualityThreshold<<
                      "), but is unflagged. Check that the system is working and not completely dazzled by RFI.");
           
           // now fill the design matrix
           const casa::uInt ant1 = itsAnt1IDs[bsln];
           ASKAPDEBUGASSERT(ant1 < dm.ncolumn()); 
           const casa::uInt ant2 = itsAnt2IDs[bsln];
           ASKAPDEBUGASSERT(ant2 < dm.ncolumn());
           ASKAPDEBUGASSERT(bsln < dm.nrow()); 
           if (ant1 != itsRefAnt) {               
               dm(bsln,ant1) = 1.;
           }
           if (ant2 != itsRefAnt) {
               dm(bsln,ant2) = -1.;
           }              
       }
  }
  ASKAPLOG_DEBUG_STR(logger, "Delays (ns) per baseline: "<<std::setprecision(9)<<delays*1e9);
  ASKAPLOG_DEBUG_STR(logger, "Quality of delay estimate: "<<std::setprecision(3)<<quality);
  // add conditions for flagged antennas to ensure zero delay
  if (excludedAntennas.size() == 0) {
      if (itsVerbose) {
          ASKAPLOG_INFO_STR(logger, "All available antennas have unflagged data");
      } else {
          ASKAPLOG_DEBUG_STR(logger, "All available antennas have unflagged data");
      }
  } else {

      ASKAPCHECK(rows2exclude.size() >= excludedAntennas.size(), "Number of excluded rows ("<<rows2exclude.size()<<
                 ") is less than the number of flagged antennas ("<<
                 excludedAntennas.size()<<") - this shouldn't happen");

      // excluded antennas are treated through additional equations in the design matrix. The required number of such eqautions is definitely less than or 
      // equal to the number of excluded baselines. So use the space available by excluded rows for these additional condition equations (the alternative would be
      // to append them and then either try to select the appropriate rows to avoid degeneracy or use SVD-based LSQ fit
      for (std::set<casa::uInt>::const_iterator rowIt = rows2exclude.begin(), antIt = excludedAntennas.begin(); antIt != excludedAntennas.end(); ++antIt) {
           if (itsVerbose) {
               ASKAPLOG_WARN_STR(logger, "Antenna "<<antennaNameString(*antIt)<<" has no valid data - result will have zero delay update");
           } else {
               ASKAPLOG_DEBUG_STR(logger, "Antenna "<<antennaNameString(*antIt)<<" has no valid data - result will have zero delay update");
           }
           ASKAPDEBUGASSERT(*rowIt < dm.nrow());
           ASKAPDEBUGASSERT(*rowIt < delays.nelements());
           ASKAPDEBUGASSERT(*antIt < dm.ncolumn());
           ASKAPDEBUGASSERT(rowIt != rows2exclude.end());
           dm(*rowIt, *antIt) = 1.; 
           delays[*rowIt] = 0.;
           ++rowIt;
      }
  }

  // condition for the reference antenna (zero ref. delay is set in the last element of delays)
  ASKAPCHECK(itsRefAnt < nAnt, "Reference antenna is not present");
  dm(itsSpcBuffer.nrow(),itsRefAnt) = 1.;
  
  // just do an explicit LSQ fit. We could've used SVD invert here.
  const casa::Matrix<double> dmt = transpose(dm);
  const casa::Matrix<double> nm = product(dmt,dm);
  
  casa::Vector<double> result = product(invert(nm), product(dmt,delays));
  
  if (itsDelayApproximation.nelements() > 0) {
      ASKAPCHECK(itsDelayApproximation.nelements() == result.nelements(), "Delay approximations should be given for all antennas. nAnt="<<nAnt);
      result += itsDelayApproximation;
  }
  return result;  
}

/// @brief set initial delay approximation
/// @details The class can optionally remove some a priori known (approximate) delay in 
/// full resolution data before averaging takes place. This allows to get a coarse delay 
/// first in full resolution data and then refine it with
/// averaging. Empty array means zero initial delay for all antennas (i.e. nothing special is done,
/// this is the default). There should be one value per antenna. An exception is thrown if
/// antenna ID greater than or equal to the size of the vector is encountered
/// @param[in] delays
void DelaySolverImpl::setApproximateDelays(const casa::Vector<double> &delays)
{
  if (delays.nelements()) {
      ASKAPLOG_DEBUG_STR(logger, "The following approximate delays (ns) will be removed before averaging "<<delays * 1e9);
  }
  itsDelayApproximation.assign(delays.copy());
}
    
/// @brief initialise the accumulation
/// @details This method reverts the state of the object to that before the first call to process method.
/// This allows to repeat delay estimation with a set of approximate delays.
void DelaySolverImpl::init()
{
  itsFreqAxis.resize(0);
  itsAnt1IDs.resize(0);
  itsAnt2IDs.resize(0);
  itsNAvg = 0;  
  itsDataShapeWarningGiven = false;
}

