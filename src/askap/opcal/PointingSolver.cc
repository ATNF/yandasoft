/// @file
///
/// @brief Solver for pointing corrections
/// @details This is one of the examples of high-level solvers which 
/// refines antenna pointing. The idea behind the algorithm is similar
/// to the pointing corrections done by online ATCA software. 
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>

// local includes
#include "opcal/PointingSolver.h"
#include "measurementequation/SynthesisParamsHelper.h"

// casa includes
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/scimath/Mathematics/MatrixMathLA.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/measures/Measures/ParAngleMachine.h>

ASKAP_LOGGER(logger, ".PointingSolver");

namespace askap {

namespace synthesis {

/// @brief constructor
/// @param[in] parset parameters of the solver
PointingSolver::PointingSolver(const LOFAR::ParameterSet& parset) : 
   itsMaxPatternSize(SynthesisParamsHelper::convertQuantity(parset.getString("maxpatternsize", "5deg"),"rad"))
{
  ASKAPLOG_INFO_STR(logger, "Maximum angular size of the pointing pattern is assumed to be "<<itsMaxPatternSize / casacore::C::pi * 180.<<" deg");
}    
    
/// @brief method to be called when the results are ready
/// @details This method is called (in derived classes) when the results
/// are ready. Two parameters describe scans and contain calibration data, respectively.
/// @param[in] scans description of scans, note separate beams are present as separate scans.
///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
///            not necessarily match the scans known to online system
/// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
///            (column index is antenna ID used in the measurement set).
void PointingSolver::process(const ScanStats &scans, const casacore::Matrix<GenericCalInfo> &caldata)
{
   // need a better way to separate scans, but for now assume a fixed order:
   // +Dec, -Dec, boresight, +RA, -RA, boresight_long (the last one to be skipped)
   ASKAPASSERT(scans.size() == caldata.nrow());
   const size_t nExtraScans  = 0;
   const size_t nScansPerPattern  = 5 + nExtraScans;
   ASKAPCHECK(scans.size() >= nScansPerPattern, "Require at least "<<nScansPerPattern<<" scans, you have "<<scans.size());
   
   std::ofstream os("result.dat");
   
   for (size_t startScan = 0, pattern = 0; startScan + nScansPerPattern <= scans.size(); startScan += nScansPerPattern,++pattern) {
        const casacore::MVDirection& phaseCentre = scans[startScan].direction();
        const double time = 0.5*(scans[startScan + nScansPerPattern/2].startTime()+
                                 scans[startScan + nScansPerPattern/2].endTime());
        const casacore::MVDirection azEl = mroAzEl(phaseCentre, time);
        ASKAPLOG_INFO_STR(logger, "Pointing pattern "<<pattern + 1<<" at "<<printDirection(phaseCentre)<<" J2000, Az="<<
                          azEl.getLong()/casacore::C::pi*180.<<" deg, El="<<
                          azEl.getLat()/casacore::C::pi*180.<<" deg:");
        os<<azEl.getLong()/casacore::C::pi*180.<<" "<<azEl.getLat()/casacore::C::pi*180.<<" ";                  
        // consistency check 
        ASKAPDEBUGASSERT(scans.size() > startScan + 4);
        const double totalDuration = scans[startScan + 4].endTime() - scans[startScan].startTime(); 
        ASKAPCHECK(totalDuration<1800., "Pattern seems to be too long: "<<totalDuration / 60.<<" min");
        for (size_t index = startScan + 1; index < startScan + nScansPerPattern; ++index) {
             ASKAPCHECK(scans[index].direction().separation(phaseCentre)<1e-6, "Phase centres seem to vary throughout the pattern");
        }
                 
        for (casacore::uInt ant=0; ant<caldata.ncolumn(); ++ant) {
             casacore::Vector<GenericCalInfo> thisPatternData = caldata(casacore::Slice(startScan, 5),
                     casacore::Slice(ant));
             ASKAPASSERT(thisPatternData.nelements() + nExtraScans == nScansPerPattern);  
             // get result in ra and dec      
             const std::pair<double, double> resultEq = solveOne(thisPatternData);
             
             /*
             // for offsets in ra,dec
             casacore::MVDirection testDir = phaseCentre;
             testDir.shift(resultEq.first, resultEq.second,true);
             
             const casacore::MVDirection azElOff = mroAzEl(testDir, time);
             const double azOff = sin(azElOff.getLong() - azEl.getLong()) * cos(azElOff.getLat()) / casacore::C::pi * 180.;
             const double elOff = (sin(azElOff.getLat())*cos(azEl.getLat()) - cos(azElOff.getLat())*sin(azEl.getLat()) *
                                   cos(azElOff.getLong() - azEl.getLong())) / casacore::C::pi * 180.;
             */
             const double azOff = resultEq.first / casacore::C::pi * 180.;
             const double elOff = resultEq.second / casacore::C::pi * 180.;
             
             os<<azOff<<" "<<elOff<<" ";             
             ASKAPLOG_INFO_STR(logger, " antenna "<<ant<<":  az="<<azOff<<", el="<<elOff);
        }
        os<<std::endl;
   }
}

/// @brief solve for corrections for just one antenna
/// @param[in] caldata vector of data for all points of the pattern
/// @return pair with offsets in two coordinates
std::pair<double, double> PointingSolver::solveOne(const casacore::Vector<GenericCalInfo> &caldata) const
{
  ASKAPASSERT(caldata.nelements() == 5);
  casacore::Vector<double> amplitudes(caldata.nelements());
  for (casacore::uInt point=0; point<caldata.nelements(); ++point) {
       ASKAPASSERT(caldata[point].gainDefined());
       amplitudes[point] = casacore::abs(caldata[point].gain());
  }
  // normalise
  amplitudes /= casacore::max(amplitudes);
  //ASKAPLOG_INFO_STR(logger, "Amplitudes: "<<amplitudes);
  bool peakAtCentre = true;
  for (casacore::uInt point = 0; point<amplitudes.nelements(); ++point) {
       if (point != 2) {
           peakAtCentre &= (amplitudes[point] < 1.);
       }
  }
  if (!peakAtCentre) {
      ASKAPLOG_WARN_STR(logger, "The middle point of the pointing pattern is expected to be closest to boresight");
      ASKAPLOG_WARN_STR(logger, "   Normalised amplitudes: "<<amplitudes);
  }
  
  // need a better way to separate scans, but for now assume a fixed order
  const double fwhm = 30./0.8635/1200.*1.02; // need to pass this as a parameter (and extract from data)
  const double fractBandwidth = 0.304/0.8635;
  // compute differences w.r.t. predicted gain
  /*
  // order:
  // +Dec, -Dec, boresight, +RA, -RA
  const double decshift = 1./180.*casacore::C::pi; // one degree
  const double rashift = decshift * cos(12.39111 / 180. * casacore::C::pi); // we just did +/- 4 min of RA
  const double raOffsets[5] = {0.,0.,0., +rashift, -rashift};
  const double decOffsets[5] = {+decshift, -decshift, 0., 0., 0.};
  */
  // order:
  // for -az,+az,boresight,-el,+el
  const double shift = fwhm / 1.02 / 2.;
  const double raOffsets[5] = {-shift,+shift,0., 0.,0.};
  const double decOffsets[5] = {0., 0., 0.,-shift, +shift};
  
  // make least-square fit
  std::pair<double,double> result(0.,0.);
  const size_t nIter = 10;
  for (size_t iter=0; iter<nIter; ++iter) {
       casacore::Matrix<double> nm(2,2,0.); // normal matrix
       casacore::Vector<double> data(2.,0.); // projected right-hand side
       for (casacore::uInt point=0; point<amplitudes.nelements(); ++point) {
            const double off1 = raOffsets[point] + result.first;
            const double off2 = decOffsets[point] + result.second;
            const double offsetSq = casacore::square(off1) + casacore::square(off2);
            const double coeff = -4. * log(2.) / casacore::square(fwhm);
            const double expectedGain = 0.5*(exp(coeff * offsetSq * casacore::square(1. - fractBandwidth / 2.)) +
                                        exp(coeff * offsetSq * casacore::square(1. + fractBandwidth / 2.)));
       
            // derivatives of the gain
            const double dg_dx =  off1 * coeff * (casacore::square(1. - fractBandwidth / 2.) * exp(coeff * offsetSq * casacore::square(1. - fractBandwidth / 2.)) +
                                        casacore::square(1. + fractBandwidth / 2.) * exp(coeff * offsetSq * casacore::square(1. + fractBandwidth / 2.)));
            const double dg_dy =  off2 * coeff * (casacore::square(1. - fractBandwidth / 2.) * exp(coeff * offsetSq * casacore::square(1. - fractBandwidth / 2.)) +
                                        casacore::square(1. + fractBandwidth / 2.) * exp(coeff * offsetSq * casacore::square(1. + fractBandwidth / 2.)));
            const double measuredDifference = amplitudes[point] - expectedGain;
            nm(0,0) += casacore::square(dg_dx);
            nm(0,1) += dg_dx * dg_dy;
            nm(1,0) += dg_dy * dg_dx;
            nm(1,1) += casacore::square(dg_dy);
            data[0] += measuredDifference * dg_dx;
            data[1] += measuredDifference * dg_dy;                                               
       }    
       ASKAPCHECK(casacore::determinate(nm) > 1e-7, "Inversion failed; nm="<<nm);
       casacore::Vector<double> corrections = casacore::product(casacore::invert(nm),data);
       ASKAPASSERT(corrections.nelements() == 2); 
       result.first += corrections[0];
       result.second += corrections[1];
       //std::cout<<"Iter="<<iter<<" "<<corrections<<std::endl;
       if (iter + 1 == nIter) {
           const double misFit = sqrt(casacore::square(corrections[0])+casacore::square(corrections[1]));
           if (misFit > 1e-7) {
               ASKAPLOG_WARN_STR(logger, "LSF failed to converge after "<<nIter<<" iterations, misFit="<<misFit);
           }
       }
  }
  return result;
}


/// @brief obtain MRO reference position
/// @return position measure
casacore::MPosition PointingSolver::mroPosition()
{
   casacore::MPosition mroPos(casacore::MVPosition(casacore::Quantity(370.81, "m"), casacore::Quantity(116.6310372795, "deg"), 
                          casacore::Quantity(-26.6991531922, "deg")), casacore::MPosition::Ref(casacore::MPosition::WGS84));
   return mroPos;
}

/// @brief obtain Az and El at MRO for the given direction
/// @param[in] dir direction of interest in J2000
/// @param[in] time time in seconds since 0 MJD
/// @return direction in AzEl
casacore::MVDirection PointingSolver::mroAzEl(const casacore::MVDirection &dir, double time)
{
  const casacore::MEpoch epoch(casacore::Quantity(time/86400.,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));
  casacore::MeasFrame frame(mroPosition(), epoch);    
  return casacore::MDirection::Convert(casacore::MDirection(dir, casacore::MDirection::J2000), 
         casacore::MDirection::Ref(casacore::MDirection::AZEL,frame))().getValue();
}

/// @brief obtain parallactic angle at MRO for the given direction
/// @param[in] dir direction of interest in J2000
/// @param[in] time time in seconds since 0 MJD
/// @return parallactic angle in radians
double PointingSolver::mroPA(const casacore::MVDirection &dir, double time)
{
  const casacore::MEpoch epoch(casacore::Quantity(time/86400.,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));
  const casacore::MPosition mroPos = mroPosition();
  const casacore::Vector<double> mroPosRad = mroPos.getAngle("rad").getValue();
  ASKAPASSERT(mroPosRad.nelements() == 2);
  casacore::MeasFrame frame(mroPos, epoch);    
  /*
  casacore::ParAngleMachine pam(casacore::MDirection(dir,casacore::MDirection::Ref(casacore::MDirection::J2000)));
  pam.set(frame);
  return pam(epoch).getValue("rad");
  */
  const casacore::MVDirection hadec = casacore::MDirection::Convert(casacore::MDirection(dir, casacore::MDirection::J2000), 
         casacore::MDirection::Ref(casacore::MDirection::HADEC,frame))().getValue();
  const double sq = cos(mroPosRad[1])*sin(hadec.getLong());
  const double cq = sin(mroPosRad[1])*cos(hadec.getLat()) - cos(mroPosRad[1])*sin(hadec.getLat())*cos(hadec.getLong());
  return atan2(sq,cq);  
}



} // namespace synthesis

} // namespace askap
