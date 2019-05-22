/// @file
///
/// @brief Solver for array baselines
/// @details This is one of the examples of high-level solvers which 
/// refines antenna positions on the ground. Unlike generic algorithms
/// which just throw measurements to a least square solver, this one
/// attempts to unwrap the phase to make it more robust against large
/// errors in antenna positions (or timing). 
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
#include <askap/opcal/BaselineSolver.h>
#include <askap/IndexedCompare.h>
#include <askap/opcal/ObservationDescription.h>
#include <askap/scimath/utils/PhaseUnwrapper.h>

// casa includes
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/scimath/Functionals/Function1D.h>
#include <casacore/scimath/Fitting/NonLinearFitLM.h>

#include <casacore/measures/Measures/MCEpoch.h>


// std
#include <algorithm>
#include <fstream>
#include <vector>
#include <iomanip>

ASKAP_LOGGER(logger, ".BaselineSolver");


namespace askap {

namespace synthesis {

// helper comparator class ordering scans by hour angle
struct LesserHAorDec {
 
   /// @brief constructor, sets up position on Earth to do the calculations for
   /// @details
   /// @param[in] pos position on the ground
   /// @param[in] workWithHA if true, hour angle is compared
   LesserHAorDec(const casacore::MPosition &pos, bool workWithHA) : itsPosition(pos), itsWorkWithHA(workWithHA) {}
   

   /// @brief compares two directions 
   /// @details We want to order scans in hour angle. This method 
   /// checks two positions on the sky and determines which one
   /// has a lesser hour angle. 
   /// @param[in] scan1 first direction
   /// @param[in] scan2 second direction
   /// @return true, if the first scan corresponds to a lesser hour angle
   /// @note Time is assumed to be in seconds since 0 MJD, directions are assumed to be in J2000
   bool operator()(const ObservationDescription &scan1, const ObservationDescription &scan2) const;
   
private:
   /// @brief position on the ground to do the calculations for
   casacore::MPosition itsPosition;
   
   /// @brief true, if hour angle is compared; false, if declination is compared
   bool itsWorkWithHA;
};

/// @brief compares two directions 
/// @details We want to order scans in hour angle. This method 
/// checks two positions on the sky and determines which one
/// has a lesser hour angle. 
/// @param[in] scan1 first direction
/// @param[in] scan2 second direction
/// @return true, if the first scan corresponds to a lesser hour angle
/// @note Time is assumed to be in seconds since 0 MJD, directions are assumed to be in J2000
bool LesserHAorDec::operator()(const ObservationDescription &scan1, const ObservationDescription &scan2) const
{
  const double time1 = 0.5*(scan1.startTime() + scan1.endTime());
  const casacore::MEpoch epoch1(casacore::Quantity(time1/86400.,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));
  casacore::MeasFrame frame1(itsPosition, epoch1);    
  casacore::MVDirection hadec1 = casacore::MDirection::Convert(casacore::MDirection(scan1.direction(),casacore::MDirection::J2000), 
         casacore::MDirection::Ref(casacore::MDirection::HADEC,frame1))().getValue();

  const double time2 = 0.5*(scan2.startTime() + scan2.endTime());
  const casacore::MEpoch epoch2(casacore::Quantity(time2/86400.,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));
  casacore::MeasFrame frame2(itsPosition, epoch2);    
  casacore::MVDirection hadec2 = casacore::MDirection::Convert(casacore::MDirection(scan2.direction(),casacore::MDirection::J2000), 
         casacore::MDirection::Ref(casacore::MDirection::HADEC,frame2))().getValue();

  if (itsWorkWithHA) {
      return hadec1.getLong() < hadec2.getLong();
  } 
  return hadec1.getLat() < hadec2.getLat();
}

////////////////////

// helper functional fitting sinusoid due to errors in X and Y component of the baseline
// to be used with the casa fitter. We could've used our own, but it doesn't solve for errors
template<typename T>
class SinusoidDueToXY : public casacore::Function1D<T> {
public:
   
   /// @brief evaluate - main method
   virtual T eval(casacore::Function<double>::FunctionArg x) const {
       return this->param_p[0]*cos(x[0]+this->param_p[1])+this->param_p[2]; 
   };
   SinusoidDueToXY() : casacore::Function1D<T>(3u) {}  

   template<typename W>
   SinusoidDueToXY(const SinusoidDueToXY<W> &other) : casacore::Function1D<T>(other) {}

   // Copy it - compulsory methods required by casacore, we never use them explicitly
   virtual casacore::Function1D<T>* clone() const { return new SinusoidDueToXY<T>(*this); };   
   virtual casacore::Function1D<typename casacore::FunctionTraits<T>::DiffType>* cloneAD() const { return new SinusoidDueToXY<typename casacore::FunctionTraits<T>::DiffType>(*this); };   
   virtual casacore::Function1D<typename casacore::FunctionTraits<T>::BaseType>* cloneNonAD() const { return new SinusoidDueToXY<typename casacore::FunctionTraits<T>::BaseType>(*this); };   
};

//////////////////


/// @brief constructor
/// @param[in] parset parameters of the solver
BaselineSolver::BaselineSolver(const LOFAR::ParameterSet& parset) : 
  itsSolveXY(parset.getBool("solvexy", true)), itsSolveZ(parset.getBool("solvez",false)),
  itsRefAnt(parset.getUint("refant",1)), itsOrigFCMLayout(parset.getString("origlayout",""))
{
   ASKAPLOG_INFO_STR(logger, "Initialising antenna location solver");
   ASKAPCHECK(itsSolveXY || itsSolveZ, "Either dX,dY or dZ, or both should be chosen to solve for");
   if (itsSolveXY && !itsSolveZ) {
       ASKAPLOG_INFO_STR(logger, "  - will solve for dX and dY");
   } else if (itsSolveXY && itsSolveZ) {
       ASKAPLOG_INFO_STR(logger, "  - will solve for dX, dY and dZ");
   } else {
       ASKAPLOG_INFO_STR(logger, "  - will solve for dZ only");
   }
}
    
/// @brief method to be called when the results are ready
/// @details This method is called (in derived classes) when the results
/// are ready. Two parameters describe scans and contain calibration data, respectively.
/// @param[in] scans description of scans, note separate beams are present as separate scans.
///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
///            not necessarily match the scans known to online system
/// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
///            (column index is antenna ID used in the measurement set).
void BaselineSolver::process(const ScanStats &scans, const casacore::Matrix<GenericCalInfo> &caldata)
{
  ASKAPASSERT(scans.size() == caldata.nrow());
  itsCorrections.resize(caldata.ncolumn(), 3);
  itsCorrections.set(0.);
  itsErrors.resize(caldata.ncolumn(), 3);
  itsErrors.set(0.);
  
  if (itsSolveXY) {
      solveForXY(scans,caldata);
  }
  // doing Z after XY allows to use the fit results to achieve a greater phase stability and, potentially 
  // a better linear fit
  if (itsSolveZ) {
      solveForZ(scans,caldata);
  } 
  
  // writing results to the log
  for (casacore::uInt ant=0; ant<itsCorrections.nrow(); ++ant) {
       if (ant != itsRefAnt) {
           ASKAPLOG_INFO_STR(logger, "Antenna "<<ant);
           ASKAPLOG_INFO_STR(logger, "dX: "<<itsCorrections(ant,0)<<" +/- "<<itsErrors(ant,0)<<" metres");
           ASKAPLOG_INFO_STR(logger, "dY: "<<itsCorrections(ant,1)<<" +/- "<<itsErrors(ant,1)<<" metres");
           if (itsSolveZ) {
               ASKAPLOG_INFO_STR(logger, "dZ: "<<itsCorrections(ant,2)<<" +/- "<<itsErrors(ant,2)<<" metres");
           }
       }       
  }
  // update positions if necessary
  if (itsOrigFCMLayout != "") {
      ASKAPLOG_INFO_STR(logger, "Reading FCM configuration from "<<itsOrigFCMLayout<<", writing corrections.dat");
      LOFAR::ParameterSet fcmParset(itsOrigFCMLayout);
      std::vector<std::string> idMap = fcmParset.getStringVector("baselinemap.antennaidx");
      std::map<std::string, casacore::uInt> antmap;
      // building antenna map
      for (casacore::uInt ant = 0; ant < idMap.size(); ++ant) {
           const std::string idStr = idMap[ant];
           ASKAPCHECK(idStr.size() >= 3, "Antenna names in the baselinemap.antennaidx are supposed to be at least 3 symbols long, you have "<<idStr);
           ASKAPCHECK(idStr.find("ak") == 0, "All names in the baselinemap.antennaidx are supposed to start with ak, you have "<<idStr);
           const std::string antkey = "ant" + utility::toString<casacore::uInt>(utility::fromString<casacore::uInt>(idStr.substr(2)));
           antmap[antkey] = ant;
      }
      std::vector<std::string> antennas = fcmParset.getStringVector("antennas");
      ASKAPCHECK(itsCorrections.nrow() <= antennas.size(), "Number of antennas defined in FCM is smaller than the number of antennas solved for");
      std::ofstream os("corrections.dat");
      for (std::vector<std::string>::const_iterator ci=antennas.begin(); ci != antennas.end(); ++ci) {
           const std::string parsetKey = "antenna."+*ci+".location.itrf";
           const casacore::Vector<double> oldXYZ = fcmParset.isDefined(parsetKey) ? 
                        fcmParset.getDoubleVector(parsetKey) : fcmParset.getDoubleVector("common."+parsetKey);
           ASKAPCHECK(oldXYZ.nelements() == 3, "Expect exactly 3 elements for antenna.ant??.location.itrf key");
           std::map<std::string, casacore::uInt>::const_iterator mapIt = antmap.find(*ci);
           if (mapIt == antmap.end()) { 
               continue; // unused antenna
           }
           const casacore::uInt ant = mapIt->second;
           ASKAPLOG_INFO_STR(logger, "Antenna "<<ant<<" is "<<*ci);
           if (ant == itsRefAnt) {
               continue;
           }
           os<<"common."<<parsetKey<<" = [";
           for (casacore::uInt elem=0; elem<itsCorrections.ncolumn(); ++elem) {
                ASKAPASSERT(elem < oldXYZ.nelements());
                ASKAPDEBUGASSERT(ant < itsCorrections.nrow());
                // corrections are determined deviations from the best XYZ, need to subtract them to apply
                os<<std::setprecision(15)<<oldXYZ[elem] - itsCorrections(ant,elem);
                if (elem == 2) {
                    os << "]"<<std::endl;
                } else {
                    os <<", ";
                }
           }
      }
  }
}
 

/// @brief solve for dX and dY and populate corrections
/// @param[in] scans description of scans, note separate beams are present as separate scans.
///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
///            not necessarily match the scans known to online system
/// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
///            (column index is antenna ID used in the measurement set).
void BaselineSolver::solveForXY(const ScanStats &scans, const casacore::Matrix<GenericCalInfo> &caldata)
{
  ASKAPASSERT(scans.size()>1);
  const casacore::MPosition mroPos = mroPosition();
  std::vector<size_t> scanIndices(scans.size());
  for (size_t scan=0; scan<scanIndices.size(); ++scan) {
       scanIndices[scan]=scan;
  }
  std::sort(scanIndices.begin(),scanIndices.end(), utility::indexedCompare<size_t>(scans.begin(), LesserHAorDec(mroPos,true)));

  ASKAPDEBUGASSERT(scanIndices.size() == scans.size());
  for (casacore::uInt ant=0; ant < caldata.ncolumn(); ++ant) {
       // figure out the number of valid points
       size_t nGoodPoints = 0;
       for (casacore::uInt cnt = 0; cnt < scans.size(); ++cnt) {
            if (caldata(scanIndices[cnt],ant).gainDefined()) {
                ++nGoodPoints;
            }
       }
       const std::string spcfilename = "phase.xy.ant"+utility::toString(ant)+".dat";
       std::ofstream os(spcfilename.c_str());
       if (nGoodPoints == 0) {
           ASKAPLOG_WARN_STR(logger, "No valid data points for antenna "<<ant);
           // just to avoid accidental reuse of an old file
           os<<"No valid data"<<std::endl;
       }
       if (nGoodPoints < scans.size()) {
           ASKAPLOG_INFO_STR(logger, "Only "<<nGoodPoints<<" valid data points out of "<<scans.size()<<
                             " possible for antenna "<<ant);
       }

       // use casa fitter
       casacore::Vector<double> hangles(nGoodPoints);
       casacore::Vector<double> phases(nGoodPoints);
       casacore::Vector<double> sigma(nGoodPoints,1.);

       casacore::NonLinearFitLM<double> fitter;
       SinusoidDueToXY<double> func;
       func.parameters()[0] = 1.;
       func.parameters()[1] = 0.;
       func.parameters()[2] = 0.;
       fitter.setFunction(func);
       fitter.setMaxIter(180);
       fitter.setCriteria(0.00005);
  
       // could've cached the hour angles, but for now leave as is
       scimath::PhaseUnwrapper<double> unwrapper;
       for (casacore::uInt cnt = 0, pointCnt = 0; cnt < scans.size(); ++cnt) {
            const GenericCalInfo& cInfo = caldata(scanIndices[cnt],ant);
            if (!cInfo.gainDefined()) {
                continue;
            }
            const ObservationDescription& scan = scans[scanIndices[cnt]];
            const double time = 0.5*(scan.startTime() + scan.endTime());
            const casacore::MEpoch epoch(casacore::Quantity(time/86400.,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));
            casacore::MeasFrame frame(mroPos, epoch);    
            
            const casacore::MVDirection hadec = casacore::MDirection::Convert(casacore::MDirection(scan.direction(),casacore::MDirection::J2000), 
                                   casacore::MDirection::Ref(casacore::MDirection::HADEC,frame))().getValue();
            
            ASKAPASSERT(pointCnt < nGoodPoints);
            hangles[pointCnt] = hadec.getLong() - mroPos.getValue().getLong(); // Hour angle at latitude 0
            
            /*
            // Phase centre in the apparent topocentric frame
            casacore::MVDirection fpc = casacore::MDirection::Convert(casacore::MDirection(scan.direction(),casacore::MDirection::J2000), 
                                   casacore::MDirection::Ref(casacore::MDirection::TOPO,frame))().getValue();
            const double gastDayFrac = casacore::MEpoch::Convert(epoch,casacore::MEpoch::Ref(casacore::MEpoch::GAST))().get("d").getValue("d");
            const double gast = (gastDayFrac - casacore::Int(gastDayFrac)) * casacore::C::_2pi; // in radians
            hangles[pointCnt] = gast - fpc.getLong();
            // to ensure corresponding local hour angle is contiguous without the need to unwrap
            if (hangles[pointCnt] + mroPos.getValue().getLong() > casacore::C::pi) {
                hangles[pointCnt] -= 2.*casacore::C::pi;
            }
            if (hangles[pointCnt] + mroPos.getValue().getLong() < -casacore::C::pi) {
                hangles[pointCnt] += 2.*casacore::C::pi;
            }
            */
            /*
            ASKAPLOG_DEBUG_STR(logger, "cnt = "<<cnt<<" newH: "<<hangles[pointCnt]*180./casacore::C::pi<<
                     " oldH:"<<(hadec.getLong() - mroPos.getValue().getLong())*180./casacore::C::pi<<" diff: "<<
                      (hangles[pointCnt] - hadec.getLong() + mroPos.getValue().getLong()) / casacore::C::pi * 648000.<<" "
                     <<" decDiff: "<<(fpc.getLat() - hadec.getLat()) / casacore::C::pi * 648000.);
            */

            const double cd = cos(hadec.getLat()); 
            ASKAPCHECK(cd > 0, "Cannot work with sources at either pole");
            const casacore::MVDirection azel = casacore::MDirection::Convert(casacore::MDirection(scan.direction(),casacore::MDirection::J2000), 
                                   casacore::MDirection::Ref(casacore::MDirection::AZEL,frame))().getValue();
            phases[pointCnt] = unwrapper(arg(cInfo.gain())) / cd;
            os<<cnt<<" "<<hangles[pointCnt] / casacore::C::pi * 180<<" "<<phases[pointCnt] / casacore::C::pi * 180.<<" "<<scan.scanID()<<" "<<scan.fieldID()<<" "<<azel.getLat() / casacore::C::pi * 180.<<" "<<azel.getLong() / casacore::C::pi * 180.<<std::endl;
            ++pointCnt;
       }
       casacore::Vector<double> param = fitter.fit(hangles,phases,sigma);
       ASKAPCHECK(param.nelements() == 3, "Expect 3 parameters out of the fitter, you have size="<<param.nelements());
       casacore::Vector<double> err = fitter.errors();
       ASKAPCHECK(err.nelements() == 3, "Expect 3 uncertainties out of the fitter, you have size="<<err.nelements());
       // for our test setup of 864.5 MHz central freq
       const double wavelength = casacore::C::c / 864.5e6; // effective wavelength in metres (to do: get it from scan's frequency)
       double ampl = param[0] / 2. / casacore::C::pi * wavelength;
       // fit can converge with either sign of the first coefficient, but we like to always have a positive amplitude
       if (ampl < 0) {
           ampl = -ampl;
           param[0] = -param[0];
           param[1] += casacore::C::pi;
       }
       const double amplErr = err[0] / 2. / casacore::C::pi * wavelength;
       ASKAPLOG_DEBUG_STR(logger, "Antenna "<<ant);
       ASKAPLOG_DEBUG_STR(logger, "Amplitude of the residual "<<ampl<<" +/- "<<amplErr<<" metres");
       ASKAPLOG_DEBUG_STR(logger, "Phase of the residual "<<param[1]/casacore::C::pi*180.<<" +/- "<<err[1] / casacore::C::pi * 180.<<" deg");
       const double cphi = cos(param[1]);
       const double sphi = sin(param[1]);
       const double cphiErr = abs(sphi)*err[1];
       const double sphiErr = abs(cphi)*err[1];
       const double dX = -ampl * cphi;
       const double dXErr = sqrt(casacore::square(ampl * cphiErr) + casacore::square(cphi * amplErr));
       const double dY = -ampl * sphi;
       const double dYErr = sqrt(casacore::square(ampl * sphiErr) + casacore::square(sphi * amplErr));
       ASKAPDEBUGASSERT(itsCorrections.nrow()>ant);
       ASKAPDEBUGASSERT(itsErrors.nrow()>ant);
       ASKAPDEBUGASSERT(itsCorrections.ncolumn()>=2);
       ASKAPDEBUGASSERT(itsErrors.ncolumn()>=2);
       itsCorrections(ant,0) = dX; 
       itsCorrections(ant,1) = dY; 
       itsErrors(ant,0) = dXErr;
       itsErrors(ant,1) = dXErr;
    
       ASKAPLOG_DEBUG_STR(logger, "dX: "<<dX<<" +/- "<<dXErr<<" metres");
       ASKAPLOG_DEBUG_STR(logger, "dY: "<<dY<<" +/- "<<dYErr<<" metres");       
  }
      
}

/// @brief obtain MRO reference position
/// @return position measure
casacore::MPosition BaselineSolver::mroPosition()
{
   casacore::MPosition mroPos(casacore::MVPosition(casacore::Quantity(370.81, "m"), casacore::Quantity(116.6310372795, "deg"), 
                          casacore::Quantity(-26.6991531922, "deg")), casacore::MPosition::Ref(casacore::MPosition::WGS84));
   return mroPos;
}


/// @brief solve for dZ and populate corrections
/// @param[in] scans description of scans, note separate beams are present as separate scans.
///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
///            not necessarily match the scans known to online system
/// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
///            (column index is antenna ID used in the measurement set).
void BaselineSolver::solveForZ(const ScanStats &scans, const casacore::Matrix<GenericCalInfo> &caldata)
{
  ASKAPASSERT(scans.size()>2);
  std::vector<size_t> scanIndices(scans.size());
  for (size_t scan=0; scan<scanIndices.size(); ++scan) {
       scanIndices[scan]=scan;
  }
  const casacore::MPosition mroPos = mroPosition();
  std::sort(scanIndices.begin(),scanIndices.end(), utility::indexedCompare<size_t>(scans.begin(), LesserHAorDec(mroPos,false)));

  for (casacore::uInt ant=0; ant < caldata.ncolumn(); ++ant) {
       // do LSF into phase vs. sin(dec)
       double sx = 0., sy = 0., sx2 = 0., sy2 = 0., sxy = 0.;
       
       scimath::PhaseUnwrapper<double> unwrapper;
       const std::string spcfilename = "phase.z.ant"+utility::toString(ant)+".dat";
       std::ofstream os(spcfilename.c_str());
       casacore::uInt counter = 0;
       for (casacore::uInt cnt = 0; cnt < scans.size(); ++cnt) {
            const GenericCalInfo& cInfo = caldata(scanIndices[cnt],ant);
            if (!cInfo.gainDefined()) {
                continue;
            }
            const ObservationDescription& scan = scans[scanIndices[cnt]];
            const double time = 0.5*(scan.startTime() + scan.endTime());
            const casacore::MEpoch epoch(casacore::Quantity(time/86400.,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));
            casacore::MeasFrame frame(mroPos, epoch);    
            casacore::MVDirection hadec = casacore::MDirection::Convert(casacore::MDirection(scan.direction(),casacore::MDirection::J2000), 
                                   casacore::MDirection::Ref(casacore::MDirection::HADEC,frame))().getValue();
            const double sd = sin(hadec.getLat());
            const double phase = unwrapper(arg(cInfo.gain()));
            /*
            // the following code can be useful is unwrapper fails to do a decent job
            // due to gaps in declination coverage and large Z correction
            if ((ant == 8) && (sd > -0.5)) {
                phase += 2.*casacore::C::pi;
            }
            

            // the following code is useful to exclude sources based on sin(dec)
            if (sd > 0.2) {
                continue;
            }
            */
            const casacore::MVDirection azel = casacore::MDirection::Convert(casacore::MDirection(scan.direction(),casacore::MDirection::J2000), 
                                   casacore::MDirection::Ref(casacore::MDirection::AZEL,frame))().getValue();

            os<<cnt<<" "<<sd<<" "<<phase / casacore::C::pi * 180.<<" "<<azel.getLat() / casacore::C::pi * 180.<<" "<<
                     azel.getLong() / casacore::C::pi * 180.<<std::endl;
            // build normal equations
            sx += sd;
            sx2 += sd*sd;
            sy += phase;
            sy2 += phase*phase;
            sxy += sd*phase;            
            ++counter;
       }
       if (counter == 0) {
           os<<"No valid data for this antenna!"<<std::endl;
           ASKAPLOG_WARN_STR(logger, "No valid data for antenna "<<ant);
           continue;
       }
       ASKAPCHECK(counter > 0, "Too few points to do the fit for ant="<<ant);
       sx /= double(counter);
       sy /= double(counter);
       sx2 /= double(counter);
       sy2 /= double(counter);
       sxy /= double(counter);
       const double denominator = sx2 - sx * sx;
       ASKAPCHECK(denominator > 0, "Degenerate case has been encountered");
       const double coeff = (sxy - sx * sy) / denominator;
 
       // for our test setup of 864.5 MHz central freq
       const double wavelength = casacore::C::c / 864.5e6; // effective wavelength in metres (to do get it from scan's frequency)
       // the formula for phase has -sin(dec)*dZ therefore, we have to swap the sign here
       const double dZ = -coeff / 2. / casacore::C::pi * wavelength;
       ASKAPDEBUGASSERT(ant < itsCorrections.nrow());
       ASKAPDEBUGASSERT(itsCorrections.ncolumn() > 2); 
       itsCorrections(ant,2) = dZ; 
       const double r=(sxy-sx*sy)/sqrt(denominator * (sy2-sy*sy));
       const double coeffErr = sqrt((sy2-sy*sy)/denominator)*sqrt((double)(1.0-r*r)/(double(scans.size())-2));       
       const double dZErr = coeffErr / 2. / casacore::C::pi * wavelength;
       itsErrors(ant,2) = dZErr;
       
       ASKAPLOG_INFO_STR(logger, "Antenna "<<ant<<" dZ: "<<dZ<<" +/- "<<dZErr<<" metres; r="<<r);
  }
}

/// @brief solve for rate-dependent decorrelation
/// @details This method is intended for an experiment with calibration solution attempting to fit
/// decorrelation due to imperfect rate compensation in the BETA system. Unless some effort is spent
/// in the future, this method is not supposed to be part of the production code and may not be
/// accessible from parset.
/// @param[in] scans description of scans, note separate beams are present as separate scans.
///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
///            not necessarily match the scans known to online system
/// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
///            (column index is antenna ID used in the measurement set).
void BaselineSolver::solveForDecorrelation(const ScanStats &scans, const casacore::Matrix<GenericCalInfo> &caldata)
{
  ASKAPASSERT(scans.size()>2);
  std::vector<size_t> scanIndices(scans.size());
  for (size_t scan=0; scan<scanIndices.size(); ++scan) {
       scanIndices[scan]=scan;
  }
  const casacore::MPosition mroPos = mroPosition();
  std::sort(scanIndices.begin(),scanIndices.end(), utility::indexedCompare<size_t>(scans.begin(), LesserHAorDec(mroPos,false)));
  
  // hardcoded antenna positions - nominal ones are enough as we're only interested in rates
  const double antxyz[6][3]={{-2556227.87416317, 5097380.43738531, -2848323.37716186},
                             {-2556084.669, 5097398.337, -2848424.133},
                             {-2556118.1113922, 5097384.72442044, -2848417.24758565},
                             {-2555389.85800794, 5097664.67850522, -2848561.91777563},
                             {-2556002.4999209, 5097320.29662981, -2848637.49793128},
                             {-2555888.89502723, 5097552.29356202, -2848324.68084774}};
  const double siderealRate = casacore::C::_2pi / 86400. / (1. - 1./365.25);

  for (casacore::uInt ant=0; ant < caldata.ncolumn(); ++ant) {
       ASKAPCHECK(ant < 6, "This experimental code only supports 6 antennas as it is BETA-specific");
       // do LSF into amp vs. rate
       double sx = 0., sy = 0., sx2 = 0., sy2 = 0., sxy = 0.;
       
       const std::string spcfilename = "amp.ant"+utility::toString(ant)+".dat";
       std::ofstream os(spcfilename.c_str());
       for (casacore::uInt cnt = 0; cnt < scans.size(); ++cnt) {
            const ObservationDescription& scan = scans[scanIndices[cnt]];
            const double time = 0.5*(scan.startTime() + scan.endTime());
            const casacore::MEpoch epoch(casacore::Quantity(time/86400.,"d"), casacore::MEpoch::Ref(casacore::MEpoch::UTC));
            //casacore::MeasFrame frame(mroPos, epoch);
            casacore::MeasFrame frame(epoch);
            double gast = casacore::MEpoch::Convert(epoch, casacore::MEpoch::Ref(casacore::MEpoch::GAST))().get("d").getValue("d");
            gast = (gast - casacore::Int(gast)) * casacore::C::_2pi; // into radians
            casacore::MDirection radec = casacore::MDirection::Convert(casacore::MDirection(scan.direction(),casacore::MDirection::J2000), 
                                   casacore::MDirection::Ref(casacore::MDirection::TOPO,frame))();
            const double H0 = gast - radec.getAngle().getValue()(0);
            const double sH0 = sin(H0);
            const double cH0 = cos(H0);                                         
            const double cd = cos(radec.getAngle().getValue()(1));
            
            // hardcoded effective LO freq
            //const double rate = (cd * sH0 * (antxyz[ant][0] - antxyz[1][0]) + cd * cH0 * (antxyz[ant][1] - antxyz[1][1])) *
            //              siderealRate * casacore::C::_2pi / casacore::C::c * 672e6;
            const double rate = (cd * sH0 * (antxyz[ant][0] - antxyz[1][0]) + cd * cH0 * (antxyz[ant][1] - antxyz[1][1])) *
                          siderealRate * casacore::C::_2pi / casacore::C::c * 864.5e6;
            const double amp = abs(caldata(scanIndices[cnt],ant).gain());
            
            os<<cnt<<" "<<rate<<" "<<amp<<std::endl;
            // build normal equations
            sx += rate;
            sx2 += rate*rate;
            sy += amp;
            sy2 += amp*amp;
            sxy += rate*amp;            
       }
       sx /= double(scans.size());
       sy /= double(scans.size());
       sx2 /= double(scans.size());
       sy2 /= double(scans.size());
       sxy /= double(scans.size());
       const double denominator = sx2 - sx * sx;
       ASKAPCHECK(denominator > 0, "Degenerate case has been encountered");
       const double coeff = (sxy - sx * sy) / denominator;
 
       const double r=(sxy-sx*sy)/sqrt(denominator * (sy2-sy*sy));
       const double coeffErr = sqrt((sy2-sy*sy)/denominator)*sqrt((double)(1.0-r*r)/(double(scans.size())-2));       
       
       ASKAPLOG_INFO_STR(logger, "Antenna "<<ant<<" coeff: "<<coeff<<" +/- "<<coeffErr<<" s; r="<<r);
  }  
}  



} // namespace synthesis

} // namespace askap

