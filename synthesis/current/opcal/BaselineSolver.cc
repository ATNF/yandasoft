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
#include <opcal/BaselineSolver.h>
#include <askap/IndexedCompare.h>
#include <opcal/ObservationDescription.h>
#include <utils/PhaseUnwrapper.h>

// casa includes
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCDirection.h>
#include <scimath/Functionals/Function.h>
#include <scimath/Fitting/NonLinearFitLM.h>

// std
#include <algorithm>

ASKAP_LOGGER(logger, ".BaselineSolver");


namespace askap {

namespace synthesis {

// helper comparator class ordering scans by hour angle
struct LesserHourAngle {
 
   /// @brief constructor, sets up position on Earth to do the calculations for
   /// @details
   /// @param[in] pos position on the ground
   LesserHourAngle(const casa::MPosition &pos) : itsPosition(pos) {}
   

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
   casa::MPosition itsPosition;
};

/// @brief compares two directions 
/// @details We want to order scans in hour angle. This method 
/// checks two positions on the sky and determines which one
/// has a lesser hour angle. 
/// @param[in] scan1 first direction
/// @param[in] scan2 second direction
/// @return true, if the first scan corresponds to a lesser hour angle
/// @note Time is assumed to be in seconds since 0 MJD, directions are assumed to be in J2000
bool LesserHourAngle::operator()(const ObservationDescription &scan1, const ObservationDescription &scan2) const
{
  const double time1 = 0.5*(scan1.startTime() + scan1.endTime());
  const casa::MEpoch epoch1(casa::Quantity(time1/86400.,"d"), casa::MEpoch::Ref(casa::MEpoch::UTC));
  casa::MeasFrame frame1(itsPosition, epoch1);    
  casa::MVDirection hadec1 = casa::MDirection::Convert(casa::MDirection(scan1.direction(),casa::MDirection::J2000), 
         casa::MDirection::Ref(casa::MDirection::HADEC,frame1))().getValue();

  const double time2 = 0.5*(scan2.startTime() + scan2.endTime());
  const casa::MEpoch epoch2(casa::Quantity(time2/86400.,"d"), casa::MEpoch::Ref(casa::MEpoch::UTC));
  casa::MeasFrame frame2(itsPosition, epoch2);    
  casa::MVDirection hadec2 = casa::MDirection::Convert(casa::MDirection(scan2.direction(),casa::MDirection::J2000), 
         casa::MDirection::Ref(casa::MDirection::HADEC,frame2))().getValue();

  return hadec1.getLong() < hadec2.getLong();
}

////////////////////

// helper functional fitting sinusoid due to errors in X and Y component of the baseline
// to be used with the casa fitter. We could've used our own, but it doesn't solve for errors
class SinusoidDueToXY : public casa::Function<double> {
public:
   
   /// @brief dimensionality
   /// @return dimensionality
   virtual casa::uInt ndim() const { return 3u; };
   
   /// @brief evaluate - main method
   virtual double eval(casa::Function<double>::FunctionArg x) const {
       return param_p[0]*cos(x[0]+param_p[1])+param_p[2]; 
   };
   // Copy it
   virtual casa::Function<double>* clone() const { ASKAPTHROW(AskapError,"Not yet implemented"); };   
};

//////////////////


/// @brief constructor
/// @param[in] parset parameters of the solver
BaselineSolver::BaselineSolver(const LOFAR::ParameterSet& parset) : 
  itsSolveXY(parset.getBool("solvexy", true)), itsSolveZ(parset.getBool("solvez",false)),
  itsRefAnt(parset.getUint("refant",1))
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
    
/// @brief abstract method to be called when the results are ready
/// @details This method is called (in derived classes) when the results
/// are ready. Two parameters describe scans and contain calibration data, respectively.
/// @param[in] scans description of scans, note separate beams are present as separate scans.
///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
///            not necessarily match the scans known to online system
/// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
///            (column index is antenna ID used in the measurement set).
void BaselineSolver::process(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata)
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
  for (casa::uInt ant=0; ant<itsCorrections.nrow(); ++ant) {
       if (ant != itsRefAnt) {
           ASKAPLOG_INFO_STR(logger, "Antenna "<<ant);
           ASKAPLOG_INFO_STR(logger, "dX: "<<itsCorrections(ant,0)<<" +/- "<<itsErrors(ant,0)<<" metres");
           ASKAPLOG_INFO_STR(logger, "dY: "<<itsCorrections(ant,1)<<" +/- "<<itsErrors(ant,1)<<" metres");
           if (itsSolveZ) {
               ASKAPLOG_INFO_STR(logger, "dZ: "<<itsCorrections(ant,2)<<" +/- "<<itsErrors(ant,2)<<" metres");
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
void BaselineSolver::solveForXY(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata)
{
  ASKAPASSERT(scans.size()>1);
  casa::MPosition mroPos(casa::MVPosition(casa::Quantity(370.81, "m"), casa::Quantity(-26.6991531922, "deg"), 
           casa::Quantity(116.6310372795, "deg")),casa::MPosition::Ref(casa::MPosition::WGS84));
  std::vector<size_t> scanIndices(scans.size());
  for (size_t scan=0; scan<scanIndices.size(); ++scan) {
       scanIndices[scan]=scan;
  }
  std::sort(scanIndices.begin(),scanIndices.end(), utility::indexedCompare<size_t>(scans.begin(), LesserHourAngle(mroPos)));

  // use casa fitter
  casa::Vector<double> hangles(scans.size());
  casa::Vector<double> phases(scans.size());
  casa::Vector<double> sigma(scans.size(),1.);
  ASKAPDEBUGASSERT(scanIndices.size() == scans.size());
  
  for (casa::uInt ant=0; ant < caldata.ncolumn(); ++ant) {
       casa::NonLinearFitLM<double> fitter;
       fitter.setFunction(SinusoidDueToXY());
       fitter.setMaxIter(80);
       fitter.setCriteria(0.001);
  
       // could've cached the hour angles, but for now leave as is
       scimath::PhaseUnwrapper<double> unwrapper;
       for (casa::uInt cnt = 0; cnt < scans.size(); ++cnt) {
            const ObservationDescription& scan = scans[scanIndices[cnt]];
            const double time = 0.5*(scan.startTime() + scan.endTime());
            const casa::MEpoch epoch(casa::Quantity(time/86400.,"d"), casa::MEpoch::Ref(casa::MEpoch::UTC));
            casa::MeasFrame frame(mroPos, epoch);    
            casa::MVDirection hadec = casa::MDirection::Convert(casa::MDirection(scan.direction(),casa::MDirection::J2000), 
                                   casa::MDirection::Ref(casa::MDirection::HADEC,frame))().getValue();
            hangles[cnt] = hadec.getLong() - mroPos.getValue().getLong(); // Hour angle at latitude 0
            const double cd = cos(hadec.getLat()); 
            ASKAPCHECK(cd > 0, "Cannot work with sources at either pole");
            phases[cnt] = unwrapper(arg(caldata(scanIndices[cnt],ant).gain())) / cd;
       }
       casa::Vector<double> param = fitter.fit(hangles,phases,sigma);
       ASKAPCHECK(param.nelements() == 3, "Expect 3 parameters out of the fitter, you have size="<<param.nelements());
       casa::Vector<double> err = fitter.errors();
       ASKAPCHECK(err.nelements() == 3, "Expect 3 uncertainties out of the fitter, you have size="<<err.nelements());
       const double wavelength = casa::C::c / 672e6; // effective wavelength in metres (to do get it from scan's frequency)
       double ampl = param[0] / 2. / casa::C::pi * wavelength;
       // fit can converge with either sign of the first coefficient, but we like to always have a positive amplitude
       if (ampl < 0) {
           ampl = -ampl;
           param[0] = -param[0];
           param[1] += casa::C::pi;
       }
       const double amplErr = err[0] / 2. / casa::C::pi * wavelength;
       ASKAPLOG_INFO_STR(logger, "Antenna "<<ant);
       ASKAPLOG_INFO_STR(logger, "Amplitude of the error "<<ampl<<" +/- "<<amplErr<<" metres");
       const double cphi = cos(param[1]);
       const double sphi = sin(param[1]);
       const double cphiErr = abs(sphi)*err[1];
       const double sphiErr = abs(cphi)*err[1];
       const double dX = -ampl * cphi;
       const double dXErr = sqrt(casa::square(ampl * cphiErr) + casa::square(cphi * amplErr));
       const double dY = -ampl * sphi;
       const double dYErr = sqrt(casa::square(ampl * sphiErr) + casa::square(sphi * amplErr));
       ASKAPLOG_INFO_STR(logger, "dX: "<<dX<<" +/- "<<dXErr<<" metres");
       ASKAPLOG_INFO_STR(logger, "dY: "<<dY<<" +/- "<<dYErr<<" metres");       
  }
      
}

/// @brief solve for dZ and populate corrections
/// @param[in] scans description of scans, note separate beams are present as separate scans.
///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
///            not necessarily match the scans known to online system
/// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
///            (column index is antenna ID used in the measurement set).
void BaselineSolver::solveForZ(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata)
{
  ASKAPTHROW(AskapError, "Not yet implemented");
}


} // namespace synthesis

} // namespace askap

