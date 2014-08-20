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
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCDirection.h>

ASKAP_LOGGER(logger, ".PointingSolver");

namespace askap {

namespace synthesis {

/// @brief constructor
/// @param[in] parset parameters of the solver
PointingSolver::PointingSolver(const LOFAR::ParameterSet& parset) : 
   itsMaxPatternSize(SynthesisParamsHelper::convertQuantity(parset.getString("maxpatternsize", "5deg"),"rad"))
{
  ASKAPLOG_INFO_STR(logger, "Maximum angular size of the pointing pattern is assumed to be "<<itsMaxPatternSize / casa::C::pi * 180.<<" deg");
}    
    
/// @brief method to be called when the results are ready
/// @details This method is called (in derived classes) when the results
/// are ready. Two parameters describe scans and contain calibration data, respectively.
/// @param[in] scans description of scans, note separate beams are present as separate scans.
///            Scans here are defined by some splitting criterion used in OpCalImpl, and do
///            not necessarily match the scans known to online system
/// @param[in] caldata calibration data. A matrix with one row per scan. Columns represent antennas
///            (column index is antenna ID used in the measurement set).
void PointingSolver::process(const ScanStats &scans, const casa::Matrix<GenericCalInfo> &caldata)
{
   // need a better way to separate scans, but for now assume a fixed order:
   // +Dec, -Dec, boresight, +RA, -RA, boresight_long (the last one to be skipped)
   ASKAPASSERT(scans.size() == caldata.nrow());
   const size_t nScansPerPattern  = 6;
   ASKAPCHECK(scans.size() >= nScansPerPattern, "Require at least "<<nScansPerPattern<<" scans, you have "<<scans.size());
   for (size_t startScan = 0, pattern = 0; startScan + nScansPerPattern <= scans.size(); startScan += nScansPerPattern,++pattern) {
        const casa::MVDirection& phaseCentre = scans[startScan].direction();
        const casa::MVDirection azEl = mroAzEl(phaseCentre, 0.5*(scans[startScan + nScansPerPattern/2].startTime()+
                                               scans[startScan + nScansPerPattern/2].endTime()));
        ASKAPLOG_INFO_STR(logger, "Pointing pattern "<<pattern + 1<<" at "<<printDirection(phaseCentre)<<" J2000, Az="<<
                          azEl.getLong()/casa::C::pi*180.<<" deg, El="<<
                          azEl.getLat()/casa::C::pi*180.<<" deg:");
        // consistency check 
        ASKAPDEBUGASSERT(nScansPerPattern > 2);
        const double totalDuration = scans[startScan + nScansPerPattern - 2].endTime() - scans[startScan].startTime(); 
        ASKAPCHECK(totalDuration<1800., "Pattern seems to be too long: "<<totalDuration / 60.<<" min");
        for (size_t index = startScan + 1; index < startScan + nScansPerPattern; ++index) {
             ASKAPCHECK(scans[index].direction().separation(phaseCentre)<1e-6, "Phase centres seem to vary throughout the pattern");
        }
        // 
        for (casa::uInt ant=0; ant<caldata.ncolumn(); ++ant) {
             casa::Vector<GenericCalInfo> thisPatternData = caldata(casa::Slice(startScan, nScansPerPattern - 1),
                     casa::Slice(ant));
             ASKAPASSERT(thisPatternData.nelements() + 1 == nScansPerPattern);        
             const std::pair<double, double> result = solveOne(thisPatternData);
             ASKAPLOG_INFO_STR(logger, " antenna "<<ant<<":  x="<<result.first<<", y="<<result.second);
        }
   }
}

/// @brief solve for corrections for just one antenna
/// @param[in] caldata vector of data for all points of the pattern
/// @return pair with offsets in two coordinates
std::pair<double, double> PointingSolver::solveOne(const casa::Vector<GenericCalInfo> &caldata) const
{
  ASKAPASSERT(caldata.nelements() == 5);
  return std::pair<double,double>(0.,0.);
}


/// @brief obtain MRO reference position
/// @return position measure
casa::MPosition PointingSolver::mroPosition()
{
   casa::MPosition mroPos(casa::MVPosition(casa::Quantity(370.81, "m"), casa::Quantity(-26.6991531922, "deg"), 
           casa::Quantity(116.6310372795, "deg")),casa::MPosition::Ref(casa::MPosition::WGS84));
   return mroPos;
}

/// @brief obtain Az and El at MRO for the given direction
/// @param[in] dir direction of interest in J2000
/// @param[in] time time in seconds since 0 MJD
/// @return direction in AzEl
casa::MVDirection PointingSolver::mroAzEl(const casa::MVDirection &dir, double time)
{
  const casa::MEpoch epoch(casa::Quantity(time/86400.,"d"), casa::MEpoch::Ref(casa::MEpoch::UTC));
  casa::MeasFrame frame(mroPosition(), epoch);    
  return casa::MDirection::Convert(casa::MDirection(dir, casa::MDirection::J2000), 
         casa::MDirection::Ref(casa::MDirection::AZEL,frame))().getValue();
}



} // namespace synthesis

} // namespace askap
