/// @file
/// @brief Common functionality for all mosaicing gridders
/// @details AProjectGridderBase class encapsulates common operations for all mosaicing 
/// gridders: CF cache support and recalculation statistics, support for the buffer in the uv-space,
/// and the factory of illumination pattrns.
///
/// This is a template implementation file. We use templates to avoid duplication of initialisation
/// code for two mosaicing gridders 
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

#ifndef A_PROJECT_GRIDDER_BASE_TCC
#define A_PROJECT_GRIDDER_BASE_TCC

#include <askap/AskapLogging.h>

namespace askap {

namespace synthesis {


/// @brief actual factory of derived gridders
/// @details Gridders derived from this class use exactly the same parameters, but doing
/// the job using slightly different algorithms. This templated method allows to create both
/// of them without duplication of the code extracting parameters from the parset file.
/// @note Parset should have the gridder name already removed from all parameter names (i.e.
/// a proper subset should be formed first, see IVisGridder::createGridder as well)
/// @param[in] parset parameter set containing description of the gridder to create
/// @return shared pointer to the gridder
template<typename GridderType>
boost::shared_ptr<GridderType> AProjectGridderBase::createAProjectGridder(const LOFAR::ParameterSet &parset)
{
  // just for logging, declare private handle to avoid issues with template
  log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(askap::generateLoggerName(std::string("createAProjectGridder")));
  //
  const double pointingTol = SynthesisParamsHelper::convertQuantity(
               parset.getString("pointingtolerance", "0.0001rad"),"rad");
  const double paTol = SynthesisParamsHelper::convertQuantity(parset.getString("patolerance", "0.1rad"),"rad");
		
  // load frequency tolerance. It can be a string "infinite", which means to bypass frequency
  // axis checks or a non-negative number. We pass "undefined" by default here as it is not a
  // recognised value, which will cause a numeric default to be adopted 
  double freqTol = -1.;
  const std::string freqTolString = toLower(parset.getString("freqtolerance", "undefined"));
  if (freqTolString != "infinite") {
      freqTol = parset.getDouble("freqtolerance", 1e-6);
		    ASKAPCHECK(freqTol>=0., 
		        "Frequency tolerance parameter is supposed to be either a non-negative number or a word infinite to by pass checks");
  }
  //
  const double wmax = parset.getDouble("wmax", 10000.0);
  const int nwplanes = parset.getInt32("nwplanes", 65);
  // strictly speaking cutoff is required for AWProject gridder only, we may need to think about a more
  // tidy way of working around this situation.
  const double cutoff=parset.getDouble("cutoff", 1e-3);
  const int oversample = parset.getInt32("oversample", 8);
  const int maxSupport = parset.getInt32("maxsupport",128);
  const int limitSupport = parset.getInt32("limitsupport",0);
  const int maxFeeds = parset.getInt32("maxfeeds", 1);
  const int maxFields = parset.getInt32("maxfields", 1);
  if (parset.isDefined("maxantennas")) {
      ASKAPLOG_WARN_STR(logger,"maxantennas is no longer used! Update your parset");
  }
  const bool freqDep=parset.getBool("frequencydependent", true);
  const string tablename=parset.getString("tablename","");
  ASKAPLOG_INFO_STR(logger,"Gridding with Antenna Illumination projection");
  if (freqDep) {
      ASKAPLOG_INFO_STR(logger,	"Antenna illumination scales with frequency");
  } else {
      ASKAPLOG_INFO_STR(logger,	"Antenna illumination independent of frequency");
  }
  ASKAPLOG_INFO_STR(logger, "Parallactic angle tolerance = "<<paTol/casa::C::pi*180.<<" deg");
  ASKAPLOG_INFO_STR(logger, "Fields offset by more than "<<pointingTol/casa::C::pi*180.*3600.<<
                             " arcsec will be considered separate fields");
  if (freqTol<0) {
	  ASKAPLOG_INFO_STR(logger,"Frequency axis is assumed constant"); 
  } else {
	  ASKAPLOG_INFO_STR(logger,"Frequency axis tolerance (relative) is "<<freqTol<<
	                   " (equivalent to "<<freqTol*3e5<<" km/s)"); 
  }

  boost::shared_ptr<GridderType> gridder(new GridderType(makeIllumination(parset),
                wmax, nwplanes, cutoff, oversample, maxSupport, limitSupport, maxFeeds, maxFields,
                pointingTol, paTol, freqTol, freqDep, tablename));
  gridder->configureWSampling(parset);
  return gridder;              
}


} // namespace synthesis

} // namespace askap

#endif // #ifndef A_PROJECT_GRIDDER_BASE_TCC


