/// @file
///
/// @brief Gridder taking w-term into account
/// @details
/// This is a base class for all gridders taking w-term into account. It manages sampling
/// in w-space (which may be non-linear, if so chosen by the user)
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

#include <gridding/WDependentGridderBase.h>
#include <gridding/PowerWSampling.h>
#include <gridding/GaussianWSampling.h>

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.wdependentgridderbase");

#include <fstream>
#include <iomanip>
#include <boost/shared_ptr.hpp>

using namespace askap;
using namespace askap::synthesis;

/// @brief constructor initialising for default linear sampling   
/// @details
/// @param[in] wmax Maximum baseline (wavelengths)
/// @param[in] nwplanes Number of w planes   
WDependentGridderBase::WDependentGridderBase(const double wmax, const int nwplanes) : itsWScale(wmax),
     itsNWPlanes(nwplanes) 
{
  ASKAPCHECK(wmax>0.0, "Baseline length must be greater than zero, you have wmax="<<wmax);
  ASKAPCHECK(nwplanes>0, "Number of w planes must be greater than zero, you have nwplanes="<<nwplanes);
  ASKAPCHECK(nwplanes%2==1, "Number of w planes must be odd, you have nwplanes="<<nwplanes);
 
  if (nwplanes>1) {
      itsWScale /= double((nwplanes-1)/2);
  } else {
      itsWScale = 1.0;
  }
}

/// @brief destructor, writes w-plane hit distribution if necessary
WDependentGridderBase::~WDependentGridderBase()
{
  if (itsWPlaneStats.size()) {
      const std::string tag = isPSFGridder() ? "(PSF gridder)" : "(non-PSF gridder)";
      ASKAPLOG_INFO_STR(logger, "W-plane hit statistics "<<tag<<":  plane(nwplanes="<<itsNWPlanes<<")  w-term(wmax="<<getWMax()<<")  Number of hits");
      for (int plane = 0; plane<int(itsWPlaneStats.size()); ++plane) {
           if (plane < itsNWPlanes) {
               const int val = itsWPlaneStats[plane];
               ASKAPLOG_INFO_STR(logger, "W-plane hit statistics:     "<<std::setw(4)<<plane<<"  "<<std::setw(9)<<std::setprecision(7)<<getWTerm(plane)<<"  "<<val);
           }
      }
  }
}

/// @brief increment w-plane hit statistics
/// @details This method is called from derived classes. It increments the cache of statistics
/// every time the appropriate plane is used.
/// @param[in] plane w-plane used (should be less than the number of wplanes)
void WDependentGridderBase::notifyOfWPlaneUse(const int plane) const
{
  // update statistics if necessary
  if (itsWPlaneStats.size()) {
      ASKAPDEBUGASSERT(int(itsWPlaneStats.size()) == itsNWPlanes);
      ASKAPDEBUGASSERT(plane < itsNWPlanes);
      ++itsWPlaneStats[plane];
  }         
}

/// @brief obtain plane number
/// @details
/// @param[in] w w-term (in wavelengths) to map
/// @return plane number
/// @note an exception is thrown if the requested w-term lies outside (-wmax,wmax) range
int WDependentGridderBase::getWPlane(const double w) const
{
  int result = -1;
  if (itsWSampling && (itsNWPlanes>1)) {
      // non-linear sampling of w-space is probably used
      const int halfNPlanes = (itsNWPlanes-1)/2;
      const double scaledPlaneNumber = itsWSampling->index(w/getWMax());
      result = halfNPlanes + nint(scaledPlaneNumber*halfNPlanes);
  } else {
      result = itsNWPlanes>1 ? (itsNWPlanes-1)/2 + nint(w/itsWScale) : 0;
  }
#ifdef ASKAP_DEBUG  
  if (result < 0) {
      ASKAPLOG_DEBUG_STR(logger, w << " " << itsWScale << " "<< result );
  }
#endif // #ifdef ASKAP_DEBUG
  
  ASKAPCHECK(result < itsNWPlanes,
           "W scaling error: recommend allowing larger range of w, you have w="<< w <<" wavelengths");
  ASKAPCHECK(result > -1,
           "W scaling error: recommend allowing larger range of w, you have w="<< w <<" wavelengths");
  return result;
}
   
/// @brief obtain w-term for a given plane
/// @details This is a reverse operation to wPlaneNumber.
/// @param[in] plane plane number
/// @return w-term (in wavelengths) corresponding to the given plane
/// @note an exception is thrown (in debug mode) if the plane is outside [0.plane) range
double WDependentGridderBase::getWTerm(const int plane) const
{
  ASKAPDEBUGASSERT( (plane >=0 ) && (plane<itsNWPlanes) );
  const int halfNPlanes = (itsNWPlanes-1)/2;
  if (itsWSampling && (itsNWPlanes>1)) {
      // non-linear sampling of w-space is probably used
      const double scaledPlaneNumber = double(plane-halfNPlanes)/double(halfNPlanes);
      return itsWSampling->map(scaledPlaneNumber)*getWMax();
  }
  return double(plane - halfNPlanes)*itsWScale;
}

/// @brief enable power law sampling in the w-space
/// @details After this method is called, w-planes will be spaced non-linearly 
/// (power law with the given exponent)
/// @param[in] exponent exponent of the power law
void WDependentGridderBase::powerLawWSampling(const double exponent)
{ 
  itsWSampling.reset(new PowerWSampling(exponent));
}

/// @brief enable gaussian sampling in the w-space
/// @details After this method is called, w-planes will be spaced non-linearly
/// (gaussian distribution defined by the given parameter)
/// @param[in] nwplanes50 The number of w-planes covering 50% of the w-term range [-wmax,wmax]. The first and the last
/// w-planes always correspond to -wmax and +wmax, while the mid-plane always corresponds to zero w-term.
void WDependentGridderBase::gaussianWSampling(const double nwplanes50)
{
  itsWSampling.reset(new GaussianWSampling(nwplanes50));
}


/// @brief configure w-sampling from the parset
/// @details This method hides all details about w-sampling common for all derived gridders
/// @param[in] parset parameter set (gridder name already removed)
void WDependentGridderBase::configureWSampling(const LOFAR::ParameterSet& parset)
{
  const std::string sampling = parset.getString("wsampling","linear");
  if (sampling != "linear") {
      if (sampling == "powerlaw") {
          const double exponent = parset.getDouble("wsampling.exponent");
          ASKAPLOG_INFO_STR(logger, "Power law sampling of the w-space, exponent = "<<exponent);
          powerLawWSampling(exponent);
      } else if (sampling == "gaussian") {
           const double nwplanes50 = parset.getDouble("wsampling.nwplanes50");
           ASKAPLOG_INFO_STR(logger, "Gaussian sampling of the w-space, number of w-planes covering 50% of sampled w-range is "<<
                                      nwplanes50);
           ASKAPCHECK( (nwplanes50>0.) && (nwplanes50 < double(itsNWPlanes)/sqrt(2.)), 
               "Number of w-planes covering 50% of w-range has to be a positive number not exceeding total number of w-planes "<<
               itsNWPlanes<<" divided by sqrt(2) or "<<double(itsNWPlanes)/sqrt(2.));
           gaussianWSampling(nwplanes50/itsNWPlanes);                           
      } else {
        ASKAPTHROW(AskapError, "W-sampling "<<sampling<<" is not implemented");
      }
  } else {
     ASKAPLOG_INFO_STR(logger, "Linear sampling of the w-space");
  }
  
  // optionally show w-plane hit statistics
  if (parset.getBool("wstats", false)) {
      ASKAPLOG_INFO_STR(logger, "W-plane usage distribution will be written into log at the end of processing");
      itsWPlaneStats.resize(nWPlanes(),0);
  }
  
  // optionally export w-terms corresponding to each w-plane into dat file (largely for debugging)
  const std::string dbgExportFileName = parset.getString("wsampling.export","");
  if (dbgExportFileName != "") {
      ASKAPLOG_INFO_STR(logger, "Exporing w-term for each w-plane and reverse translation into "<<dbgExportFileName);
      boost::shared_ptr<std::ofstream> os;
      if (dbgExportFileName != "log") {
          os.reset(new std::ofstream(dbgExportFileName.c_str()));
          ASKAPCHECK(*os, "Unable to open file "<<dbgExportFileName<<" for writing");
          *os<<"# plane  w-term  plane(w-term)"<<std::endl;
      }
      for (int plane = 0; plane<nWPlanes(); ++plane) {
           const double wTerm = getWTerm(plane);
           if (os) {
               *os<<plane<<" "<<wTerm<<" "<<getWPlane(wTerm)<<std::endl;
           } else {
              // exporting to the log instead
              ASKAPLOG_INFO_STR(logger, "w-plane="<<plane<<" -> w-term="<<wTerm<<" wavelengths -> plane="<<getWPlane(wTerm));
           }
      }
  }
}

