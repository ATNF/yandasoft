/// @file
/// 
/// @brief Estimation of the restoring beam
/// @details This method encapsulate operations with restoring beam parameters.
/// We support both the explicit definition and the PSF fitting. However, in some
/// cases PSF may not present at the time of initialisation, so we cannot fit it
/// straight away. This class allows us to delay this fitting until the PSF is
/// calculated and keep all related information together.
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

#include <measurementequation/RestoringBeamHelper.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <askap/AskapError.h>

namespace askap {

namespace synthesis {

/// @brief default constructor - uninitialised class
/// @details An exception is thrown if one attempts to access beam parameters
RestoringBeamHelper::RestoringBeamHelper() : itsCutoff(-1.) {}

/// @brief construct with explicitly defined beam parameters
/// @param[in] beam beam parameters (should be 3 elements)
RestoringBeamHelper::RestoringBeamHelper(const casa::Vector<casa::Quantum<double> > &beam) : itsBeam(beam.copy()), itsCutoff(1.) 
{
  ASKAPCHECK(beam.nelements() == 3, "Bean parameters should be given in a vector with 3 elements, you have "<<beam.nelements());
}

/// @brief copy constructor
/// @param[in] other other instance
RestoringBeamHelper::RestoringBeamHelper(const RestoringBeamHelper &other) : itsBeam(other.itsBeam.copy()),
       itsCutoff(other.itsCutoff) {}

/// @brief assignment operator
/// @param[in] other other instance
RestoringBeamHelper& RestoringBeamHelper::operator=(const RestoringBeamHelper &other)
{
  if (&other != this) {
      itsBeam.assign(other.itsBeam.copy());
      itsCutoff = other.itsCutoff;
  }
  return *this;
}
   
/// @brief construct for a delayed fit
/// @param[in] cutoff relative cutoff to determine which pixels are included in the fit
RestoringBeamHelper::RestoringBeamHelper(const double cutoff) : itsCutoff(cutoff) 
{
  ASKAPCHECK(cutoff >= 0., "RestoringBeamHelper::configureFit - negative cutoff is not allowed, you have cutoff="<<cutoff);
}
   
/// @brief initialise with explicitly defined beam parameters
/// @param[in] beam beam parameters (should be 3 elements)
void RestoringBeamHelper::assign(const casa::Vector<casa::Quantum<double> > &beam)
{
  ASKAPCHECK(beam.nelements() == 3, "Bean parameters should be given in a vector with 3 elements, you have "<<beam.nelements());
  itsBeam.assign(beam.copy());
  itsCutoff = 1.; // just a flag that the object is now initialised
}
   
/// @brief initialise for a delayed fit
/// @param[in] cutoff relative cutoff to determine which pixels are included in the fit
void RestoringBeamHelper::configureFit(const double cutoff) {
  ASKAPCHECK(cutoff >= 0., "RestoringBeamHelper::configureFit - negative cutoff is not allowed, you have cutoff="<<cutoff);
  itsCutoff = cutoff;
  itsBeam.resize(0);
}
   
/// @return true, if the class is initialised
bool RestoringBeamHelper::valid() const 
{
  return itsCutoff >= 0.;
}
   
/// @return true, if PSF fit is required
bool RestoringBeamHelper::fitRequired() const
{
  return !valid() || (itsBeam.nelements() != 3);
}
   
/// @brief perform the fit
/// @details This method performs the fit, it should be called if fitRequired() returns
/// true before any attempt to access the result
/// @param[in] ip parameters (the first ecountered PSF parameter is fitted)
void RestoringBeamHelper::fitBeam(const scimath::Params &ip)
{
   ASKAPCHECK(valid(), "RestoringBeamHelper::fitBeam is called before the fit is properly configured");
   // we could also move fitBeam into this class from SynthesisParamsHelper
   itsBeam.assign(SynthesisParamsHelper::fitBeam(ip,itsCutoff).copy());
   ASKAPDEBUGASSERT(itsBeam.nelements() == 3);
}
   
/// @brief access the result
/// @return parameters of the restoring beam (3-element vector)
const casa::Vector<casa::Quantum<double> >& RestoringBeamHelper::value() const
{
   ASKAPCHECK(itsBeam.nelements() == 3, "call to RestoringBeamHelper::value() is made for uninitialised object or prior to the PSF fit");
   return itsBeam;
}


} // namespace synthesis

} // namespace askap

