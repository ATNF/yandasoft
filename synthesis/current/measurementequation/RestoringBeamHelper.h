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

#ifndef RESTORING_BEAM_HELPER_H
#define RESTORING_BEAM_HELPER_H

// casa includes
#include <casa/Arrays/Vector.h>
#include <casa/Quanta.h>

// own includes
#include <fitting/Params.h>

// std vector
#include <string>
#include <vector>

namespace askap {

namespace synthesis {

/// @brief Estimation of the restoring beam
/// @details This method encapsulate operations with restoring beam parameters.
/// We support both the explicit definition and the PSF fitting. However, in some
/// cases PSF may not present at the time of initialisation, so we cannot fit it
/// straight away. This class allows us to delay this fitting until the PSF is
/// calculated and keep all related information together.
/// @ingroup measurementequation
class RestoringBeamHelper {
public:
   /// @brief default constructor - uninitialised class
   /// @details An exception is thrown if one attempts to access beam parameters
   RestoringBeamHelper();

   /// @brief construct with explicitly defined beam parameters
   /// @param[in] beam beam parameters (should be 3 elements)
   explicit RestoringBeamHelper(const casa::Vector<casa::Quantum<double> > &beam);
   
   /// @brief construct for a delayed fit
   /// @param[in] cutoff relative cutoff to determine which pixels are included in the fit
   explicit RestoringBeamHelper(const double cutoff);
   
   /// @brief copy constructor
   /// @param[in] other other instance
   explicit RestoringBeamHelper(const RestoringBeamHelper &other);

   /// @brief assignment operator
   /// @param[in] other other instance
   RestoringBeamHelper& operator=(const RestoringBeamHelper &other);
   
   
   /// @brief initialise with explicitly defined beam parameters
   /// @param[in] beam beam parameters (should be 3 elements)
   void assign(const casa::Vector<casa::Quantum<double> > &beam);
   
   /// @brief initialise for a delayed fit
   /// @param[in] cutoff relative cutoff to determine which pixels are included in the fit
   void configureFit(const double cutoff);
   
   /// @return true, if the class is initialised
   bool valid() const;
   
   /// @return true, if PSF fit is required
   bool fitRequired() const;
   
   /// @brief perform the fit
   /// @details This method performs the fit, it should be called if fitRequired() returns
   /// true before any attempt to access the result
   /// @param[in] ip parameters (the first ecountered PSF parameter is fitted)
   void fitBeam(const scimath::Params &ip);
   
   /// @brief access the result
   /// @return parameters of the restoring beam (3-element vector)
   const casa::Vector<casa::Quantum<double> >& value() const;
   
private:
   /// @brief parameters of the restoring beam
   /// @details This vector should always contain 3 elements. Otherwise, it is assumed that
   /// a fit is required.
   casa::Vector<casa::Quantum<double> > itsBeam;
   
   /// @brief relative cutoff for the pixels used for the PSF fit
   /// @details This data field is negative for an uninitialised object
   double itsCutoff;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef RESTORING_BEAM_HELPER_H

