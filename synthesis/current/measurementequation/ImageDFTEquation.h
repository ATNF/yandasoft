/// @file
///
/// ImageDFTEquation: Equation for discrete Fourier transform of an image
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///
#ifndef SYNIMAGEDFTEQUATION_H_
#define SYNIMAGEDFTEQUATION_H_

#include <fitting/Params.h>
#include <fitting/GenericEquation.h>

#include <dataaccess/SharedIter.h>
#include <dataaccess/IDataIterator.h>

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>

namespace askap
{
  namespace synthesis
  {

    /// @brief Discrete Fourier Transform-based image equations
    ///
    /// @details This class does predictions and calculates normal equations
    /// images. Parameter names are image.{i,q,u,v}.*
    /// @ingroup measurementequation
    class ImageDFTEquation : public askap::scimath::GenericEquation
    {
      public:
        /// Standard constructor
        /// @param ip Parameters
        /// @param idi Data iterator
        ImageDFTEquation(const askap::scimath::Params& ip,
          accessors::IDataSharedIter& idi);
        
        /// Constructor with default parameters
        /// @param idi Data iterator
        ImageDFTEquation(accessors::IDataSharedIter& idi);

        /// Copy constructor
        ImageDFTEquation(const ImageDFTEquation& other);
        
        /// Assignment operator
        ImageDFTEquation& operator=(const ImageDFTEquation& other);

        virtual ~ImageDFTEquation();

        /// Return the default parameters
        static askap::scimath::Params defaultParameters();


/// Predict model visibility
        virtual void predict() const;

/// Calculate the normal equations
/// @param ne Normal equations
        virtual void calcGenericEquations(askap::scimath::GenericNormalEquations& ne) const;

        /// Clone this into a shared pointer
        /// @return shared pointer to a copy
        virtual ImageDFTEquation::ShPtr clone() const;
      

      private:
      /// Iterator giving access to the visibility data
        accessors::IDataSharedIter itsIdi;

        void init();
        /// Calculate visibility, and optionally the derivatives.
        /// @param imagepixels Image pixels
        /// @param raStart Start of the RA axis (rad)
        /// @param raEnd End of the RA axis (rad)
        /// @param raCells Number of cells on RA axis
        /// @param decStart Start of the Dec axis (rad)
        /// @param decEnd End of the Dec axis (rad)
        /// @param decCells Number of cells on Dec axis
        /// @param freq Observing frequency (Hz)
        /// @param uvw UVW in a vector of rigid vectors
        /// @param vis Output visibility
        /// @param doDeriv Do we want to calculate the derivative of the 
        /// visibility with respect to each pixel?
        /// @param imageDeriv Image derivative with respect to each pixel
        static void calcVisDFT(const casa::Array<double>& imagePixels,
          const double raStart, const double raEnd, const int raCells,
          const double decStart, const double decEnd, const int decCells,
          const casa::Vector<double>& freq,
          const casa::Vector<casa::RigidVector<double, 3> >& uvw,
          casa::Matrix<double>& vis, bool doderiv, casa::Matrix<double>& imageDeriv);
    };

  }

}
#endif
