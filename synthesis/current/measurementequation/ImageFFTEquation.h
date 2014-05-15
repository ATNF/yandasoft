/// @file
///
/// ImageDFTEquation: Equation for discrete Fourier transform of an image
/// using gridding and FFTs.
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
#ifndef SYNIMAGEFFTEQUATION_H_
#define SYNIMAGEFFTEQUATION_H_

#include <fitting/Params.h>
#include <fitting/ImagingEquation.h>
#include <utils/ChangeMonitor.h>

#include <gridding/IVisGridder.h>
#include <dataaccess/SharedIter.h>
#include <dataaccess/IDataIterator.h>
#include <measurementequation/IVisCubeUpdate.h>

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>

#include <map>

#include <boost/shared_ptr.hpp>

namespace askap
{
  namespace synthesis
  {

    /// @brief FFT-based image equations
    ///
    /// @details This class does predictions and calculates normal equations
    /// images. Parameter names are image.{i,q,u,v}.*
    /// The transforms are done using gridding and FFTs.
    /// @ingroup measurementequation
    class ImageFFTEquation : public askap::scimath::ImagingEquation
    {
      public:

        /// Standard constructor
        /// @param ip Parameters
        /// @param idi Data iterator
        ImageFFTEquation(const askap::scimath::Params& ip,
          accessors::IDataSharedIter& idi);
        
        /// Constructor with default parameters
        /// @param idi Data iterator
        ImageFFTEquation(accessors::IDataSharedIter& idi);

        /// Standard constructor with specified gridder
        /// @param ip Parameters
        /// @param idi Data iterator
        /// @param gridder Shared pointer to a gridder
        ImageFFTEquation(const askap::scimath::Params& ip,
          accessors::IDataSharedIter& idi, IVisGridder::ShPtr gridder);
        
        /// Constructor with default parameters with specified gridder
        /// @param idi Data iterator
        /// @param gridder Shared pointer to a gridder
        ImageFFTEquation(accessors::IDataSharedIter& idi, IVisGridder::ShPtr gridder);

        /// Copy constructor
        ImageFFTEquation(const ImageFFTEquation& other);
        
        /// Assignment operator
        ImageFFTEquation& operator=(const ImageFFTEquation& other);


        virtual ~ImageFFTEquation();
        
        /// Return the default parameters
        static askap::scimath::Params defaultParameters();

/// Predict model visibility
        virtual void predict() const;

/// Calculate the normal equations
/// @param ne Normal equations
        virtual void calcImagingEquations(askap::scimath::ImagingNormalEquations& ne) const;

        /// Clone this into a shared pointer
        /// @return shared pointer to a copy
        virtual ImageFFTEquation::ShPtr clone() const;

        /// @brief assign a different iterator
        /// @details This is a temporary method to assign a different iterator.
        /// All this business is a bit ugly, but should go away when all
        /// measurement equations are converted to work with accessors.
        /// @param idi shared pointer to a new iterator
        void setIterator(accessors::IDataSharedIter& idi);


        /// @brief define whether the default spheroidal function gridder is used for PSF
        /// @details We have an option to build PSF using the default spheriodal function
        /// gridder, i.e. no w-term and no primary beam is simulated, as an alternative 
        /// to the same user-defined gridder as used for the model. Apart from the speed, 
        /// it probably makes the overall approximation better (i.e. removes some factors 
        /// of spatial dependence).
        /// @param[in] useSphFunc true, if spheroidal function gridder is to be used for PSF
        void useSphFuncForPSF(bool useSphFunc);
        
        /// @brief setup object function to update degridded visibilities
        /// @details For the parallel implementation of the measurement equation we need
        /// inter-rank communication. To avoid introducing cross-dependency of the measurement
        /// equation and the MPI one can use polymorphic object function to sum degridded visibilities 
        /// across all required ranks in the distributed case and do nothing otherwise.
        /// By default, this class doesn't alter degridded visibilities.
        /// @param[in] obj new object function (or an empty shared pointer to turn this option off)
        void setVisUpdateObject(const boost::shared_ptr<IVisCubeUpdate> &obj);
        
      private:
      
      /// Pointer to prototype gridder
        IVisGridder::ShPtr itsGridder;
        
        /// Map of gridders for the model
        mutable std::map<string, IVisGridder::ShPtr> itsModelGridders;
        
        /// Map of gridders for the residuals
        mutable std::map<string, IVisGridder::ShPtr> itsResidualGridders;
        
        /// Map of PSF gridders
        mutable std::map<string, IVisGridder::ShPtr> itsPSFGridders;

        /// Iterator giving access to the data
        mutable accessors::IDataSharedIter itsIdi;

        /// @brief change monitors per image parameter
        /// @details This objects are used to determine whether a new 
        /// initialise degrid is necessary (i.e. image or coordinate system
        /// has been updated since the last call).
        mutable std::map<std::string, scimath::ChangeMonitor> itsImageChangeMonitors; 
        
        /// @brief helper method to verify whether a parameter had been changed 
        /// @details This method checks whether a particular parameter is tracked. If 
        /// yes, its change monitor is used to verify the status since the last call of
        /// the method, otherwise new tracking begins and true is returned (i.e. to 
        /// update all dependent cache).
        /// @param[in] name name of the parameter
        /// @return true if parameter has been updated since the previous call
        bool notYetDegridded(const std::string &name) const;
        
        void init();
        
        /// @brief true, if the PSF is built using the default spheroidal function gridder
        /// @details We have an option to build PSF using the default spheriodal function
        /// gridder, i.e. no w-term and no primary beam is simulated. Apart from speed, 
        /// it probably makes the approximation better (i.e. removes some factors of spatial 
        /// dependence).
        bool itsSphFuncPSFGridder;

        /// @brief if set, visibility cube will be passed through this object function
        /// @details For the parallel implementation of the measurement equation we need
        /// inter-rank communication. To avoid introducing cross-dependency of the measurement
        /// equation and the MPI one can use polymorphic object function to sum degridded visibilities 
        /// across all required ranks in the distributed case and do nothing otherwise.
        boost::shared_ptr<IVisCubeUpdate> itsVisUpdateObject;
    };

  }

}
#endif
