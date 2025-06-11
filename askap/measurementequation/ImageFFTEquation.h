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

#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/ImagingEquation.h>
#include <askap/scimath/utils/ChangeMonitor.h>

#include <askap/gridding/IVisGridder.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/IDataIterator.h>
#include <askap/measurementequation/IVisCubeUpdate.h>
#include <askap/gridding/IUVWeightAccessor.h>
#include <askap/gridding/UVWeightCollection.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>

#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>

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

          /// Standard constructor with specified gridder
          /// @param ip Parameters - shared
          /// @param idi Data iterator
          /// @param gridder Shared pointer to a gridder
        ImageFFTEquation(const askap::scimath::Params::ShPtr& ip,
          accessors::IDataSharedIter& idi, IVisGridder::ShPtr gridder);

        /// Standard constructor with specified gridder
        /// @param ip Parameters -shared
        /// @param idi Data iterator
        /// @param gridder Shared pointer to a gridder
        /// @param parset parameter set to check for PSF/PCF options.
        ImageFFTEquation(const askap::scimath::Params::ShPtr& ip,
          accessors::IDataSharedIter& idi, IVisGridder::ShPtr gridder,
          const LOFAR::ParameterSet& parset);

        /// Constructor with default parameters with specified gridder
        /// @param idi Data iterator
        /// @param gridder Shared pointer to a gridder
        ImageFFTEquation(accessors::IDataSharedIter& idi, IVisGridder::ShPtr gridder);

        /// Copy constructor
        ImageFFTEquation(const ImageFFTEquation& other);

        /// Assignment operator
        ImageFFTEquation& operator=(const ImageFFTEquation& other);

        virtual ~ImageFFTEquation();

        /// @brief configure from parset
        /// @details set all configuration options from a parset
        /// @param parset parameter set
        void configure(const LOFAR::ParameterSet& parset);

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

        /// @brief setup object function to update degridded visibilities
        /// @details For the parallel implementation of the measurement equation we need
        /// inter-rank communication. To avoid introducing cross-dependency of the measurement
        /// equation and the MPI one can use polymorphic object function to sum degridded visibilities
        /// across all required ranks in the distributed case and do nothing otherwise.
        /// By default, this class doesn't alter degridded visibilities.
        /// @param[in] obj new object function (or an empty shared pointer to turn this option off)
        void setVisUpdateObject(const boost::shared_ptr<IVisCubeUpdate> &obj);

        /// I want access to the gridders but I dont want to change them

        IVisGridder::ShPtr getResidualGridder(std::string name) const {
            return itsResidualGridders[name];
        };

        IVisGridder::ShPtr getPreconGridder(std::string name) const {
            return itsPreconGridders[name];
        };

        IVisGridder::ShPtr getPSFGridder(std::string name) const {
            return itsPSFGridders[name];
        };

        /// DDCALTAG
        /// @brief set the number of DD calibration directions
        void setNDir(casacore::uInt nDir) const { itsNDir = nDir; }

        /// @brief set calibration direction map
        void setCalDirMap(const std::map<std::string, int>& map)
        { itsCalDirMap = map; itsNDir = (map.size() == 0 ? 1 : map.size());}

// DAM TRADITIONAL
        /// @brief set robustness
        void setRobustness(float value) const { itsRobustness = value; }

// DAM TRADITIONAL
    /// @brief get robustness
        boost::optional<float> getRobustness() const { return itsRobustness; }

      protected:
      /// @brief define whether to use an alternative gridder for the PSF
      /// and/or the preconditioner function.
      /// @details We have an option to build the PSF using the default
      /// spheriodal function gridder or the box (nearest neighbour) gridder,
      /// i.e. no w-term and no primary beam is simulated, as an alternative
      /// to the same user-defined gridder as used for the model. Apart from
      /// speed, it probably makes the approximation better (i.e. removes
      /// some factors of spatial  dependence) and may lead to cleaner
      /// preconditioning. However it can filter out features of the W and AW
      /// kernels during preconditioning when robustness approaches uniform
      /// weighting. A separate preconditioner function can also be selected.
      /// @param[in] parset imager parameter set to check for PSF options.
      /// Current options: sphfuncforpsf, boxforpsf and preconditioner.preservecf.
      void useAlternativePSF(const LOFAR::ParameterSet& parset);


        /// @brief helper method to make accessor for the given image parameter
        /// @details It encapsulates handling the Taylor terms the right way
        /// (same uv-weight for all Taylor terms) and translation the image name
        /// into parameter name in the model (via the appropriate Params helper class).
        /// Also index translation is encapsulated.
        /// @param[in] name image parameter name (the full one with "image" prefix -
        /// we always deal with the full name makes the code more readable, although we
        /// could've cut down some operations if we take the name without the leading
        /// "image").
        /// @return shared pointer to the uv-weight accessor object accepted by gridders
        boost::shared_ptr<IUVWeightAccessor> makeUVWeightAccessor(const std::string &name) const;

        /// @brief helper method to assign uv-weight accessor to the given gridder
        /// @param[in] gridder gridder to work with
        /// @param[in] acc uv-weight accessor to assign
        /// @note if the accessor is empty nothing is done. Otherwise, if the gridder is of a wrong type which
        /// doesn't support setting of an accessor, an exception is thrown
        static void assignUVWeightAccessorIfNecessary(const boost::shared_ptr<IVisGridder> &gridder, const boost::shared_ptr<IUVWeightAccessor const> &acc);

      private:

      /// Pointer to prototype gridders
        IVisGridder::ShPtr itsGridder, itsAltGridder;

        /// Map of gridders for the model
        mutable std::map<string, IVisGridder::ShPtr> itsModelGridders;

        /// Map of gridders for the residuals
        mutable std::map<string, IVisGridder::ShPtr> itsResidualGridders;

        /// Map of PSF gridders
        mutable std::map<string, IVisGridder::ShPtr> itsPSFGridders;

        /// Map of Wiener preconditioning gridders
        mutable std::map<string, IVisGridder::ShPtr> itsPreconGridders;

        /// Map of coordinate systems
        mutable std::map<string, casacore::CoordinateSystem> itsCoordSystems;

        /// reuse PSF and PCF grids?
        mutable bool itsReuseGrids;

        /// Map of bools indicating we can reuse gridding results
        mutable std::map<string, bool> itsReuseGrid;

        /// Iterator giving access to the data
        mutable accessors::IDataSharedIter itsIdi;

        /// DDCALTAG
        /// The number of separate directions (equations)
        mutable casacore::uInt itsNDir;

        /// The map of calibration directions
        mutable std::map<string, int> itsCalDirMap;

// DAM TRADITIONAL
        /// @brief store the desired robustness if given as a parameter
        mutable boost::optional<float> itsRobustness;

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
        /// dependence). However it can filter out features of the W and AW kernels during
        /// preconditioning when robustness approaches uniform weighting.
        bool itsSphFuncPSFGridder;

        /// @brief true, if the PSF is built using the box (nearest neighbour) gridder
        /// @details We have an option to build PSF using the box (nearest neighbour)
        /// gridder, i.e. no w-term and no primary beam is simulated. Pros and cons
        /// are similar to the itsSphFuncPSFGridder option, however this gridder can lead
        /// to cleaner preconditioning (esp. as robustness approaches uniform weighting).
        bool itsBoxPSFGridder;

        bool itsUsePreconGridder;

        /// @brief true, if we grid the offset fields using the default spheroidal function gridder
        bool itsSphFuncOffsetFields;

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
