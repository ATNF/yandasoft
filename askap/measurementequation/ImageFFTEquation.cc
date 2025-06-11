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

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.imagefftequation");

#include <askap/askap/AskapError.h>

#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/MemBufferDataAccessor.h>
#include <askap/dataaccess/DDCalBufferDataAccessor.h>
#include <askap/dataaccess/DDCalOnDemandNoiseAndFlagDA.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/measurementequation/CalibrationIterator.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/gridding/BoxVisGridder.h>
#include <askap/gridding/SphFuncVisGridder.h>
#include <askap/scimath/fitting/ImagingNormalEquations.h>
#include <askap/scimath/fitting/DesignMatrix.h>
#include <askap/scimath/fitting/Axes.h>
#include <askap/profile/AskapProfiler.h>
#include <askap/gridding/GenericUVWeightAccessor.h>
#include <askap/gridding/UVWeightParamsHelper.h>
#include <askap/gridding/VisGridderFactory.h>
#include <askap/measurementequation/ImageParamsHelper.h>

#include <casacore/scimath/Mathematics/RigidVector.h>

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <stdexcept>

using askap::scimath::Params;
using askap::scimath::Axes;
using askap::scimath::ImagingNormalEquations;
using askap::scimath::DesignMatrix;
using namespace askap::accessors;

namespace askap
{
  namespace synthesis
  {

    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params& ip,
        IDataSharedIter& idi) : scimath::Equation(ip),
      askap::scimath::ImagingEquation(ip), itsIdi(idi),
      itsSphFuncPSFGridder(false), itsBoxPSFGridder(false),
      itsUsePreconGridder(false), itsSphFuncOffsetFields(false), itsNDir(1), itsReuseGrids(false)
    {
      itsGridder = IVisGridder::ShPtr(new SphFuncVisGridder());
      init();
    }


    ImageFFTEquation::ImageFFTEquation(IDataSharedIter& idi) :
      itsIdi(idi), itsSphFuncPSFGridder(false), itsBoxPSFGridder(false),
      itsUsePreconGridder(false), itsSphFuncOffsetFields(false), itsNDir(1), itsReuseGrids(false)
    {
      itsGridder = IVisGridder::ShPtr(new SphFuncVisGridder());
      reference(defaultParameters().clone());
      init();
    }

    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params& ip,
        IDataSharedIter& idi, IVisGridder::ShPtr gridder) :
      scimath::Equation(ip), askap::scimath::ImagingEquation(ip),
      itsGridder(gridder), itsIdi(idi), itsSphFuncPSFGridder(false),
      itsBoxPSFGridder(false), itsSphFuncOffsetFields(false), itsUsePreconGridder(false), itsNDir(1),
      itsReuseGrids(false)
    {
      init();
    }

    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params::ShPtr& ip,
      accessors::IDataSharedIter& idi, IVisGridder::ShPtr gridder):
    scimath::Equation(ip), askap::scimath::ImagingEquation(ip),
    itsGridder(gridder), itsIdi(idi), itsSphFuncPSFGridder(false),
    itsBoxPSFGridder(false), itsSphFuncOffsetFields(false), itsUsePreconGridder(false), itsNDir(1),
    itsReuseGrids(false)
    {
      init();
    }


    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params::ShPtr& ip,
        IDataSharedIter& idi, IVisGridder::ShPtr gridder,
        const LOFAR::ParameterSet& parset) : scimath::Equation(ip),
      askap::scimath::ImagingEquation(ip), itsGridder(gridder), itsIdi(idi),
      itsSphFuncPSFGridder(false), itsBoxPSFGridder(false),
      itsUsePreconGridder(false), itsSphFuncOffsetFields(false), itsNDir(1), itsReuseGrids(false)
    {
      configure(parset);
      init();
    }

    ImageFFTEquation::ImageFFTEquation(IDataSharedIter& idi,
        IVisGridder::ShPtr gridder) :
      itsGridder(gridder), itsIdi(idi), itsSphFuncPSFGridder(false),
      itsBoxPSFGridder(false), itsSphFuncOffsetFields(false), itsUsePreconGridder(false), itsNDir(1), itsReuseGrids(false)
    {
      reference(defaultParameters().clone());
      init();
    }

    ImageFFTEquation::~ImageFFTEquation()
    {
    }

    void ImageFFTEquation::configure(const LOFAR::ParameterSet& parset)
    {
    // DAM TRADITIONAL
       if (parset.isDefined("gridder.robustness")) {
           const float robustness = parset.getFloat("gridder.robustness");
           ASKAPCHECK((robustness>=-2) && (robustness<=2), "gridder.robustness should be in the range [-2,2]");
           // also check that it is spectral line? Or do that earlier
           //  - won't work in continuum imaging, unless combo is done with combinechannels on a single worker
           //  - won't work with Taylor terms
           setRobustness(robustness);
       }
       useAlternativePSF(parset);
       itsReuseGrids = parset.getBool("reusegrids",false);
       if (itsReuseGrids) {
           ASKAPLOG_INFO_STR(logger, "Will reuse the PSF/PCF grids each major cycle");
       }
       itsSphFuncOffsetFields = parset.getBool("sphfuncforoffset", false);
       if (itsSphFuncOffsetFields  && !itsAltGridder) {
          const string gridder = parset.getString("gridder","");
          LOFAR::ParameterSet altParset = parset.makeSubset("");
          altParset.subtractSubset("gridder."+gridder);
          altParset.replace("gridder","SphFunc");
          altParset.adoptCollection(parset.makeSubset("gridder."+gridder),"gridder.SphFunc");
          itsAltGridder = VisGridderFactory::make(altParset);
       }
    }

    /// @brief define whether to use an alternative gridder for the PSF
    /// and/or the preconditioner function
    void ImageFFTEquation::useAlternativePSF(const LOFAR::ParameterSet& parset)
    {
      const bool useGentlePCF = parset.getBool("preconditioner.preservecf", parset.isDefined("preconditioner.Names"));
      const bool useBoxPSF = parset.getBool("boxforpsf", false);
      const bool useSphPSF = parset.getBool("sphfuncforpsf", false);

      if (useGentlePCF) {
         itsUsePreconGridder = true;
         ASKAPLOG_INFO_STR(logger,
             "A separate tophat-style gridder will be used to calculate the preconditioner function");
      }

      if (useBoxPSF) {
         itsBoxPSFGridder = true;
         ASKAPLOG_INFO_STR(logger,
             "The box (nearest neighbour) gridder will be used to calculate the PSF");
         if (useSphPSF) {
             ASKAPLOG_WARN_STR(logger,
                 "The spheroidal function will not be used to calculate the PSF");
         }
      } else if (useSphPSF) {
         itsSphFuncPSFGridder = true;
         ASKAPLOG_INFO_STR(logger,
             "The default spheroidal function gridder will be used to calculate the PSF");
      } else {
         ASKAPLOG_INFO_STR(logger,
             "The PSF will be calculated by the same gridder type as used for the model and residuals");
      }
    }

    askap::scimath::Params ImageFFTEquation::defaultParameters()
    {
      Params ip(true);
      ip.add("image");
      return ip;
    }

    ImageFFTEquation::ImageFFTEquation(const ImageFFTEquation& other) :
          Equation(), ImagingEquation()
    {
      operator=(other);
    }

    ImageFFTEquation& ImageFFTEquation::operator=(const ImageFFTEquation& other)
    {
      if(this!=&other)
      {
        static_cast<askap::scimath::Equation*>(this)->operator=(other);
        itsIdi=other.itsIdi;
        itsGridder = other.itsGridder;
        itsAltGridder = other.itsAltGridder;
        itsSphFuncPSFGridder = other.itsSphFuncPSFGridder;
        itsBoxPSFGridder = other.itsBoxPSFGridder;
        itsUsePreconGridder = other.itsUsePreconGridder;
        itsVisUpdateObject = other.itsVisUpdateObject;
      }
      return *this;
    }

    void ImageFFTEquation::init()
    {
    }

    /// Clone this into a shared pointer
    /// @return shared pointer to a copy
    ImageFFTEquation::ShPtr ImageFFTEquation::clone() const
    {
      return ImageFFTEquation::ShPtr(new ImageFFTEquation(*this));
    }

    /// @brief setup object function to update degridded visibilities
    /// @details For the parallel implementation of the measurement equation we need
    /// inter-rank communication. To avoid introducing cross-dependency of the measurement
    /// equation and the MPI one can use polymorphic object function to sum degridded visibilities
    /// across all required ranks in the distributed case and do nothing otherwise.
    /// By default, this class doesn't alter degridded visibilities.
    /// @param[in] obj new object function (or an empty shared pointer to turn this option off)
    void ImageFFTEquation::setVisUpdateObject(const boost::shared_ptr<IVisCubeUpdate> &obj)
    {
      itsVisUpdateObject = obj;
    }

    /// @brief helper method to verify whether a parameter had been changed
    /// @details This method checks whether a particular parameter is tracked. If
    /// yes, its change monitor is used to verify the status since the last call of
    /// the method, otherwise new tracking begins and true is returned (i.e. to
    /// update all dependent cache).
    /// @param[in] name name of the parameter
    /// @return true if parameter has been updated since the previous call
    bool ImageFFTEquation::notYetDegridded(const std::string &name) const
    {
      std::map<std::string,scimath::ChangeMonitor>::iterator it = itsImageChangeMonitors.find(name);
      bool result = true;
      if (it != itsImageChangeMonitors.end()) {
          result = parameters().isChanged(name,it->second);
          it->second = parameters().monitorChanges(name);
      } else {
          itsImageChangeMonitors[name] = parameters().monitorChanges(name);
      }
      return result;
    }


    void ImageFFTEquation::predict() const
    {
      ASKAPTRACE("ImageFFTEquation::predict");
      const std::vector<std::string> completions(parameters().completions("image"));

      // To minimize the number of data passes, we keep copies of the gridders in memory, and
      // switch between these. This optimization may not be sufficient in the long run.

      itsIdi.chooseOriginal();

      // DDCALTAG -- set increased buffer size
      if (itsNDir > 1) {
         try {
             DDCalBufferDataAccessor& accBuffer = dynamic_cast<DDCalBufferDataAccessor&>(*itsIdi);
             ASKAPLOG_DEBUG_STR(logger, "calling accBuffer.setNDir("<<itsNDir<<")");
             accBuffer.setNDir(itsNDir);
         } catch (const std::bad_cast&) {
             ASKAPTHROW(AskapError, "Wrong accessor type for DDCalibration in ImageFFTEquation::predict()");
         }
      }
      int dirIndex = itsNDir > 1 ? -1 : 0;

      ASKAPLOG_DEBUG_STR(logger, "Initialising for model degridding");
      for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        string imageName("image"+(*it));
        SynthesisParamsHelper::clipImage(parameters(),imageName);

        if(itsModelGridders.count(imageName)==0) {
          itsModelGridders[imageName]=itsGridder->clone();
        }
	    itsModelGridders[imageName]->customiseForContext(*it);

        if (notYetDegridded(imageName)) {
            ASKAPLOG_DEBUG_STR(logger, "Degridding image "<<imageName);
            const Axes axes(parameters().axes(imageName));
            casacore::Array<imtype> imagePixels;
            #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                imagePixels = parameters().valueF(imageName);
            #else
                imagePixels = parameters().value(imageName);
            #endif
            const casacore::IPosition imageShape(imagePixels.shape());
            itsModelGridders[imageName]->initialiseDegrid(axes, imagePixels);
        }

        // DDCALTAG -- set parameter for increased buffer size
        // could alternatively pass dirIndex * itsIdi->nRow(). See if one is better than the other
        if (itsNDir > 1) {
            if ( imageName.find(".taylor.") != std::string::npos ) {
                // Taylor terms are used, so only increment dirIndex if this is taylor.0
                if ( imageName.find(".taylor.0") != std::string::npos ) {
                    ASKAPLOG_DEBUG_STR(logger, imageName<<" is taylor 0, so incrementing DD-cal buffer index");
                    dirIndex++;
                } else {
                    ASKAPLOG_DEBUG_STR(logger, imageName<<" is not taylor 0 -- not incrementing DD-cal buffer index");
                }
            } else {
                // Taylor terms are not used, so increment dirIndex
                dirIndex++;
            }
            ASKAPLOG_DEBUG_STR(logger, "degridding "<<imageName<<" into buffer "<<dirIndex);
            const boost::shared_ptr<TableVisGridder const> &tGridder =
                   boost::dynamic_pointer_cast<TableVisGridder const>(itsModelGridders[imageName]);
            if (!tGridder) {
                ASKAPTHROW(AskapError,
                    "Wrong gridder type for DDCalibration in ImageFFTEquation::predict()");
            } else {
                tGridder->setSourceIndex(dirIndex);
            }
        }

      }
      // Loop through degridding the data
      ASKAPLOG_DEBUG_STR(logger, "Starting to degrid model" );

      // report every 5000000 degridded rows into log in the debug mode
      #ifdef ASKAP_DEBUG
      unsigned long report_every = 5000000;
      unsigned long total_rows = 0;
      unsigned long current_rows = 0;
      #endif // #ifdef ASKAP_DEBUG

      for (itsIdi.init();itsIdi.hasMore();itsIdi.next())
      {
        itsIdi->rwVisibility().set(0.0);
        for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
        {
            string imageName("image"+(*it));
            itsModelGridders[imageName]->degrid(*itsIdi);
        }
        #ifdef ASKAP_DEBUG
        const casacore::uInt nRow = itsIdi->nRow();
        total_rows += nRow;
        current_rows += nRow;
        if (current_rows > report_every) {
            current_rows = 0;
            ASKAPLOG_DEBUG_STR(logger, "Degridded "<<total_rows<<" rows of data");
        }
        #endif // #ifdef ASKAP_DEBUG
      }
      ASKAPLOG_DEBUG_STR(logger, "Finished degridding model" );
    };

    /// @brief assign a different iterator
    /// @details This is a temporary method to assign a different iterator.
    /// All this business is a bit ugly, but should go away when all
    /// measurement equations are converted to work with accessors.
    /// @param idi shared pointer to a new iterator
    void ImageFFTEquation::setIterator(IDataSharedIter& idi)
    {
      itsIdi = idi;
    }

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
    boost::shared_ptr<IUVWeightAccessor> ImageFFTEquation::makeUVWeightAccessor(const std::string &name) const
    {
      ASKAPTRACE("ImageFFTEquation::makeUVWeightAccessor");
      ImageParamsHelper iph(ImageParamsHelper::replaceLeadingWordWith(name, "image.",""));
      const std::string parName = iph.facetName();
      // a bit of technical debt - we don't actually need to write parameters here, but this is the only way to get shared pointer and
      // the interface is always read/write, so we can't easily implement a version of the constructor accepting const reference here
      // without splitting the helper class into two classes (const and non-const)
      const UVWeightParamsHelper hlp(rwParameters());
      if (hlp.exists(parName)) {
          const boost::shared_ptr<IUVWeightIndexTranslator> ttor = hlp.getIndexTranslator(parName);
          const boost::shared_ptr<UVWeightCollection> wts = hlp.getUVWeights(parName);
          boost::shared_ptr<GenericUVWeightAccessor> wtAcc(new GenericUVWeightAccessor(wts, ttor));
          return wtAcc;
      }
      return boost::shared_ptr<IUVWeightAccessor>();
    }

    /// @brief helper method to assign uv-weight accessor to the given gridder
    /// @param[in] gridder gridder to work with
    /// @param[in] acc uv-weight accessor to assign
    /// @note if the accessor is empty nothing is done. Otherwise, if the gridder is of a wrong type which
    /// doesn't support setting of an accessor, an exception is thrown
    void ImageFFTEquation::assignUVWeightAccessorIfNecessary(const boost::shared_ptr<IVisGridder> &gridder, const boost::shared_ptr<IUVWeightAccessor const> &acc)
    {
       if (acc) {
           const boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(gridder);
           ASKAPCHECK(tvg, "Gridder is either not setup or of a wrong type which doesn't support setting of the UV weight accessor");
           tvg->setUVWeightAccessor(acc);
       }
    }

    // Calculate the residual visibility and image. We transform the model on the fly
    // so that we only have to read (and write) the data once. This uses more memory
    // but cuts down on IO
    void ImageFFTEquation::calcImagingEquations(askap::scimath::ImagingNormalEquations& ne) const
    {
      ASKAPTRACE("ImageFFTEquation::calcImagingEquations");

      // We will need to loop over all completions i.e. all sources
      const std::vector<std::string> completions(parameters().completions("image"));
      const bool ddCal = itsCalDirMap.size() > 0;

      // To minimize the number of data passes, we keep copies of the gridders in memory, and
      // switch between these. This optimization may not be sufficient in the long run.
      // Set up initial gridders for model and for the residuals. This enables us to
      // do both at the same time.

      // we use the first flag to optionally change gridder after the first image
      string firstName;
      for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        const string imageName("image"+(*it));
        // remove taylor or facet parts from the name
        const string baseName(ImageParamsHelper(imageName).name());
        if (firstName.size() == 0) {
          firstName = baseName;
        }
        const bool first = (baseName == firstName);
        SynthesisParamsHelper::clipImage(parameters(),imageName);
        if(itsModelGridders.count(imageName)==0) {
          if (first || !itsSphFuncOffsetFields) {
            itsModelGridders[imageName]=itsGridder->clone();
          } else {
            itsModelGridders[imageName]= itsAltGridder->clone();
            ASKAPLOG_INFO_STR(logger, "Using Spheroidal gridder for "<<imageName);
          }
        }
        // obtain uv-weights accessor if the appropriate details are present in the model
        // (otherwise an empty shared pointer is returned). The logic inside makeUVWeightAccessor
        // ensures Taylor terms are handled appropriately
        const boost::shared_ptr<IUVWeightAccessor const> wtAcc = makeUVWeightAccessor(imageName);
        if (wtAcc) {
            ASKAPLOG_DEBUG_STR(logger, "UV Weight will be applied during gridding for "<<imageName);
        } else {
            ASKAPLOG_DEBUG_STR(logger, "UV Weight will not be applied during gridding for "<<imageName);
        }

        if(itsResidualGridders.count(imageName)==0) {
          if (first || !itsSphFuncOffsetFields) {
            itsResidualGridders[imageName]=itsGridder->clone();
          } else {
            itsResidualGridders[imageName]= itsAltGridder->clone();
          }
          assignUVWeightAccessorIfNecessary(itsResidualGridders[imageName], wtAcc);
        }

        if(itsPSFGridders.count(imageName)==0) {
          if (itsBoxPSFGridder) {
             boost::shared_ptr<BoxVisGridder> psfGridder(new BoxVisGridder);
             itsPSFGridders[imageName] = psfGridder;
          } else if (itsSphFuncPSFGridder) {
             boost::shared_ptr<SphFuncVisGridder> psfGridder(new SphFuncVisGridder);
             itsPSFGridders[imageName] = psfGridder;
          } else if (!first && itsSphFuncOffsetFields) {
             itsPSFGridders[imageName] = itsAltGridder->clone();
          } else {
             itsPSFGridders[imageName] = itsGridder->clone();
          }
          assignUVWeightAccessorIfNecessary(itsPSFGridders[imageName], wtAcc);
          itsReuseGrid[imageName] = false;
          boost::shared_ptr<TableVisGridder> tvgPSF = boost::dynamic_pointer_cast<TableVisGridder>(itsPSFGridders[imageName]);
          if (tvgPSF && itsReuseGrids) {
              tvgPSF->doClearGrid(false);
              ASKAPLOG_DEBUG_STR(logger, "Setting clear grid to false, reuse PSF grid for "<<imageName);
          }
        } else {
          // reuse the grid if gridder has been configured to allow this
          if (itsReuseGrids) {
              itsReuseGrid[imageName] = true;
              ASKAPLOG_DEBUG_STR(logger, "Will reuse PSF grid for "<<imageName);
          }

        }

        if(itsUsePreconGridder && itsPreconGridders.count(imageName)==0) {
           // preconditioning of higher order terms is set from term 0
           bool isMFS = (imageName.find(".taylor.") != std::string::npos);
           bool isTT0 = (imageName.find(".taylor.0") != std::string::npos);
           if (isTT0 || !isMFS) {
             // Should this be a clone of the psf or the image?
             //itsPreconGridders[imageName] = itsGridder->clone();
             itsPreconGridders[imageName] = itsPSFGridders[imageName]->clone();
             // technically, cloning of the PSF gridder should copy the weight accessor (by reference) if set
             boost::shared_ptr<TableVisGridder> tvgPrecon = boost::dynamic_pointer_cast<TableVisGridder>(itsPreconGridders[imageName]);
             if (tvgPrecon && itsReuseGrids) {
                 tvgPrecon->doClearGrid(false);
                 ASKAPLOG_DEBUG_STR(logger, "Setting clear grid to false, reuse PCF grid for "<<imageName);
             }
           }
        } else {
            // reuse the grid if gridder has been configured to allow this
            if (itsReuseGrids) {
                ASKAPLOG_DEBUG_STR(logger, "Will reuse PCF grid for "<<imageName);
            }
        }

        if (itsCoordSystems.count(imageName) == 0) {
          itsCoordSystems[imageName] = SynthesisParamsHelper::coordinateSystem(parameters(),imageName);
        }
      }
      // Now we initialise appropriately
      ASKAPLOG_DEBUG_STR(logger, "Initialising for model degridding and residual gridding" );
      if (completions.size() == 0) {
          ASKAPLOG_WARN_STR(logger,
              "Found no free image parameters, this rank will not contribute usefully to normal equations");
      }
      bool somethingHasToBeDegridded = false;
      for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        string imageName("image"+(*it));
        ASKAPLOG_DEBUG_STR(logger, "Initialising for " << imageName);
        const Axes axes(parameters().axes(imageName));
        casacore::Array<imtype> imagePixels;
        #ifdef ASKAP_FLOAT_IMAGE_PARAMS
            imagePixels = parameters().valueF(imageName);
        #else
            imagePixels = parameters().valueT(imageName);
        #endif
        const casacore::IPosition imageShape(imagePixels.shape());
        /// First the model
        itsModelGridders[imageName]->customiseForContext(*it);
        itsModelGridders[imageName]->initialiseDegrid(axes, imagePixels);
        if (!itsModelGridders[imageName]->isModelEmpty()) {
            somethingHasToBeDegridded = true;
        }
        /// Now the residual images, dopsf=false, dopcf=false
        itsResidualGridders[imageName]->customiseForContext(*it);
        itsResidualGridders[imageName]->initialiseGrid(axes, imageShape, false);
        // DDCALTAG
        if (ddCal) {
            const boost::shared_ptr<TableVisGridder const> &tmGridder =
                boost::dynamic_pointer_cast<TableVisGridder const>(itsModelGridders[imageName]);
            const boost::shared_ptr<TableVisGridder const> &trGridder =
                boost::dynamic_pointer_cast<TableVisGridder const>(itsResidualGridders[imageName]);
            ASKAPCHECK(tmGridder!=0 && trGridder!=0,
                "Wrong gridder type for DDCalibration in ImageFFTEquation::calcImagingEquations()");
            const int index = itsCalDirMap[ImageParamsHelper(imageName).name()];
            ASKAPCHECK(index >= 0, "Invalid source index for DDCAL");
            ASKAPLOG_DEBUG_STR(logger, "Setting DDCal source index for "<<imageName<<" gridders to " << index);
            tmGridder->setSourceIndex(index);
            trGridder->setSourceIndex(index);
        }
        // and PSF gridders, dopsf=true, dopcf=false
        if (!itsReuseGrid[imageName]) {
            itsPSFGridders[imageName]->customiseForContext(*it);
            itsPSFGridders[imageName]->initialiseGrid(axes, imageShape, true);
            // and PCF gridders, dopsf=false, dopcf=true
            if (itsUsePreconGridder && (itsPreconGridders.count(imageName)>0)) {
                itsPreconGridders[imageName]->customiseForContext(*it);
                itsPreconGridders[imageName]->initialiseGrid(axes, imageShape, false, true);
            }
        }
      }
      // synchronise emtpy flag across multiple ranks if necessary
      if (itsVisUpdateObject) {
          itsVisUpdateObject->aggregateFlag(somethingHasToBeDegridded);
      }

// DAM TRADITIONAL
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// in the case that different completions are different Taylor terms, really only want the first one.
//  - so remove the loop and set imageName to "image"+completions.begin()?
//  - are there other things that this will break? Multiple output images with different dir, cellsize, etc.
//  - won't work in continuum imaging, unless combo is done with combinechannels on a single worker
//  - won't work with Taylor terms
//
// I think these are regenerated each major cycles. Need to stop that.
//
      boost::shared_ptr<BoxVisGridder> uvSamplingGridder(new BoxVisGridder);
      if (getRobustness()) {
          // do an initial pass over the entire dataset to generate the uv sampling function
          ASKAPLOG_INFO_STR(logger, "Traditional Briggs weighting with robustness = "<<*getRobustness());

          // basing the weights on the first image name...
          string imageName("image"+(*completions.begin()));
          const Axes axes(parameters().axes(imageName));
          // does this need to be a copy?
          casacore::Array<imtype> imagePixels;
          #ifdef ASKAP_FLOAT_IMAGE_PARAMS
              imagePixels = parameters().valueF(imageName);
          #else
              imagePixels = parameters().value(imageName).copy();
          #endif
          const casa::IPosition imageShape(imagePixels.shape());
          uvSamplingGridder->customiseForContext(*completions.begin());
          uvSamplingGridder->initialiseGrid(axes, imageShape, true);

          // Do an initial loop through all the data to set up weights
          ASKAPLOG_DEBUG_STR(logger, "Initial gridding of uv sampling function" );
          for (itsIdi.init();itsIdi.hasMore();itsIdi.next()) {
              MemBufferDataAccessor accBuffer(*itsIdi);
              uvSamplingGridder->grid(accBuffer);
          }

          // add conjugate points to ensure weights aren't split between positive and negative frequencies
          ASKAPLOG_DEBUG_STR(logger, "Adding conjugate visibilties before calculating robust weights" );
          uvSamplingGridder->addConjugates();

          // convert the uv sampling function to robust weights
          ASKAPLOG_DEBUG_STR(logger, "Convert to robust weights" );
          uvSamplingGridder->setRobustness(*getRobustness());
      }
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      // Get access to explicit calibration operations (predict/correct)
      boost::shared_ptr<ICalibrationApplicator> calME;
      if (ddCal) {
        auto calIter = itsIdi.dynamicCast<CalibrationIterator>();
        if (calIter) {
            calME = calIter->calApplicator();
        }
      }

      // Now we loop through all the data
      ASKAPLOG_INFO_STR(logger, "Starting degridding model and gridding residuals" );
      size_t counterGrid = 0, counterDegrid = 0;
      bool once = true;
      for (itsIdi.init();itsIdi.hasMore();itsIdi.next())
      {
        // buffer-accessor, used as a replacement for proper buffers held in the subtable
        // effectively, an array with the same shape as the visibility cube is held by this class

        // We need the buffer to have both DDCal and NoiseAndFlag features
        DDCalOnDemandNoiseAndFlagDA accBuffer(*itsIdi);
        if (ddCal) {
            // we need a buffer with space for multiple directions
            accBuffer.setNDir(itsNDir);
        }

        // Accumulate model visibility for all models
        // different directions go into separate parts of the buffer
        accBuffer.rwVisibility().set(0.0);
        if (somethingHasToBeDegridded) {
            for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();++it) {
                 const std::string imageName("image"+(*it));
                 const std::map<std::string, IVisGridder::ShPtr>::iterator grdIt = itsModelGridders.find(imageName);
                 ASKAPDEBUGASSERT(grdIt != itsModelGridders.end());
                 const IVisGridder::ShPtr degridder = grdIt->second;
                 ASKAPDEBUGASSERT(degridder);
                 if (!degridder->isModelEmpty()) {
                     degridder->degrid(accBuffer);
                     counterDegrid+=accBuffer.nRow();
                 }
            }
            if (calME) {
                // corrupt the model in ddcal case - each DD section with its own calibration
                calME->predict(accBuffer);
            }
            // optional aggregation of visibilities in the case of distributed model
            // somethingHasToBeDegridded is supposed to have consistent value across all participating ranks
            if (itsVisUpdateObject) {
                itsVisUpdateObject->update(accBuffer.rwVisibility());
            }
            //
        }
        if (!ddCal) {
            // calculate residual visibilities
            accBuffer.rwVisibility() -= itsIdi->visibility();
            accBuffer.rwVisibility() *= float(-1.);
        } else {
            // calculate residual visibilities - subtract all models and calibrate for each direction
            DDCalBufferDataAccessor residBuffer(*itsIdi);
            residBuffer.setNDir(itsNDir);
            casacore::Cube<casacore::Complex> & vis = residBuffer.rwVisibility();
            const casacore::Cube<casacore::Complex> & model = accBuffer.visibility();
            ASKAPCHECK(vis.shape()==model.shape(),"Mismatch in shape between residual and model visibilities");
            const casacore::Slice all;
            const auto nrow = residBuffer.nRow();
            const casacore::Slice row0Slice(0, nrow);
            casacore::Cube<casacore::Complex> vis0(vis(all,all,row0Slice));
            vis0 = itsIdi->visibility();
            // loop over directions - subtract all models from first block of visibilities
            if (somethingHasToBeDegridded) {
                for (int dir = 0; dir < itsNDir; dir++) {
                    const casacore::Slice rowSlice(dir * nrow, nrow);
                    vis0 -= model(all,all,rowSlice);
                }
            }
            // replicate result to other directions
            for (int dir = 1; dir < itsNDir; dir++) {
                const casacore::Slice rowSlice(dir * nrow, nrow);
                vis(all,all,rowSlice) = vis0;
            }
            // Now assign to output buffer
            accBuffer.rwVisibility() = vis;
            // and calibrate for each direction if calibrating
            if (calME) {
                calME->correct(accBuffer);
            }
        }

        /// Now we can calculate the residual visibility and image
        size_t tempCounter = 0;

// DAM TRADITIONAL
// setWeights changes the visibility noise weights (or, rather, those cached in accBuffer),
// This needs to be redone each major cycle, because accBuffer is reset each time (new residual vis)
        if (getRobustness()) {
            uvSamplingGridder->setWeights(accBuffer);
        }

        for (size_t i = 0; i<completions.size(); ++i) {
            const string imageName("image"+completions[i]);
            if (parameters().isFree(imageName)) {
                itsResidualGridders[imageName]->grid(accBuffer);
                if (!itsReuseGrid[imageName]) {
                    itsPSFGridders[imageName]->grid(accBuffer);
                    if (itsUsePreconGridder && (itsPreconGridders.count(imageName)>0)) {
                        itsPreconGridders[imageName]->grid(accBuffer);
                    }
                } else if (once) {
                    ASKAPLOG_DEBUG_STR(logger, "Skipped gridding for PSF/PCF: reuse grid for "<<imageName);
                    once = false;
                }
                tempCounter += accBuffer.nRow();
            }
        }
        counterGrid += tempCounter;
      }
      ASKAPLOG_DEBUG_STR(logger, "Finished degridding model and gridding residuals" );
      ASKAPLOG_DEBUG_STR(logger, "Number of accessor rows iterated through is "<<counterGrid<<" (gridding) and "<<
                        counterDegrid<<" (degridding)");

      // We have looped over all the data, so now we have to complete the
      // transforms and fill in the normal equations with the results from the
      // residual gridders
      if (itsUsePreconGridder) {
        ASKAPLOG_DEBUG_STR(logger,
            "Adding residual image, PSF, preconditioner function and weights image to the normal equations" );
      } else {
        ASKAPLOG_DEBUG_STR(logger,
            "Adding residual image, PSF and weights image to the normal equations" );
      }
      for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();++it)
      {
        const string imageName("image"+(*it));
        const casacore::IPosition imageShape(parameters().shape(imageName));
        ASKAPLOG_INFO_STR(logger,"Name: " << imageName << " Shape: " << imageShape);
        casacore::Array<imtype> imageDeriv(imageShape);
        casacore::Array<imtype> imagePSF(imageShape);
        casacore::Array<imtype> imageWeight(imageShape);
        itsResidualGridders[imageName]->finaliseGrid(imageDeriv);

        // for debugging/research, store grid prior to FFT
        // boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(itsResidualGridders[imageName]);
        // if (tvg) {
        //    tvg->storeGrid("uvcoverage"+(*it),0);
        //}
        // end debugging code

        itsPSFGridders[imageName]->finaliseGrid(imagePSF);
        itsResidualGridders[imageName]->finaliseWeights(imageWeight);

        itsModelGridders[imageName]->finaliseDegrid();

        /*{
          casacore::Array<imtype> imagePSFWeight(imageShape);
          itsPSFGridders[imageName]->finaliseWeights(imagePSFWeight);
          const double maxPSFWeight = casacore::max(imagePSFWeight);
          if (maxPSFWeight > 0) {
              const double psfScalingFactor = casacore::max(imageWeight)/maxPSFWeight;
              //std::cout<<"psf peak = "<<casacore::max(imagePSF)<<" maxPSFWeight = "<<maxPSFWeight<<" factor="<<
              //     psfScalingFactor<<" maxImageWeight="<<casacore::max(imageWeight)<<std::endl;
              imagePSF *= psfScalingFactor;
              // now psf has the same peak as the weight image
          } else {
             if (maxPSFWeight < 0) {
                 ASKAPTHROW(AskapError, "PSF weight is supposed to be non-negative, you have "<<maxPSFWeight);
             }
             // do nothing for zero weight, it just means that this part of normal equations is empty. However,
             // we may still be able to have data after summing all parts of the NE
          }
        }*/

        casacore::IPosition vecShape(1, imagePSF.nelements());

        casacore::Vector<imtype> imagePreconVec;
        if (itsUsePreconGridder && (itsPreconGridders.count(imageName)>0)) {
          casacore::Array<imtype> imagePrecon(imageShape);
          itsPreconGridders[imageName]->finaliseGrid(imagePrecon);
          imagePreconVec.reference(imagePrecon.reform(vecShape));
        }

        {
          casacore::IPosition reference(4, imageShape(0)/2, imageShape(1)/2, 0, 0);
          casacore::Vector<imtype> imagePSFVec(imagePSF.reform(vecShape));
          casacore::Vector<imtype> imageWeightVec(imageWeight.reform(vecShape));
          casacore::Vector<imtype> imageDerivVec(imageDeriv.reform(vecShape));
          ne.addSlice(imageName, imagePSFVec, imageWeightVec, imagePreconVec,
              imageDerivVec, imageShape, reference,itsCoordSystems[imageName]);
        }
      }
    }

  }

}
