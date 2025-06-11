/// @file
///
/// ImagerParallel: Support for parallel applications using the measurement equation
/// classes.
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
#ifndef ASKAP_SYNTHESIS_IMAGERPARALLEL_H_
#define ASKAP_SYNTHESIS_IMAGERPARALLEL_H_

// ASKAPsoft includes
#include <askap/askapparallel/AskapParallel.h>
#include <askap/scimath/fitting/Solver.h>
#include <Common/ParameterSet.h>

// Local package includes
#include <askap/parallel/AdviseParallel.h>
#include <askap/parallel/MEParallelApp.h>
#include <askap/measurementequation/IMeasurementEquation.h>
#include <askap/calibaccess/ICalSolutionConstSource.h>
#include <askap/dataaccess/IDataSource.h>
#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/gridding/IUVWeightCalculator.h>
#include <askap/gridding/GenericUVWeightBuilder.h>
#include <askap/dataaccess/IConstDataIterator.h>

namespace askap
{
  namespace synthesis
  {
    /// @brief Support for parallel algorithms implementing imaging
    ///
    /// @details Provides imaging using the measurement equation framework.
    ///
    /// The control parameters are specified in a parset file. For example:
    /// @code
    ///  	Cimager.datacolumnset           = DATACOL     # default is DATA
    ///  	Cimager.dataset                 = [data/spw_1/sim.ms]
    ///  	#Feed                           = 5
    ///
    ///  	Cimager.Images.Names                    = [image.i.sim]
    ///  	Cimager.Images.image.i.sim.shape                = [4096,4096]
    ///  	Cimager.Images.image.i.sim.cellsize     = [6arcsec, 6arcsec]
    ///  	Cimager.Images.image.i.sim.frequency    = [1.164e9, 1.420e9]
    ///  	Cimager.Images.image.i.sim.nchan                = 1
    ///  	Cimager.Images.image.i.sim.direction    = [12h30m00.00, -45.00.00.00, J2000]
    ///
    ///  	#Cimager.gridder                          = WProject
    ///  	Cimager.gridder.WProject.wmax            = 8000
    ///  	Cimager.gridder.WProject.nwplanes        = 64
    ///  	Cimager.gridder.WProject.oversample     = 1
    ///  	Cimager.gridder.WProject.cutoff         = 0.001
    ///
    ///  	Cimager.gridder                                 = AntennaIllum
    ///  	Cimager.gridder.AntennaIllum.wmax               = 8000
    ///  	Cimager.gridder.AntennaIllum.nwplanes           = 64
    ///  	Cimager.gridder.AntennaIllum.oversample         = 1
    ///  	Cimager.gridder.AntennaIllum.diameter           = 12m
    ///  	Cimager.gridder.AntennaIllum.blockage           = 1m
    ///  	Cimager.gridder.AntennaIllum.maxfeeds           = 18
    ///  	Cimager.gridder.AntennaIllum.maxsupport         = 256
    ///
    ///  	Cimager.solver                                  = Dirty
    ///  	Cimager.solver.Dirty.cycles                     = 1
    ///  	Cimager.solver.Dirty.threshold                  = 10%
    ///  	Cimager.solver.Dirty.verbose                    = True
    ///
    ///  	Cimager.solver.Clean.algorithm                  = MultiScale
    ///  	Cimager.solver.Clean.niter                      = 1000
    ///  	Cimager.solver.Clean.gain                       = 0.7
    ///  	Cimager.solver.Clean.cycles                     = 5
    ///  	Cimager.solver.Clean.verbose                    = True
    ///
    ///  	Cimager.restore                                 = True
    ///  	Cimager.restore.beam                            = [30arcsec, 30arcsec, 0deg]
    /// @endcode
    /// @ingroup parallel
    class ImagerParallel : public MEParallelApp
    {
  public:

      /// @brief Constructor from ParameterSet
      /// @details The parset is used to construct the internal state. We could
      /// also support construction from a python dictionary (for example).
      /// The command line inputs are needed solely for MPI - currently no
      /// application specific information is passed on the command line.
      /// @param comms communication object
      /// @param parset ParameterSet for inputs
      ImagerParallel(askap::askapparallel::AskapParallel& comms,
          const LOFAR::ParameterSet& parset);

      /// @brief Estimate any appropriate parameters that were not specified in the parset.
      /// @details The Cadvise application is used to fill in missing parameters.
      /// @param comms communication object
      /// @param parset initial ParameterSet
      /// @return updated ParameterSet
      static LOFAR::ParameterSet autoSetParameters(askap::askapparallel::AskapParallel& comms,
          const LOFAR::ParameterSet &parset);

      /// @brief Calculate the normalequations (runs in the prediffers)
      /// @details ImageFFTEquation and the specified gridder (set in the parset
      /// file) are used to calculate the normal equations. The image parameters
      /// are defined in the parset file.
      virtual void calcNE();

      /// @brief Solve the normal equations (runs in the solver)
      /// @details Either a dirty image can be constructed or the
      /// multiscale clean can be used, as specified in the parset file.
      virtual void solveNE();

      /// @brief Write the results (runs in the solver)
      /// @details The model images are written as AIPS++ images. In addition,
      /// the images may be restored using the specified beam.
      /// @param[in] postfix this string is added to the end of each name
      /// (used to separate images at different iterations)
      virtual void writeModel(const std::string &postfix = std::string());

      /// @brief make sensitivity image
      /// @details This is a helper method intended to be called from writeModel. It
      /// converts the given weights image into a sensitivity image and exports it.
      /// This method is intended to be called if itsExportSensitivityImage is true.
      /// @param[in] wtImage weight image parameter name
      void makeSensitivityImage(const std::string &wtImage) const;

      /// @brief Helper method to zero all model images
      /// @details We need this for dirty solver only, as otherwise restored image
      /// (which is crucial for faceting) will be wrong.
      void zeroAllModelImages() const;

      /// @brief test whether to advise on wmax, which can require an extra pass over the data.
      /// @param parset ParameterSet to be updated
      /// @return if advice is needed, returns the name of the gridder. Otherwise, returns an empty string.
      static string wMaxAdviceNeeded(LOFAR::ParameterSet &parset);

      // methods to setup traditional weighting. We probably should have a separate class to manage this, but I (MV) feel that not
      // all the required details are fleshed out yet. Leave this here for now, although it is a bit of the technical debt

      /// @brief factory method creating uv weight calculator based on the parset
      /// @details The main parameter controlling the mode of traditional weighting is
      /// Cimager.uvweight which either can take a keyword describing some special method
      /// of getting the weights (which doesn't require iteration over data), e.g. reading from disk
      /// or a list of "effects" which should be applied to the density of uv samples obtained via
      /// iteration over data. This method acts as a factory for weight calculators (i.e. the second
      /// case with the list of effects) or returns an empty pointer if no iteration over data is required
      /// (i.e. either some special algorithm is in use or there is no uv-weighting)
      /// @param[in] parset configuration parset to use
      /// @return shared pointer to the uv-weight calculator object
      static boost::shared_ptr<IUVWeightCalculator> createUVWeightCalculator(const LOFAR::ParameterSet &parset);

      /// @brief factory method creating uv weight calculator based on the current parset
      /// @details Unlike the static method which gets the parset as the parameter and returns the shared pointer to
      /// the calculator object, this method uses the current parset passed to the imager in the constructor and assigns
      /// the result to itsUVWeightCalculator. There is a bit of the technical debt here and it would probably be better to
      /// factor out this factory into a separate class (and remove this code from the imager). For now, it seems to be the 
      /// quicker way to be able to quiery the traditional weighting setup without configuring the imager.
      /// @note This method updates itsUVWeightCalculator which will be either non-zero shared pointer to the weight
      /// calculator object to be applied to the density of uv samples, or an empty shared pointer which implies that
      /// there is no need obtaining the density because either no traditional weighting is done or
      /// we're using some special algorithm which does not require iteration over data
      void createUVWeightCalculator();

      /// @brief set uv-weight calculator object explicitly
      /// @details This method may be useful, if the weight calculator is created via the static createUVWeightCalculator method.
      /// This ambiguity w.r.t. dealing with weight calculator objects is hopefully temporary (and is indicative of the technical debt
      /// because we probably shouldn't have this functionality as part of the imager and/or have approach with two imagers during accumulation).
      /// @param[in] calc shared pointer to the uv-weight calculator object to set
      void setUVWeightCalculator(const boost::shared_ptr<IUVWeightCalculator> &calc) { itsUVWeightCalculator = calc; }

      /// @brief check if sample density grid needs to be built
      /// @details For now, use the shared pointer carrying uv weight calculator as a flag that we need to
      /// build grid of weights (this is what uv weight calculator works with). It is a bit of the technical debt
      /// to do it this way (as we don't need calculator on the ranks which don't compute weights), but allows us
      /// to hide all parset interpretation to the createUVWeightCalculator factory method (for the price that it
      /// needs to be called for all ranks for which we need the consistent picture). But these calculator objects are
      /// very lightweight, so it is a lesser of possible evils. More pure design would involve interpretation of the parset in
      /// one place handled by a separate class (which we either run on all ranks or just on a master and then distribute).
      /// Therefore, it is expected that createUVWeightCalculator should be called before this method.
      /// @return true if sample density grid needs to be built
      bool inline isSampleDensityGridNeeded() const { return static_cast<bool>(itsUVWeightCalculator); }

      /// @brief replace normal equations with an adapter handling uv weights
      /// @details We use normal equations infrastructure to merge uv weights built in parallel
      /// (via EstimatorAdapter). This method sets up the appropriate builder class based on the
      /// parameters in the parset, wraps it in the adapter and assigns to itsNE in the base class.
      /// @note An instance of appropriate normal equations class should be setup before normal major
      /// cycles can resume. It can be done with a call to recereateNormalEquations or calcNE
      void setupUVWeightBuilder();

      /// @brief compute uv weights using data stored via adapter in the normal equations
      /// @details This method gets access to the uv-weight builder handled via EstimatorAdapter and
      /// stored instead of normal equations by shared pointer. It then runs the finalisation step using
      /// the stored shared pointer to the uv-weight calculator object. The resulting weight is added as
      /// a parameter to the existing model.
      /// @note This method assumes that it is called from the right place (i.e. on the correct rank) and
      /// all merging/reduction has already been done. In other words, it is agnostic of the parallelism.
      void computeUVWeights() const;


      /// @brief iterate over given data and accumulate samples for uv weights
      /// @details This method is used to build the sample density in the uv-plane via the appropriate gridder
      /// and weight builder class. It expects the builder already setup and accessible via the normal equations
      /// shared pointer. The data iterator to work with is passed as a parameter. The image details are extracted from
      /// the model (to initialise sample grid).
      /// @param[in] iter shared pointer to the iterator to use (note it is advanced by this method to iterate over
      /// all available data)
      void accumulateUVWeights(const boost::shared_ptr<accessors::IConstDataIterator> &iter) const;

      /// @brief recreate imaging normal equations object
      /// @details If sample density grid is built, normal equations are setup with
      /// an adapter which is an incompatible type. This method resets the object to a pristine
      /// state similar to that before the first major cycle if no traditional weighting is done.
      void recreateNormalEquations();

  protected:

      /// @brief make calibration iterator if necessary, otherwise return unchanged interator
      /// @details This method wraps the iterator passed as the input into into a calibration iterator adapter
      /// if calibration is to be performed (i.e. if solution source is defined).
      /// @param[in] origIt original iterator to uncalibrated data
      /// @return shared pointer to the data iterator with on-the-fly calibration application, if necessary
      /// or the original iterator otherwise
      accessors::IDataSharedIter makeCalibratedDataIteratorIfNeeded(const accessors::IDataSharedIter &origIt) const;

      /// @brief helper method to extract weight builder object out of normal equations
      /// @details For traditional weighting with distributed data we use normal equation merging
      /// mechanism to do the gather operation (and EstimatorAdapter). This method does required
      /// casts and checks to get the required shared pointer (which is guaranteed to be non-empty)
      /// @return shared pointer to the uv-weight builder object stored in the current normal equations
      boost::shared_ptr<GenericUVWeightBuilder> getUVWeightBuilder() const;

      /// @brief a helper method to extract peak residual
      /// @details This object actually manipulates with the normal equations. We need
      /// to be able to stop iterations on the basis of maximum residual, which is a
      /// data vector of the normal equations. This helper method is designed to extract
      /// peak residual. It is then added to a model as a parameter (the model is
      /// shipped around).
      /// @return absolute value peak of the residuals corresponding to the current normal
      /// equations
      double getPeakResidual() const;

      /// @brief helper method to indentify model parameters to broadcast
      /// @details We use itsModel to buffer some derived images like psf, weights, etc
      /// which are not required for prediffers. It just wastes memory and CPU time if
      /// we broadcast them. At the same time, some auxilliary parameters like peak
      /// residual value need to be broadcast (so the major cycle can terminate in workers).
      /// This method returns the vector with all parameters to be broadcast. By default
      /// it returns all parameter names, so it is overridden here to broadcast only
      /// model images and the peak_residual metadata.
      /// @return a vector with parameters to broadcast
      virtual std::vector<std::string> parametersToBroadcast() const;

      /// @brief helper method to access the solution source
      /// @details This object is initialised by workers. It knows how to
      /// retrieve calibration solutions (from a parset file, casa table or a database).
      /// Uninitialised shared pointer means that no calibration is done and the data are
      /// imaged as they are.
      /// @return the shared_ptr
      boost::shared_ptr<accessors::ICalSolutionConstSource> getSolutionSource() const {
          return itsSolutionSource; }

      /// @brief obtain measurement equation cast to ImageFFTEquation
      /// @details This helper method encapsulates operations common to a number of methods of this and derived classes to obtain the
      /// current measurement equation with the original type as created (i.e. ImageFFTEquation) and
      /// does the appropriate checks (so the return is guaranteed to be a non-null shared pointer).
      /// @return shared pointer of the appropriate type to the current measurement equation
      boost::shared_ptr<ImageFFTEquation> getMeasurementEquation() const;

      /// @brief make data iterator
      /// @details This helper method makes an iterator based on the configuration in the current parset and
      /// given data source object
      /// @param[in] ds datasource object to use
      /// @return shared pointer to the iterator over data
      accessors::IDataSharedIter makeDataIterator(const accessors::IDataSource &ds) const;

      /// @brief map of image names to calibration directions
      /// @return map  the map of image names to calibration directions
      std::map<std::string, int> calDirMap() {return itsCalDirMap;}

  private:

      /// @brief check whether any preProcess advice is needed.
      /// @param parset initial ParameterSet
      /// @return true if there are parameters that need to be set.
      static bool checkForMissingParameters(const LOFAR::ParameterSet &parset);

      /// @brief copy any required Cimager parameters to AdviseParallel parameters.
      /// @param parset ParameterSet to be updated
      static void addAdviseParameters(LOFAR::ParameterSet &parset);

      /// @brief add missing parameters based on advice from AdviseParallel.
      /// @param advice AdviseParallel statistics
      /// @param parset ParameterSet to be updated
      static void addMissingParameters(const VisMetaDataStats &advice, LOFAR::ParameterSet &parset);

      /// @brief remove any added AdviseParallel parameters.
      /// @param parset updated ParameterSet
      static void cleanUpAdviseParameters(LOFAR::ParameterSet &parset);

      /// Calculate normal equations for one data set
      /// @param ms Name of data set
      /// @param discard Discard old equation?
      void calcOne(const string& dataset, bool discard=false);

      /// @brief Do we want a restored image?
      bool itsRestore;

      /// @brief write the first pass restored images if there is a second pass
      bool itsWriteFirstRestore;

      /// @brief Do we want a residual image
      bool itsWriteResidual;

      /// @brief write 'raw', unnormalised, natural weight psf
      bool itsWritePsfRaw;

      /// @brief write normalised, preconditioned psf
      bool itsWritePsfImage;

      /// @brief write weights image
      bool itsWriteWtImage;

      /// @brief write a weights log
      bool itsWriteWtLog;

      /// @brief write out the mask image
      bool itsWriteMaskImage;

      /// @brief write out the scalemask image
      bool itsWriteScaleMask;

      /// @brief write out the (clean) model image
      bool itsWriteModelImage;

      /// @brief export sensitivity image?
      /// @details If true, the weight image is converted to sensitivity image and
      /// exported in addition to the original weight image.
      bool itsWriteSensitivityImage;

      /// @brief fraction of peak sensitivity to crop on
      /// @details We normally don't want to divide by small numbers when calculating
      /// sensitivity images. This field gives the fraction of the maximum weight
      /// below which the sensitivity image will be set to 0.
      double itsSensitivityCutoff;

      /// @brief write out the gridded data, pcf and psf
      bool itsWriteGrids;

      /// @brief write out multiple images if there are more than 1
      bool itsWriteMultiple;

      /// @brief write out multiple (clean) model images if there are more than 1
      bool itsWriteMultipleModels;

      /// @brief number of directions to use for calibration
      int itsNDir;

      // @brief map of image names to calibration directions
      std::map<std::string, int> itsCalDirMap;

      /// @brief name of first image (after removing image prefix)
      std::string itsFirstImageName;

      /// @brief solution source to get calibration data from
      /// @details This object is initialised by workers. It knows how to
      /// retrieve calibration solutions (from a parset file, casa table or a database).
      /// Uninitialised shared pointer means that no calibration is done and the data are
      /// imaged as they are.
      boost::shared_ptr<accessors::ICalSolutionConstSource> itsSolutionSource;

      /// @brief void measurement equation
      /// @details Does nothing, just returns calls to predict and
      /// calcNormalEquations. This shared pointer is initialized at the
      /// first use and then passed as the perfect MeasurementEquation, if
      /// required. We could have created a brand new object each time.
      boost::shared_ptr<IMeasurementEquation> itsVoidME;

      /// @brief object calculating uv weights from sample density
      /// @details This data member stores shared pointer to the calculator object
      /// returned by the factory method. Note, technically this object only needs
      /// to be created if we're going to use it (i.e. run computeUVWeights method on the
      /// given rank), but for now it is handy to use the presence of this object as a flag
      /// that we're doing traditional weighting and need to compute sample density. It is
      /// a bit of the technical debt, but a small overhead (setting up calculators is cheap).
      /// On the positive side, it allows us to avoid code duplication in parset interpretation.
      /// (in the future, we can have a specialised class containing all the required operations)
      boost::shared_ptr<IUVWeightCalculator> itsUVWeightCalculator;

    };

  }
}
#endif
