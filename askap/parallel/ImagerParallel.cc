/// @file ImagerParallel.cc
///
/// Performs synthesis imaging from a data source, using any of a number of
/// image solvers. Can run in serial or parallel (MPI) mode.
///
/// The data are accessed from the DataSource. This is and will probably remain
/// disk based. The images are kept purely in memory until the end.
///
/// Control parameters are passed in from a LOFAR ParameterSet file.
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

#include <askap/parallel/ImagerParallel.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

#include <askap/askap/AskapError.h>
// need it just for null deleter
#include <askap/askap/AskapUtil.h>

#include <askap/askapparallel/AskapParallel.h>
#include <askap/dataaccess/DataAccessError.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/imageaccess/WeightsLog.h>


#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/ImageRestoreSolver.h>
#include <askap/measurementequation/MEParsetInterface.h>
#include <askap/measurementequation/CalibrationIterator.h>
#include <askap/measurementequation/NoXPolGain.h>
#include <askap/measurementequation/ImageParamsHelper.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>
#include <askap/scimath/fft/FFTUtils.h>

#include <askap/measurementequation/ImageSolverFactory.h>
#include <askap/measurementequation/ImageCleaningSolver.h>
#include <askap/calibaccess/CalibAccessFactory.h>
#include <askap/measurementequation/CalibrationApplicatorME.h>
#include <askap/profile/AskapProfiler.h>
#include <askap/parallel/GroupVisAggregator.h>
#include <askap/parallel/AdviseParallel.h>
#include <askap/utils/StatsAndMask.h>

// this is used for uv-weight calculators factory, to be removed
// when/if we factor this out into a separate class
#include <askap/gridding/CompositeUVWeightCalculator.h>
#include <askap/gridding/RobustUVWeightCalculator.h>
#include <askap/gridding/ReciprocalUVWeightCalculator.h>
#include <askap/gridding/ConjugatesAdderFFT.h>

#include <askap/gridding/GenericUVWeightBuilder.h>
#include <askap/scimath/utils/EstimatorAdapter.h>
#include <askap/gridding/UVWeightParamsHelper.h>
#include <askap/gridding/UVWeightGridder.h>


#include <casacore/casa/aips.h>
#include <casacore/casa/OS/Timer.h>

#include <Common/ParameterSet.h>

#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>

#include <boost/optional.hpp>

using namespace askap;
using namespace askap::scimath;
using namespace askap::synthesis;
using namespace askap::accessors;
using namespace askap::askapparallel;

namespace askap
{
  namespace synthesis
  {
    using utility::toString;
    using askap::operator<<;

    ImagerParallel::ImagerParallel(askap::askapparallel::AskapParallel& comms,
        const LOFAR::ParameterSet& parset) :
      MEParallelApp(comms,parset,true)
    {
      if (itsComms.isMaster())
      {
        itsRestore=parset.getBool("restore", false); // do restore and write restored image
        itsWriteFirstRestore=parset.getBool("write.firstrestore",false); // write first restore products if alt specified
        itsWriteResidual=parset.getBool("residuals",false); // write residual image
        itsWriteResidual=parset.getBool("write.residualimage",itsWriteResidual); // alternative param name
        itsWritePsfRaw = parset.getBool("write.psfrawimage", false); // write unnormalised, natural wt psf
        itsWritePsfImage = parset.getBool("write.psfimage", true); // write normalised, preconditioned psf
        itsWriteWtLog = parset.getBool("write.weightslog", false); // write weights log file
        itsWriteWtImage = parset.getBool("write.weightsimage", false); // write weights image
        itsWriteMaskImage = parset.getBool("write.maskimage", false); // write mask image
        itsWriteModelImage = parset.getBool("write.modelimage", !itsRestore); // write clean model
        itsWriteSensitivityImage = parset.getBool("write.sensitivityimage", false);
        itsWriteGrids = parset.getBool("dumpgrids", false); // write (dump) the gridded data, psf and pcf
        itsWriteGrids = parset.getBool("write.grids",itsWriteGrids); // new name
        itsWriteMultiple = parset.getBool("write.multiple", false); // Only write out the first image by default
        itsWriteMultipleModels = parset.getBool("write.multiplemodels",true); // write all model images
        itsWriteScaleMask = parset.getBool("write.scalemask",false); // write out the scale mask(s)
        itsSensitivityCutoff = parset.getDouble("sensitivityimage.cutoff", 0.01);
        bool reuseModel = parset.getBool("Images.reuse", false); // continue solving from existing model

        if (itsWriteSensitivityImage) {
            ASKAPLOG_INFO_STR(logger,
               "Theoretical sensitivity images will be generated in addition to weights images, cutoff="<<
                itsSensitivityCutoff);
        }

        ASKAPCHECK(itsModel, "itsModel is supposed to be initialized at this stage");

        if (reuseModel) {
            ASKAPLOG_INFO_STR(logger, "Reusing model images stored on disk");
            SynthesisParamsHelper::loadImages(itsModel,parset.makeSubset("Images."));
        } else {
            ASKAPLOG_INFO_STR(logger, "Initializing the model images");

            /// Create the specified images from the definition in the
            /// parameter set. We can solve for any number of images
            /// at once (but you may/will run out of memory!)
            SynthesisParamsHelper::setUpImages(itsModel,
                                      parset.makeSubset("Images."));
        }
        ASKAPCHECK(parset.getStringVector("Images.Names").size()>0,"No image names specified");
        itsFirstImageName = parset.getStringVector("Images.Names")[0];
        ASKAPCHECK(itsFirstImageName.find("image")==0,"Image names should start with 'image'");
        itsFirstImageName = itsFirstImageName.substr(5);

        /// Create the solver from the parameterset definition
        LOFAR::ParameterSet tmpset = parset.makeSubset("");
        if (itsWriteScaleMask) {
            tmpset.replace(LOFAR::KVpair("Solver.Clean.writescalemask",true));
        }
        itsSolver = ImageSolverFactory::make(tmpset);
        ASKAPCHECK(itsSolver, "Solver not defined correctly");
      }
      if (itsComms.isWorker())
      {
        const bool doCalib = parset.getBool("calibrate",false);
        if (doCalib) {
            ASKAPCHECK(!parset.isDefined("gainsfile"), "Deprecated 'gainsfile' keyword is found together with calibrate=true, please remove it");
            // setup solution source from the parset directly using the factory
            itsSolutionSource = CalibAccessFactory::roCalSolutionSource(parset);
            ASKAPASSERT(itsSolutionSource);
        } else if (parset.isDefined("gainsfile")) {
            // temporary code to support deprecated gainsfile parameter (used in various tests)
            const std::string gainsFile = parset.getString("gainsfile");
            ASKAPLOG_WARN_STR(logger, "The parset has deprecated 'gainsfile = "<<gainsFile<<
                       "' parameter defined. Use calibrate=true and calibaccess.parset=filename instead");
            ASKAPCHECK(!parset.isDefined("calibaccess") && !parset.isDefined("calibaccess.parset"),
                       "Detected a mix of deprecated and new parameters defining calibration solution access, remove gainsfile!");
            LOFAR::ParameterSet tmpParset(parset);
            tmpParset.add("calibaccess.parset", gainsFile);
            tmpParset.add("calibaccess", "parset");
            // setup solution source from the temporary parset
            itsSolutionSource = CalibAccessFactory::roCalSolutionSource(tmpParset);
            ASKAPASSERT(itsSolutionSource);
        }
        if (itsSolutionSource) {
            ASKAPLOG_INFO_STR(logger, "Data will be calibrated before imaging");
            const std::vector<std::string> names = parset.getStringVector("Images.Names");
            const bool doDDCal = parset.getBool("calibrate.directiondependent",false) &&
                (names.size() > 1);
            if (doDDCal) {
                // Now map names to directions
                for (int i=0; i<names.size(); i++) {
                    itsCalDirMap[names[i]] = i;
                }
                ASKAPLOG_INFO_STR(logger, "Using "<<names.size()<< " directions for DDCAL");
            }
        } else {
            ASKAPLOG_INFO_STR(logger, "No calibration will be performed");
        }
      }
    }

    /// Estimate any appropriate parameters that were not specified in the parset
    LOFAR::ParameterSet ImagerParallel::autoSetParameters(askap::askapparallel::AskapParallel& comms,
        const LOFAR::ParameterSet &parset) {

      // make a copy of the initial parset to fill with extra parameters and return.
      LOFAR::ParameterSet fullset = parset;

      // if advice is needed, use the AdviseParallel class to estimate missing parameters.
      if ( checkForMissingParameters(fullset) ) {
          ASKAPLOG_INFO_STR(logger, "Using the advise application to fill unset parameters.");

          addAdviseParameters(fullset); // copy any required Cimager parameters to AdviseParallel parameters.

          AdviseParallel cadvise(comms, fullset);
          cadvise.estimate(); // generate the statistics

          const VisMetaDataStats &advice = cadvise.estimator();
          addMissingParameters(advice, fullset); /// add missing parameters based on advice from AdviseParallel.

          cleanUpAdviseParameters(fullset); // remove any added AdviseParallel parameters.
      }

      return fullset;
    }

    /// check whether any preProcess advice is needed.
    bool ImagerParallel::checkForMissingParameters(const LOFAR::ParameterSet &parset) {

      // set to true if any parameters are missing
      bool paramTest = false;

      // these parameters can be set globally or individually
      bool cellsizeNeeded = false;
      bool shapeNeeded = false;

      // test for missing image-specific parameters:
      const std::vector<std::string> imageNames = parset.getStringVector("Images.Names", false);
      for (size_t img = 0; img < imageNames.size(); ++img) {
          if ( !parset.isDefined("Images."+imageNames[img]+".cellsize") ) cellsizeNeeded = true;
          if ( !parset.isDefined("Images."+imageNames[img]+".shape") ) shapeNeeded = true;
          if ( !parset.isDefined("Images."+imageNames[img]+".nchan") ) {
              paramTest = true;
              break;
          } else if ( !parset.isDefined("Images."+imageNames[img]+".frequency") ) {
              paramTest = true;
              break;
          } else if ( !parset.isDefined("Images."+imageNames[img]+".direction") ) {
              paramTest = true;
              break;
          }
      }
      if (paramTest) return paramTest;

      // test for general missing parameters:
      if ( !parset.isDefined("nUVWMachines") ) {
          paramTest = true;
      } else if ( cellsizeNeeded && !parset.isDefined("Images.cellsize") ) {
          paramTest = true;
      } else if ( shapeNeeded && !parset.isDefined("Images.shape") ) {
          paramTest = true;
      }

      return paramTest;
    }

    /// copy any required Cimager parameters to AdviseParallel parameters.
    void ImagerParallel::addAdviseParameters(LOFAR::ParameterSet &parset) {

      // Add Cimager parameters in a form that AdviseParallel is expecting.
      // In case any of these parameters ever become 'Cimager.' parameters, add an extra parameter
      // for each with suffix ".advised" so they can be correctly identified during clean up.
      // Also, use variable 'tangentDefined' to avoid such a problem here.

      string param;

      // Use the tangent of the first image as the advise tangent (defaults to direction if no tangent is given).
      const std::vector<std::string> imageNames = parset.getStringVector("Images.Names", false);
      if (imageNames.size() > 1) {
          ASKAPLOG_WARN_STR(logger, "  Multiple images. Only the first will be inspected for tangent information.");
      }
      bool tangentDefined = parset.isDefined("tangent");
      if (parset.isDefined("Images."+imageNames[0]+".tangent") && !tangentDefined) {
          param = "Images."+imageNames[0]+".tangent";
          string tangent = parset.getString(param);
          if (tangent.find("J2000")==string::npos) { // Advise will exit if epoch is not J2000
              ASKAPLOG_WARN_STR(logger, "  AdviseParallel only accepts J2000 at the moment. Not using " <<
                    param<<" in Advise.");
          } else {
              ASKAPLOG_INFO_STR(logger, "  Adding tangent for advise: " << tangent);
              parset.add("tangent", tangent);
              parset.add("tangent.advised", tangent);
              tangentDefined = true;
          }
      }
      else if ( parset.isDefined("Images."+imageNames[0]+".direction") && !tangentDefined) {
          param = "Images."+imageNames[0]+".direction";
          string tangent = parset.getString(param);
          if (tangent.find("J2000")==string::npos) { // Advise will exit if epoch is not J2000
              ASKAPLOG_WARN_STR(logger, "  AdviseParallel only accepts J2000 at the moment. Not using " <<
                    param<<" in Advise.");
          } else {
              ASKAPLOG_INFO_STR(logger, "  Adding tangent for advise: " << tangent);
              parset.add("tangent", tangent);
              parset.add("tangent.advised", tangent);
              tangentDefined = true;
          }
      }

      const string wMaxGridder = wMaxAdviceNeeded(parset);
      // Use the snapshotimaging wtolerance as the advise wtolerance. But only if wmax is needed.
      bool snapShot = parset.getBool("gridder.snapshotimaging",false);
      param = "gridder.snapshotimaging.wtolerance";
      if (snapShot && parset.isDefined(param) && !parset.isDefined("wtolerance") && (wMaxGridder!="")) {
          string wtolerance = parset.getString(param);
          ASKAPLOG_INFO_STR(logger, "  Adding wtolerance for advise: " << wtolerance);
          parset.add("wtolerance", wtolerance);
          parset.add("wtolerance.advised", wtolerance);
      }
      if (wMaxGridder!="") {
          param = "gridder."+wMaxGridder+".wpercentile";
          if (parset.isDefined(param)) {
              string wpercentile = parset.get(param);
              parset.add("wpercentile", wpercentile);
              parset.add("wpercentile.advised", wpercentile);
          }
        }
    }

    /// test whether to advise on wmax, which can require an extra pass over the data.
    string ImagerParallel::wMaxAdviceNeeded(LOFAR::ParameterSet &parset) {

      // make a vector containing all gridders that require the wmax parameter
      std::vector<std::string> wGridders;
      wGridders.push_back("MPIWProject");
      wGridders.push_back("WProject");
      wGridders.push_back("WStack");
      wGridders.push_back("AWProject");
      wGridders.push_back("AProjectWStack");
      for(std::vector<std::string>::const_iterator i = wGridders.begin(); i != wGridders.end(); ++i) {
          if ((parset.getString("gridder")==*i) && !parset.isDefined("gridder."+*i+".wmax") ) {
              return *i;
          }
      }

      return "";

    }

    /// remove any added AdviseParallel parameters.
    void ImagerParallel::cleanUpAdviseParameters(LOFAR::ParameterSet &parset) {

      if (parset.isDefined("tangent.advised")) {
          parset.remove("tangent");
          parset.remove("tangent.advised");
      }

      if (parset.isDefined("wtolerance.advised")) {
          parset.remove("wtolerance");
          parset.remove("wtolerance.advised");
      }
      if (parset.isDefined("wpercentile.advised")) {
          parset.remove("wpercentile");
          parset.remove("wpercentile.advised");
      }

    }

    /// add missing parameters based on advice from AdviseParallel.
    /// if adding any more parameters here, also add them to checkforMissingParameters(). Parameters are:
    /// - advice.nAntennas()
    /// - advice.nVis()
    /// - advice.maxU() // [wavelengths]
    /// - advice.maxV() // [wavelengths]
    /// - advice.maxW() // [wavelengths]
    /// - if (wTolerance >= 0.) advice.maxResidualW() // [wavelengths]
    /// - advice.nBeams()
    /// - advice.maxOffsets() // Largest beam offset: maxOffsets().first, maxOffsets().second [rad]
    /// - advice.squareCellSize(); // [arcsec]
    /// - advice.squareFieldSize(0); // side length out to about the first sidelobe (about the 'average' pointing) [deg]
    /// - advice.squareFieldSize(1); // side length out to about the first sidelobe (about the tangent point) [deg]
    /// - advice.minFreq() // [Hz]
    /// - advice.maxFreq() // [Hz]
    void ImagerParallel::addMissingParameters(const VisMetaDataStats &advice, LOFAR::ParameterSet &parset) {

      string param;
      // these parameters can be set globally or individually
      bool cellsizeNeeded = false;
      bool shapeNeeded = false;
      int nTerms = 1;

      // for image specific parameters, cycle through the images and check whether any advised parameters are undefined.
      const std::vector<std::string> imageNames = parset.getStringVector("Images.Names", false);
      for (size_t img = 0; img < imageNames.size(); ++img) {
          // could set this up to add shape if it's missing but cellsize is set, like below for the global parameters.
          if (!parset.isDefined("Images."+imageNames[img]+".cellsize")) cellsizeNeeded = true;
          if (!parset.isDefined("Images."+imageNames[img]+".shape")) shapeNeeded = true;
          int nChan = 1;
          param = "Images."+imageNames[img]+".nchan"; // if the number of image frequency channels undefined, set to 1.
          if (parset.isDefined(param)) {
              nChan = parset.getInt(param);
          } else {
              std::ostringstream pstr;
              pstr<<nChan;
              ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << ": " << pstr.str().c_str());
              parset.add(param, pstr.str().c_str());
          }
          param = "Images."+imageNames[img]+".frequency"; // if freq is undefined, use the info from the ms.
          if (!parset.isDefined(param)) {
              std::ostringstream pstr;
              if (nChan==1) {
                  const double aveFreq = 0.5*(advice.minFreq()+advice.maxFreq());
                  pstr<<"["<<aveFreq<<","<<aveFreq<<"]";
              } else {
                  pstr<<"["<<advice.minFreq()<<","<<advice.maxFreq()<<"]";
              }
              ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << ": " << pstr.str().c_str());
              parset.add(param, pstr.str().c_str());
          }
          param = "Images."+imageNames[img]+".direction"; // if the image centre is undefined, use the ms phase centre.
          if (!parset.isDefined(param)) {
              std::ostringstream pstr;
              // Only J2000 is implemented at the moment.
              pstr<<"["<<printLon(advice.centre())<<", "<<printLat(advice.centre())<<", J2000]";
              ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param << ": " << pstr.str().c_str());
              parset.add(param, pstr.str().c_str());
          }
          param = "Images."+imageNames[img]+".nterms"; // if nterms is set, store it for later
          if (parset.isDefined(param)) {
              if ((nTerms>1) && (nTerms!=parset.getInt(param))) {
                  ASKAPLOG_WARN_STR(logger, "  Imaging with different nterms may not work");
              }
              nTerms = parset.getInt(param);
          }
      }

      if (nTerms > 1) { // check required MFS parameters
         param = "visweights"; // set to "MFS" if unset and nTerms > 1
         if (!parset.isDefined(param)) {
             std::ostringstream pstr;
             pstr<<"MFS";
             ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
             parset.add(param, pstr.str().c_str());
         }
         param = "visweights.MFS.reffreq"; // set to average frequency if unset and nTerms > 1
         if ((parset.getString("visweights")=="MFS") && !parset.isDefined(param)) {
             std::ostringstream pstr;
             const double aveFreq = 0.5*(advice.minFreq()+advice.maxFreq());
             pstr<<aveFreq;
             ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
             parset.add(param, pstr.str().c_str());
         }
      }

      param = "nUVWMachines"; // if the number of uvw machines is undefined, set it to the number of beams.
      if (!parset.isDefined(param)) {
          std::ostringstream pstr;
          pstr<<advice.nBeams();
          ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
          parset.add(param, pstr.str().c_str());
      }
      std::vector<double> cellSize(2, advice.squareCellSize());
      param = "Images.cellsize"; // if cellsize is undefined, use the advice. Otherwise update cellSize for later.
      if (parset.isDefined(param)) {
          cellSize = SynthesisParamsHelper::convertQuantity(parset.getStringVector(param),"arcsec");
      } else if (cellsizeNeeded) {
          std::ostringstream pstr;
          pstr<<"["<<cellSize[0]<<"arcsec,"<<cellSize[1]<<"arcsec]";
          ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
          parset.add(param, pstr.str().c_str());
      }
      param = "Images.shape"; // if image shape is undefined, use the advice.
      if (shapeNeeded && !parset.isDefined(param)) {
          std::ostringstream pstr;
          const double fieldSize = advice.squareFieldSize(1); // in deg
          const int lSize = scimath::goodFFTSize(int(fieldSize * 3600 / cellSize[0]) + 1);
          const int mSize = scimath::goodFFTSize(int(fieldSize * 3600 / cellSize[1]) + 1);
          pstr<<"["<<lSize<<","<<mSize<<"]";
          ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
          parset.add(param, pstr.str().c_str());
      }
      string gridder = wMaxAdviceNeeded(parset); // returns empty string if wmax is not required.
      if (gridder!="") {
          param = "gridder."+gridder+".wmax"; // if wmax is undefined but needed, use the advice.
          if (!parset.isDefined(param)) {
              string pstr;
              // both should be set if either is, but include the latter to ensure that maxResidualW is generated.
              if (parset.isDefined("gridder.snapshotimaging.wtolerance") && parset.isDefined("wtolerance") ) {
                  pstr = toString(advice.maxResidualW()); // could use parset.getString("gridder.snapshotimaging.wtolerance");
                  ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr);
                  parset.add(param, pstr);
              } else if (parset.isDefined("gridder."+gridder+".wpercentile")){
                  pstr = toString(advice.wPercentile());
                  ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr << " with specified percentile value");
                  parset.add(param, pstr);
                  param = "gridder."+gridder+".wmaxclip";
                  if (!parset.isDefined(param)) parset.add(param, "true");
              } else {
                  pstr = toString(advice.maxW());
                  ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr);
                  parset.add(param, pstr);
              }
          }
      }

      // add Nyquist gridding parameters if needed. Wait until after doing others requiring VisMetaDataStats.
      SynthesisParamsHelper::setNyquistSampling(advice, parset);

    }

    /// @brief make calibration iterator if necessary, otherwise return unchanged interator
    /// @details This method wraps the iterator passed as the input into into a calibration iterator adapter
    /// if calibration is to be performed (i.e. if solution source is defined).
    /// @param[in] origIt original iterator to uncalibrated data
    /// @return shared pointer to the data iterator with on-the-fly calibration application, if necessary
    /// or the original iterator otherwise
    accessors::IDataSharedIter ImagerParallel::makeCalibratedDataIteratorIfNeeded(const accessors::IDataSharedIter &origIt) const
    {
       if (getSolutionSource()) {
            ASKAPLOG_DEBUG_STR(logger, "Calibration will be performed using solution source");
            const boost::shared_ptr<ICalibrationApplicator> calME(new CalibrationApplicatorME(getSolutionSource()));
            // fine tune parameters
            ASKAPDEBUGASSERT(calME);
            calME->scaleNoise(parset().getBool("calibrate.scalenoise",false));
            calME->allowFlag(parset().getBool("calibrate.allowflag",false));
            calME->beamIndependent(parset().getBool("calibrate.ignorebeam", false));
            calME->interpolateTime(parset().getBool("calibrate.interpolatetime",false));

            // calibration iterator to replace the original one for the purpose of measurement equation creation
            const IDataSharedIter calIter(new CalibrationIterator(origIt,calME,itsCalDirMap.size()>0));
            return calIter;
       }
       ASKAPLOG_DEBUG_STR(logger,"Not applying calibration");
       return origIt;
    }

    /// @brief make data iterator
    /// @details This helper method makes an iterator based on the configuration in the current parset and
    /// given data source object
    /// @param[in] ds datasource object to use
    /// @return shared pointer to the iterator over data
    accessors::IDataSharedIter ImagerParallel::makeDataIterator(const accessors::IDataSource &ds) const
    {
        IDataSelectorPtr sel=ds.createSelector();
        sel->chooseCrossCorrelations();
        sel << parset();
        IDataConverterPtr conv=ds.createConverter();
        conv->setFrequencyFrame(getFreqRefFrame(), "Hz");
        conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));
        // ensure that time is counted in seconds since 0 MJD
        conv->setEpochFrame();
        return ds.createIterator(sel, conv);
    }

    void ImagerParallel::calcOne(const string& ms, bool discard)
    {
      ASKAPDEBUGTRACE("ImagerParallel::calcOne");
      casacore::Timer timer;
      timer.mark();
      ASKAPLOG_INFO_STR(logger, "Calculating normal equations for " << ms );
      // First time around we need to generate the equation
      if ((!itsEquation)||discard)
      {
        ASKAPLOG_INFO_STR(logger, "Creating measurement equation" );

        // MEMORY_BUFFERS mode opens the MS readonly
        TableDataSource ds(ms, TableDataSource::MEMORY_BUFFERS, dataColumn());
        ds.configureUVWMachineCache(uvwMachineCacheSize(),uvwMachineCacheTolerance());
        IDataSharedIter it=makeCalibratedDataIteratorIfNeeded(makeDataIterator(ds));

        ASKAPCHECK(itsModel, "Model not defined");
        ASKAPCHECK(gridder(), "Gridder not defined");
        boost::shared_ptr<ImageFFTEquation> fftEquation(new ImageFFTEquation (itsModel, it, gridder()));
        ASKAPDEBUGASSERT(fftEquation);
        fftEquation->configure(parset());
        fftEquation->setVisUpdateObject(GroupVisAggregator::create(itsComms));
        fftEquation->setCalDirMap(itsCalDirMap);
        itsEquation = fftEquation;
      }
      else {
        ASKAPLOG_INFO_STR(logger, "Reusing measurement equation and updating with latest model images" );
        itsEquation->setParameters(*itsModel);
      }
      ASKAPCHECK(itsEquation, "Equation not defined");
      ASKAPCHECK(itsNe, "NormalEquations not defined");
      itsEquation->calcEquations(*itsNe);
      ASKAPLOG_INFO_STR(logger, "Calculated normal equations for "<< ms << " in "<< timer.real()
                         << " seconds ");
    }

    /// @brief obtain measurement equation cast to ImageFFTEquation
    /// @details This helper method encapsulates operations common to a number of methods of this and derived classes to obtain the
    /// current measurement equation with the original type as created (i.e. ImageFFTEquation) and
    /// does the appropriate checks (so the return is guaranteed to be a non-null shared pointer).
    /// @return shared pointer of the appropriate type to the current measurement equation
    boost::shared_ptr<ImageFFTEquation> ImagerParallel::getMeasurementEquation() const
    {
      ASKAPCHECK(itsEquation, "Equation not defined");
      const boost::shared_ptr<ImageFFTEquation> fftEquation = boost::dynamic_pointer_cast<ImageFFTEquation>(itsEquation);
      ASKAPCHECK(fftEquation, "Incompatible type of the measurement equation is in use (this shouldn't happen - logic error suspected).");
      return fftEquation;
    }

    /// Calculate the normal equations for a given measurement set
    void ImagerParallel::calcNE()
    {
      ASKAPTRACE("ImagerParallel::calcNE");
      /// Now we need to recreate the normal equations
      itsNe=ImagingNormalEquations::ShPtr(new ImagingNormalEquations(*itsModel));

      if (itsComms.isWorker())
      {
        ASKAPCHECK(gridder(), "Gridder not defined");
        ASKAPCHECK(itsModel, "Model not defined");
        //				ASKAPCHECK(measurementSets().size()>0, "Data sets not defined");

        ASKAPCHECK(itsNe, "NormalEquations not defined");

        if (itsComms.isParallel())
        {
          calcOne(measurementSets()[itsComms.rank()-1]);
          sendNE();
        }
        else
        {
          ASKAPCHECK(itsSolver, "Solver not defined correctly");
          itsSolver->init();
          for (size_t iMs=0; iMs<measurementSets().size(); ++iMs)
          {
            calcOne(measurementSets()[iMs],true);
            itsSolver->addNormalEquations(*itsNe);
          }
        }
      }
    }

    /// @brief helper method to indentify model parameters to broadcast
    /// @details We use itsModel to buffer some derived images like psf, weights, etc
    /// which are not required for prediffers. It just wastes memory and CPU time if
    /// we broadcast them. At the same time, some auxilliary parameters like peak
    /// residual value need to be broadcast (so the major cycle can terminate in workers).
    /// This method returns the vector with all parameters to be broadcast. By default
    /// it returns all parameter names, so it is overridden here to broadcast only
    /// model images and the peak_residual metadata.
    /// @return a vector with parameters to broadcast
    std::vector<std::string> ImagerParallel::parametersToBroadcast() const
    {
       ASKAPDEBUGASSERT(itsModel);
       const std::vector<std::string> names = itsModel->names();
       std::vector<std::string> result;
       result.reserve(names.size());
       for (std::vector<std::string>::const_iterator ci=names.begin(); ci!=names.end(); ++ci) {
            if ((ci->find("image") == 0) || (ci->find("peak_residual") == 0) ||
                (ci->find("uvweight") == 0) || (ci->find("noise_threshold_reached") == 0)) {
                result.push_back(*ci);
            }
       }
       return result;
    }


    void ImagerParallel::solveNE()
    {
      ASKAPTRACE("ImagerParallel::solveNE");
      if (itsComms.isMaster())
      {
        // Receive the normal equations
        if (itsComms.isParallel())
        {
          receiveNE();
        }
        ASKAPLOG_INFO_STR(logger, "Solving normal equations");
        casacore::Timer timer;
        timer.mark();
        Quality q;
        ASKAPDEBUGASSERT(itsModel);
        itsSolver->solveNormalEquations(*itsModel,q);
        ASKAPLOG_INFO_STR(logger, "Solved normal equations in "<< timer.real() << " seconds "
                           );

        // we will probably send all of them out in the future, but for now
        // let's extract the largest residual
        const std::vector<std::string> peakParams = itsModel->completions("peak_residual.",true);

        double peak = peakParams.size() == 0 ? getPeakResidual() : -1.e-10;

        // note we use a negative peak val to signal deconvolution divergence and pass that on here
        for (std::vector<std::string>::const_iterator peakParIt = peakParams.begin();
             peakParIt != peakParams.end(); ++peakParIt) {
             const double tempval = itsModel->scalarValue("peak_residual."+*peakParIt);
             ASKAPLOG_INFO_STR(logger, "Peak residual for "<< *peakParIt << " is "<<std::abs(tempval));
             if (std::abs(tempval) > std::abs(peak)) {
                 peak = tempval;
             }
        }

        if (itsModel->has("peak_residual")) {
            itsModel->update("peak_residual",peak);
        } else {
            itsModel->add("peak_residual",peak);
        }
        itsModel->fix("peak_residual");
      }

      // check if all images have reached the noise threshold
      bool allDone = true;
      const std::vector<std::string> noiseParams = itsModel->completions("noise_threshold_reached.",true);
      for (const std::string& name : noiseParams) {
          if (itsModel->scalarValue("noise_threshold_reached."+name) < 0.0) {
              allDone = false;
          }
      }
      if (itsModel->has("noise_threshold_reached")) {
          itsModel->update("noise_threshold_reached", allDone ? 1.0 : -1.0);
      } else {
          itsModel->add("noise_threshold_reached", allDone ? 1.0 : -1.0);
      }
      itsModel->fix("noise_threshold_reached");

    }

    /// @brief Helper method to zero all model images
    /// @details We need this for dirty solver only, as otherwise restored image
    /// (which is crucial for faceting) will be wrong.
    void ImagerParallel::zeroAllModelImages() const
    {
      ASKAPCHECK(itsModel, "Model should not be empty at this stage!");
      ASKAPLOG_INFO_STR(logger, "Dirty solver mode, setting all model images to 0.");
      SynthesisParamsHelper::zeroAllModelImages(itsModel);
    }


    /// @brief a helper method to extract peak residual
    /// @details This object actually manipulates with the normal equations. We need
    /// to be able to stop iterations on the basis of maximum residual, which is a
    /// data vector of the normal equations. This helper method is designed to extract
    /// peak residual. It is then added to a model as a parameter (the model is
    /// shipped around).
    /// @return absolute value peak of the residuals corresponding to the current normal
    /// equations
    double ImagerParallel::getPeakResidual() const
    {
      ASKAPDEBUGASSERT(itsNe);
      // we need a specialized method of the imaging normal equations to get the peak
      // for all images. Multiple images can be represented by a single normal equations class.
      // We could also use the dataVector method of the interface (INormalEquations). However,
      // it is a bit cumbersome to iterate over all parameters. It is probably better to
      // leave this full case for a future as there is no immediate use case.
      boost::shared_ptr<ImagingNormalEquations> ine =
                    boost::dynamic_pointer_cast<ImagingNormalEquations>(itsNe);
      // we could have returned some special value (e.g. negative), but throw exception for now
      ASKAPCHECK(ine, "Current code to calculate peak residuals works for imaging-specific normal equations only");
      double peak = -1.;
      const std::map<string, casacore::Vector<imtype> >& dataVector = ine->dataVector();
      const std::map<string, casacore::Vector<imtype> >& diag = ine->normalMatrixDiagonal();
      for (std::map<string, casacore::Vector<imtype> >::const_iterator ci = dataVector.begin();
           ci!=dataVector.end(); ++ci) {
           if (ci->first.find("image") == 0) {
               // this is an image
               ASKAPASSERT(ci->second.nelements() != 0);
               std::map<std::string, casacore::Vector<imtype> >::const_iterator diagIt =
                            diag.find(ci->first);
               ASKAPDEBUGASSERT(diagIt != diag.end());
               const double maxDiag = casacore::max(diagIt->second);
               // hard coded at this stage
               const double cutoff=1e-2*maxDiag;
               ASKAPDEBUGASSERT(diagIt->second.nelements() == ci->second.nelements());
               for (casacore::uInt elem = 0; elem<diagIt->second.nelements(); ++elem) {
                    const double thisDiagElement = std::abs(diagIt->second[elem]);
                    if (thisDiagElement > cutoff) {
                        const double tempPeak =  ci->second[elem]/thisDiagElement;
                        if (tempPeak > peak) {
                            peak = tempPeak;
                        }
                    }
               }
           }
      }
      return peak;
    }

    /// @brief make sensitivity image
    /// @details This is a helper method intended to be called from writeModel. It
    /// converts the given weights image into a sensitivity image and exports it.
    /// This method is intended to be called if itsWriteSensitivityImage is true.
    /// @param[in] wtImage weight image parameter name
    void ImagerParallel::makeSensitivityImage(const std::string &wtImage) const
    {
      ASKAPTRACE("ImagerParallel::makeSensitivityImage");
      ASKAPLOG_INFO_STR(logger, "Making sensitivity image from weights image "<<wtImage);
      scimath::Params tempPar;
      ASKAPDEBUGASSERT(itsModel);
      ASKAPDEBUGASSERT(itsModel->has(wtImage));
      ASKAPASSERT(wtImage.find("weights") == 0);
      ASKAPCHECK(wtImage.size()>7, "Weights image parameter name should be longer, you have "<<wtImage);
      const std::string outParName = "sensitivity" + wtImage.substr(7);
      scimath::Axes axes = itsModel->axes(wtImage);
      //
      casacore::Array<imtype> wtArr = itsModel->valueT(wtImage);
      casacore::Array<imtype> sensitivityArr(wtArr.shape());
      const double cutoff = casacore::max(wtArr) * itsSensitivityCutoff;

      for (imagemath::MultiDimArrayPlaneIter iter(wtArr.shape()); iter.hasMore(); iter.next()) {
           const casacore::Vector<imtype> wtPlane = iter.getPlaneVector(wtArr);
           casacore::Vector<imtype> sensitivityPlane = iter.getPlaneVector(sensitivityArr);
           for (casacore::uInt elem = 0; elem < wtPlane.nelements(); ++elem) {
                const double wt = wtPlane[elem];
                if (wt > cutoff) {
                    // at this stage - just reciprocal. Still need to work on the normalisation
                    sensitivityPlane[elem] = 1./sqrt(wt);
                } else {
                    sensitivityPlane[elem] = 0.;
                }
           }
      }

      tempPar.add(outParName,sensitivityArr,axes);
      ASKAPLOG_INFO_STR(logger, "Saving " << outParName);
      const LOFAR::ParameterSet keywords = parset().makeSubset("header.");

      if (parset().isDefined("Images.extraoversampling")) {
          // should only be defined if it is set to a legitimate value. Check anyway.
          const float extraOSfactor = parset().getFloat("Images.extraoversampling");
          ASKAPDEBUGASSERT(extraOSfactor > 1.);
          SynthesisParamsHelper::saveImageParameter(tempPar, outParName, outParName, extraOSfactor, keywords);
      } else {
          SynthesisParamsHelper::saveImageParameter(tempPar, outParName, outParName, boost::none, keywords);
      }
    }

    /// @brief replace normal equations with an adapter handling uv weights
    /// @details We use normal equations infrastructure to merge uv weights built in parallel
    /// (via EstimatorAdapter). This method sets up the appropriate builder class based on the
    /// parameters in the parset, wraps it in the adapter and assigns to itsNE in the base class.
    /// @note An instance of appropriate normal equations class should be setup before normal major
    /// cycles can resume. It can be done with a call to recereateNormalEquations or calcNE
    void ImagerParallel::setupUVWeightBuilder()
    {
        // need to read config from parset regarding index translation and setup generic builder
        const bool perBeamHandling = parset().getBool("uvweight.perbeam", false);
        // by default, set this parameter to something larger than or equal to the number of beams we're likely to encounter. This way, the user
        // wouldn't need to worry about setting it
        const casacore::uInt fieldCoeff = perBeamHandling ? parset().getUint32("uvweight.maxbeams", 36u) : 0u;
        const casacore::uInt beamCoeff = perBeamHandling ? 1u : 0u;

        // don't bother with the 3rd coefficient for now, setting it to zero means it would be ignored
        const boost::shared_ptr<GenericUVWeightBuilder> builder(new GenericUVWeightBuilder(beamCoeff, fieldCoeff, 0u));
        const boost::shared_ptr<scimath::EstimatorAdapter<GenericUVWeightBuilder> > adapter(new scimath::EstimatorAdapter<GenericUVWeightBuilder>(builder));
        // finally, assign the adapter as "NormalEquations"
        setNE(adapter);
    }

    /// @brief helper method to extract weight builder object out of normal equations
    /// @details For traditional weighting with distributed data we use normal equation merging
    /// mechanism to do the gather operation (and EstimatorAdapter). This method does required
    /// casts and checks to get the required shared pointer (which is guaranteed to be non-empty)
    /// @return shared pointer to the uv-weight builder object stored in the current normal equations
    boost::shared_ptr<GenericUVWeightBuilder> ImagerParallel::getUVWeightBuilder() const {
       const boost::shared_ptr<scimath::EstimatorAdapter<GenericUVWeightBuilder> > adapter =
             boost::dynamic_pointer_cast<scimath::EstimatorAdapter<GenericUVWeightBuilder> >(getNE());
       ASKAPCHECK(adapter, "Incompatible type of normal equations detected while trying to access UV Weight builder");
       const boost::shared_ptr<GenericUVWeightBuilder> builder = adapter->get();
       ASKAPCHECK(builder, "Empty shared pointer found in the EstimatorAdapter while trying to access UV Weight builder");
       return builder;
    }

    /// @brief compute uv weights using data stored via adapter in the normal equations
    /// @details This method gets access to the uv-weight builder handled via EstimatorAdapter and
    /// stored instead of normal equations by shared pointer. It then runs the finalisation step using
    /// the stored shared pointer to the uv-weight calculator object. The resulting weight is added as
    /// a parameter to the existing model.
    /// @note This method assumes that it is called from the right place (i.e. on the correct rank) and
    /// all merging/reduction has already been done. In other words, it is agnostic of the parallelism.
    void ImagerParallel::computeUVWeights() const
    {
       ASKAPDEBUGASSERT(itsUVWeightCalculator);
       const boost::shared_ptr<GenericUVWeightBuilder> builder = getUVWeightBuilder();
       UVWeightCollection& wts = builder->finalise(*itsUVWeightCalculator);
       // now need to store results into model. Note, technically the framework can support multiple image parameters for the same data
       // in this case, weighting can be performed separately (or done with the same weight). This possible complexity is not handled by the
       // high-level code at the moment. Some changes would be necessary here to handle such cases (and probably in other places as well)
       ASKAPDEBUGASSERT(itsModel);
       // MV: some part of the code below could probably go into UVWeightParamsHelper
       const std::vector<std::string> completions = itsModel->completions("image");
       // although we currently allow only one simultaneous image there could still be many parameters sharing the same weight (e.g. Taylor terms)
       // need to loop over all completions and parse the parameter name anyway
       std::set<std::string> currentParamNames;
       for (std::vector<std::string>::const_iterator ci = completions.begin(); ci != completions.end(); ++ci) {
            ImageParamsHelper iph(ImageParamsHelper::replaceLeadingWordWith("image"+*ci, "image.",""));
            currentParamNames.insert(iph.facetName());
       }
       ASKAPCHECK(currentParamNames.size() > 0, "Unable to find any suitable image parameter name to do traditional weighting for, check that the model has been setup");
       // for now abort if there is more than one free image parameter, although it would be very straightforward to setup e.g. the same weighting for all of them
       // this may be needed for facets to work with traditional weighting!
       ASKAPCHECK(currentParamNames.size() == 1, "We currently support only one free image parameter with traditional weighting you have "<<currentParamNames.size());
       const std::string paramName = *currentParamNames.begin();
       // for now, figure out and copy index translation details from the builder and pass it on as is. However, here we can setup more logic
       // to do non-trivial stuff, e.g. select a particular weight grid and apply to other data, etc
       const casacore::uInt coeffBeam = builder->indexOf(1u, 0u, 0u);
       const casacore::uInt coeffField = builder->indexOf(0u, 1u, 0u);
       const casacore::uInt coeffSource = builder->indexOf(0u, 0u, 1u);
       const boost::shared_ptr<GenericUVWeightIndexTranslator> translator(new GenericUVWeightIndexTranslator(coeffBeam, coeffField, coeffSource));
       UVWeightParamsHelper hlp(itsModel);
       hlp.addUVWeights(paramName, wts, translator);
    }

    /// @brief iterate over given data and accumulate samples for uv weights
    /// @details This method is used to build the sample density in the uv-plane via the appropriate gridder
    /// and weight builder class. It expects the builder already setup and accessible via the normal equations
    /// shared pointer. The data iterator to work with is passed as a parameter. The image details are extracted from
    /// the model (to initialise sample grid).
    /// @param[in] iter shared pointer to the iterator to use (note it is advanced by this method to iterate over
    /// all available data)
    void ImagerParallel::accumulateUVWeights(const boost::shared_ptr<accessors::IConstDataIterator> &iter) const
    {
       // first, figure out the image parameter to work with as we need the cell size and image shape to setup the weight grid
       // Again similar code to that in ImageFFTEquation but we can't easily reuse it (although, perhaps, some refactoring is possible). 
       // There is a complication that we can have more than one free image parameter (even after removing Taylor terms). For now, only 
       // support a single image, which would break faceting. It is straightforward to extend the code to support multiple facets if we 
       // build just one weight for all facets, otherwise more thoughts are needed (but, perhaps, using the 3rd index we could have facets 
       // factored in and even build separate weight for each facet).

       // again some code duplication with what we have above, but here we want to retain the full parameter name
       // it would be nice to think about some refactoring
       ASKAPDEBUGASSERT(itsModel);
       const std::vector<std::string> completions = itsModel->completions("image");
       std::set<std::string> currentParamNames;
       for (std::vector<std::string>::const_iterator ci = completions.begin(); ci != completions.end(); ++ci) {
            ImageParamsHelper iph("image"+*ci);
            if (iph.isTaylorTerm()) {
                // in the case of multiple Taylor terms work with order 0, also it may be absent in distributed Taylor terms case - see below
                iph.makeTaylorTerm(0);
            }
            // MV: we could've selected the required facet here too, but don't bother with it at this stage
            // passing all facets as they are, if present
            currentParamNames.insert(iph.paramName());
       }
       ASKAPCHECK(currentParamNames.size() > 0, "Unable to find any free image parameter in the current model, there is nothing to build the uv-weight for!");
       // the next check can be commented out, if we always want to use the parameters extracted from the first image (or some logic is necessary to choose the
       // one we want). Leave the check in place for now, as it can alert us to some unexpected use cases
       ASKAPCHECK(currentParamNames.size() == 1, "Only one image parameter is currently supported for traditional weighting");
       const std::string imageParamName = *currentParamNames.begin();
       // for distributed Taylor terms, term 0 may actually be on a different rank. But we can skip weight calculation for other Taylor terms
       // (the merge logic would merge the empty builder with whatever it will try to merge in; more thoughts would be needed if there are multiple images)
       if (!itsModel->has(imageParamName)) {
           ASKAPLOG_DEBUG_STR(logger, "Distributed Taylor terms case, skip weight accumulation for non zero Taylor terms");
           return;
       }

       // now get the builder (should already be setup by the time this method is called) and setup the weight gridder
       const boost::shared_ptr<GenericUVWeightBuilder> builder = getUVWeightBuilder();
       UVWeightGridder gridder(builder);

       // setup additional parameters of the specialised gridder, more complex logic can be added here later on
       // we can also probably factor out the code configuring the gridder into a separate method to have better structure of the code
       const bool perBeamHandling = parset().getBool("uvweight.perbeam", false);
       gridder.doBeamAndFieldSelection(!parset().getBool("uvweight.alldatawt", perBeamHandling));

       // there is a bit of the technical debt here with code duplication with VisGridderFactory, but this is an unrelated class
       // which accepts the same parameters (and it is handy to force it to use the gridder parameters)
       if (parset().isDefined("gridder.MaxPointingSeparation")) {
           const double threshold = SynthesisParamsHelper::convertQuantity(
                 parset().getString("gridder.MaxPointingSeparation","-1rad"),"rad");
           ASKAPLOG_INFO_STR(logger,"MaxPointingSeparation is used, data from pointing centres further than "<<
                  threshold*180./casacore::C::pi<<" deg from the image centre will be rejected in uv-weight calculation");
           gridder.maxPointingSeparation(threshold);
       } else {
           ASKAPLOG_INFO_STR(logger,"MaxPointingSeparation is not used for uv-weight gridder, all data will be used");
       }

       if (parset().isDefined("gridder.padding")) {
           const float padding = parset().getFloat("gridder.padding");
           ASKAPLOG_INFO_STR(logger, "Use padding at the weight gridder level, padding factor = " << padding);
           gridder.setPaddingFactor(padding);
       } else {
           ASKAPLOG_INFO_STR(logger,"No padding at the weight gridder level");
       }

       // technical debt!
       // this is hopefully a temporary hack - due to the current way to place oversampling planes, the weight gridder needs to know
       // the oversampling factor used by the actual data gridder. The following code gets it from the parset and sets to the weight gridder.
       // And unfortunately, the defaults are different for different gridders which adds to technical debt.
       // Note, the weight gridder doesn't oversample its grid, it just needs to sample the same way for the first oversampling plane
       const std::string dataGridderName = parset().getString("gridder");
       int defaultOversample = dataGridderName == "Box" ? 1 : 8;
       if (dataGridderName == "SphFunc" || dataGridderName == "WStack") {
           defaultOversample = 128;
       } 
       const int oversample = parset().getInt32("gridder."+dataGridderName+".oversample", defaultOversample);
       ASKAPCHECK(oversample > 0, "Oversampling factor is supposed to be positive, you have "<<oversample);
       ASKAPLOG_INFO_STR(logger, "The weight gridder will assume that the data gridder ("<<dataGridderName<<") is using the oversampling factor of "<<oversample);
       ASKAPLOG_DEBUG_STR(logger, "data gridder parameters: "<<parset().makeSubset("gridder."+dataGridderName));
       gridder.setOversampleFactor(oversample);
       //

       // now figure out the cell size and image shape for the image parameter selected above and pass them to the gridder
       const scimath::Axes axes(itsModel->axes(imageParamName));
       const casacore::IPosition imageShape = itsModel->shape(imageParamName);
       gridder.initialise(axes, imageShape);

       // now the setup is done, so we can iterate over data and accumulate
       ASKAPLOG_INFO_STR(logger, "Iterating over data to compute the density of samples for uv-weight construction");
       IConstDataSharedIter it(iter);
       for (it.init(); it.hasMore(); it.next()) {
            gridder.accumulate(*it);
       }
    }

    /// @brief recreate imaging normal equations object
    /// @details If sample density grid is built, normal equations are setup with
    /// an adapter which is an incompatible type. This method resets the object to a pristine
    /// state similar to that before the first major cycle if no traditional weighting is done.
    void ImagerParallel::recreateNormalEquations()
    {
       ASKAPDEBUGASSERT(itsModel);
       ASKAPLOG_DEBUG_STR(logger,"Recreating NE from model "<<*itsModel);
       ImagingNormalEquations::ShPtr newNE(new ImagingNormalEquations(*itsModel));
       setNE(newNE);
    }

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
    void ImagerParallel::createUVWeightCalculator()
    {
       // this is hopefully a legacy method, need some cleanup
       itsUVWeightCalculator = createUVWeightCalculator(parset());
    }

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
    boost::shared_ptr<IUVWeightCalculator> ImagerParallel::createUVWeightCalculator(const LOFAR::ParameterSet &parset) 
    {
       const std::string keyword = "uvweight";
       if (!parset.isDefined(keyword)) {
           return boost::shared_ptr<IUVWeightCalculator>();
       }
       const std::vector<std::string> wtCalcList = parset.getStringVector(keyword);
       ASKAPCHECK(wtCalcList.size() > 0u, "Cimager.uvweight should contain either a single keyword describing how the weight is obtained or a vector with procedure names to apply these to measured uv-density");

       // also need to check here later on that Cimager.uvweight parameter is set to one of the resereved keywords and return an empty
       // shared pointer if this is the case. Or factor out this method into a separate class where the same logic would be done
       // some other way

       std::vector<boost::shared_ptr<IUVWeightCalculator> >  calculators(wtCalcList.size());
       ASKAPLOG_INFO_STR(logger, "Setting up uv-weight calculation:");
       // the same could've been achieved with iterators but performance gain is very small, if any and this code is cleaner
       for (size_t index = 0; index < wtCalcList.size(); ++index) {
            // the actual factory code, we can factor this out into a separate method
            const std::string name = wtCalcList[index];
            if (name == "Robust") {
                // we could've had parameters in the form Cimager.uvweight.Robust.robustness to follow a more structured apporach - can be changed if we want it
                const float robustness = parset.getFloat(keyword + ".robustness");
                ASKAPLOG_INFO_STR(logger, "        + "<<name<<": robust weighting with robustness = "<<robustness);
                const boost::shared_ptr<RobustUVWeightCalculator> calc(new RobustUVWeightCalculator(robustness));
                calculators[index]  = calc;
            } else {
                if (name == "ConjugatesAdderFFT") {
                    ASKAPLOG_INFO_STR(logger, "        + "<<name<<": ensuring conjugate symmetry via FFT");
                    const boost::shared_ptr<ConjugatesAdderFFT> calc(new ConjugatesAdderFFT());
                    calculators[index] = calc;
                } else {
                    if (name == "Reciprocal") {
                        const float threshold = parset.getFloat(keyword + ".recipthreshold", 1e-5);
                        ASKAPLOG_INFO_STR(logger, "        + "<<name<<": calculating reciprocal for weight application (threshold = "<<threshold<<")");
                        const boost::shared_ptr<ReciprocalUVWeightCalculator> calc(new ReciprocalUVWeightCalculator(threshold));
                        calculators[index] = calc;
                    } else {
                       // taper class comes here when we have it
                       ASKAPTHROW(AskapError, "Unknown type of the uv-weight calculator: "<<name);
                    }
                }
            }
       }

       if (calculators.size() == 1) {
           return calculators[0];
       } 

       // there are several effects which need to be applied one by one, create composite calculator to achieve this
       const boost::shared_ptr<CompositeUVWeightCalculator> result(new CompositeUVWeightCalculator(calculators.begin(), calculators.end()));
       return result;
    }


    /// Write the results out
    /// @param[in] postfix this string is added to the end of each name
    /// (used to separate images at different iterations)
    void ImagerParallel::writeModel(const std::string &postfix)
    {
      ASKAPTRACE("ImagerParallel::writeModel");
      if (itsComms.isMaster())
      {
        ASKAPLOG_INFO_STR(logger, "Writing out results as images");
        ASKAPDEBUGASSERT(itsModel);
        std::vector<std::string> resultimages=itsModel->names();
        bool hasWeights = false;
        for (std::vector<std::string>::const_iterator it=resultimages.begin(); it
            !=resultimages.end(); it++) {
            if (it->find("weights") == 0) {
                hasWeights = true;
            }
        }
        if (!hasWeights && (itsWriteWtImage || itsWriteSensitivityImage || itsWriteWtLog)) {
            ASKAPDEBUGASSERT(itsSolver);
            boost::shared_ptr<ImageSolver> image_solver = boost::dynamic_pointer_cast<ImageSolver>(itsSolver);
            ASKAPDEBUGASSERT(image_solver);
            image_solver->saveWeights(*itsModel);
            resultimages=itsModel->names();
        }

        // Check whether or not the model has been stored with a higher resolution
        boost::optional<float> extraOSfactor;
        if (parset().isDefined("Images.extraoversampling")) {
            extraOSfactor = parset().getFloat("Images.extraoversampling");
            // The parameter should only be defined if has a legitimate value (is set by the code). Check anyway.
            ASKAPDEBUGASSERT(*extraOSfactor > 1.);
        }

        // Get any header keywords to add to the images from the parset
        LOFAR::ParameterSet keywords = parset().makeSubset("header.");

        if (hasWeights && !itsWriteWtImage) {
            ASKAPLOG_INFO_STR(logger,"Writing weights "<< (itsWriteWtLog ? "log":"keyword"));
            askap::accessors::WeightsLog weightslog;
            string name;
            // look for name of weights image
            for (std::vector<std::string>::const_iterator it=resultimages.begin();
                it !=resultimages.end(); it++) {
                if (it->find("weights") == 0) {
                    name = *it;
                    // could be taylor2 - change to taylor0
                    int n = name.find(".taylor.");
                    if (n != string::npos) {
                        name.replace(n+8,1,"0");
                    }
                    break;
                }
            }

            casacore::Array<float> wts = itsModel->valueF(name);
            float wt = wts.data()[0];
            if (allEQ(wts,wt)) {
                if (itsWriteWtLog) {
                    weightslog.weightslist()[0] = wt;
                    weightslog.setFilename(name + postfix + ".txt");
                    weightslog.write();
                }
                // always write the keyword if possible?
                keywords.replace("IMWEIGHT","["+std::to_string(wt)+",Imaging Weight]");
            } else {
                ASKAPLOG_WARN_STR(logger,"Weights are not identical across image, disabling weight log");
                itsWriteWtImage = true;
            }
        }

        // get the imageHistory keyword in the parset
        std::vector<std::string> historyLines;
        if ( parset().isDefined("imageHistory") ) {
            historyLines = parset().getStringVector("imageHistory");
        } else {
            ASKAPLOG_DEBUG_STR(logger, "---> imageHistory is not defined");
        }

        for (std::vector<std::string>::const_iterator it=resultimages.begin(); it
            !=resultimages.end(); it++) {
            const ImageParamsHelper iph(*it);
            const bool doWrite = itsWriteMultiple || it->find(itsFirstImageName) != std::string::npos;

            // if true, "image.*" is retained for degridding but "fullres.*" is used for cleaning & restoring
            // always write model if writeAtMajorCycle is true (sets postfix)
            if (extraOSfactor) {
                if ((it->find("fullres") == 0) && (itsWriteModelImage || postfix!="")) {
                    // change "fullres" back to "image" for output
                    string tmpname = *it;
                    tmpname.replace(0,7,"image");
                    if (doWrite || itsWriteMultipleModels) {
                        ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << tmpname+postfix );
                        accessors::IImageAccess<>& imageAccessor = SynthesisParamsHelper::imageHandler();
                        SynthesisParamsHelper::saveImageParameter(*itsModel, *it, tmpname+postfix, boost::none, keywords, historyLines);
                        // write the image stats to the image table
                        askap::utils::StatsAndMask::writeStatsToImageTable(itsComms,imageAccessor,tmpname+postfix,parset());
                    }
                }
            } else {
                if ((it->find("image") == 0) && (itsWriteModelImage || postfix!="")) {
                    if (doWrite || itsWriteMultipleModels) {
                        ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                        accessors::IImageAccess<>& imageAccessor = SynthesisParamsHelper::imageHandler();
                        SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, boost::none, keywords, historyLines);
                        // write the image stats to the image table
                        askap::utils::StatsAndMask::writeStatsToImageTable(itsComms,imageAccessor,*it+postfix,parset());
                    }
                }
            }
            if ((it->find("weights") == 0) && (itsWriteWtImage||itsWriteSensitivityImage))  {
                if (itsWriteWtImage && doWrite) {
                    ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                    SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, extraOSfactor, keywords, historyLines);
                }
                if (itsWriteSensitivityImage && (it->find("weights") == 0)
                    && (postfix == "") && doWrite) {
                        makeSensitivityImage(*it);
                }
            }
            if ((it->find("mask") == 0) && itsWriteMaskImage && doWrite)  {
                ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, extraOSfactor, keywords, historyLines);
            }
            if ((it->find("residual") == 0) && itsWriteResidual) {
                if (!iph.isFacet()) {
                    if (!itsRestore && itsWriteResidual && doWrite) {
                        ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                        SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, extraOSfactor, keywords, historyLines);
                        // write the image stats to the image table
                        accessors::IImageAccess<>& imageAccessor = SynthesisParamsHelper::imageHandler();
                        askap::utils::StatsAndMask::writeStatsToImageTable(itsComms,imageAccessor,*it+postfix,parset());
                    }
                }
                else {
                    if (itsWriteResidual && doWrite) {
                        ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                        SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, extraOSfactor, keywords, historyLines);
                        // write the image stats to the image table
                        accessors::IImageAccess<>& imageAccessor = SynthesisParamsHelper::imageHandler();
                        askap::utils::StatsAndMask::writeStatsToImageTable(itsComms,imageAccessor,*it+postfix,parset());
                    }
                }

            }
            if ((it->find("psf") == 0) && itsWritePsfRaw && doWrite) {
                ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, extraOSfactor, keywords, historyLines);
            }
            if ((it->find("scalemask") == 0) && itsWriteScaleMask && doWrite) {
                ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, boost::none, keywords, historyLines);
            }
        }

        if (itsRestore && postfix == "") {
            // add a second pass if a separate restore preconditioner is defined
            // make a deep copy of the parset
            LOFAR::ParameterSet tmpset = parset().makeSubset("");
            string restore_suffix;
            uint n_passes = 1;
            // if extra passes are required, save anything that is changed so it can be reset
            std::map<string, boost::shared_ptr<casacore::ImageInterface<float> > > saved_models;

            // add a second pass if a separate restore preconditioner is defined
            if (tmpset.isDefined("restore.preconditioner.Names")) {
                n_passes += 1;

                // save the initial values of image model free parameters
                std::vector<std::string> names(itsModel->completions("image"));
                map<string,int> facetmap;
                SynthesisParamsHelper::listFacets(names, facetmap);
                for (map<string,int>::const_iterator ci=facetmap.begin();ci!=facetmap.end();++ci) {
                    string name;
                    if (extraOSfactor) {
                     name = "fullres"+ci->first;
                    }
                    else {
                     name = "image"+ci->first;
                    }
                    if ((ci->second != 1) && !itsModel->has(name)) {
                        // this is a multi-facet image, add a fixed parameter representing the whole image
                        ASKAPLOG_INFO_STR(logger, "Adding a fixed parameter "<<name<<
                        " representing faceted image with "<<ci->second<<" facets");
                        SynthesisParamsHelper::add(*itsModel,name,ci->second);
                        // this isn't a free parameter, it is made from free parameters (i.e. the facets)
                        itsModel->fix(name);
                    }
                    saved_models[name] = SynthesisParamsHelper::tempImage(*itsModel, name);
                }


            }

            for (uint pass=0; pass<n_passes; ++pass) {
                if (pass == 0) {
                    ASKAPLOG_INFO_STR(logger, "Restoring images" <<
                      (n_passes == 1 || itsWriteFirstRestore ? " and writing them to disk" :""));
                    restore_suffix = "";
                }
                else {
                    ASKAPLOG_INFO_STR(logger, "Restoring again with a second preconditioner");
                    // replace any existing preconditioner params with the restore.* set
                    tmpset.subtractSubset("preconditioner.");
                    tmpset.adoptCollection(parset().makeSubset("restore.preconditioner.","preconditioner."));
                    restore_suffix = "."+parset().getString("restore.preconditioner.suffix","alt");
                    // reset image models to be free parameters with their initial values
                    for (std::map<string, boost::shared_ptr<casacore::ImageInterface<float> > >::iterator
                        it=saved_models.begin(); it!=saved_models.end(); ++it)
                    {
                        SynthesisParamsHelper::update(*itsModel, it->first, *(it->second));
                    }

                    std::vector<std::string> names(itsModel->completions("image"));
                    map<string,int> facetmap;
                    SynthesisParamsHelper::listFacets(names, facetmap);
                    for (map<string,int>::const_iterator ci=facetmap.begin();ci!=facetmap.end();++ci) {
                       if (ci->second != 1) {
                           // this isn't a free parameter, it is merged from free parameters (i.e. the facets)
                           itsModel->fix("image"+ci->first);
                       }
                    }

                }
                // Set parset parameters for ImageRestoreSolver to avoid saving things we don't need
                tmpset.replace(LOFAR::KVpair("restore.updateresiduals",itsWriteResidual));
                tmpset.replace(LOFAR::KVpair("restore.savepsfimage",itsWritePsfImage));
                boost::shared_ptr<ImageRestoreSolver>
                ir = ImageRestoreSolver::createSolver(tmpset.makeSubset("restore."));

                ASKAPDEBUGASSERT(ir);
                ASKAPDEBUGASSERT(itsSolver);

                // configure restore solver
                if (extraOSfactor) {
                    ASKAPLOG_INFO_STR(logger,"Configuring restore solver with an extra oversampling factor of "<<
                        *extraOSfactor);
                    ir->setExtraOversampling(*extraOSfactor);
                }

                // configure restore solver
                boost::shared_ptr<ImageSolver> template_solver = boost::dynamic_pointer_cast<ImageSolver>(itsSolver);
                ASKAPDEBUGASSERT(template_solver);
                // use existing preconditioners for pass 0, special restore ones for pass 1
                if (pass>0) {
                    ImageSolverFactory::configurePreconditioners(tmpset,ir);
                }
                ir->configureSolver(*template_solver,pass==0);
                ir->copyNormalEquations(*template_solver);
                Quality q;

                // this solves for residual images and convolves and adds each model image
                // FIX: if itsWriteMultiple==false we do restore work we throw away
                ir->solveNormalEquations(*itsModel,q);
                // merged image should be a fixed parameter without facet suffixes
                if (n_passes == 1 || itsWriteFirstRestore || pass > 0) {
                    std::vector<std::string> resultimages2=itsModel->names();
                    for (std::vector<std::string>::const_iterator
                            ci=resultimages2.begin(); ci!=resultimages2.end(); ++ci) {
                        if (itsWriteMultiple || ci->find(itsFirstImageName)!=std::string::npos) {
                            const ImageParamsHelper iph(*ci);
                            // if true, "image.*" is retained for degridding but "fullres.*" is used for restoring
                            if (extraOSfactor) {
                                if (!iph.isFacet() && ((ci->find("fullres") == 0)))  {
                                    string tmpname = *ci;
                                    tmpname.replace(0,7,"image");
                                    ASKAPLOG_INFO_STR(logger, "Saving restored image " << *ci << " with name "
                                            << tmpname+restore_suffix+".restored" );
                                    SynthesisParamsHelper::saveImageParameter(*itsModel, *ci, tmpname+restore_suffix+".restored", boost::none, keywords, historyLines);
                                    // write the image stats to the image table
                                    accessors::IImageAccess<>& imageAccessor = SynthesisParamsHelper::imageHandler();
                                    askap::utils::StatsAndMask::writeStatsToImageTable(itsComms,imageAccessor,tmpname+restore_suffix+".restored",parset());
                                }
                            } else {
                                if (!iph.isFacet() && ((ci->find("image") == 0)))  {
                                    ASKAPLOG_INFO_STR(logger, "Saving restored image " << *ci << " with name "
                                            << *ci+restore_suffix+".restored" );
                                    SynthesisParamsHelper::saveImageParameter(*itsModel, *ci, *ci+restore_suffix+".restored", boost::none, keywords, historyLines);
                                    // write the image stats to the image table
                                    accessors::IImageAccess<>& imageAccessor = SynthesisParamsHelper::imageHandler();
                                    askap::utils::StatsAndMask::writeStatsToImageTable(itsComms,imageAccessor,*ci+restore_suffix+".restored",parset());
                                }
                            }
                            if (!iph.isFacet() && ((ci->find("psf.image") == 0) && itsWritePsfImage))  {
                                ASKAPLOG_INFO_STR(logger, "Saving psf image " << *ci << " with name "
                                        << *ci+restore_suffix );
                                SynthesisParamsHelper::saveImageParameter(*itsModel, *ci, *ci+restore_suffix, extraOSfactor, keywords, historyLines);
                            }
                            if (!iph.isFacet() && ((ci->find("residual") == 0) && itsWriteResidual))  {
                                ASKAPLOG_INFO_STR(logger, "Saving residual image " << *ci << " with name "
                                        << *ci+restore_suffix );
                                SynthesisParamsHelper::saveImageParameter(*itsModel, *ci, *ci+restore_suffix, extraOSfactor, keywords, historyLines);
                                // write the image stats to the image table
                                accessors::IImageAccess<>& imageAccessor = SynthesisParamsHelper::imageHandler();
                                askap::utils::StatsAndMask::writeStatsToImageTable(itsComms,imageAccessor,*ci+restore_suffix,parset());
                            }
                        }
                    }
                }
            }

            // remove parts of each faceted image
            std::vector<std::string> names(itsModel->completions("image"));
            for (std::vector<std::string>::const_iterator ci=names.begin(); ci !=names.end(); ++ci) {
                const string name="image"+*ci;
                ImageParamsHelper iph(name);
                if (iph.isFacet()) {
                    ASKAPLOG_INFO_STR(logger, "Remove facet patch "<<name<<" from the parameters");
                    itsModel->remove(name);
                }
            }
            ASKAPLOG_DEBUG_STR(logger, "Writing out additional parameters made by restore solver as images");
            // (MHW) Not sure what is suppposed to be written here, turned off psf.image as that is already written above
            std::vector<std::string> resultimages2=itsModel->names();
            for (std::vector<std::string>::const_iterator it=resultimages2.begin(); it !=resultimages2.end(); it++) {
                ASKAPLOG_DEBUG_STR(logger, "Checking "<<*it);
                if ((it->find("psf") == 0)) {
                    ASKAPLOG_DEBUG_STR(logger, "Found " <<*it);

                    if (std::find(resultimages.begin(),resultimages.end(),*it) == resultimages.end()) {
                        if (it->find("psf.image")!=0) {
                            ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                            SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, extraOSfactor, keywords);
                        }
                    }
                    else {
                        ASKAPLOG_DEBUG_STR(logger, "Not Saving as " << *it << " is in the original params list");
                    }
                }
            }
        }
      }
    }
  }
}
