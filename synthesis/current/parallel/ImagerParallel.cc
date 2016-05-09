///
/// @file 
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

#include <parallel/ImagerParallel.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

#include <askap/AskapError.h>
// need it just for null deleter
#include <askap/AskapUtil.h>

#include <askapparallel/AskapParallel.h>
#include <dataaccess/DataAccessError.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>

#include <measurementequation/ImageFFTEquation.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <measurementequation/ImageRestoreSolver.h>
#include <measurementequation/MEParsetInterface.h>
#include <measurementequation/CalibrationIterator.h>
#include <measurementequation/NoXPolGain.h>
#include <measurementequation/ImageParamsHelper.h>
#include <fitting/Params.h>
#include <utils/MultiDimArrayPlaneIter.h>

#include <measurementequation/ImageSolverFactory.h>
#include <calibaccess/CalibAccessFactory.h>
#include <measurementequation/CalibrationApplicatorME.h>
#include <profile/AskapProfiler.h>
#include <parallel/GroupVisAggregator.h>
#include <parallel/AdviseParallel.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/OS/Timer.h>

#include <Common/ParameterSet.h>

#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>

using namespace askap;
using namespace askap::scimath;
using namespace askap::synthesis;
using namespace askap::accessors;
using namespace askap::askapparallel;

namespace askap
{
  namespace synthesis
  {
    using askap::operator<<;

    ImagerParallel::ImagerParallel(askap::askapparallel::AskapParallel& comms,
        const LOFAR::ParameterSet& parset) :
      MEParallelApp(comms,parset),
      itsExportSensitivityImage(false), itsExpSensitivityCutoff(0.)
    {
      if (itsComms.isMaster())
      {      
        itsRestore=parset.getBool("restore", false);
        
        bool reuseModel = parset.getBool("Images.reuse", false);
        
        itsUseMemoryBuffers = parset.getBool("memorybuffers", false);
        
        itsExportSensitivityImage = parset.getBool("sensitivityimage", true);
        
        itsExpSensitivityCutoff = parset.getDouble("sensitivityimage.cutoff", 0.01);
        
        if (itsExportSensitivityImage) {
            ASKAPLOG_INFO_STR(logger, 
               "Theoretical sensitivity images will be generated in addition to weights images, cutoff="<<
                itsExpSensitivityCutoff);
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

        /// Create the solver from the parameterset definition
        itsSolver = ImageSolverFactory::make(parset);
        ASKAPCHECK(itsSolver, "Solver not defined correctly");
      }
      if (itsComms.isWorker())
      {
        bool doCalib = parset.getBool("calibrate",false);
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
      const vector<string> imageNames = parset.getStringVector("Images.Names", false);
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
      if ( cellsizeNeeded && !parset.isDefined("nUVWMachines") ) {
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
      const vector<string> imageNames = parset.getStringVector("Images.Names", false);
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

      // Use the snapshotimaging wtolerance as the advise wtolerance. But only if wmax is needed.
      param = "gridder.snapshotimaging.wtolerance";
      if (parset.isDefined(param) && !parset.isDefined("wtolerance") && (wMaxAdviceNeeded(parset)!="")) {
          string wtolerance = parset.getString(param);
          ASKAPLOG_INFO_STR(logger, "  Adding wtolerance for advise: " << wtolerance);
          parset.add("wtolerance", wtolerance);
          parset.add("wtolerance.advised", wtolerance);
      }

    }

    /// test whether to advise on wmax, which can require an extra pass over the data.
    string ImagerParallel::wMaxAdviceNeeded(LOFAR::ParameterSet &parset) {

      // make a vector containing all gridders that require the wmax parameter
      vector<string> wGridders;
      wGridders.push_back("WProject"); 
      wGridders.push_back("WStack"); 
      wGridders.push_back("AWProject"); 
      wGridders.push_back("AProjectWStack"); 
      for(vector<string>::const_iterator i = wGridders.begin(); i != wGridders.end(); ++i) {
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
      const vector<string> imageNames = parset.getStringVector("Images.Names", false);
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
          const long lSize = long(fieldSize * 3600 / cellSize[0]) + 1;
          const long mSize = long(fieldSize * 3600 / cellSize[1]) + 1;
          pstr<<"["<<lSize<<","<<mSize<<"]";
          ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
          parset.add(param, pstr.str().c_str());
      }
      string gridder = wMaxAdviceNeeded(parset); // returns empty string if wmax is not required.
      if (gridder!="") {
          param = "gridder."+gridder+".wmax"; // if wmax is undefined but needed, use the advice.
          if (!parset.isDefined(param)) {
              std::ostringstream pstr;
              // both should be set if either is, but include the latter to ensure that maxResidualW is generated.
              if (parset.isDefined("gridder.snapshotimaging.wtolerance") && parset.isDefined("wtolerance") ) {
                  pstr<<advice.maxResidualW(); // could use parset.getString("gridder.snapshotimaging.wtolerance");
                  ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
                  parset.add(param, pstr.str().c_str());
              } else {
                  pstr<<advice.maxW();
                  ASKAPLOG_INFO_STR(logger, "  Advising on parameter " << param <<": " << pstr.str().c_str());
                  parset.add(param, pstr.str().c_str());
              }
          }
      }

    }

    void ImagerParallel::calcOne(const string& ms, bool discard)
    {
      ASKAPDEBUGTRACE("ImagerParallel::calcOne");
      casa::Timer timer;
      timer.mark();
      ASKAPLOG_INFO_STR(logger, "Calculating normal equations for " << ms );
      // First time around we need to generate the equation 
      if ((!itsEquation)||discard)
      {
        ASKAPLOG_INFO_STR(logger, "Creating measurement equation" );

        // just to print the current mode to the log
        if (itsUseMemoryBuffers) {
            ASKAPLOG_INFO_STR(logger, "Scratch data will be held in memory" );
        } else {
            ASKAPLOG_INFO_STR(logger, "Scratch data will be written to the subtable of the original dataset" );
        }
        
        TableDataSource ds(ms, (itsUseMemoryBuffers ? TableDataSource::MEMORY_BUFFERS : TableDataSource::DEFAULT), 
                           dataColumn());
        ds.configureUVWMachineCache(uvwMachineCacheSize(),uvwMachineCacheTolerance());                   
        IDataSelectorPtr sel=ds.createSelector();
        sel->chooseCrossCorrelations();
        sel << parset();
        IDataConverterPtr conv=ds.createConverter();
        conv->setFrequencyFrame(getFreqRefFrame(), "Hz");
        conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
        // ensure that time is counted in seconds since 0 MJD
        conv->setEpochFrame();
        
        IDataSharedIter it=ds.createIterator(sel, conv);
        ASKAPCHECK(itsModel, "Model not defined");
        ASKAPCHECK(gridder(), "Gridder not defined");
        if (!itsSolutionSource) {
            ASKAPLOG_INFO_STR(logger, "No calibration is applied" );
            boost::shared_ptr<ImageFFTEquation> fftEquation(new ImageFFTEquation (*itsModel, it, gridder()));
            ASKAPDEBUGASSERT(fftEquation);
            fftEquation->useAlternativePSF(parset());
            fftEquation->setVisUpdateObject(GroupVisAggregator::create(itsComms));
            itsEquation = fftEquation;
        } else {
            ASKAPLOG_INFO_STR(logger, "Calibration will be performed using solution source");
            boost::shared_ptr<ICalibrationApplicator> calME(new CalibrationApplicatorME(itsSolutionSource));
            // fine tune parameters
            ASKAPDEBUGASSERT(calME);
            calME->scaleNoise(parset().getBool("calibrate.scalenoise",false));
            calME->allowFlag(parset().getBool("calibrate.allowflag",false));
            calME->beamIndependent(parset().getBool("calibrate.ignorebeam", false));
            //
            IDataSharedIter calIter(new CalibrationIterator(it,calME));
            boost::shared_ptr<ImageFFTEquation> fftEquation(
                          new ImageFFTEquation (*itsModel, calIter, gridder()));
            ASKAPDEBUGASSERT(fftEquation);
            fftEquation->useAlternativePSF(parset());
            fftEquation->setVisUpdateObject(GroupVisAggregator::create(itsComms));
            itsEquation = fftEquation;
        }
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
            if ((ci->find("image") == 0) || (ci->find("peak_residual") == 0)) {
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
        casa::Timer timer;
        timer.mark();
        Quality q;
        ASKAPDEBUGASSERT(itsModel);
        itsSolver->solveNormalEquations(*itsModel,q);
        ASKAPLOG_INFO_STR(logger, "Solved normal equations in "<< timer.real() << " seconds "
                           );
        
        // we will probably send all of them out in the future, but for now
        // let's extract the largest residual
        const std::vector<std::string> peakParams = itsModel->completions("peak_residual.",true);
        
        double peak = peakParams.size() == 0 ? getPeakResidual() : -1.;
        for (std::vector<std::string>::const_iterator peakParIt = peakParams.begin();
             peakParIt != peakParams.end(); ++peakParIt) {
             const double tempval = std::abs(itsModel->scalarValue("peak_residual."+*peakParIt));
             ASKAPLOG_INFO_STR(logger, "Peak residual for "<< *peakParIt << " is "<<tempval);
             if (tempval > peak) {
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
      const std::map<string, casa::Vector<double> >& dataVector = ine->dataVector();
      const std::map<string, casa::Vector<double> >& diag = ine->normalMatrixDiagonal();
      for (std::map<string, casa::Vector<double> >::const_iterator ci = dataVector.begin();
           ci!=dataVector.end(); ++ci) {
           if (ci->first.find("image") == 0) {
               // this is an image
               ASKAPASSERT(ci->second.nelements() != 0);               
               std::map<std::string, casa::Vector<double> >::const_iterator diagIt = 
                            diag.find(ci->first);
               ASKAPDEBUGASSERT(diagIt != diag.end());
               const double maxDiag = casa::max(diagIt->second);
               // hard coded at this stage
               const double cutoff=1e-2*maxDiag;
               ASKAPDEBUGASSERT(diagIt->second.nelements() == ci->second.nelements());               
               for (casa::uInt elem = 0; elem<diagIt->second.nelements(); ++elem) {
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
    /// This method is intended to be called if itsExportSensitivityImage is true.
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
      casa::Array<double> wtArr = itsModel->value(wtImage);
      const double cutoff = casa::max(wtArr) * itsExpSensitivityCutoff;
      casa::Array<double> sensitivityArr(wtArr.shape());
      
      for (scimath::MultiDimArrayPlaneIter iter(wtArr.shape()); iter.hasMore(); iter.next()) {
           const casa::Vector<double> wtPlane = iter.getPlaneVector(wtArr);
           casa::Vector<double> sensitivityPlane = iter.getPlaneVector(sensitivityArr);
           for (casa::uInt elem = 0; elem < wtPlane.nelements(); ++elem) {
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
      SynthesisParamsHelper::saveImageParameter(tempPar, outParName, outParName);      
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
        vector<string> resultimages=itsModel->names();
        bool hasWeights = false;
        for (vector<string>::const_iterator it=resultimages.begin(); it
            !=resultimages.end(); it++) {
            if (it->find("weights") == 0) {
                hasWeights = true;
            }
        }
        if (!hasWeights) {
            ASKAPDEBUGASSERT(itsSolver);            
            boost::shared_ptr<ImageSolver> image_solver = boost::dynamic_pointer_cast<ImageSolver>(itsSolver);
            ASKAPDEBUGASSERT(image_solver);
            image_solver->saveWeights(*itsModel);
            resultimages=itsModel->names();
        }
        for (vector<string>::const_iterator it=resultimages.begin(); it
            !=resultimages.end(); it++) {
            if ((it->find("image") == 0) || (it->find("psf") == 0) ||
                (it->find("weights") == 0) || (it->find("mask") == 0) ||
                (it->find("residual")==0)) {
                ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix);
                if (itsExportSensitivityImage && (it->find("weights") == 0) && (postfix == "")) {
                    makeSensitivityImage(*it); 
                }
            }
        }
        if (itsRestore && postfix == "")
        {
          ASKAPLOG_INFO_STR(logger, "Restore images and writing them to disk");
          boost::shared_ptr<ImageRestoreSolver> ir = ImageRestoreSolver::createSolver(parset().makeSubset("restore."));
          ASKAPDEBUGASSERT(ir);
          ASKAPDEBUGASSERT(itsSolver);
          // configure restore solver the same way as normal imaging solver
          boost::shared_ptr<ImageSolver> template_solver = boost::dynamic_pointer_cast<ImageSolver>(itsSolver);
          ASKAPDEBUGASSERT(template_solver);
          ImageSolverFactory::configurePreconditioners(parset(),ir);
          ir->configureSolver(*template_solver);
          ir->copyNormalEquations(*template_solver);
          Quality q;
          ir->solveNormalEquations(*itsModel,q);
          // merged image should be a fixed parameter without facet suffixes
          resultimages=itsModel->fixedNames();
          for (vector<string>::const_iterator ci=resultimages.begin(); ci!=resultimages.end(); ++ci) {
               const ImageParamsHelper iph(*ci);
               if (!iph.isFacet() && (ci->find("image") == 0)) {
                   ASKAPLOG_INFO_STR(logger, "Saving restored image " << *ci << " with name "
                               << *ci+string(".restored") );
                   SynthesisParamsHelper::saveImageParameter(*itsModel, *ci,*ci+string(".restored"));
               }
          }
        }
        ASKAPLOG_INFO_STR(logger, "Writing out additional parameters made by restore solver as images");
        vector<string> resultimages2=itsModel->names();
        for (vector<string>::const_iterator it=resultimages2.begin(); it
            !=resultimages2.end(); it++) {
            // ASKAPLOG_INFO_STR(logger, "Checking "<<*it);
            if ((it->find("psf") == 0) && (std::find(resultimages.begin(), 
                 resultimages.end(),*it) == resultimages.end())) {
                ASKAPLOG_INFO_STR(logger, "Saving " << *it << " with name " << *it+postfix );
                SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix);
            }
        }
      }
    }

  }
}
