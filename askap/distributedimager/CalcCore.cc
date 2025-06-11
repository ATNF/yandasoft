/// @file CalcCore.cc
///
/// @copyright (c) 2009 CSIRO
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>

// Include own header file first
#include "CalcCore.h"

// System includes
#include <string>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>
#include <Common/ParameterSet.h>
#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/ImagingNormalEquations.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/Solver.h>
#include <askap/scimath/fitting/Quality.h>
#include <askap/measurementequation/ImageSolverFactory.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/ImageRestoreSolver.h>
#include <askap/measurementequation/IImagePreconditioner.h>
#include <askap/measurementequation/WienerPreconditioner.h>
#include <askap/measurementequation/GaussianTaperPreconditioner.h>
#include <askap/measurementequation/ImageMultiScaleSolver.h>
#include <askap/measurementequation/ImageParamsHelper.h>
#include <askap/measurementequation/CalibrationApplicatorME.h>
#include <askap/measurementequation/CalibrationIterator.h>
#include <askap/calibaccess/CalibAccessFactory.h>
#include <casacore/casa/OS/Timer.h>
#include <askap/dataaccess/IDataSource.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/parallel/GroupVisAggregator.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>
#include <askap/gridding/IVisGridder.h>
#include <askap/gridding/TableVisGridder.h>
#include <askap/gridding/VisGridderFactory.h>
#include <askap/imagemath/linmos/LinmosAccumulator.h>

#include <boost/optional.hpp>

// Local includes


// Using
using namespace askap;
using namespace askap::accessors;
using namespace askap::cp;
using namespace askap::scimath;
using namespace askap::imagemath;
using namespace askap::synthesis;

ASKAP_LOGGER(logger, ".CalcCore");

CalcCore::CalcCore(LOFAR::ParameterSet& parset,
                       askap::askapparallel::AskapParallel& comms,
                       accessors::IDataSource& ds, int localChannel, double frequency, bool initialiseSolver)
    : ImagerParallel(comms,parset), itsComms(comms),itsDataSource(ds),itsChannel(localChannel),itsFrequency(frequency)
{
    // MV: it is untidy to update the parset and get defaults logic implemented that way
    // in particular, the second constructor doesn't do it which could lead to bugs. Leave as is for now.

    // tell gridder it can throw the grids away if we don't need to write them out
    bool writeGrids = parset.getBool("dumpgrids",false);
    writeGrids = parset.getBool("write.grids",writeGrids); // new name
    parset.replace(LOFAR::KVpair("gridder.cleargrids",!writeGrids));
    // tell restore solver to save the raw (unnormalised, unpreconditioned) psf
    parset.replace(LOFAR::KVpair("restore.saverawpsf",writeGrids));
    // only switch on updateResiduals if we want the residuals written out
    const bool writeResiduals = parset.getBool("write.residualimage",false);
    parset.replace(LOFAR::KVpair("restore.updateresiduals",writeResiduals));
    // only switch on savepsfimage if we want the preconditioned psf written out
    const bool writePsfImage = parset.getBool("write.psfimage",false);
    parset.replace(LOFAR::KVpair("restore.savepsfimage",writePsfImage));
    //
    // MV: it's very hacky and untidy to shadow data members in the base class(es). Although using the data member in
    // the base class is also ugly. Leave as is for now
    if (initialiseSolver) {
        itsSolver = ImageSolverFactory::make(parset);
    }
    itsGridder = VisGridderFactory::make(parset); // this is private to an inherited class so have to make a new one
    itsRestore = parset.getBool("restore", false);
}
CalcCore::CalcCore(LOFAR::ParameterSet& parset,
                       askap::askapparallel::AskapParallel& comms,
                       accessors::IDataSource& ds, askap::synthesis::IVisGridder::ShPtr gdr,
                       int localChannel, double frequency, bool initialiseSolver)
    : ImagerParallel(comms,parset), itsComms(comms),itsDataSource(ds),itsGridder(gdr), itsChannel(localChannel),
      itsFrequency(frequency)
{
  if (initialiseSolver) {
      itsSolver = ImageSolverFactory::make(parset);
  }
  itsRestore = parset.getBool("restore", false);
}

/// @brief reset measurement equation
/// @details We create measurement equation (i.e. ImageFFTEquation) on demand. However, it
/// has grids which are heavy objects. This method resets the appropriate shared pointer which
/// should free up the memory.
void CalcCore::resetMeasurementEquation()
{
   // MV: it breaks encapsulation accessing the data member of the base class but the whole code was written this way
   // One day we should clean up this technical debt.
   // In principle, I could've put this method with the appropriate base class. But leave it here for now as it is only used
   // in the new imager at the moment
   itsEquation.reset();
}

/// @brief make data iterator
/// @details This helper method makes an iterator based on the configuration in the current parset and
/// data fields of this class such as itsChannel and itsFrequency
/// @return shared pointer to the iterator over original data
accessors::IDataSharedIter CalcCore::makeDataIterator() const
{
   IDataSelectorPtr sel = itsDataSource.createSelector();

   sel->chooseCrossCorrelations();
   sel << parset();

   // This is the logic that switches on the combination of channels.
   // Earlier logic has updated the Channels parameter in the parset ....
   const bool combineChannels = parset().getBool("combinechannels",false);
   const bool dopplerTracking = parset().getBool("dopplertracking",false);

   if (!combineChannels) {
       if (dopplerTracking) {
           // To allow a doppler tracking reference position to be specified we need
           // a chooseFrequencies function that takes a MFrequency with reference frame
           const std::vector<string> direction = parset().getStringVector("dopplertracking.direction",{},false);
           casacore::MFrequency::Ref freqRef = getFreqRefFrame();
           if (direction.size() == 3) {
               const casacore::MeasFrame frame(asMDirection(direction));
               freqRef.set(frame);
           }
           sel->chooseFrequencies(1, casacore::MFrequency(casacore::MVFrequency(itsFrequency),freqRef), casacore::MVFrequency(0));
       } else {
           sel->chooseChannels(1, itsChannel);
       }
   }

   IDataConverterPtr conv = itsDataSource.createConverter();
   conv->setFrequencyFrame(casacore::MFrequency::Ref(casacore::MFrequency::TOPO), "Hz");
   conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));
   conv->setEpochFrame();

   return itsDataSource.createIterator(sel, conv);
}

/// @brief iterate over data and accumulate samples for uv weights
/// @details This method is used to build the sample density in the uv-plane via the appropriate gridder
/// and weight builder class. It expects the builder already setup and accessible via the normal equations
/// shared pointer. Unlike the variant from the base class which works with the iterator supplied as a parameter,
/// this version uses the iterator returned by makeDataIterator (wrapped into the calibration adapter, if needed)
void CalcCore::accumulateUVWeights() const
{
   // technically, calibration application can alter the flags, so we have to apply calibration
   // but it is a valid short-cut to skip calibration application (and have less than ideal weights) and may be even to ignore flags completely
   // (and we can have an accessor adapter which ignores flags, it would also speed up the first iteration)
   //
   // In principle, we can build different iterators here depending on the parset and make this behaviour configurable
   boost::shared_ptr<accessors::IDataIterator> it = makeCalibratedDataIteratorIfNeeded(makeDataIterator());
   // call version of the base class
   accumulateUVWeights(it);
}

/// @brief configure normal equation for linear mosaicing
/// @details When linmos is expected to happen during merge of normal equations we need to configure
/// NEs appropriately to interpret weight correctly. This helper method does it.
/// @note Normal equations should already be setup (although could be empty) before this method is called.
/// Otherwise, an exception will be thrown. Also, we could've do this setup automatically based on the
/// gridder type. But at the moment the same approach is followed as we had prior to refactoring.
void CalcCore::configureNormalEquationsForMosaicing() const
{
   const boost::shared_ptr<ImagingNormalEquations> ine = boost::dynamic_pointer_cast<ImagingNormalEquations>(getNE());
   ASKAPCHECK(ine, "Logic error - normal equations are of the wrong type or uninitialised in CalcCore::configureNormalEquationsForMosaicing");
   // enums used below are defined in LinmosAccumulator
   ine->weightType(FROM_WEIGHT_IMAGES);
   ine->weightState(WEIGHTED);
}

/// @brief merge normal equations from another CalcCore
/// @details This is a convenience method to merge in normal equations held by other CalcCore
/// object. In principle, we can have this method in one of the base classes (and require
/// broader type rather than CalcCore as the input) because all of the required functionality is
/// in the base classes. But we only use it with CalcCore, so keep it in this class as well.
/// @note Normal equations should be initialised (and with the consistent type) in both
/// this and other CalcCore instances, but could be empty. The method is const because it doesn't change
/// this class (only changes normal equations held by pointer).
/// @param[in] other an instance of CalcCore to merge from
void CalcCore::mergeNormalEquations(const CalcCore &other) const
{
   const boost::shared_ptr<INormalEquations> thisNE = getNE();
   const boost::shared_ptr<INormalEquations> otherNE = other.getNE();
   ASKAPCHECK(thisNE && otherNE, "Logic error - either this or other CalcCore object has empty normal equations in merge");
   thisNE->merge(*otherNE);
}

/// @brief create measurement equation
/// @details This method creates measurement equation as appropriate (with calibration application or without) using
/// internal state of this class and the parset
void CalcCore::createMeasurementEquation()
{
   // Setup data iterator
   accessors::IDataSharedIter it = makeCalibratedDataIteratorIfNeeded(makeDataIterator());

   ASKAPCHECK(itsModel, "Model not defined");
   ASKAPCHECK(gridder(), "Prototype gridder not defined");

   ASKAPLOG_DEBUG_STR(logger, "building FFT/measurement equation" );
   // the ImageFFTEquation actually clones the gridders and stores them internally
   // which is good - but you do not get the expected behaviour here. You would think that
   // this gridder is the one that is being used - unfortunately it is not.
   // You therefore get no benefit from initialising the gridder.
   // Also this is why you cannot get at the grid from outside FFT equation
   // Changed itsModel argument to reference (like in doCalc)
   const boost::shared_ptr<ImageFFTEquation> fftEquation(new ImageFFTEquation (itsModel, it, gridder(), parset()));
   ASKAPDEBUGASSERT(fftEquation);

   fftEquation->configure(parset());
   fftEquation->setVisUpdateObject(GroupVisAggregator::create(itsComms));
   if (calDirMap().size()>0) {
       fftEquation->setCalDirMap(calDirMap());
   }
   // MV: it is not great that the code breaks encapsulation here by changing the data member of a base class, leave it as is for now
   itsEquation = fftEquation;
}

void CalcCore::doCalc()
{
    ASKAPTRACE("CalcCore::doCalc");

    casacore::Timer timer;
    timer.mark();

    ASKAPLOG_DEBUG_STR(logger, "Calculating NE .... for channel " << itsChannel);
    if (!itsEquation) {
        createMeasurementEquation();
    } else {
        ASKAPLOG_INFO_STR(logger, "Reusing measurement equation and updating with latest model images" );
        itsEquation->reference(itsModel);
    }
    ASKAPCHECK(itsEquation, "Equation not defined");
    ASKAPCHECK(itsNe, "NormalEquations not defined");

    itsEquation->calcEquations(*itsNe);

    ASKAPLOG_INFO_STR(logger,"Calculated normal equations in "<< timer.real()
                      << " seconds ");
}

/// @brief first image name in the model
/// @details This is a helper method to obtain the name of the first encountered image parameter in the model.
/// @note It is written as part of the refactoring of various getGrid methods. However, in principle we could have multiple
/// image parameters simultaneously. The original approach getting the first one won't work in this case.
/// @return name of the first encountered image parameter in the model
std::string CalcCore::getFirstImageName() const
{
   ASKAPASSERT(itsModel);
   const std::vector<std::string> completions(itsModel->completions("image"));
   const std::vector<std::string>::const_iterator it=completions.begin();
   ASKAPCHECK(it != completions.end(), "There are no images in the current model!");
   return "image"+(*it);
}

casacore::Array<casacore::Complex> CalcCore::getGrid() const
{
   ASKAPLOG_INFO_STR(logger,"Extracting vis grid for channel " << itsChannel);
   const boost::shared_ptr<ImageFFTEquation> fftEquation = getMeasurementEquation();
   const string imageName = getFirstImageName();
   // note, it's ok to pass null pointer to dynamic cast, no need to check it separately beforehand
   const boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(fftEquation->getResidualGridder(imageName));
   ASKAPCHECK(tvg, "Incompatible type of residual gridder is used in FFTEquation");
   return tvg->getGrid();
}

casacore::Array<casacore::Complex> CalcCore::getPCFGrid() const
{
   ASKAPLOG_INFO_STR(logger,"Extracting pcf grid for channel " << itsChannel);
   const boost::shared_ptr<ImageFFTEquation> fftEquation = getMeasurementEquation();
   const string imageName = getFirstImageName();
   // in principle, we can pass the shared pointer on interface straight to dynamic cast and test the result only
   // (it will be null pointer if getPreconGridder returns the null pointer). But doing the separate check allows us to
   // distinguish the situations when preconditioning gridder is not defined (possible use case) vs. when
   // the preconditioning gridder is of an unsupported type (logic error somewhere).
   const boost::shared_ptr<IVisGridder> gridder = fftEquation->getPreconGridder(imageName);
   if (gridder) {
       const boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(gridder);
       ASKAPCHECK(tvg, "Incompatible type of PCF gridder is used in FFTEquation");
       return tvg->getGrid();
   }

   ASKAPLOG_WARN_STR(logger,"PreconGridder not defined, make sure preservecf is set to true");
   return casacore::Array<casacore::Complex>();
}

casacore::Array<casacore::Complex> CalcCore::getPSFGrid() const
{
   ASKAPLOG_INFO_STR(logger,"Extracting psf grid for channel " << itsChannel);
   const boost::shared_ptr<ImageFFTEquation> fftEquation = getMeasurementEquation();
   const string imageName = getFirstImageName();
   // note, it's ok to pass null pointer to dynamic cast, no need to check it separately beforehand
   const boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(fftEquation->getPSFGridder(imageName));
   ASKAPCHECK(tvg, "Incompatible type of PSF gridder is used in FFTEquation");
   return tvg->getGrid();
}

/// @brief store all complex grids in the model object for future writing
/// @details This method calls getGrid, getPCFGrid and getPSFGrid and stores
/// returned arrays in the model so they can be exported later. If the model
/// object already has grids, the new values are added. Shape must conform.
/// @param[in] storage shared pointer to the model where grids will be stored
void CalcCore::addGridsToModel(const boost::shared_ptr<scimath::Params> &storage)
{
   ASKAPLOG_INFO_STR(logger,"Adding grid.slice");
   ASKAPASSERT(storage);
   casacore::Array<casacore::Complex> garr = getGrid();
   casacore::Vector<casacore::Complex> garrVec(garr.reform(casacore::IPosition(1,garr.nelements())));
   if (storage->has("grid.slice")) {
       // MV: probably unnecessary complex <-> two floats conversion underneath. Leave as is for now
       // but in general we seem to be doing unnecessary copy a lot with the current Params class
       garrVec += storage->complexVectorValue("grid.slice");
       storage->updateComplexVector("grid.slice",garrVec);
   } else {
      storage->addComplexVector("grid.slice",garrVec);
   }
   ASKAPLOG_INFO_STR(logger,"Adding pcf.slice");
   casacore::Array<casacore::Complex> pcfarr = getPCFGrid();
   if (pcfarr.nelements()) {
       ASKAPLOG_INFO_STR(logger,"Adding pcf.slice");
       casacore::Vector<casacore::Complex> pcfVec(pcfarr.reform(casacore::IPosition(1,pcfarr.nelements())));
       if (storage->has("pcf.slice")) {
           pcfVec += storage->complexVectorValue("pcf.slice");
           storage->updateComplexVector("pcf.slice",pcfVec);
       } else {
         storage->addComplexVector("pcf.slice",pcfVec);
       }
   }
   ASKAPLOG_INFO_STR(logger,"Adding psfgrid.slice");
   casacore::Array<casacore::Complex> psfarr = getPSFGrid();
   casacore::Vector<casacore::Complex> psfVec(psfarr.reform(casacore::IPosition(1,psfarr.nelements())));
   if (storage->has("psfgrid.slice")) {
      psfVec += storage->complexVectorValue("psfgrid.slice");
      storage->updateComplexVector("psfgrid.slice",psfVec);
   } else {
      storage->addComplexVector("psfgrid.slice",psfVec);
   }
}

void CalcCore::calcNE()
{

    init();

    doCalc();

}

void CalcCore::zero() const {

  ASKAPCHECK(itsNe, "Normal equations are not setup");
  // the following would throw bad_cast exception if the wrong type is used
  ImagingNormalEquations &zeroRef =
  dynamic_cast<ImagingNormalEquations&>(*itsNe);

  zeroRef.zero(*itsModel);
}

void CalcCore::updateSolver() const {
  ASKAPLOG_INFO_STR(logger,"Updating the Ne in the solver with the current NE set");
  ASKAPCHECK(itsSolver, "Solver uninitialised inside CalcCore::updateSolver");
  itsSolver->init();
  itsSolver->addNormalEquations(*itsNe);
}

void CalcCore::init()
{
  // MV - leave as is for now, but reset would throw an exception if itsNe is uninitialised, so the following if-statement is redundant

  reset();

  if (!itsNe) {
      recreateNormalEquations();
  }
  ASKAPCHECK(gridder(), "Gridder not defined");
  ASKAPCHECK(itsModel, "Model not defined");
  ASKAPCHECK(itsNe, "NormalEquations not defined");

}

void CalcCore::reset() const
{
    ASKAPLOG_DEBUG_STR(logger,"Reset NE");
    ASKAPCHECK(itsNe, "Normal equations are not setup inside CalcCore::reset");
    itsNe->reset();
    ASKAPLOG_DEBUG_STR(logger,"Reset NE - done");
}

void CalcCore::check() const
{
    ASKAPCHECK(itsNe, "Normal equations are not defined inside CalcCore::check");
    std::vector<std::string> names = itsNe->unknowns();
    const ImagingNormalEquations &checkRef =
    dynamic_cast<const ImagingNormalEquations&>(*itsNe);
    ASKAPCHECK(names.size() > 0, "CalcCore::check has been called for an empty NE");

    casacore::Vector<imtype> diag(checkRef.normalMatrixDiagonal(names[0]));
    casacore::Vector<imtype> dv = checkRef.dataVectorT(names[0]);
    casacore::Vector<imtype> slice(checkRef.normalMatrixSlice(names[0]));
    casacore::Vector<imtype> pcf(checkRef.preconditionerSlice(names[0]));

    ASKAPLOG_DEBUG_STR(logger, "Max data: " << max(dv) << " Max PSF: " << max(slice) << " Normalised: " << max(dv)/max(slice)<<" ("<<names[0]<<")");
}
void CalcCore::solveNE()
{
    casacore::Timer timer;
    timer.mark();

    ASKAPCHECK(itsSolver, "Solver is not defined in solveNE!");
    itsSolver->init();
    itsSolver->addNormalEquations(*itsNe);

    ASKAPLOG_DEBUG_STR(logger, "Solving Normal Equations");
    askap::scimath::Quality q;

    ASKAPDEBUGASSERT(itsModel);
    itsSolver->solveNormalEquations(*itsModel, q);
    ASKAPLOG_INFO_STR(logger, "Solved normal equations in " << timer.real()
                       << " seconds ");

    // Extract the largest residual
    const std::vector<std::string> peakParams = itsModel->completions("peak_residual.",true);

    // note we use a negative peak val to signal deconvolution divergence and pass that on here
    double peak = peakParams.size() == 0 ? getPeakResidual() : -1.e-10;
    for (std::vector<std::string>::const_iterator peakParIt = peakParams.begin();
            peakParIt != peakParams.end(); ++peakParIt) {
        const double tempval = itsModel->scalarValue("peak_residual." + *peakParIt);
        if (std::abs(tempval) > std::abs(peak)) {
            peak = tempval;
        }
    }

    if (itsModel->has("peak_residual")) {
        itsModel->update("peak_residual", peak);
    } else {
        itsModel->add("peak_residual", peak);
    }
    itsModel->fix("peak_residual");

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

// This code is not called from anywhere at present
void CalcCore::writeLocalModel(const std::string &postfix) const {

    ASKAPLOG_DEBUG_STR(logger, "Writing out results as images");
    ASKAPDEBUGASSERT(itsModel);
    std::vector<std::string> resultimages=itsModel->names();
    bool hasWeights = false;
    for (std::vector<std::string>::const_iterator it=resultimages.begin(); it
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

    // Check whether or not the model has been stored at a higher resolution
    boost::optional<float> extraOSfactor;
    if (parset().isDefined("Images.extraoversampling")) {
        extraOSfactor = parset().getFloat("Images.extraoversampling");
        ASKAPDEBUGASSERT(*extraOSfactor > 1.);
    }

    if (itsRestore && postfix == "")
    {
        ASKAPLOG_DEBUG_STR(logger, "Restore images and writing them to disk");
        boost::shared_ptr<ImageRestoreSolver> ir = ImageRestoreSolver::createSolver(parset().makeSubset("restore."));
        ASKAPDEBUGASSERT(ir);
        ASKAPDEBUGASSERT(itsSolver);
        // configure restore solver the same way as normal imaging solver
        if (extraOSfactor) {
            ASKAPLOG_INFO_STR(logger,
                "Configuring restore solver with an extra oversampling factor of "<<*extraOSfactor);
            ir->setExtraOversampling(*extraOSfactor);
        }
        boost::shared_ptr<ImageSolver> template_solver = boost::dynamic_pointer_cast<ImageSolver>(itsSolver);
        ASKAPDEBUGASSERT(template_solver);
        ir->configureSolver(*template_solver);
        ir->copyNormalEquations(*template_solver);
        Quality q;
        ir->solveNormalEquations(*itsModel,q);
        // merged image should be a fixed parameter without facet suffixes
        resultimages=itsModel->fixedNames();
        for (std::vector<std::string>::const_iterator ci=resultimages.begin(); ci!=resultimages.end(); ++ci) {
            const ImageParamsHelper iph(*ci);
            if (extraOSfactor) {
                if (!iph.isFacet() && (ci->find("fullres") == 0)) {
                    string tmpname = *ci;
                    tmpname.replace(0,7,"image");
                    ASKAPLOG_DEBUG_STR(logger, "Saving restored image " << *ci << " with name "
                                  << tmpname+string(".restored") );
                    SynthesisParamsHelper::saveImageParameter(*itsModel, *ci, tmpname+string(".restored"));
                }
            } else {
                if (!iph.isFacet() && (ci->find("image") == 0)) {
                    ASKAPLOG_DEBUG_STR(logger, "Saving restored image " << *ci << " with name "
                                  << *ci+string(".restored") );
                    SynthesisParamsHelper::saveImageParameter(*itsModel, *ci, *ci+string(".restored"));
                }
            }
        }
    }
    ASKAPLOG_DEBUG_STR(logger, "Writing out additional parameters made by restore solver as images");
    std::vector<std::string> resultimages2=itsModel->names();
    for (std::vector<std::string>::const_iterator it=resultimages2.begin(); it
        !=resultimages2.end(); it++) {
        ASKAPLOG_DEBUG_STR(logger, "Checking "<<*it);
        if ((it->find("psf") == 0) && (std::find(resultimages.begin(),
            resultimages.end(),*it) == resultimages.end())) {
            ASKAPLOG_DEBUG_STR(logger, "Saving " << *it << " with name " << *it+postfix );
            SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix, extraOSfactor);
        }
    }

}
void CalcCore::restoreImage() const
{
    ASKAPDEBUGASSERT(itsModel);
    boost::shared_ptr<ImageRestoreSolver>
    ir = ImageRestoreSolver::createSolver(parset().makeSubset("restore."));
    ASKAPDEBUGASSERT(ir);
    ASKAPDEBUGASSERT(itsSolver);

    if (parset().isDefined("Images.extraoversampling")) {
        const float extraOSfactor = parset().getFloat("Images.extraoversampling");
        // The parameter should only be defined if has a legitimate value (is set by the code). Check anyway.
        ASKAPDEBUGASSERT(extraOSfactor > 1.);
        ASKAPLOG_INFO_STR(logger, "Configuring restore solver with an extra oversampling factor of "<<extraOSfactor);
        ir->setExtraOversampling(extraOSfactor);
    }
    // configure restore solver the same way as normal imaging solver
    boost::shared_ptr<ImageSolver>
    template_solver = boost::dynamic_pointer_cast<ImageSolver>(itsSolver);
    ASKAPDEBUGASSERT(template_solver);

    // Can we copy the preconditioners from itsSolver to avoid some work & memory?
    // Both Wiener and Gaussian keep a cache we should try to reuse
    // added code to configureSolver to do this.
    ir->configureSolver(*template_solver);

    try {
      ir->copyNormalEquations(*template_solver);
    }
    catch (...) {
      ASKAPLOG_WARN_STR(logger, "Adding missing normal equations for restore");
      template_solver->addNormalEquations(*itsNe);
      try {
          ir->copyNormalEquations(*template_solver);
      }
      catch (...) {
        throw;
      }
    }

    Quality q;
    ir->solveNormalEquations(*itsModel, q);
    std::vector<std::string> resultimages=itsModel->completions("image",true);

    for (std::vector<std::string>::const_iterator ci=resultimages.begin(); ci!=resultimages.end(); ++ci) {

        ASKAPLOG_INFO_STR(logger, "Restored image " << "image"+*ci);

    }
    ASKAPDEBUGASSERT(itsModel);

}

/// @brief stash current normal equations in the buffer
/// @details It simply copies shared pointer to the normal equations into itsSavedNE. Note, an exception is
/// thrown if the buffer is not empty (cross check as we currently plan to have a single-element stack).
void CalcCore::stashNormalEquations() 
{
   ASKAPCHECK(!itsSavedNE, "Logic error - attempting to stash normal equations while the buffer is not empty!");
   itsSavedNE = getNE();
}

/// @brief pop normal equations from the buffer
/// @details This method restores previously stashed normal equations. An exception is thrown if the buffer is
/// empty.
void CalcCore::popNormalEquations()
{
   ASKAPCHECK(itsSavedNE, "Logic error - attempting to pop normal equations from an empty buffer!");
   setNE(itsSavedNE);
   itsSavedNE.reset();
}
