/// @file SolverCore.cc
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
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/parallel/GroupVisAggregator.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>
#include <askap/gridding/IVisGridder.h>
#include <askap/gridding/TableVisGridder.h>
#include <askap/gridding/VisGridderFactory.h>

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
                       accessors::TableDataSource ds, int localChannel, double frequency)
    : ImagerParallel(comms,parset), itsComms(comms),itsData(ds),itsChannel(localChannel),itsFrequency(frequency)
{
    /// We need to set the calibration info here
    /// the ImagerParallel constructor will do the work to
    /// obtain the itsSolutionSource - but that is a private member of
    /// the parent class.
    /// Not sure whether to use it directly or copy it.
    const std::string solver_par = parset.getString("solver");
    const std::string algorithm_par = parset.getString("solver.Clean.algorithm", "BasisfunctionMFS");
    // tell gridder it can throw the grids away if we don't need to write them out
    bool writeGrids = parset.getBool("dumpgrids",false);
    writeGrids = parset.getBool("write.grids",writeGrids); // new name
    parset.replace(LOFAR::KVpair("gridder.cleargrids",!writeGrids));
    // tell restore solver to save the raw (unnormalised, unpreconditioned) psf
    parset.replace(LOFAR::KVpair("restore.saverawpsf",writeGrids));
    // only switch on updateResiduals if we want the residuals written out
    bool writeResiduals = parset.getBool("write.residualimage",false);
    parset.replace(LOFAR::KVpair("restore.updateresiduals",writeResiduals));
    // only switch on savepsfimage if we want the preconditioned psf written out
    bool writePsfImage = parset.getBool("write.psfimage",false);
    parset.replace(LOFAR::KVpair("restore.savepsfimage",writePsfImage));
    itsSolver = ImageSolverFactory::make(parset);
    itsGridder_p = VisGridderFactory::make(parset); // this is private to an inherited class so have to make a new one
    itsRestore = parset.getBool("restore", false);
}
CalcCore::CalcCore(LOFAR::ParameterSet& parset,
                       askap::askapparallel::AskapParallel& comms,
                       accessors::TableDataSource ds, askap::synthesis::IVisGridder::ShPtr gdr,
                       int localChannel, double frequency)
    : ImagerParallel(comms,parset), itsComms(comms),itsData(ds),itsGridder_p(gdr), itsChannel(localChannel),
      itsFrequency(frequency)
{
  const std::string solver_par = parset.getString("solver");
  const std::string algorithm_par = parset.getString("solver.Clean.algorithm", "BasisfunctionMFS");
  itsSolver = ImageSolverFactory::make(parset);
  itsRestore = parset.getBool("restore", false);
}

CalcCore::~CalcCore()
{

}
void CalcCore::doCalc()
{

    casacore::Timer timer;
    timer.mark();

    ASKAPLOG_DEBUG_STR(logger, "Calculating NE .... for channel " << itsChannel);
    if (!itsEquation) {

        accessors::TableDataSource ds = itsData;

        // Setup data iterator

        IDataSelectorPtr sel = ds.createSelector();

        sel->chooseCrossCorrelations();
        sel << parset();

        // This is the logic that switches on the combination of channels.
        // Earlier logic has updated the Channels parameter in the parset ....
        bool combineChannels = parset().getBool("combinechannels",false);
        bool dopplerTracking = parset().getBool("dopplertracking",false);

        if (!combineChannels) {
            if (dopplerTracking) {
                sel->chooseFrequencies(1, casacore::MVFrequency(itsFrequency), casacore::MVFrequency(0),
                    casacore::MFrequency::castType(getFreqRefFrame().getType()));
            } else {
                sel->chooseChannels(1, itsChannel);
            }
        }

        IDataConverterPtr conv = ds.createConverter();
        conv->setFrequencyFrame(casacore::MFrequency::Ref(casacore::MFrequency::TOPO), "Hz");
        conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));
        conv->setEpochFrame();

        IDataSharedIter it = ds.createIterator(sel, conv);


        ASKAPCHECK(itsModel, "Model not defined");
        ASKAPCHECK(gridder(), "Prototype gridder not defined");
        // calibration can go below if required

        if (!getSolutionSource()) {
            ASKAPLOG_DEBUG_STR(logger,"Not applying calibration");
            ASKAPLOG_DEBUG_STR(logger, "building FFT/measurement equation" );
            // the ImageFFTEquation actually clones the gridders and stores them internally
            // which is good - but you do not get the expected behaviour here. You would think that
            // this gridder is the one that is being used - unfortunately it is not.
            // You therefore get no benefit from initialising the gridder.
            // Also this is why you cannot get at the grid from outside FFT equation
            boost::shared_ptr<ImageFFTEquation> fftEquation(new ImageFFTEquation (*itsModel, it, gridder()));
            ASKAPDEBUGASSERT(fftEquation);
// DAM TRADITIONAL
            if (parset().isDefined("gridder.robustness")) {
                const float robustness = parset().getFloat("gridder.robustness");
                ASKAPCHECK((robustness>=-2) && (robustness<=2), "gridder.robustness should be in the range [-2,2]");
                // also check that it is spectral line? Or do that earlier
                //  - won't work in continuum imaging, unless combo is done with combinechannels on a single worker
                //  - won't work with Taylor terms
                fftEquation->setRobustness(parset().getFloat("gridder.robustness"));
            }
            fftEquation->useAlternativePSF(parset());
            fftEquation->setVisUpdateObject(GroupVisAggregator::create(itsComms));
            itsEquation = fftEquation;
        } else {
            ASKAPLOG_DEBUG_STR(logger, "Calibration will be performed using solution source");
            boost::shared_ptr<ICalibrationApplicator> calME(\
            new CalibrationApplicatorME(getSolutionSource()));
            // fine tune parameters
            ASKAPDEBUGASSERT(calME);
            calME->scaleNoise(parset().getBool("calibrate.scalenoise",false));
            calME->allowFlag(parset().getBool("calibrate.allowflag",false));
            calME->beamIndependent(parset().getBool("calibrate.ignorebeam", false));
            calME->interpolateTime(parset().getBool("calibrate.interpolatetime",false));

            //
            IDataSharedIter calIter(new CalibrationIterator(it,calME));
            boost::shared_ptr<ImageFFTEquation> fftEquation( \
            new ImageFFTEquation (*itsModel, calIter, gridder()));
            ASKAPDEBUGASSERT(fftEquation);
// DAM TRADITIONAL
            if (parset().isDefined("gridder.robustness")) {
                const float robustness = parset().getFloat("gridder.robustness");
                ASKAPCHECK((robustness>=-2) && (robustness<=2), "gridder.robustness should be in the range [-2,2]");
                fftEquation->setRobustness(parset().getFloat("gridder.robustness"));
            }
            fftEquation->useAlternativePSF(parset());
            fftEquation->setVisUpdateObject(GroupVisAggregator::create(itsComms));
            itsEquation = fftEquation;
        }

    }
    else {
        ASKAPLOG_INFO_STR(logger, "Reusing measurement equation and updating with latest model images" );
        // Try changing this to reference instead of copy - passes tests
        //itsEquation->setParameters(*itsModel);
        itsEquation->reference(itsModel);
    }
    ASKAPCHECK(itsEquation, "Equation not defined");
    ASKAPCHECK(itsNe, "NormalEquations not defined");
    itsEquation->calcEquations(*itsNe);

    ASKAPLOG_INFO_STR(logger,"Calculated normal equations in "<< timer.real()
                      << " seconds ");

}

casacore::Array<casacore::Complex> CalcCore::getGrid() {

    ASKAPCHECK(itsEquation, "Equation not defined");
    ASKAPLOG_INFO_STR(logger,"Dumping vis grid for channel " << itsChannel);
    boost::shared_ptr<ImageFFTEquation> fftEquation = boost::dynamic_pointer_cast<ImageFFTEquation>(itsEquation);

    // We will need to loop over all completions i.e. all sources
    const std::vector<std::string> completions(itsModel->completions("image"));

    std::vector<std::string>::const_iterator it=completions.begin();
    const string imageName("image"+(*it));
    boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(fftEquation->getResidualGridder(imageName));
    return tvg->getGrid();
}
casacore::Array<casacore::Complex> CalcCore::getPCFGrid() {

    ASKAPCHECK(itsEquation, "Equation not defined");
    ASKAPLOG_INFO_STR(logger,"Dumping pcf grid for channel " << itsChannel);
    boost::shared_ptr<ImageFFTEquation> fftEquation = boost::dynamic_pointer_cast<ImageFFTEquation>(itsEquation);
    // We will need to loop over all completions i.e. all sources
    const std::vector<std::string> completions(itsModel->completions("image"));

    std::vector<std::string>::const_iterator it=completions.begin();
    const string imageName("image"+(*it));
    boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(fftEquation->getPreconGridder(imageName));
    ASKAPCHECK(tvg,"PreconGridder not defined, make sure preservecf is set to true")

    return tvg->getGrid();
}
casacore::Array<casacore::Complex> CalcCore::getPSFGrid() {

    ASKAPCHECK(itsEquation, "Equation not defined");
    ASKAPLOG_INFO_STR(logger,"Dumping psf grid for channel " << itsChannel);
    boost::shared_ptr<ImageFFTEquation> fftEquation = boost::dynamic_pointer_cast<ImageFFTEquation>(itsEquation);
    // We will need to loop over all completions i.e. all sources
    const std::vector<std::string> completions(itsModel->completions("image"));

    std::vector<std::string>::const_iterator it=completions.begin();
    const string imageName("image"+(*it));
    boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(fftEquation->getPSFGridder(imageName));
    ASKAPCHECK(tvg,"PSFGridder not defined")

    return tvg->getGrid();
}
void CalcCore::calcNE()
{



    init();

    doCalc();




}
void CalcCore::zero() {

  ImagingNormalEquations &zeroRef =
  dynamic_cast<ImagingNormalEquations&>(*itsNe);

  zeroRef.zero(*itsModel);
}

void CalcCore::updateSolver() {
  ASKAPLOG_INFO_STR(logger,"Updating the Ne in the solver with the current NE set");
  itsSolver->init();
  itsSolver->addNormalEquations(*itsNe);
}
void CalcCore::init()
{

  reset();

  if (!itsNe) {
      ASKAPLOG_DEBUG_STR(logger,"Recreating NE from model");
      itsNe=ImagingNormalEquations::ShPtr(new ImagingNormalEquations(*itsModel));
      ASKAPLOG_DEBUG_STR(logger,"Done recreating model");
  }
  ASKAPCHECK(gridder(), "Gridder not defined");
  ASKAPCHECK(itsModel, "Model not defined");
  ASKAPCHECK(itsNe, "NormalEquations not defined");

}

void CalcCore::reset()
{

    ASKAPLOG_DEBUG_STR(logger,"Reset NE");
    itsNe->reset();
    ASKAPLOG_DEBUG_STR(logger,"Reset NE - done");
}

void CalcCore::check()
{
    std::vector<std::string> names = itsNe->unknowns();
    const ImagingNormalEquations &checkRef =
    dynamic_cast<const ImagingNormalEquations&>(*itsNe);

    casacore::Vector<imtype> diag(checkRef.normalMatrixDiagonal(names[0]));
    casacore::Vector<imtype> dv = checkRef.dataVectorT(names[0]);
    casacore::Vector<imtype> slice(checkRef.normalMatrixSlice(names[0]));
    casacore::Vector<imtype> pcf(checkRef.preconditionerSlice(names[0]));

    ASKAPLOG_DEBUG_STR(logger, "Max data: " << max(dv) << " Max PSF: " << max(slice) << " Normalised: " << max(dv)/max(slice));

}
void CalcCore::solveNE()
{


    casacore::Timer timer;
    timer.mark();

    itsSolver->init();
    itsSolver->addNormalEquations(*itsNe);

    ASKAPLOG_DEBUG_STR(logger, "Solving Normal Equations");
    askap::scimath::Quality q;

    ASKAPDEBUGASSERT(itsModel);
    itsSolver->solveNormalEquations(*itsModel, q);
    ASKAPLOG_DEBUG_STR(logger, "Solved normal equations in " << timer.real()
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

}

// This code is not called from anywhere at present
void CalcCore::writeLocalModel(const std::string &postfix) {

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
void CalcCore::restoreImage()
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
