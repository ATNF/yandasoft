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
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <Common/ParameterSet.h>
#include <fitting/INormalEquations.h>
#include <fitting/ImagingNormalEquations.h>
#include <fitting/Params.h>
#include <fitting/Solver.h>
#include <fitting/Quality.h>
#include <measurementequation/ImageSolverFactory.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <measurementequation/ImageRestoreSolver.h>
#include <measurementequation/IImagePreconditioner.h>
#include <measurementequation/WienerPreconditioner.h>
#include <measurementequation/GaussianTaperPreconditioner.h>
#include <measurementequation/ImageMultiScaleSolver.h>
#include <measurementequation/ImageParamsHelper.h>
#include <measurementequation/CalibrationApplicatorME.h>
#include <measurementequation/CalibrationIterator.h>
#include <calibaccess/CalibAccessFactory.h>
#include <casacore/casa/OS/Timer.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>
#include <measurementequation/ImageFFTEquation.h>
#include <parallel/GroupVisAggregator.h>
#include <utils/MultiDimArrayPlaneIter.h>
#include <gridding/IVisGridder.h>
#include <gridding/VisGridderFactory.h>

// Local includes


// Using
using namespace askap;
using namespace askap::accessors;
using namespace askap::cp;
using namespace askap::scimath;
using namespace askap::synthesis;

ASKAP_LOGGER(logger, ".CalcCore");

CalcCore::CalcCore(LOFAR::ParameterSet& parset,
                       askap::askapparallel::AskapParallel& comms,
                       accessors::TableDataSource ds, int localChannel)
    : ImagerParallel(comms,parset), itsParset(parset), itsComms(comms),itsData(ds),itsChannel(localChannel)
{
    /// We need to set the calibration info here
    /// the ImagerParallel constructor will do the work to
    /// obtain the itsSolutionSource - but that is a provate member of
    /// the parent class.
    /// Not sure whether to use it directly or copy it.
    const std::string solver_par = parset.getString("solver");
    const std::string algorithm_par = parset.getString("solver.Clean.algorithm", "MultiScale");
    itsSolver = ImageSolverFactory::make(parset);
    itsGridder_p = VisGridderFactory::make(parset); // this is private to an inherited class so have to make a new one
    itsRestore = itsParset.getBool("restore", false);
}
CalcCore::CalcCore(LOFAR::ParameterSet& parset,
                       askap::askapparallel::AskapParallel& comms,
                       accessors::TableDataSource ds, askap::synthesis::IVisGridder::ShPtr gdr,
                       int localChannel)
    : ImagerParallel(comms,parset), itsParset(parset), itsComms(comms),itsData(ds),itsGridder_p(gdr), itsChannel(localChannel)
{
  const std::string solver_par = parset.getString("solver");
  const std::string algorithm_par = parset.getString("solver.Clean.algorithm", "MultiScale");
  itsSolver = ImageSolverFactory::make(parset);
  itsRestore = parset.getBool("restore", false);
}

CalcCore::~CalcCore()
{

}
void CalcCore::doCalc()
{

    casa::Timer timer;
    timer.mark();

    ASKAPLOG_DEBUG_STR(logger, "Calculating NE .... for channel " << itsChannel);
    if (!itsEquation) {

        accessors::TableDataSource ds = itsData;

        // Setup data iterator

        IDataSelectorPtr sel = ds.createSelector();

        sel->chooseCrossCorrelations();
        sel << itsParset;
        sel->chooseChannels(1, itsChannel);

        IDataConverterPtr conv = ds.createConverter();
        conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO), "Hz");
        conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
        conv->setEpochFrame();

        IDataSharedIter it = ds.createIterator(sel, conv);


        ASKAPCHECK(itsModel, "Model not defined");
        ASKAPCHECK(gridder(), "Gridder not defined");
        // calibration can go below if required

        if (!getSolutionSource()) {
            ASKAPLOG_DEBUG_STR(logger,"Not applying calibration");
            ASKAPLOG_DEBUG_STR(logger, "building FFT/measurement equation" );
            boost::shared_ptr<ImageFFTEquation> fftEquation(new ImageFFTEquation (*itsModel, it, gridder()));
            ASKAPDEBUGASSERT(fftEquation);
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
            //
            IDataSharedIter calIter(new CalibrationIterator(it,calME));
            boost::shared_ptr<ImageFFTEquation> fftEquation( \
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

    ASKAPLOG_INFO_STR(logger,"Calculated normal equations in "<< timer.real()
                      << " seconds ");

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

    casa::Vector<double> diag(checkRef.normalMatrixDiagonal(names[0]));
    casa::Vector<double> dv = checkRef.dataVector(names[0]);
    casa::Vector<double> slice(checkRef.normalMatrixSlice(names[0]));
    casa::Vector<double> pcf(checkRef.preconditionerSlice(names[0]));

    ASKAPLOG_DEBUG_STR(logger, "Max data: " << max(dv) << " Max PSF: " << max(slice) << " Normalised: " << max(dv)/max(slice));

}
void CalcCore::solveNE()
{


    casa::Timer timer;
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

    double peak = peakParams.size() == 0 ? getPeakResidual() : -1.;
    for (std::vector<std::string>::const_iterator peakParIt = peakParams.begin();
            peakParIt != peakParams.end(); ++peakParIt) {
        const double tempval = std::abs(itsModel->scalarValue("peak_residual." + *peakParIt));
        if (tempval > peak) {
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
void CalcCore::writeLocalModel(const std::string &postfix) {

    ASKAPLOG_DEBUG_STR(logger, "Writing out results as images");
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

    if (itsRestore && postfix == "")
    {
        ASKAPLOG_DEBUG_STR(logger, "Restore images and writing them to disk");
        boost::shared_ptr<ImageRestoreSolver> ir = ImageRestoreSolver::createSolver(parset().makeSubset("restore."));
        ASKAPDEBUGASSERT(ir);
        ASKAPDEBUGASSERT(itsSolver);
        // configure restore solver the same way as normal imaging solver
        boost::shared_ptr<ImageSolver> template_solver = boost::dynamic_pointer_cast<ImageSolver>(itsSolver);
        ASKAPDEBUGASSERT(template_solver);
        ImageSolverFactory::configurePreconditioners(itsParset,ir);
        ir->configureSolver(*template_solver);
        ir->copyNormalEquations(*template_solver);
        Quality q;
        ir->solveNormalEquations(*itsModel,q);
        // merged image should be a fixed parameter without facet suffixes
        resultimages=itsModel->fixedNames();
        for (vector<string>::const_iterator ci=resultimages.begin(); ci!=resultimages.end(); ++ci) {
            const ImageParamsHelper iph(*ci);
            if (!iph.isFacet() && (ci->find("image") == 0)) {
                ASKAPLOG_DEBUG_STR(logger, "Saving restored image " << *ci << " with name "
                              << *ci+string(".restored") );
                SynthesisParamsHelper::saveImageParameter(*itsModel, *ci,*ci+string(".restored"));
            }
        }
    }
    ASKAPLOG_DEBUG_STR(logger, "Writing out additional parameters made by restore solver as images");
    vector<string> resultimages2=itsModel->names();
    for (vector<string>::const_iterator it=resultimages2.begin(); it
        !=resultimages2.end(); it++) {
        ASKAPLOG_DEBUG_STR(logger, "Checking "<<*it);
        if ((it->find("psf") == 0) && (std::find(resultimages.begin(),
            resultimages.end(),*it) == resultimages.end())) {
            ASKAPLOG_DEBUG_STR(logger, "Saving " << *it << " with name " << *it+postfix );
            SynthesisParamsHelper::saveImageParameter(*itsModel, *it, *it+postfix);
        }
    }

}
void CalcCore::restoreImage()
{
    ASKAPDEBUGASSERT(itsModel);
    boost::shared_ptr<ImageRestoreSolver>
    ir = ImageRestoreSolver::createSolver(itsParset.makeSubset("restore."));
    ASKAPDEBUGASSERT(ir);
    ASKAPDEBUGASSERT(itsSolver);
    // configure restore solver the same way as normal imaging solver
    boost::shared_ptr<ImageSolver>
    template_solver = boost::dynamic_pointer_cast<ImageSolver>(itsSolver);
    ASKAPDEBUGASSERT(template_solver);
    ImageSolverFactory::configurePreconditioners(itsParset, ir);
    ir->configureSolver(*template_solver);

    try {
      ir->copyNormalEquations(*template_solver);
    }
    catch (...) {
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
    vector<string> resultimages=itsModel->names();

    for (vector<string>::const_iterator ci=resultimages.begin(); ci!=resultimages.end(); ++ci) {

        ASKAPLOG_INFO_STR(logger, "Restored image " << *ci);

    }
    ASKAPDEBUGASSERT(itsModel);

}
