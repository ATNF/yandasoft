/// @file
///
/// ContSubtractParallel: Support for parallel continuum subtraction using model
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>
///

#include <askap/parallel/ContSubtractParallel.h>
#include <askap/parallel/SimParallel.h>
#include <askap/scimath/fitting/NormalEquationsStub.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/dataaccess/TableDataIterator.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/dataaccess/DDCalBufferDataAccessor.h>
#include <askap/dataaccess/DataIteratorStub.h>
#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/scimath/fitting/Equation.h>
#include <askap/measurementequation/ComponentEquation.h>
#include <askap/measurementequation/IMeasurementEquation.h>
#include <askap/measurementequation/ImagingEquationAdapter.h>
#include <askap/measurementequation/CalibrationME.h>
#include <askap/measurementequation/CalibParamsMEAdapter.h>
#include <askap/measurementequation/CalibrationApplicatorME.h>
#include <askap/calibaccess/CalibAccessFactory.h>
#include <askap/measurementequation/NoXPolGain.h>

#include <askap/askap/AskapError.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/utils/TilingUtils.h>
#include <casacore/scimath/Fitting/LinearFitSVD.h>

#include <casacore/casa/Arrays/ArrayMath.h>

// logging stuff
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".contsubtract");

#include <casacore/casa/OS/Timer.h>

#include <vector>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

/// @brief Constructor from ParameterSet
/// @details The parset is used to construct the internal state. We could
/// also support construction from a python dictionary (for example).
/// The command line inputs are needed solely for MPI - currently no
/// application specific information is passed on the command line.
/// @param comms communication object
/// @param parset ParameterSet for inputs
ContSubtractParallel::ContSubtractParallel(askapparallel::AskapParallel& comms,
      const LOFAR::ParameterSet& parset) : MEParallelApp(comms,parset,true)
{
  // the stub allows to reuse MEParallelApp code although we're not solving
  // for the normal equations here
  itsNe.reset(new scimath::NormalEquationsStub);

  itsModelReadByMaster = parset.getBool("modelReadByMaster", true);
  itsDoUVlin = parset.getBool("doUVlin", false);
  itsOrder = parset.getInt("uvlin.order", 1);
  itsHarmonic = parset.getInt("uvlin.harmonic", 1);
  itsWidth = parset.getInt("uvlin.width",0); // 0 = whole spectrum
  itsOffset = min(max(0,parset.getInt("uvlin.offset",0)),itsWidth);
  itsThreshold = max(0.0f,parset.getFloat("uvlin.threshold",2.5));
  const std::vector<std::string> dir = parset.getStringVector("uvlin.direction",{},true);
  itsRotate = itsDoUVlin && (dir.size() > 0);
  if (itsRotate) {
      if (dir.size()==1) {
        // assume we have something like "SUN" that doesn't need lat,long
        casacore::MDirection::Types type;
        casacore::MDirection::getType(type, dir[0]);
        itsUVlinDirection = casacore::MDirection(type);
      } else if (dir.size()==3) {
        itsUVlinDirection = asMDirection(dir);
    } else {
        ASKAPLOG_WARN_STR(logger,"uvlin.direction specified incorrectly - phase rotation disabled");
        itsRotate = false;
    }
  }
  if (itsWidth > 0) ASKAPCHECK(itsOffset < itsWidth,"The offset needs to be less than the width");
  if (itsDoUVlin) ASKAPLOG_INFO_STR(logger, "Doing uvlin operation with order = "
    << itsOrder<<", harmonic = "<<itsHarmonic<< ", width = "<< itsWidth
    <<" channels, offset = "<< itsOffset << " and threshold = "<< itsThreshold);
  itsDoSubtraction = parset.getBool("doSubtraction", true);
  itsDoReplaceByModel = parset.getBool("doReplaceByModel", false);
}

/// @brief Initialise continuum subtractor
/// @details The parameters are taken from the parset file supplied in the constructor.
/// This method does initialisation which may involve communications in the parallel case
/// (i.e. distribution of the models between workers). Technically, we could've done this in
/// the constructor.
void ContSubtractParallel::init()
{
    if (!itsDoSubtraction) {
        return;
    }
    if (itsModelReadByMaster) {
        if (itsComms.isMaster()) {
            readModels();
            broadcastModel();
        }
        if (itsComms.isWorker()) {
            receiveModel();
        }
    } else if (doWork()) {
      readModels();
    }
    bool doDDCal = parset().getBool("calibrate.directiondependent",false);
    // Get number of sources with components
    itsNDir = params()->completions("sourceID").size();
    if (itsNDir > 0) {
        ASKAPLOG_INFO_STR(logger, "Found "<<itsNDir<< " component models");
    } else {
        //no component models, try image model
        itsNDir = params()->completions("image").size();
        if (itsNDir > 0) {
            ASKAPLOG_INFO_STR(logger, "Found "<<itsNDir<< " image models "<<params()->names());
        }
    }
    if (!doDDCal || itsNDir < 1) {
        itsNDir = 1;
    }
    if (doDDCal) {
        ASKAPLOG_INFO_STR(logger, "Found "<<itsNDir<< " directions for DDCAL");
    }
}

/// @brief initialise measurement equation
/// @details This method initialises measurement equation
void ContSubtractParallel::initMeasurementEquation(IDataSharedIter& it)
{
   ASKAPLOG_INFO_STR(logger, "Creating measurement equation" );

   // it doesn't matter which iterator to pass to the measurement equations
   // as we're only using accessor-based interface
   IDataSharedIter stubIter(new DataIteratorStub(1));

   ASKAPCHECK(params(), "Model is not defined");

   // a part of the equation defined via image
   askap::scimath::Equation::ShPtr imgEquation;

   if (SynthesisParamsHelper::hasImage(params())) {
       ASKAPLOG_INFO_STR(logger, "Sky model contains at least one image, building an image-specific equation");
       // it should ignore parameters which are not applicable (e.g. components)
       ASKAPCHECK(gridder(), "Gridder is not defined");
       imgEquation.reset(new ImageFFTEquation(*params(), stubIter, gridder()));
   }

   // a part of the equation defined via components
   boost::shared_ptr<ComponentEquation> compEquation;

   if (SynthesisParamsHelper::hasComponent(params())) {
       // model is a number of components
       ASKAPLOG_INFO_STR(logger, "Sky model contains at least one component, building a component-specific equation");
       // it doesn't matter which iterator is passed below. It is not used
       // it should ignore parameters which are not applicable (e.g. images)
       compEquation.reset(new ComponentEquation(*params(), stubIter));
       if (itsNDir > 1 && !SynthesisParamsHelper::hasImage(params())) {
           compEquation->setNDir(itsNDir);
       }

   }

   if (imgEquation && !compEquation) {
       ASKAPLOG_INFO_STR(logger, "Pure image-based model (no components defined)");
       itsEquation = imgEquation;
   } else if (compEquation && !imgEquation) {
       ASKAPLOG_INFO_STR(logger, "Pure component-based model (no images defined)");
       itsEquation = compEquation;
   } else if (imgEquation && compEquation) {
       ASKAPLOG_INFO_STR(logger, "Making a sum of image-based and component-based equations");
       itsEquation = imgEquation;
       SimParallel::addEquation(itsEquation, compEquation, stubIter);
   } else {
       ASKAPTHROW(AskapError, "No sky models are defined");
   }

   // we need accessor-based equation for the actual iteration. Make it now, if necessary.
   // Note, that component-based and composite equations are already accessor-based, so no
   // additional fiddling  is needed

   boost::shared_ptr<IMeasurementEquation> accessorBasedEquation =
        boost::dynamic_pointer_cast<IMeasurementEquation>(itsEquation);

   if (!accessorBasedEquation) {
       ASKAPLOG_INFO_STR(logger,"Creating accessor based equation");
        // form a replacement equation first
        const boost::shared_ptr<ImagingEquationAdapter> new_equation(new ImagingEquationAdapter);
        // the actual equation (from itsEquation) will be locked inside ImagingEquationAdapter
        // in a shared pointer. We can change itsEquation after the following line
        new_equation->assign(itsEquation);
        if (itsNDir > 1) {
            new_equation->setNDir(itsNDir);
        }

        // replacing the original equation with an accessor-based adapter
        itsEquation = new_equation;
        // this should work now
        accessorBasedEquation = boost::dynamic_pointer_cast<IMeasurementEquation>(itsEquation);
   }
   ASKAPDEBUGASSERT(accessorBasedEquation);
   if (parset().getBool("calibrate",false)) {
       itsSolutionSource = CalibAccessFactory::roCalSolutionSource(parset());
   }
   // corrupt model with calibration
   if (itsSolutionSource) {
       if (parset().getBool("calibrate.usecalapplicator",true)) {
           ASKAPLOG_INFO_STR(logger, "Using CalibrationApplicator for gains");
           itsCalApplicator.reset(new CalibrationApplicatorME(itsSolutionSource));
           itsCalApplicator->scaleNoise(parset().getBool("calibrate.scalenoise",false));
           itsCalApplicator->allowFlag(parset().getBool("calibrate.allowflag",false));
           itsCalApplicator->beamIndependent(parset().getBool("calibrate.ignorebeam", false));
           itsCalApplicator->channelIndependent(parset().getBool("calibrate.ignorechannel", false));
           itsCalApplicator->interpolateTime(parset().getBool("calibrate.interpolatetime",false));
       } else {
           ASKAPCHECK(itsNDir==1,"Use calibrate.usecalapplicator=true for DD calibration");
           boost::shared_ptr<CalibrationMEBase> calME;
           // Using the CalibrationME class is more general, but slower than the
           // CalibratorApplicatorME class
           // Not sure how to fit the DDCal classes into this
           ASKAPLOG_INFO_STR(logger, "Only parallel-hand gains will be applied. Polarisation leakage will not be applied.");
           calME.reset(new CalibrationME<NoXPolGain>(scimath::Params(), stubIter, accessorBasedEquation));
           ASKAPDEBUGASSERT(calME);
           // set up the adapter
           itsEquation.reset(new CalibParamsMEAdapter(calME, itsSolutionSource, stubIter));
       }
   }

}

void ContSubtractParallel::modelSpectrum(casacore::Vector<casacore::Float> & model,
        const casacore::Vector<casacore::Float>& spec, const casacore::Vector<casacore::Bool>& mask)
{
    const int nChan= spec.size();
    const int nParams = itsOrder+1+itsHarmonic*2;
    casacore::LSQaips fitter(nParams);
    casacore::Vector<casacore::Double> xx(nParams);
    casacore::Vector<casacore::Double> solution(nParams);
    casacore::VectorSTLIterator<casacore::Double> it(xx);
    casacore::Vector<casacore::Bool> tmask(nChan);

    // If we are doing outlier rejection against the model iterate a few times
    const int niter = (itsThreshold > 0 ? 3 : 1);
    // initial mask is the data flags
    tmask = mask;
    // initial model is zero
    model = 0;

    // we are doing the fitting in channel bins
    if (itsWidth == 0) itsWidth = nChan;
    casacore::Vector<casacore::Float> y(itsWidth);
    std::vector<casacore::Float> tmp;
    uint lastDof = nParams;
    for (int binStart=-itsOffset; binStart<nChan; binStart+=itsWidth) {
        int start = max(0,binStart);
        int end = min(binStart+itsWidth,nChan);
        int binWidth = end - start;
        for (int iter = 0; iter<niter; iter++) {
            // do thresholding of values before fitting?
            if (itsThreshold > 0) {
                // count how many valid values
                int n = 0;
                // collect valid values
                for (int i=start; i < end; i++) {
                    if (tmask(i)) {
                        y(n++) = spec(i) - model(i);
                    }
                }
                // work out robust sigma and median
                if (n>0) {
                    const casacore::Float q25 = casacore::fractile(y(casacore::Slice(0,n)), tmp, 0.25f, casacore::False, casacore::True);
                    const casacore::Float q50 = casacore::fractile(y(casacore::Slice(0,n)), tmp, 0.50f, casacore::False, casacore::True);
                    const casacore::Float q75 = casacore::fractile(y(casacore::Slice(0,n)), tmp, 0.75f, casacore::False, casacore::True);
                    const casacore::Float sigma = (q75-q25)/1.35; // robust sigma estimate
                    int count = 0;
                    // flag outliers
                    for (int i=start; i < end; i++) {
                        if (mask(i) && (abs((spec(i)-model(i)) - q50) > itsThreshold * sigma)) {
                            tmask(i) = false;
                            count++;
                        } else {
                            tmask(i) = mask(i);
                        }
                    }
                    // extend the mask by a few channels - disabled for now
                    const int nextend = 0;
                    if (nextend > 0) {
                        bool clip = false;
                        int first = end, last = 0;
                        int i = start;
                        while (i < end) {
                            if (mask(i) && !tmask(i)) {
                                if (!clip) first = i;
                                clip = true;
                                i++;
                            } else {
                                if (clip) last = i - 1;
                                clip = false;
                                if (last - first > 2) {
                                    for (int k=0; k<nextend; k++) {
                                        if (first-1-k >= start) tmask(first-1-k) = false;
                                        if (i+k < end) tmask(i+k) = false;
                                    }
                                    i+=nextend;
                                    first = end;
                                } else {
                                    i++;
                                }
                            }
                        }
                    }

                    int count2 = 0;
                    for (int i=start; i<end; i++) {
                        if (mask(i) && !tmask(i)) {
                            count2++;
                        }
                    }
                    if (iter > 0 && count2 == 0) {
                        break;
                        // no further change expected
                    }
                }
            }

            // apply external mask - would need to read this from file
            //for (int i=start; i<end; i++) if (i > 40 && i < 60)  tmask(i) = false;

            // We may need to limit #degrees of freedom (like miriad uvlin) if there are large gaps.
            // Higher orders blow up quicker for given gap size; order<=1 is safest with large gaps.
            // If valid channels < 60% -> reduce #dof, if valid channels < 5*(#dof) -> reduce dof
            int valid = 0;
            for (int i=start; i< end; i++) if (tmask(i)) valid++;
            int order = itsOrder;
            int harm = itsHarmonic;
            if (float(valid) < 0.6*binWidth) {
                if (order > harm) {
                    order=max(0,order-1);
                } else {
                    harm=max(0,harm-1);
                }
            }
            uint dof = order + 1 + 2 * harm;
            while (valid < 5 * dof) {
                if (order > harm) {
                    order=max(0,order-1);
                } else {
                    harm=max(0,harm-1);
                }
                dof = order + 1 + 2 * harm;
                if (dof == 1) break;
            }
            // make sure fitter is ready for new fit
            if (lastDof != dof) {
                fitter.set(dof);
                lastDof = dof;
            } else {
                fitter.reset();
            }
            // set the basis functions: polynomial and sin, cos terms
            for (int i=start; i < end; i++) {
                // we could use itsWidth instead of binWidth to keep 'frequency' the same for sine
                const float x = (i-start) / float(binWidth);
                if (tmask(i)) {
                    xx(0) = 1;
                    for (int j=1; j<order+1; j++) {
                        xx(j) = xx(j-1) * x;
                    }
                    for (int j=0; j<harm; j++) {
                        xx(order+1+2*j)   = sin((j+1)*casacore::C::pi*x);
                        xx(order+1+2*j+1) = cos((j+1)*casacore::C::pi*x);
                    }
                    fitter.makeNorm(it,1.0,casacore::Double(spec(i)));
                }
            }

            casacore::uInt nr1;
            const casacore::Bool ok = fitter.invert(nr1);
            if (ok) {
                fitter.solve(solution.data());
                //casacore::Float chisq = fitter.getChi();
                //casacore::Float sd1 = fitter.getSD();
                //casacore::cerr << "Fit="<< ok <<", rank="<<nr1<<", chisq="<<chisq<<", sd="<<sd1<<", sol"<<solution<<casacore::endl;

                // evaluate the solution to generate the model
                for (int i=start; i < end; i++) {
                    const float x = (i-start) / float(binWidth);
                    model(i) = solution(order);
                    for (int j=order-1; j>=0; j--) {
                        model(i) = x * model(i) + solution(j);
                    }
                    for (int j=0; j<harm; j++) {
                        model(i) += solution(order+1+2*j)   * sin((j+1)*casacore::C::pi*x);
                        model(i) += solution(order+1+2*j+1) * cos((j+1)*casacore::C::pi*x);
                    }
                }
            } else {
                break; // no point iterating if the fit failed, keep last model or zero
            }
        }
    }
}

void ContSubtractParallel::subtractContFit(casacore::Cube<casacore::Complex>& vis,
        const casacore::Cube<casacore::Bool>& flag, const casacore::Matrix<casacore::Complex>& phasor) {
    const int nPol = vis.shape()(0);
    const int nChan = vis.shape()(1);
    const int nRow = vis.shape()(2);
    casacore::Vector<casacore::Float> visreal(nChan), visimag(nChan), modelreal(nChan), modelimag(nChan);
    casacore::Vector<casacore::Bool> mask(nChan);
    const bool rotate = phasor.nelements() > 0;
    ASKAPDEBUGASSERT(!rotate || (phasor.ncolumn()==nRow && phasor.nrow()==nChan));
    for (int row=0; row<nRow; row++) {
        for (int pol=0; pol<nPol; pol++) {
            for (int chan=0; chan<nChan; chan++) {
                casacore::Complex v = vis(pol,chan,row);
                if (rotate) {
                    v *= phasor(chan,row);
                }
                visreal(chan) = casacore::real(v);
                visimag(chan) = casacore::imag(v);
                mask(chan) = !flag(pol,chan,row);
            }
            modelSpectrum(modelreal,visreal,mask);
            modelSpectrum(modelimag,visimag,mask);
            for (int chan=0; chan<nChan; chan++) {
                vis(pol,chan,row) -= casacore::Complex(modelreal(chan),modelimag(chan));
                if (rotate) {
                    vis(pol,chan,row) *= conj(phasor(row,chan));
                }
            }
        }
    }
}

void ContSubtractParallel::computePhasor(const accessors::IDataSharedIter& it,
    casacore::Matrix<casacore::Complex>& phasor)
{
    casacore::MDirection newDir = itsUVlinDirection;
    if (itsUVlinDirection.getRef().getType() != casacore::MDirection::J2000) {
        // need to do a conversion since uvwRotationDelay doesn't deal with non J2000
        const casacore::MEpoch epoch(casacore::MVEpoch(it->time()/casacore::C::day),
            casacore::MEpoch::Ref(casacore::MEpoch::UTC));
        // where to get position? only avalable inside table iterator classes
        const boost::shared_ptr<accessors::TableConstDataIterator> tableIt =
            it.dynamicCast<accessors::TableConstDataIterator>();
        if (!it) {
            ASKAPTHROW(AskapError, "Bad cast in ContSubtractParallel::computePhasor, most likely this means "
                   "there is a logical error");
        }
        const casacore::MPosition mroPos = tableIt->subtableInfo().getAntenna().getPosition(0);
        const casacore::MeasFrame frame(epoch, mroPos);
        const casacore::MDirection j2000dir = casacore::MDirection::Convert(newDir,
            casacore::MDirection::Ref(casacore::MDirection::J2000, frame))();
        static bool once = true;
        if (once) {
            ASKAPLOG_INFO_STR(logger, "Input direction: "<<newDir.toString()<<" -> J2000 : "
                <<j2000dir.toString());
            once = false;
        }

        newDir = j2000dir;
    }
    const casacore::Vector<double> &delay = it->uvwRotationDelay(newDir, newDir);
    const uint nRow = it->nRow();
    const uint nChan = it->nChannel();
    const casacore::Vector<casacore::Double>& freq = it->frequency();
    phasor.resize(nChan, nRow);
    for (uint row = 0; row < nRow; row++) {
        for (uint chan = 0; chan < nChan; chan++) {
            // Calculate the delay phasor - note delay is in meters
            const double phase = casacore::C::_2pi / casacore::C::c * freq[chan] * delay(row);
            phasor(chan, row) = casacore::Complex(cos(phase), sin(phase));
        }
    }
}

/// @brief perform the subtraction for the given dataset
/// @details This method iterates over the given dataset, predicts visibilities according to the
/// model and subtracts these model visibilities from the original visibilities in the dataset.
/// This is the core operation of the doSubtraction method, which manages the parallel aspect of it.
/// All actual calculations are done inside this helper method.
/// @param[in] ms measurement set name
/// @param[in] distributeByTile do distribution by tile if possible
void ContSubtractParallel::calcOne(const std::string &ms, bool distributeByTile)
{
    casacore::Timer timer;
    timer.mark();
    boost::shared_ptr<IMeasurementEquation> accessorBasedEquation;
    ASKAPLOG_INFO_STR(logger, "Performing continuum model subtraction for " << ms );

    // Open readonly, accessor will reopen table r/w when needed
    TableDataSource ds(ms, TableDataSource::MEMORY_BUFFERS | TableDataSource::WRITE_DATA_ONLY,
         dataColumn());
    ds.configureUVWMachineCache(uvwMachineCacheSize(),uvwMachineCacheTolerance());
    IDataSelectorPtr sel=ds.createSelector();
    if (distributeByTile) {
        utils::distributeByTile(sel, dataColumn(),nWorkers(),workerRank());
    }
    sel << parset();
    IDataConverterPtr conv=ds.createConverter();
    conv->setFrequencyFrame(getFreqRefFrame(), "Hz");
    conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));
    IDataSharedIter it=ds.createIterator(sel, conv);

    if (itsDoSubtraction) {
        if (!itsEquation) {
            initMeasurementEquation(it);
        } else {
            ASKAPLOG_INFO_STR(logger, "Reusing measurement equation" );
        }
        accessorBasedEquation = boost::dynamic_pointer_cast<IMeasurementEquation>(itsEquation);
        ASKAPDEBUGASSERT(accessorBasedEquation);
    }

    uint niter = 0;
    casacore::Matrix<casacore::Complex> phasor;

    for (; it.hasMore(); it.next()) {
        // iteration over the dataset
        niter++;
        casacore::Cube<casacore::Complex>& vis = it->rwVisibility();
        if (itsDoSubtraction) {
            DDCalBufferDataAccessor acc(*it);
            // acc.rwVisibility().set(0.); // already done in predict
            accessorBasedEquation->predict(acc);
            const casacore::Cube<casacore::Complex>& model = acc.visibility();

            // If using the cal applicator, we apply calibration separately
            // from predicting the visibilities, otherwise it is already done
            if (itsCalApplicator) {
                itsCalApplicator->predict(acc);
            }
            ASKAPDEBUGASSERT(model.nrow() == vis.nrow());
            ASKAPDEBUGASSERT(model.ncolumn() == vis.ncolumn());
            if (itsNDir == 1) {
                if (itsDoReplaceByModel) {
                    vis = model;
                } else {
                    vis -= model;
                }
            } else {
                const casacore::Slice all;
                // loop over directions
                for (int dir = 0; dir < itsNDir; dir++) {
                    const auto nrow = vis.nplane();
                    const casacore::Slice rowSlice(nrow * dir, nrow);
                    vis -= model(all,all,rowSlice);
                }
            }
        }
        // Do the uvlin fit/subtract on the residuals after the model subtraction
        if (itsDoUVlin) {
            if (itsRotate) {
                computePhasor(it, phasor);
            }
            subtractContFit(vis, it->flag(), phasor);
        }
    }

    ASKAPLOG_INFO_STR(logger, "Finished continuum subtraction for "<< ms << " in "<< timer.real()
    << " seconds, using "<<niter <<" iterations");
}

/// @brief perform the subtraction
/// @details This method iterates over one or more datasets, predicts visibilities according to
/// the model and subtracts these model visibilities from the original visibilities in the
/// dataset. In parallel mode with a single dataset we can optionally distribute work over tiles.
void ContSubtractParallel::doSubtraction()
{
    if (itsComms.isParallel()) {
        if (doWork()) {
            const uint rank = workerRank();
            ASKAPLOG_INFO_STR(logger, "Worker "<<rank<< " is processing "<<measurementSets()[rank]);
            // do automatic distribution over tiles if requested and
            //   if all measurementset names are the same
            const std::vector<std::string>& v = measurementSets();
            const bool distributeByTile = parset().isDefined("Tiles") &&
                parset().getString("Tiles")=="auto" &&
                (std::adjacent_find(v.begin(), v.end(),
                    std::not_equal_to<std::string>()) == v.end());
            calcOne(measurementSets()[rank], distributeByTile);
        }
    } else {
        for (size_t iMs=0; iMs<measurementSets().size(); ++iMs) {
            calcOne(measurementSets()[iMs]);
        }
    }
}
