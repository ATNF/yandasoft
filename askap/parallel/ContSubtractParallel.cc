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
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/dataaccess/MemBufferDataAccessor.h>
#include <askap/dataaccess/DataIteratorStub.h>
#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/scimath/fitting/Equation.h>
#include <askap/measurementequation/ComponentEquation.h>
#include <askap/measurementequation/IMeasurementEquation.h>
#include <askap/measurementequation/ImagingEquationAdapter.h>
#include <askap/askap/AskapError.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <casacore/scimath/Fitting/LinearFitSVD.h>

#include <casacore/casa/Arrays/ArrayMath.h>

// logging stuff
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

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
ContSubtractParallel::ContSubtractParallel(askap::askapparallel::AskapParallel& comms,
      const LOFAR::ParameterSet& parset) : MEParallelApp(comms,parset,true)
{
  // the stub allows to reuse MEParallelApp code although we're not solving
  // for the normal equations here
  itsNe.reset(new scimath::NormalEquationsStub);

  itsModelReadByMaster = parset.getBool("modelReadByMaster", true);
  itsDoUVlin = parset.getBool("doUVlin", false);
  itsOrder = parset.getInt("order", 1);
  itsHarmonic = parset.getInt("harmonic", 1);
  itsWidth = parset.getInt("width",0); // 0 = whole spectrum
  itsOffset = min(max(0,parset.getInt("offset",0)),itsWidth);
  itsThreshold = max(0.0f,parset.getFloat("threshold",2.5));
  if (itsWidth > 0) ASKAPCHECK(itsOffset < itsWidth,"The offset needs to be less than the width");
  if (itsDoUVlin) ASKAPLOG_INFO_STR(logger, "Doing uvlin operation with order = "
    << itsOrder<<", harmonic = "<<itsHarmonic<< ", width = "<< itsWidth
    <<" channels, offset = "<< itsOffset << " and threshold = "<< itsThreshold);
}

/// @brief Initialise continuum subtractor
/// @details The parameters are taken from the parset file supplied in the constructor.
/// This method does initialisation which may involve communications in the parallel case
/// (i.e. distribution of the models between workers). Technically, we could've done this in
/// the constructor.
void ContSubtractParallel::init()
{
  if (itsComms.isMaster()) {
      if (itsModelReadByMaster) {
          readModels();
          broadcastModel();
      }
  }
  if (itsComms.isWorker()) {
      if (itsModelReadByMaster) {
          receiveModel();
      } else {
          readModels();
      }
  }
}

/// @brief initialise measurement equation
/// @details This method initialises measurement equation
void ContSubtractParallel::initMeasurementEquation()
{
   ASKAPLOG_INFO_STR(logger, "Creating measurement equation" );

   // it doesn't matter which iterator to pass to the measurement equations
   // as we're only using accessor-based interface
   IDataSharedIter stubIter(new DataIteratorStub(1));

   ASKAPCHECK(itsModel, "Model is not defined");
   ASKAPCHECK(gridder(), "Gridder is not defined");

   // a part of the equation defined via image
   askap::scimath::Equation::ShPtr imgEquation;

   if (SynthesisParamsHelper::hasImage(itsModel)) {
       ASKAPLOG_INFO_STR(logger, "Sky model contains at least one image, building an image-specific equation");
       // it should ignore parameters which are not applicable (e.g. components)
       imgEquation.reset(new ImageFFTEquation(*itsModel, stubIter, gridder()));
   }

   // a part of the equation defined via components
   boost::shared_ptr<ComponentEquation> compEquation;

   if (SynthesisParamsHelper::hasComponent(itsModel)) {
       // model is a number of components
       ASKAPLOG_INFO_STR(logger, "Sky model contains at least one component, building a component-specific equation");
       // it doesn't matter which iterator is passed below. It is not used
       // it should ignore parameters which are not applicable (e.g. images)
       compEquation.reset(new ComponentEquation(*itsModel, stubIter));
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
        // form a replacement equation first
        const boost::shared_ptr<ImagingEquationAdapter> new_equation(new ImagingEquationAdapter);
        // the actual equation (from itsEquation) will be locked inside ImagingEquationAdapter
        // in a shared pointer. We can change itsEquation after the following line
        new_equation->assign(itsEquation);
        // replacing the original equation with an accessor-based adapter
        itsEquation = new_equation;
   }

}

void ContSubtractParallel::modelSpectrum(casa::Vector<casa::Float> & model,
        const casacore::Vector<casa::Float>& spec, const casa::Vector<casa::Bool>& mask)
{
    int nChan= spec.size();
    int nParams = itsOrder+1+itsHarmonic*2;
    casa::LSQaips fitter(nParams);
    casa::Vector<casa::Double> xx(nParams);
    casa::Vector<casa::Double> solution(nParams);
    casa::VectorSTLIterator<casa::Double> it(xx);
    casa::Vector<casa::Bool> tmask(nChan);

    // If we are doing outlier rejection against the model iterate a few times
    const int niter = (itsThreshold > 0 ? 3 : 1);
    // initial mask is the data flags
    tmask = mask;
    // initial model is zero
    model = 0;

    // we are doing the fitting in channel bins
    if (itsWidth == 0) itsWidth = nChan;
    casa::Vector<casa::Float> y(itsWidth);
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
                    casa::Float q25 = casa::fractile(y(casa::Slice(0,n)), tmp, 0.25f, casa::False, casa::True);
                    casa::Float q50 = casa::fractile(y(casa::Slice(0,n)), tmp, 0.50f, casa::False, casa::True);
                    casa::Float q75 = casa::fractile(y(casa::Slice(0,n)), tmp, 0.75f, casa::False, casa::True);
                    casa::Float sigma = (q75-q25)/1.35; // robust sigma estimate
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
                    for (int i=start; i<end; i++) if (mask(i) && !tmask(i)) count2++;
                    //casa::cerr<<"iter="<<iter<<", n="<<n<<", start="<<start<<", median="<<q50<<", rsigma="<<sigma<<", outliers flagged: "<<count<<", extended to "<<count2<<casa::endl;
                    //casa::cerr<<"model("<<start<<")="<<model(start)<<", model("<<end-1<<")="<<model(end-1)<<
                    //    ", y(0)="<<y(0)<<", y("<<n-1<<")="<<y(n-1)<<casa::endl;
                    if (iter > 0 && count2 == 0) break; // no further change expected
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
                float x = (i-start) / float(binWidth);
                if (tmask(i)) {
                    xx(0) = 1;
                    for (int j=1; j<order+1; j++) {
                        xx(j) = xx(j-1) * x;
                    }
                    for (int j=0; j<harm; j++) {
                        xx(order+1+2*j)   = sin((j+1)*casa::C::pi*x);
                        xx(order+1+2*j+1) = cos((j+1)*casa::C::pi*x);
                    }
                    fitter.makeNorm(it,1.0,casa::Double(spec(i)));
                }
            }

            casa::uInt nr1;
            casa::Bool ok = fitter.invert(nr1);
            if (ok) {
                fitter.solve(solution.data());
                //casa::Float chisq = fitter.getChi();
                //casa::Float sd1 = fitter.getSD();
                //casa::cerr << "Fit="<< ok <<", rank="<<nr1<<", chisq="<<chisq<<", sd="<<sd1<<", sol"<<solution<<casa::endl;

                // evaluate the solution to generate the model
                for (int i=start; i < end; i++) {
                    float x = (i-start) / float(binWidth);
                    model(i) = solution(order);
                    for (int j=order-1; j>=0; j--) {
                        model(i) = x * model(i) + solution(j);
                    }
                    for (int j=0; j<harm; j++) {
                        model(i) += solution(order+1+2*j)   * sin((j+1)*casa::C::pi*x);
                        model(i) += solution(order+1+2*j+1) * cos((j+1)*casa::C::pi*x);
                    }
                }
            } else {
                break; // no point iterating if the fit failed, keep last model or zero
            }
        }
    }
}

void ContSubtractParallel::subtractContFit(casacore::Cube<casacore::Complex>& vis,
        const casacore::Cube<casacore::Bool>& flag) {
    int nPol = vis.shape()(0);
    int nChan = vis.shape()(1);
    int nRow = vis.shape()(2);
    casa::Vector<casa::Float> visreal(nChan), visimag(nChan), modelreal(nChan), modelimag(nChan);
    casa::Vector<casa::Bool> mask(nChan);
    for (int row=0; row<nRow; row++) {
        for (int pol=0; pol<nPol; pol++) {
            for (int chan=0; chan<nChan; chan++) {
                casa::Complex v = vis(pol,chan,row);
                visreal(chan) = casa::real(v);
                visimag(chan) = casa::imag(v);
                mask(chan) = !flag(pol,chan,row);
            }
            modelSpectrum(modelreal,visreal,mask);
            modelSpectrum(modelimag,visimag,mask);
            for (int chan=0; chan<nChan; chan++) {
                vis(pol,chan,row) -= casa::Complex(modelreal(chan),modelimag(chan));
            }
        }
    }
}

/// @brief perform the subtraction for the given dataset
/// @details This method iterates over the given dataset, predicts visibilities according to the
/// model and subtracts these model visibilities from the original visibilities in the dataset.
/// This is the core operation of the doSubtraction method, which manages the parallel aspect of it.
/// All actual calculations are done inside this helper method.
/// @param[in] ms measurement set name
void ContSubtractParallel::calcOne(const std::string &ms)
{
   casacore::Timer timer;
   timer.mark();
   ASKAPLOG_INFO_STR(logger, "Performing continuum model subtraction for " << ms );

   if (!itsEquation) {
       initMeasurementEquation();
   } else {
      ASKAPLOG_INFO_STR(logger, "Reusing measurement equation" );
   }


   boost::shared_ptr<IMeasurementEquation> accessorBasedEquation =
        boost::dynamic_pointer_cast<IMeasurementEquation>(itsEquation);
   ASKAPDEBUGASSERT(accessorBasedEquation);

   TableDataSource ds(ms, TableDataSource::WRITE_PERMITTED, dataColumn());
   ds.configureUVWMachineCache(uvwMachineCacheSize(),uvwMachineCacheTolerance());
   IDataSelectorPtr sel=ds.createSelector();
   sel << parset();
   IDataConverterPtr conv=ds.createConverter();
   conv->setFrequencyFrame(getFreqRefFrame(), "Hz");
   conv->setDirectionFrame(casacore::MDirection::Ref(casacore::MDirection::J2000));
   IDataSharedIter it=ds.createIterator(sel, conv);
   for (; it.hasMore(); it.next()) {
        // iteration over the dataset
        MemBufferDataAccessor acc(*it);
        acc.rwVisibility().set(0.);
        accessorBasedEquation->predict(acc);
        const casacore::Cube<casacore::Complex>& model = acc.visibility();
        casacore::Cube<casacore::Complex>& vis = it->rwVisibility();
        ASKAPDEBUGASSERT(model.nrow() == vis.nrow());
        ASKAPDEBUGASSERT(model.ncolumn() == vis.ncolumn());
        ASKAPDEBUGASSERT(model.nplane() == vis.nplane());
        vis -= model;
        if (itsDoUVlin) {
            subtractContFit(vis,acc.flag());
        }
   }

   ASKAPLOG_INFO_STR(logger, "Finished continuum subtraction for "<< ms << " in "<< timer.real()
                   << " seconds ");
}


/// @brief perform the subtraction
/// @details This method iterates over one or more datasets, predicts visibilities according to
/// the model and subtracts these model visibilities from the original visibilities in the
/// dataset. The intention is to call this method in a worker.
void ContSubtractParallel::doSubtraction()
{
  if (itsComms.isWorker()) {
      if (itsComms.isParallel()) {
          calcOne(measurementSets()[itsComms.rank()-1]);
      } else {
          for (size_t iMs=0; iMs<measurementSets().size(); ++iMs) {
               calcOne(measurementSets()[iMs]);
          }
      }
  }
}
