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
#include <askap/AskapError.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <casacore/scimath/Fitting/LinearFitSVD.h>


// logging stuff
#include <askap/askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

#include <casacore/casa/OS/Timer.h>


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
      const LOFAR::ParameterSet& parset) : MEParallelApp(comms,parset)
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
  if (itsWidth > 0) ASKAPCHECK(itsOffset < itsWidth,"The offset needs to be less than the width");
  if (itsDoUVlin) ASKAPLOG_INFO_STR(logger, "Doing uvlin operation with order = "
    << itsOrder<<", harmonic = "<<itsHarmonic<< ", width = "<< itsWidth <<" channels, offset = "<< itsOffset);
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
        const casacore::Vector<casa::Float>& spec, const casa::Vector<casa::Bool>& mask) {
    int nChan= spec.size();
    int nParams = itsOrder+1+itsHarmonic*2;
    casa::LSQaips fitter(nParams);
    casa::Vector<casa::Double> xx(nParams);
    casa::VectorSTLIterator<casa::Double> it(xx);

     // we are doing the fitting in channel bins
     if (itsWidth == 0) itsWidth = nChan;
     for (int binStart=-itsOffset; binStart<nChan; binStart+=itsWidth) {
         int start = max(0,binStart);
         int end = min(binStart+itsWidth,nChan);
         int binWidth = end - start;

         for (int i=start; i < end; i++) {
             if (mask(i)) {
                 xx(0) = 1;
                 for (int j=1; j<itsOrder+1; j++) {
                     xx(j) = xx(j-1) * (i-start);
                 }
                 for (int j=0; j<itsHarmonic; j++) {
                     xx(itsOrder+1+2*j)   = sin((j+1)*casa::C::pi*(i-start)/binWidth);
                     xx(itsOrder+1+2*j+1) = cos((j+1)*casa::C::pi*(i-start)/binWidth);
                     // we could use itsWidth instead of binWidth to keep 'frequency' the same
                 }
                 fitter.makeNorm(it,1.0,casa::Double(spec(i)));
             }
         }

         casa::uInt nr1;
         casa::Vector<casa::Double> solution(nParams);
         fitter.invert(nr1);
         fitter.solve(solution);

         for (int i=start; i < end; i++) {
             model(i) = solution(itsOrder);
             for (int j=itsOrder-1; j>=0; j--) {
                 model(i) = (i-start) * model(i) + solution(j);
             }
             for (int j=0; j<itsHarmonic; j++) {
                  model(i) += solution(itsOrder+1+2*j)   * sin((j+1)*casa::C::pi*(i-start)/binWidth);
                  model(i) += solution(itsOrder+1+2*j+1) * cos((j+1)*casa::C::pi*(i-start)/binWidth);
             }
         }
         fitter.reset();
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
