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

#ifndef CONT_SUBTRACT_PARALLEL_H
#define CONT_SUBTRACT_PARALLEL_H

// ASKAPsoft includes
#include <Common/ParameterSet.h>
#include <askap/parallel/MEParallelApp.h>
#include <askap/calibaccess/ICalSolutionConstSource.h>
#include <askap/measurementequation/ICalibrationApplicator.h>

namespace askap {

namespace synthesis {

/// @brief parallel helper for continuum subtraction
/// @details This class does the core operation to subtract continuum by doing visibility
/// prediction from the given model in parallel.
/// @ingroup parallel
class ContSubtractParallel : public MEParallelApp
{
public:
   /// @brief Constructor from ParameterSet
   /// @details The parset is used to construct the internal state. We could
   /// also support construction from a python dictionary (for example).
   /// The command line inputs are needed solely for MPI - currently no
   /// application specific information is passed on the command line.
   /// @param comms communication object
   /// @param parset ParameterSet for inputs
   ContSubtractParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset);

   /// @brief Initialise continuum subtractor
   /// @details The parameters are taken from the parset file supplied in the constructor.
   /// This method does initialisation which may involve communications in the parallel case
   /// (i.e. distribution of the models between workers). Technically, we could've done this in
   /// the constructor.
   void init();

   /// @brief perform the subtraction
   /// @details This method iterates over one or more datasets, predicts visibilities according to
   /// the model and subtracts these model visibilities from the original visibilities in the
   /// dataset. The intention is to call this method in a worker.
   void doSubtraction();

 protected:
   /// @brief read the models from the parset file
   inline void readModels() const { SynParallel::readModels(itsModel); }

   /// @brief initialise measurement equation
   /// @param[in] it the data iterator
   /// @details This method initialises measurement equation
   void initMeasurementEquation(accessors::IDataSharedIter& it);

   /// @brief perform the subtraction for the given dataset
   /// @details This method iterates over the given dataset, predicts visibilities according to the
   /// model and subtracts these model visibilities from the original visibilities in the dataset.
   /// This is the core operation of the doSubtraction method, which manages the parallel aspect of it.
   /// All actual calculations are done inside this helper method.
   /// @param[in] ms measurement set name
   /// @param[in] distributeByTile set to true to do distribution by tiles if possible
   void calcOne(const std::string &ms, bool distributeByTile = false);

   // stubs for pure virtual methods which we don't use
   /// @brief calculate normal equations
   inline void calcNE() {}

   /// @brief solve normal equations
   inline void solveNE() {}

   /// @brief write results
   inline void writeModel(const std::string&) {}

 private:
   /// @brief generate a model spectrum by fitting the input spectrum where mask is true
   /// @details Use polynomial fitting to derive a model spectrum from the input spectrum
   /// @param[in] spec input spectrum
   /// @param[in] mask input mask
   /// @param[out] model output model spectrum
   void modelSpectrum(casa::Vector<casa::Float> & model,
            const casacore::Vector<casa::Float>& spec, const casa::Vector<casa::Bool>& mask);

   /// @brief subtract continuum fit from the input visibilities
   /// @details process input visibilities spectrum by spectrum, fitting and subtracting continuum
   /// @param[in/out] vis the visibilities, will be modified on output
   /// @param[in] flag the associated flags
   void subtractContFit(casacore::Cube<casacore::Complex>& vis,
            const casacore::Cube<casacore::Bool>& flag,
            const casacore::Matrix<casacore::Complex>& phasor);

   /// @brief Compute the phase rotation phasor for rotation of visibilities to new direction
   void computePhasor(const accessors::IDataSharedIter& it,//const accessors::IConstDataAccessor& acc,
            casacore::Matrix<casacore::Complex>& phasor);

   /// @brief Load the subtables into memory
   /// @details We load the subtables into memory before accessing the data in r/w mode
   /// to avoid parallel IO errors
   void loadSubtables(const accessors::IDataSharedIter& it);

   /// @brief model is read by the master and distributed?
   /// @details Depending on the model file name (containing %w or not), the model
   /// can either be read in the master and distributed across the workers or read
   /// by workers directly. If true, the model is read by the master
   bool itsModelReadByMaster;
   /// @brief do continuum subtraction based on a component or image model
   bool itsDoSubtraction;
   /// @brief Instead of subtracting, replace the data by the model
   bool itsDoReplaceByModel;
   /// @brief do 'uvlin' like fit and subtract of residual continuum emission
   bool itsDoUVlin;
   /// @brief order of 'uvlin' like fit and subtract of residual continuum emission
   int itsOrder;
   /// @brief number of 'harmonic', i.e., sin and cos terms used in the fit
   int itsHarmonic;
   /// @brief number of channels to use for the fit, 0 means all
   int itsWidth;
   /// @brief offset of channel bins, first bin will have (width-offset) channels
   /// @details the askap beamforming intervals can be a multiple of the basic 1 MHz
   /// interval, but they don't necessarily start at channel 0 in the recorded spectrum
   int itsOffset;
   /// @brief threshold for rejection of channel from the fit, 0 mean no rejection
   /// @details reject channels if |value - median|>threshold*sigma_IQR (robust sigma estimate)
   float itsThreshold;
   /// @brief rotate visibilities to new direction before uvlin fit/subtraction
   casacore::MDirection itsUVlinDirection;
   /// @brief are we rotating visibilities?
   bool itsRotate;

   /// @brief solution source to get calibration data from
   /// @details This object is initialised by workers. It knows how to
   /// retrieve calibration solutions (from a parset file, casa table or a database).
   /// Uninitialised shared pointer means no calibration
   boost::shared_ptr<accessors::ICalSolutionConstSource> itsSolutionSource;

   /// @brief calibration applicator (if needed)
   /// @details If the applicator method of calibration is selected
   /// this stores the CalibrationApplicatorME object.
   /// If uninitialised, we are using the CalibrationMEBase classes
   boost::shared_ptr<ICalibrationApplicator> itsCalApplicator;

   /// @brief number of beams/calibration directions
   int itsNDir;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef CONT_SUBTRACT_PARALLEL_H
