/// @file
/// 
/// @brief Measurement equation adapter filling calibration parameters on demand
/// @details This is an adapter intended to be used when calibration effects
/// are simulated. It implements IMeasurementEquation interface, but only
/// predict method is expected to be used. An exception is thrown if 
/// one requests normal equations to be computed with this class. The predict
/// method checks whether new antenna/beam combinations are present in the 
/// current visibility chunk. It then creates new parameters (this class follows
/// the name convension enforced by accessors::CalParamNameHelper, but this 
/// behavior can be overridden in derived classes, if necessary) and updates
/// the existing ones using calibration solution source supplied at construction.
/// After parameters are created, the predict method of the wrapped measurement
/// equation is called to predict visibilities. 
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

// own includes
#include <askap_synthesis.h>
#include <measurementequation/CalibParamsMEAdapter.h>
#include <askap/AskapError.h>
#include <fitting/Equation.h>
#include <calibaccess/JonesJTerm.h>
#include <calibaccess/JonesDTerm.h>
#include <askap/AskapLogging.h>
#include <calibaccess/CalParamNameHelper.h>
ASKAP_LOGGER(logger, ".measurementequation");


namespace askap {

namespace synthesis {

/// @brief standard constructor
/// @details initialises the class with the given solution source and iterator
/// (if necessary). It also initialises the shared pointer to the slave measurement 
/// equation which does the actual prediction of visibilities. Only accessor-based
/// functionality of this slave measurement equation is used.
/// @param[in] ime shared pointer to the slave measurement equation
/// @param[in] css shared pointer to solution source
/// @param[in] idi data iterator (if iterator-based interface is required)
CalibParamsMEAdapter::CalibParamsMEAdapter(const boost::shared_ptr<IMeasurementEquation> &ime,
                const boost::shared_ptr<accessors::ICalSolutionConstSource> &css, 
                const accessors::IDataSharedIter& idi) :                 
                MultiChunkEquation(idi), CalibrationSolutionHandler(css), itsSlaveME(ime) 
{
  ASKAPCHECK(itsSlaveME, "An attempt to initialise CalibParamsMEAdapter with a void shared pointer to the slave ME");
  // deliberately reuse shared pointer between the slaved measurement equation and this adapter to
  // ensure reference semantics for parameters if we could. It is not clear what to do if we can't extract the
  // parameters. Throwing an exception for now because it is not our intended use case.
  boost::shared_ptr<scimath::Equation> eqn = boost::dynamic_pointer_cast<scimath::Equation>(itsSlaveME);
  ASKAPCHECK(eqn, "Attempt to initialise CalibParamsMEAdapter with an incompatible type of slave measurement equation");
  reference(eqn->rwParameters());
  ASKAPASSERT(rwParameters());
}

/// @brief copy constructor
/// @details We need it because some data members have non-trivial types (shared pointers)
/// @param[in] other other object to copy from
CalibParamsMEAdapter::CalibParamsMEAdapter(const CalibParamsMEAdapter &other) : Equation(other), MultiChunkEquation(other), 
       CalibrationSolutionHandler(other), 
       itsCalSolutionCM(other.itsCalSolutionCM), itsAntBeamPairs(other.itsAntBeamPairs) 
{
  ASKAPCHECK(other.itsSlaveME, "Slave ME is uninitialised in a call to the copy constructor");
  // deliberately reuse shared pointer between the slaved measurement equation and this adapter to
  // ensure reference semantics for parameters if we could. It is not clear what to do if we can't extract the
  // parameters. Throwing an exception for now because it is not our intended use case.
  boost::shared_ptr<scimath::Equation> eqn = boost::dynamic_pointer_cast<scimath::Equation>(other.itsSlaveME);
  ASKAPCHECK(eqn, "Attempt to initialise CalibParamsMEAdapter with an incompatible type of slave measurement equation");
  itsSlaveME = boost::dynamic_pointer_cast<IMeasurementEquation>(eqn->clone());
  ASKAPDEBUGASSERT(itsSlaveME);
  eqn = boost::dynamic_pointer_cast<scimath::Equation>(itsSlaveME);
  ASKAPDEBUGASSERT(eqn);
  reference(eqn->rwParameters());
  ASKAPASSERT(rwParameters());
  // 
  ASKAPTHROW(AskapError, "This code has not been tested and technically we don't need it");
}

/// @brief assignment operator
/// @details We need it because some data members have non-trivial types (shared pointers)
/// @param[in] other other object to copy from
const CalibParamsMEAdapter& CalibParamsMEAdapter::operator=(const CalibParamsMEAdapter &other)
{
  if (this != &other) {
      MultiChunkEquation::operator=(other);
      CalibrationSolutionHandler::operator=(other);
      itsCalSolutionCM = other.itsCalSolutionCM;
      itsAntBeamPairs = other.itsAntBeamPairs;
      ASKAPCHECK(other.itsSlaveME, "Slave ME is uninitialised in a call to the copy constructor");
      // deliberately reuse shared pointer between the slaved measurement equation and this adapter to
      // ensure reference semantics for parameters if we could. It is not clear what to do if we can't extract the
      // parameters. Throwing an exception for now because it is not our intended use case.
      boost::shared_ptr<scimath::Equation> eqn = boost::dynamic_pointer_cast<scimath::Equation>(other.itsSlaveME);
      ASKAPCHECK(eqn, "Attempt to initialise CalibParamsMEAdapter with an incompatible type of slave measurement equation");
      itsSlaveME = boost::dynamic_pointer_cast<IMeasurementEquation>(eqn->clone());
      ASKAPDEBUGASSERT(itsSlaveME);
      eqn = boost::dynamic_pointer_cast<scimath::Equation>(itsSlaveME);
      ASKAPDEBUGASSERT(eqn);
      reference(eqn->rwParameters());
      ASKAPASSERT(rwParameters());
      // 
      ASKAPTHROW(AskapError, "This code has not been tested and technically we don't need it");      
  }
  return *this;
}

   
/// @brief Predict model visibilities for one accessor (chunk).
/// @details This prediction is done for single chunk of data only.
/// This method overrides pure virtual method of the base class. 
/// @param[in] chunk a read-write accessor to work with
void CalibParamsMEAdapter::predict(accessors::IDataAccessor &chunk) const
{
  updateAccessor(chunk.time());
  if (itsCalSolutionCM != changeMonitor()) {
      itsCalSolutionCM = changeMonitor();
      // there was an update, we need to rebuild parameters before calling predict
      // method of the slave ME.
      ASKAPLOG_INFO_STR(logger, "CalibParamsMEAdapter - new calibration solution (time="<<chunk.time()<<
                        " seconds since 0 MJD), regenerating calibration parameters");
      itsAntBeamPairs.clear();      
  }
  // iterate over all rows as even if there is no change in calibration solution the new
  // visibility chunk may have antennas/beams not previously seen in the data.
  const casa::Vector<casa::uInt>& ant1IDs = chunk.antenna1();
  const casa::Vector<casa::uInt>& ant2IDs = chunk.antenna2();
  const casa::Vector<casa::uInt>& beam1IDs = chunk.feed1();
  const casa::Vector<casa::uInt>& beam2IDs = chunk.feed2();
  const casa::uInt nRow = chunk.nRow();
  ASKAPDEBUGASSERT(nRow == ant1IDs.nelements());
  ASKAPDEBUGASSERT(nRow == ant2IDs.nelements());
  ASKAPDEBUGASSERT(nRow == beam1IDs.nelements());
  ASKAPDEBUGASSERT(nRow == beam2IDs.nelements());
  for (casa::uInt row = 0; row < nRow; ++row) {
       processAntBeamPair(ant1IDs[row], beam1IDs[row]);
       processAntBeamPair(ant2IDs[row], beam2IDs[row]);
  }
  // parameters are now ready, predict visibilities using the slave measurement equation
  itsSlaveME->predict(chunk);
}
   
/// @brief Calculate the normal equation for one accessor (chunk).
/// @details This method is not supposed to be used. It overrides
/// pure virtual method of the base class and throws an exception if called 
void CalibParamsMEAdapter::calcEquations(const accessors::IConstDataAccessor &,
                          askap::scimath::INormalEquations&) const
{
  ASKAPTHROW(AskapError, "CalibParamsMEAdapter::calcEquations is not supposed to be used");
}  

/// @brief manage update for a one antenna/beam
/// @details This method manages cache and calls virtual 
/// updateSingleAntBeam method if an update is necessary
/// @param[in] ant antenna index   
/// @param[in] beam beam index
void CalibParamsMEAdapter::processAntBeamPair(const casa::uInt ant, const casa::uInt beam) const
{
   const accessors::JonesIndex index(ant, beam);
   if (itsAntBeamPairs.find(index) == itsAntBeamPairs.end()) {
       itsAntBeamPairs.insert(index);
       updateSingleAntBeam(index);    
   }   
}
                        

/// @brief process parameters for a given antenna/beam   
/// @details This method encapsulates update of the parameters 
/// corresponding to the given antenna/beam pair according to 
/// the current calibration solution accessor. Override this
/// method for non-standard parameter names (use updateParameter
/// for the actual update).
/// @param[in] index ant/beam index
void CalibParamsMEAdapter::updateSingleAntBeam(const accessors::JonesIndex &index) const
{
  // it is not clear how to deal with the bandpass, ignore it for now
  // we could've also used the validity flag for calibration information to adjust the
  // flags for visibilities, but don't bother for now for simplicity (this code
  // is intended to be used with the simulator only anyway).
  
  const accessors::JonesJTerm gain = calSolution().gain(index);  
  updateParameter(accessors::CalParamNameHelper::paramName(index, casa::Stokes::XX), 
                  gain.g1IsValid() ? gain.g1() : casa::Complex(1.,0.));
  updateParameter(accessors::CalParamNameHelper::paramName(index, casa::Stokes::YY), 
                  gain.g2IsValid() ? gain.g2() : casa::Complex(1.,0.));
                  
  const accessors::JonesDTerm leakage = calSolution().leakage(index);    
  updateParameter(accessors::CalParamNameHelper::paramName(index, casa::Stokes::XY), 
                  leakage.d12IsValid() ? leakage.d12() : casa::Complex(0.,0.));
  updateParameter(accessors::CalParamNameHelper::paramName(index, casa::Stokes::YX),  
                  leakage.d21IsValid() ? leakage.d21() : casa::Complex(0.,0.));
}
   
/// @brief helper method to update a given parameter if necessary
/// @details This method checks whether the parameter is new and
/// adds or updates as required in the parameter class held by the
/// slave measurement equation.
/// @param[in] name parameter name
/// @param[in] val new value
void CalibParamsMEAdapter::updateParameter(const std::string &name, const casa::Complex &val) const
{
  ASKAPDEBUGASSERT(rwParameters());
  if (rwParameters()->has(name)) {
      rwParameters()->update(name,val);
  } else {
      rwParameters()->add(name,val);
  }
}


} // namespace synthesis

} // namespace askap


