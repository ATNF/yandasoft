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

#ifndef CALIB_PARAMS_ME_ADAPTER_H
#define CALIB_PARAMS_ME_ADAPTER_H

// own includes
#include <measurementequation/IMeasurementEquation.h>
#include <measurementequation/MultiChunkEquation.h>
#include <measurementequation/CalibrationSolutionHandler.h>
#include <dataaccess/SharedIter.h>
#include <calibaccess/ICalSolutionConstSource.h>
#include <calibaccess/JonesIndex.h>
#include <utils/ChangeMonitor.h>

// boost includes
#include <boost/shared_ptr.hpp>

// std includes
#include <string>
#include <set>


namespace askap {

namespace synthesis {

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
/// @note The class is derived from MultiChunkEquation rather than IMeasurementEquation
/// to provide immediate an appropriate interface for our immediate needs. It is safe
/// to initialise the adapter with an empty iterator if
/// only accessor-based version of the predict method is used.
/// @ingroup measurementequation
class CalibParamsMEAdapter : public MultiChunkEquation,
                             protected CalibrationSolutionHandler {
public:
   /// @brief standard constructor
   /// @details initialises the class with the given solution source and iterator
   /// (if necessary). It also initialises the shared pointer to the slave measurement 
   /// equation which does the actual prediction of visibilities. Only accessor-based
   /// functionality of this slave measurement equation is used.
   /// @param[in] ime shared pointer to the slave measurement equation
   /// @param[in] css shared pointer to solution source
   /// @param[in] idi data iterator (if iterator-based interface is required)
   explicit CalibParamsMEAdapter(const boost::shared_ptr<IMeasurementEquation> &ime,
                const boost::shared_ptr<accessors::ICalSolutionConstSource> &css = boost::shared_ptr<accessors::ICalSolutionConstSource>(), 
                const accessors::IDataSharedIter& idi = accessors::IDataSharedIter());
   
   /// @brief copy constructor
   /// @details We need it because some data members have non-trivial types (shared pointers)
   /// @param[in] other other object to copy from
   CalibParamsMEAdapter(const CalibParamsMEAdapter &other);
   
   /// @brief assignment operator
   /// @details We need it because some data members have non-trivial types (shared pointers)
   /// @param[in] other other object to copy from
   const CalibParamsMEAdapter& operator=(const CalibParamsMEAdapter &other);
   
   /// @brief clone method
   /// @details
   /// @return shared pointer to a clone
   virtual boost::shared_ptr<scimath::Equation> clone() const 
      { return boost::shared_ptr<scimath::Equation>(new CalibParamsMEAdapter(*this));}
   
   /// @brief Predict model visibilities for one accessor (chunk).
   /// @details This prediction is done for single chunk of data only.
   /// This method overrides pure virtual method of the base class. 
   /// @param[in] chunk a read-write accessor to work with
   virtual void predict(accessors::IDataAccessor &chunk) const;
   
   /// @brief Calculate the normal equation for one accessor (chunk).
   /// @details This method is not supposed to be used. It overrides
   /// pure virtual method of the base class and throws an exception if called 
   /// @param[in] chunk a read-write accessor to work with
   /// @param[in] ne Normal equations
   virtual void calcEquations(const accessors::IConstDataAccessor &chunk,
                          askap::scimath::INormalEquations& ne) const;
   
   // to make this method public
   using CalibrationSolutionHandler::setCalSolutionSource;
   // make iterator-based methods visible here
   using MultiChunkEquation::predict;
   using MultiChunkEquation::calcEquations;

protected:   
   /// @brief process parameters for a given antenna/beam   
   /// @details This method encapsulates update of the parameters 
   /// corresponding to the given antenna/beam pair according to 
   /// the current calibration solution accessor. Override this
   /// method for non-standard parameter names (use updateParameter
   /// for the actual update).
   /// @param[in] index ant/beam index   
   virtual void updateSingleAntBeam(const accessors::JonesIndex &index) const;
   
   /// @brief manage update for a one antenna/beam
   /// @details This method manages cache and calls virtual 
   /// updateSingleAntBeam method if an update is necessary
   /// @param[in] ant antenna index   
   /// @param[in] beam beam index
   void processAntBeamPair(const casa::uInt ant, const casa::uInt beam) const;
   
   
   /// @brief helper method to update a given parameter if necessary
   /// @details This method checks whether the parameter is new and
   /// adds or updates as required in the parameter class held by the
   /// slave measurement equation.
   /// @param[in] name parameter name
   /// @param[in] val new value
   void updateParameter(const std::string &name, const casa::Complex &val) const;
private:
    
    /// @brief slave measurement equation
    boost::shared_ptr<IMeasurementEquation> itsSlaveME;  
    
    /// @brief change monitor corresponding to the current parameters held by this class
    /// @details It is used for caching, parameters are rebuilt if change of calibration
    /// solution did take place.
    mutable scimath::ChangeMonitor itsCalSolutionCM;
    
    /// @brief antenna/beam pairs for which calibration parameters have been generated
    /// @details We skip generation of parameters if it has been done already. This set
    /// keeps track of the antenna/beam pairs for which parameters have been added into
    /// scimath::Params class held by this adapter (and accessible via rwParameters).
    mutable std::set<accessors::JonesIndex> itsAntBeamPairs;
    
}; // class CalibParamsMEAdapter                            

} // namespace synthesis

} // namespace askap

#endif // #ifndef CALIB_PARAMS_ME_ADAPTER_H

