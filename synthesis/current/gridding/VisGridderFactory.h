/// @file
///
/// VisGridderFactory: Factory class for visibility
//  gridders.
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
///
#ifndef ASKAP_SYNTHESIS_VISGRIDDERFACTORY_H_
#define ASKAP_SYNTHESIS_VISGRIDDERFACTORY_H_

// System includes
#include <map>

// ASKAPsoft includes
#include <Common/ParameterSet.h>
#include <scimath/Mathematics/Interpolate2D.h>
#include <boost/shared_ptr.hpp>

// Loacl package includes
#include <gridding/IVisGridder.h>
#include <gridding/IBasicIllumination.h>

namespace askap
{
  namespace synthesis
  {
    /// @brief Factory class for visibility gridders
    /// @ingroup gridding
    class VisGridderFactory
    {
    public:
      /// @brief Signature of a function creating a Gridder object.
      /// All functions creating a IVisGridder object must have
      /// this signature. Preferably such a function is a static
      /// function in that gridder class.
      typedef IVisGridder::ShPtr GridderCreator
      (const LOFAR::ParameterSet&);

      /// @brief Register a function creating a gridder object.
      /// @param name The name of the gridder.
      /// @param creatorFunc pointer to creator function.
      static void registerGridder (const std::string& name,
                                   GridderCreator* creatorFunc);

      /// @brief Try to create a non-standard gridder.
      /// Its name is looked up in the creator function registry.
      /// If the gridder name is unknown, a shared library with that name
      /// (in lowercase) is loaded and it executes its register<name>
      /// function which must register its creator function in the registry
      /// using function registerGridder.
      /// @param name The name of the gridder.
      /// @param parset ParameterSet containing description of
      /// gridder to be constructed
      static IVisGridder::ShPtr createGridder (const std::string& name,
                                               const LOFAR::ParameterSet& parset);


      /// @brief Factory class for all gridders.
      /// @todo Python version of factory 
      VisGridderFactory();

      /// @brief Make a shared pointer for a visibility gridder
      /// @param parset ParameterSet containing description of
      /// gridder to be constructed.
      /// If needed, the gridder code is loaded dynamically.
      static IVisGridder::ShPtr make(const LOFAR::ParameterSet& parset);

    protected:
      /// @brief helper template method to add pre-defined gridders
      template<typename GridderType> static inline void addPreDefinedGridder()  {
          registerGridder(GridderType::gridderName(), GridderType::createGridder);
      }
      
      
    private:
      static casa::Interpolate2D::Method interpMethod(casa::String str);

      static std::map<std::string, GridderCreator*> theirRegistry;
    };

  }
}
#endif 
