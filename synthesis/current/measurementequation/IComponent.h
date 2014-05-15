/// @file
///
/// @brief Abstract component
/// @details
/// IComponent is a base class for components working with ComponentEquation
/// examples of components include, e.g. Gaussian or point sources.
/// 
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

#ifndef I_COMPONENT_H
#define I_COMPONENT_H

// own includes
#include <askap/AskapError.h>

// std includes
#include <vector>

// casa includes
#include <scimath/Mathematics/AutoDiff.h>
#include <scimath/Mathematics/AutoDiffMath.h>
#include <casa/Arrays/Vector.h>
#include <scimath/Mathematics/RigidVector.h>
#include <measures/Measures/Stokes.h>


namespace askap {

namespace synthesis {

/// @brief Abstract component
/// @details
/// IComponent is a base class for components working with ComponentEquation
/// examples of components include, e.g. Gaussian or point sources.
/// @note Overloaded virtual function calculate most likely will call
/// templated method for casa::Double and casa::AutoDiff<casa::Double>.
/// We can't have templated method & polymorphism together. 
/// @ingroup measurementequation  
struct IComponent {
  /// virtual destructor to keep the compiler happy
  virtual ~IComponent(); 
  
  /// @brief get number of parameters
  /// @return a number of parameters  this component depends upon
  virtual size_t nParameters() const throw() = 0;
  
  /// @brief calculate visibilities for this component
  /// @details This variant of the method calculates just the visibilities
  /// (without derivatives) for a number of frequencies. The result is stored 
  /// to the provided buffer, which is resized to twice the given 
  /// number of spectral points. Complex values are stored as two consequtive 
  /// double values. The first one is a real part, the second is imaginary part.
  /// @param[in] uvw  baseline spacings (in metres)
  /// @param[in] freq vector of frequencies to do calculations for
  /// @param[in] pol required polarization 
  /// @param[out] result an output buffer used to store values
  virtual void calculate(const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    casa::Stokes::StokesTypes pol,
                    std::vector<double> &result) const = 0;
  
  /// @brief calculate visibilities and derivatives for this component
  /// @details This variant of the method does simultaneous calculations of
  /// the values and derivatives. The result is written to the provided buffer.
  /// See the another version of this method for sizes/description of the buffer
  /// structure.
  /// @param[in] uvw  baseline spacings (in metres)
  /// @param[in] freq vector of frequencies to do calculations for
  /// @param[in] pol required polarization 
  /// @param[out] result an output buffer used to store values
  virtual void calculate(const casa::RigidVector<casa::Double, 3> &uvw,
                    const casa::Vector<casa::Double> &freq,
                    casa::Stokes::StokesTypes pol,
                    std::vector<casa::AutoDiff<double> > &result) const = 0;                    

  /// @brief convert StokesTypes into an index 0..3
  /// @details It is decided that all components have to be defined in
  /// terms of IQUV stokes parameters. It is not prohibited that the 
  /// constructors of actual components accept other stokes parameters like
  /// XX, etc. However, in the latter case, these parameters should be converted
  /// to IQUV at the time of the object construction. Most likely actual 
  /// components will hold an array of fluxes for each stokes parameter. Therefore,
  /// it is necessary to convert quickly from StokesTypes to the index.
  /// This method gives a mapping of I to 0, Q to 1, U to 2 and V to 3. For
  /// other values an exception is thrown.
  /// @param[in] pol required polarization
  /// @return an index (I: 0, Q: 1, U: 2 and V: 3)
  static size_t stokesIndex(casa::Stokes::StokesTypes pol) throw(AskapError);
  
}; 

} // namespace synthesis

} // namespace askap


#endif // #ifndef I_COMPONENT_H

