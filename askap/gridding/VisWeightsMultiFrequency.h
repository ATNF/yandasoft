///
/// VisWeightsMultiFrequency: 
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
/// @author Urvashi Rau <rurvashi@aoc.nrao.edu>
///
#ifndef VISWEIGHTSMULTIFREQUENCY_H_
#define VISWEIGHTSMULTIFREQUENCY_H_

#include <askap/gridding/IVisWeights.h>

#include <string>

#include <casacore/casa/BasicSL/Complex.h>


namespace askap
{
  namespace synthesis
  {

    /// @brief Class to calculate visibility weights for Multi-Frequency Synthesis
    ///
    /// Blah blah...
    ///
    /// @ingroup gridding
    class VisWeightsMultiFrequency final : virtual public IVisWeights
    {
  public:

      /// @brief default constructor
      VisWeightsMultiFrequency();
      
      /// @brief constructor 
      /// @param[in] reffreq reference frequency
      explicit VisWeightsMultiFrequency(const casacore::Double & reffreq);

      /// @brief copy constructor
      /// @param[in] other input object
      VisWeightsMultiFrequency(const VisWeightsMultiFrequency &other);

      VisWeightsMultiFrequency& operator=(const VisWeightsMultiFrequency &other) = delete;

      /// @brief clone the object
      /// @return shared pointer to the copy
      virtual IVisWeights::ShPtr clone() const override;
     
      /// @brief Set the context
      /// @param order The order of the Taylor term
      virtual void setParameters(int order) override;

      /// @brief Calculate the visibility weight.
      /// @param freq frequency (same units as reffreq)
      virtual float getWeight(double freq) const override;
      
  protected:


  private:
      // reference frequency.
      casacore::Double itsRefFreq;
      // Taylor term 'order'
      casacore::Int itsOrder;


    };
  }
}
#endif
