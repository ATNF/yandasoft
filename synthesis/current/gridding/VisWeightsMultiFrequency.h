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

#include <gridding/IVisWeights.h>

#include <string>

#include <casa/BasicSL/Complex.h>


namespace askap
{
  namespace synthesis
  {

    /// @brief Class to calculate visibility weights for Multi-Frequency Synthesis
    ///
    /// Blah blah...
    ///
    /// @ingroup gridding
    class VisWeightsMultiFrequency : public IVisWeights
    {
  public:

      /// @brief default constructor
      VisWeightsMultiFrequency();
      
      /// @brief constructor 
      /// @param[in] reffreq reference frequency
      VisWeightsMultiFrequency(casa::Double & reffreq);

      /// @brief copy constructor
      /// @param[in] other input object
      VisWeightsMultiFrequency(const VisWeightsMultiFrequency &other);

      
      ~VisWeightsMultiFrequency();

      virtual IVisWeights::ShPtr clone();
     
      /// @brief Set the context
      /// @param order The order of the Taylor term
      void setParameters(int order);

      /// @brief Calculate the visibility weight.
      /// @param i Sample Index
      /// @param freq frequency
      /// @param pol Polarization index
      float getWeight(int i, double freq, int pol);
      
  protected:


  private:
      // reference frequency.
      casa::Double itsRefFreq;
      // Taylor term 'order'
      casa::Int itsOrder;


    };
  }
}
#endif
