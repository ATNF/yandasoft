/// @file
///
/// TestLoadGridder: Test visibility gridder.
///
/// @copyright (c) 2009 CSIRO
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
/// @author Ger van Diepen <diepen at astron dot nl>
///
#ifndef TESTLOADGRIDDER_H_
#define TESTLOADGRIDDER_H_

#include <gridding/TableVisGridder.h>
#include <Common/ParameterSet.h>

namespace askap
{
  namespace synthesis
  {
    /// @brief Gridder to test dynamic loading of gridders
    ///
    /// @ingroup testloadgridder
    class TestLoadGridder : public TableVisGridder
    {
    public:

      // Standard two dimensional box gridding
      TestLoadGridder();

      /// Clone this Gridder
      virtual IVisGridder::ShPtr clone();

      virtual ~TestLoadGridder();

      /// @brief Function to create the gridder from a parset.
      /// This function will be registered in the gridder registry.
      static IVisGridder::ShPtr makeGridder (const LOFAR::ParameterSet&);
      /// @brief Return the (unique) name of the gridder.
      static const std::string& gridderName();

      /// @brief Register the gridder create function with its name.
      static void registerGridder();

      /// @brief Initialise the indices
      /// @param[in] acc const data accessor to work with
      virtual void initIndices(const accessors::IConstDataAccessor& acc);

      /// @brief Correct for gridding convolution function
      /// @param image image to be corrected
      virtual void correctConvolution(casa::Array<double>& image);
				
    protected:
      /// Initialize convolution function
      /// @param[in] acc const data accessor to work with
      virtual void initConvolutionFunction(const accessors::IConstDataAccessor& acc);
    };

  }
}
#endif
