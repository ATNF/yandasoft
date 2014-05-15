/// @file
///
/// BoxVisGridder: Box-based visibility gridder.
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
#ifndef BOXVISGRIDDER_H_
#define BOXVISGRIDDER_H_

#include <gridding/TableVisGridder.h>
#include <dataaccess/IConstDataAccessor.h>

namespace askap
{
	namespace synthesis
	{
		/// @brief Minimal box-car convolution (aka nearest neighbour) gridder.
		///
		/// @details It doesn't work well but it is fast and simple.
		/// @ingroup gridding
		class BoxVisGridder : public TableVisGridder
		{
			public:

				// Standard two dimensional box gridding
				BoxVisGridder();

				/// Clone a copy of this Gridder
				virtual IVisGridder::ShPtr clone();
				
				/// @brief static method to get the name of the gridder
				/// @details We specify parameters per gridder type in the parset file.
				/// This method returns the gridder name which should be used to extract
				/// a subset of parameters for createGridder method.
				static inline std::string gridderName() { return "Box";}
				
				/// @brief static method to create gridder
			    /// @details Each gridder should have a static factory method, which is
			    /// able to create a particular type of the gridder and initialise it with
			    /// the parameters taken form the given parset. It is assumed that the 
		    	/// method receives a subset of parameters where the gridder name is already
                /// taken out. 
			    /// @param[in] parset input parset file
			    /// @return a shared pointer to the gridder instance					 
			    static IVisGridder::ShPtr createGridder(const LOFAR::ParameterSet& parset);
				

				virtual ~BoxVisGridder();

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
