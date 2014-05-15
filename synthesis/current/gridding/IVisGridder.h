/// @file
///
/// IVisGridder: Interface definition for visibility gridders
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
#ifndef ASKAP_SYNTHESIS_IVISGRIDDER_H_
#define ASKAP_SYNTHESIS_IVISGRIDDER_H_

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>

#include <fitting/Axes.h>

#include <boost/shared_ptr.hpp>
#include <dataaccess/IDataAccessor.h>
#include <dataaccess/IConstDataAccessor.h>
#include <dataaccess/SharedIter.h>

#include <boost/shared_ptr.hpp>
#include <Common/ParameterSet.h>


#include <gridding/IVisWeights.h>
#include <string>

namespace askap
{
	namespace synthesis
	{

		/// @brief Abstract Base Class for all gridders.
		/// A gridder puts the synthesis data onto a grid and transforms
		/// as necessary. To allow all the important possibilities, the
		/// Fourier transforms are performed here rather than externally.
		/// 
		/// There is a separate class for degridding.
		/// @ingroup gridding
		class IVisGridder
		{
			public:

			/// Shared pointer definition
			typedef boost::shared_ptr<IVisGridder> ShPtr;

			IVisGridder();
			virtual ~IVisGridder();
			
			/// Clone a copy of this Gridder
			virtual ShPtr clone() = 0;

			/// @brief Initialise the gridding
			/// @param axes axes specifications
			/// @param shape Shape of output image: cube: u,v,pol,chan
			/// @param dopsf Make the psf?
			virtual void initialiseGrid(const scimath::Axes& axes,
					const casa::IPosition& shape, const bool dopsf=true) = 0;

                        /// @brief Grid the visibility data.
                        /// @param acc const data accessor to work with
                        virtual void grid(accessors::IConstDataAccessor& acc) = 0;

			/// Form the final output image
			/// @param out Output double precision image or PSF
			virtual void finaliseGrid(casa::Array<double>& out) = 0;

			/// Form the sum of the convolution function squared, multiplied by the weights for each
			/// different convolution function. This is used in the evaluation of the second derivative.
			/// @param out Output double precision sum of weights images
			virtual void finaliseWeights(casa::Array<double>& out) = 0;

			/// @brief Initialise the degridding
			/// @param axes axes specifications
			/// @param image Input image: cube: u,v,pol,chan
			virtual void initialiseDegrid(const scimath::Axes& axes,
					const casa::Array<double>& image) = 0;

			/// @brief Make context-dependant changes to the gridder behaviour
			/// @param[in] context context description
			virtual void customiseForContext(const std::string &context) = 0;
			
			/// @brief set visibility weights
			/// @param[in] viswt shared pointer to visibility weights
			virtual void initVisWeights(const IVisWeights::ShPtr &viswt) = 0;

                        /// @brief Degrid the visibility data.
                        /// @param[in] acc non-const data accessor to work with  
                        virtual void degrid(accessors::IDataAccessor& acc) = 0;

			/// @brief Finalise
			virtual void finaliseDegrid() = 0;
			
			/// @brief static method to create gridder
			/// @details Each gridder should have a static factory method, which is
			/// able to create a particular type of the gridder and initialise it with
			/// the parameters taken form the given parset. It is assumed that the 
			/// method receives a subset of parameters where the gridder name is already
			/// taken out. 
			/// @param[in] parset input parset file
			/// @return a shared pointer to the gridder instance					 
			/// @note This method just throws an exception in this basic interface. It is 
			/// added to ensure that all derived classes have this method defined. We have 
			/// to use a static method as opposed to pure virtual function because we plan to use it
			/// to create a brand new instance of the gridder (and hence no object would
			/// exist at that stage)	
			static ShPtr createGridder(const LOFAR::ParameterSet& parset);
			
			/// @brief check whether the model is empty
			/// @details A simple check allows us to bypass heavy calculations if the input model
			/// is empty (all pixels are zero). This makes sense for degridding only.
			/// @brief true, if the model is empty
			virtual bool isModelEmpty() const = 0; 
		};
	}
}
#endif                                            /*IVISGRIDDER_H_*/
