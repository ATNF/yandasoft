/// @file
///
/// WStackVisGridder: W projection gridding

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
#ifndef ASKAP_SYNTHESIS_WSTACKVISGRIDDER_H_
#define ASKAP_SYNTHESIS_WSTACKVISGRIDDER_H_

#include <gridding/WDependentGridderBase.h>
#include <dataaccess/IConstDataAccessor.h>

namespace askap
{
	namespace synthesis
	{
		/// @brief Visibility gridder using W stacking
		/// @details The visibilities are gridded using a convolution
		/// function of compact support - actually a spheroidal function.
		/// To correct for the w term in the full synthesis measurement equation 
		/// the data are first partitioned in w and then gridded onto separate
		/// planes. At the end, all planes are Fourier transformed and stacked
		/// after multiplication by the w-dependent complex phasor image.
		///
		/// The scaling is fast in data points, slow in w planes.
		///
		/// @ingroup gridding
		class WStackVisGridder : public WDependentGridderBase
		{
			public:

				/// @brief Construct a gridder for W stacking
				/// @param wmax Maximum baseline (wavelengths)
				/// @param nwplanes Number of w planes
				WStackVisGridder(const double wmax, const int nwplanes);

				virtual ~WStackVisGridder();
				
				/// @brief copy constructor
				/// @details It is required to decouple internal arrays between
				/// input object and the copy
				/// @param[in] other input object
				WStackVisGridder(const WStackVisGridder &other);

				/// @brief Initialise the gridding
				/// @param axes axes specifications
				/// @param shape Shape of output image: u,v,pol,chan
				/// @param dopsf Make the psf?
				virtual void initialiseGrid(const scimath::Axes& axes,
				    const casa::IPosition& shape, const bool dopsf=true);
				
				/// Form the final output image
				/// @param out Output double precision image or PSF
				virtual void finaliseGrid(casa::Array<double>& out);

				/// @brief Initialise the degridding
				/// @param axes axes specifications
				/// @param image Input image: cube: u,v,pol,chan
				virtual void initialiseDegrid(const scimath::Axes& axes,
				    const casa::Array<double>& image);

				/// Clone a copy of this Gridder
				virtual IVisGridder::ShPtr clone();

                /// @brief static method to get the name of the gridder
                /// @details We specify parameters per gridder type in the parset file.
                /// This method returns the gridder name which should be used to extract
                /// a subset of parameters for createGridder method.
                static inline std::string gridderName() { return "WStack";}

                /// @brief static method to create gridder
                /// @details Each gridder should have a static factory method, which is
                /// able to create a particular type of the gridder and initialise it with
                /// the parameters taken form the given parset. It is assumed that the 
                /// method receives a subset of parameters where the gridder name is already
                /// taken out. 
                /// @param[in] parset input parset file
                /// @return a shared pointer to the gridder instance					 
                static IVisGridder::ShPtr createGridder(const LOFAR::ParameterSet& parset);

			protected:
				/// @brief Initialise the indices
				/// @param[in] acc const accessor to work with
				virtual void initIndices(const accessors::IConstDataAccessor& acc);

				/// Offset into grid
				/// @param row Row number
				/// @param pol Polarisation
				/// @param chan Channel number
				virtual int gIndex(int row, int pol, int chan);

				/// Multiply by the phase screen
				/// @param scratch To be multiplied
				/// @param i Index
				void multiply(casa::Array<casa::DComplex>& scratch, int i);
				
				/// Mapping from row, pol, and channel to planes of grid
				casa::Cube<int> itsGMap;
            private:
    	        /// @brief assignment operator
				/// @details It is required as private to avoid being called
				/// @param[in] other input object
				/// @return reference to itself
				WStackVisGridder& operator=(const WStackVisGridder &other);    
		};
	}
}
#endif
