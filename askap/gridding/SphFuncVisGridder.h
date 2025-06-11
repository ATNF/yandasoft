/// @file SphFuncVisGridder.h
///
/// @copyright (c) 2007,2015,2016 CSIRO
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
#ifndef SPHVISGRIDDER_H_
#define SPHVISGRIDDER_H_

#include <askap/gridding/TableVisGridder.h>
#include <askap/dataaccess/IConstDataAccessor.h>
#include <askap/scimath/utils/SpheroidalFunction.h>

namespace askap
{
	namespace synthesis
	{

		/// @brief SphFuncVisGridder: Spheroidal function-based visibility gridder.
		/// @details The gridding function is a prolate spheroidal function identical to the
		/// one used in AIPS, AIPS++, and probably other packages. At some point
		/// we should revisit the tradeoffs since the choice to use this was made
		/// about twenty years ago and computers are quite different now.
		///
		/// The spheroidal function has m = 6, alpha = 1 using the rational
		/// approximations discussed by fred schwab in 'indirect imaging'.
		/// The gridding function is (1-nu**2)*grdsf(nu) where nu is the distance
		/// to the edge. the grid correction function is just 1/grdsf(nu) where nu
		/// is now the distance to the edge of the image.
		/// @ingroup gridding
		class SphFuncVisGridder : public TableVisGridder
		{
                      public:

				/// @brief Standard two dimensional gridding
				/// @param[in] alpha spheroidal function alpha value
				/// @param[in] support support size in pixels (spheroidal
				/// function with m=2*support will be generated)
				/// @param[in] oversample number of oversampling planes
				explicit SphFuncVisGridder(const float alpha = 1.,
                                           const int support = 3,
                                           const int oversample = 128);

				virtual ~SphFuncVisGridder();

				/// Clone a copy of this Gridder
				virtual IVisGridder::ShPtr clone();

				/// @brief static method to get the name of the gridder
				/// @details We specify parameters per gridder type in the parset file.
				/// This method returns the gridder name which should be used to extract
				/// a subset of parameters for createGridder method.
				static inline std::string gridderName() { return "SphFunc";}

				/// @brief static method to create gridder
			    /// @details Each gridder should have a static factory method, which is
			    /// able to create a particular type of the gridder and initialise it with
			    /// the parameters taken form the given parset. It is assumed that the
			    /// method receives a subset of parameters where the gridder name is already
			    /// taken out.
			    /// @param[in] parset input parset file
			    /// @return a shared pointer to the gridder instance
			    static IVisGridder::ShPtr createGridder(const LOFAR::ParameterSet& parset);

				/// @brief Correct for gridding convolution function
				/// @details Doing the Spheroidal grid correction is used in
				/// various places, this static function makes it more widely available
				/// @param image image to be corrected
				/// @param[in] sf spheroidal function to use
				/// @param[in] support support size in pixels (spheroidal
				/// function with m=2*support will be generated)
				/// @param[in] interpolate if true, interpolate the edge values
				static void correctConvolution(casacore::Array<imtype>& image,
					scimath::SpheroidalFunction& sf, int support = 3,
					bool interpolate = true);

			protected:
				/// @brief Initialize the convolution function
				/// @param[in] acc const data accessor to work with
				virtual void initConvolutionFunction(const accessors::IConstDataAccessor& acc);

				/// @brief Initialise the indices
				/// @param[in] acc const data accessor to work with
				virtual void initIndices(const accessors::IConstDataAccessor& acc);

				/// Correct for gridding convolution function
				/// @param image image to be corrected
				virtual void correctConvolution(casacore::Array<imtype>& image);

				/// Calculate prolate spheroidal function
				/// @param nu Argument for spheroidal function
				inline double grdsf(double nu) const { return itsSphFunc(nu); }

				//double grdsf1(double nu) const;

				/// @brief calculator of spheroidal function
				scimath::SpheroidalFunction itsSphFunc;


				/// @brief iterpolate the spheroidal function at nu=1
                /// @details The function is undefined and set to zero at nu=1,
                /// but that is not the numerical limit. Setting itsInterp true
                /// will use neighbouring values to estimate it (to 2nd order).
                /// @param[in] func vector to be interpolated
                template<typename T>
                static void interpolateEdgeValues(casacore::Vector<T> &func);

                /// @brief check whether to interpolate spheroidal function
                /// @return the value of itsInterp
                /// @note there is a bit of the technical debt around how the interpolation
                /// is controlled. To limit the spread of the issue encapsulate all access in
                /// this method to avoid accessing data member from derived classes.
                /// The setter method can be written if needed, although I (MV) sense that
                /// the logic should be moved upstream into scimath.
                inline bool doInterpolation() const { return itsInterp;}

          private:
                /// @brief whether to iterpolate the spheroidal function at nu=1
                /// @details The function is undefined and set to zero at nu=1,
                /// but that is not the numerical limit. Setting itsInterp true
                /// will use neighbouring values to estimate it (to 2nd order).
                /// @note We probably can get rid of this data member as it is only set to true by the looks of it.
                /// I (MV) will make it const to highlight the fact that we don't change it and channel all access to it
                /// via doInterpolation method to make sure no important use case has been overlooked. We need to clear this
                /// technical debt at some point.
                const bool itsInterp;

                /// @brief prolate spheroidal alpha parameter
                double itsAlpha;

          };

	}
}

#include <askap/gridding/SphFuncVisGridder.tcc>

#endif
