/// @file 
/// @brief An array of data in the uv-domain
/// @details This class is designed to hold either an illumination pattern or
/// a convolution function generated from it. The common thing between 
/// illumination patterns and convolution functions is that they all are
/// represented by an array in the uv-domain and have a limited support. This
/// class binds together an array of values itself and an optimized coordinate
/// system: uv-cell sizes and an oversampling factor. The centre is always
/// assumed to be in the middle of the interval. 
/// @note It doesn't seem to be necessary to have an hierarchy of classes based
/// on this class. Therefore, this class serves as an interface without an
/// abstract base class defined (and this should deliver a better performance as
/// virtual functions are not used). However, if it is found necessary later to
/// have a tree of classes representing different kind of patterns (e.g.
/// frequency dependent patterns held in a Cube), a proper interface should be
/// created.
///
/// @copyright (c) 2008 CSIRO
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

#ifndef UV_PATTERN_H
#define UV_PATTERN_H

// casa includes
#include <casa/Arrays/Matrix.h>
#include <casa/BasicSL/Complex.h>

namespace askap {

namespace synthesis {

/// @brief An array of data in the uv-domain
/// @details This class is designed to hold either an illumination pattern or
/// a convolution function generated from it. The common thing between 
/// illumination patterns and convolution functions is that they all are
/// represented by an array in the uv-domain and have a limited support. This
/// class binds together an array of values itself and an optimized coordinate
/// system: uv-cell sizes and an oversampling factor. The centre is always
/// assumed to be in the middle of the interval. 
/// @note It doesn't seem to be necessary to have an hierarchy of classes based
/// on this class. Therefore, this class serves as an interface without an
/// abstract base class defined (and this should deliver a better performance as
/// virtual functions are not used). However, if it is found necessary later to
/// have a tree of classes representing different kind of patterns (e.g.
/// frequency dependent patterns held in a Cube), a proper interface should be
/// created.
/// @ingroup gridding
struct UVPattern {
   /// @brief construct a pattern array
   /// @details This constructor initializes the data members, array is
   /// not initialized, just resized
   /// @param[in] uSize size in the direction of u-coordinate
   /// @param[in] vSize size in the direction of v-coordinate
   /// @param[in] uCellSize size of the uv-cell in the direction of 
   ///            u-coordinate (in wavelengths)
   /// @param[in] vCellSize size of the uv-cell in the direction of 
   ///            v-coordinate (in wavelengths)
   /// @param[in] overSample oversampling factor (default is 1)
   UVPattern(casa::uInt uSize, casa::uInt vSize, double uCellSize,
             double vCellSize, casa::uInt overSample = 1) :
      itsArray(uSize,vSize), itsUCellSize(uCellSize),
      itsVCellSize(vCellSize), itsOverSample(overSample),
      itsMaxSupport(uSize > vSize ? uSize : vSize) {} 
   
   /// @brief construct a default pattern array
   /// @details This constructor initializes the data members with their
   /// default constructors. There is no much use of such an object except that
   /// the whole class can be stored in the stl containers.
   UVPattern() {}
   
   /// @brief copy constructor
   /// @details It is necessary as we have non-tivial type of one data member
   /// (casa matrix)
   /// @param[in] other an object to copy from
   UVPattern(const UVPattern &other) : itsArray(other.itsArray.copy()), itsUCellSize(other.itsUCellSize),
             itsVCellSize(other.itsVCellSize), itsOverSample(other.itsOverSample), itsMaxSupport(other.itsMaxSupport) {}

   /// @brief read-only access to a pattern
   /// @details This method allows a direct access to the pattern array
   /// @return const reference to the pattern
   inline const casa::Matrix<casa::DComplex>& pattern() const 
   { return itsArray; }

   /// @brief read-write access to a pattern
   /// @details This method allows a direct access to the pattern array
   /// @return non-const reference to the pattern
   inline casa::Matrix<casa::DComplex>& pattern()
   { return itsArray; }
   
   /// @brief read-only access to an individual element
   /// @details This method allows a direct access to the given element
   /// of the pattern array
   /// @param[in] iu index in the u-coordinate
   /// @param[in] iv index in the v-coordinate
   /// @return const reference to the requested element
   inline const casa::DComplex& operator()(casa::uInt iu, casa::uInt iv) const
   { return itsArray(iu,iv); }
   
   /// @brief read-write access to an individual element
   /// @details This method allows a direct access to the given element
   /// of the pattern array
   /// @param[in] iu index in the u-coordinate
   /// @param[in] iv index in the v-coordinate
   /// @return non-const reference to the requested element
   inline casa::DComplex& operator()(casa::uInt iu, casa::uInt iv)
   { return itsArray(iu,iv); }
   
   /// @brief obtain a size of the uv-cell in the u-direction 
   /// @return a size of the uv-cell in wavelengths
   inline double uCellSize() const { return itsUCellSize;}

   /// @brief obtain a size of the uv-cell in the v-direction 
   /// @return a size of the uv-cell in wavelengths
   inline double vCellSize() const { return itsVCellSize;}
   
   /// @brief obtain an oversampling factor
   /// @return an oversampling factor
   inline casa::uInt overSample() const { return itsOverSample;}
   
   /// @brief obtain a size of the grid in the u-direction
   /// @return a size of the array in the u-direction
   inline casa::uInt uSize() const { return itsArray.nrow(); }
   
   /// @brief obtain a size of the grid in the v-direction
   /// @return a size of the array in the v-direction
   inline casa::uInt vSize() const { return itsArray.ncolumn(); }
   
   /// @brief obtain the upper limit on the support
   /// @return maximum possible support size (known a priori, e.g. dish size)
   inline casa::uInt maxSupport() const { return itsMaxSupport; }
   
   /// @brief assign a new value to the upper limit of the support
   /// @param[in] supp support size (always assume that it is centred)
   inline void setMaxSupport(casa::uInt supp) { itsMaxSupport = supp; }
   
private:
   /// @brief array of values describing the pattern
   casa::Matrix<casa::DComplex> itsArray;
   
   /// @brief a size of the uv-cell in the direction of u-coordinate (in wavelengths)
   double itsUCellSize;

   /// @brief a size of the uv-cell in the direction of v-coordinate (in wavelengths)
   double itsVCellSize;   
   
   /// @brief oversampling factor
   casa::uInt itsOverSample;
   
   /// @brief max support
   /// @details The upper limit can often be placed on the support a priori
   /// All routines filling this buffer class with actual numbers, can
   /// amend this field to speed up calculations of the actual support.
   /// By default it is the largest size of the buffer.
   /// @note Do we need two values, one for each axis?
   casa::uInt itsMaxSupport;
};
		
} // namespace synthesis

} // namespace askap

#endif // #ifndef UV_PATTERN_H
