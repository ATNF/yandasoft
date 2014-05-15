/// @file
/// 
/// @brief Generic template to iterate inner contour of a 2D image
/// @details Several applications across the synthesis require estimation of
/// some statistics along the inner contour around the peak in a 2D image.
/// This generic template acts as an iterator over the points of the given contour
/// enclosing the peak. Each value point referenced by this iterator is an IPosition of
/// a contour point. The points may appear at an arbitrary order, so some kind of sorting
/// is necessary if one wants to join the nearest neighbours. The contour is defined by
/// a predicate (template argument). It is a locus of points where the predicate becomes true,
/// which is closest to the maximum. The predicate can be an arbitrary type with the operator()
/// defined, which receives the value of the array and returns true or false. This class
/// is generic enough to be moved to an upper level (e.g. at least to Base, but may be even
/// to casacore). If it is ever going to become a part of casacore, it would be good to
/// generalize the class to be able to get several contours at once. In this case
/// iterator can return IPosition and the contour index.
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
///

#ifndef CONTOUR_FINDER_H
#define CONTOUR_FINDER_H

// askap includes - just the exception class
#include <askap/AskapError.h>

// casa includes
#include <casa/Arrays/Array.h>
#include <casa/Arrays/IPosition.h>

namespace askap {

namespace synthesis {

/// @brief Generic template to iterate inner contour of a 2D image
/// @details Several applications across the synthesis require estimation of
/// some statistics along the inner contour around the peak in a 2D image.
/// This generic template acts as an iterator over the points of the given contour
/// enclosing the peak. Each value point referenced by this iterator is an IPosition of
/// a contour point. The points may appear at an arbitrary order, so some kind of sorting
/// is necessary if one wants to join the nearest neighbours. The contour is defined by
/// a predicate (template argument). It is a locus of points where the predicate becomes true,
/// which is closest to the maximum. The predicate can be an arbitrary type with the operator()
/// defined, which receives the value of the array and returns true or false. This class
/// is generic enough to be moved to an upper level (e.g. at least to Base, but may be even
/// to casacore). If it is ever going to become a part of casacore, it would be good to
/// generalize the class to be able to get several contours at once. In this case
/// iterator can return IPosition and the contour index.
///
/// template types:
/// T data type of the array
/// P predicate type
/// @note Reference semantics used for the array.
/// @ingroup measurementequation
template<typename T, typename P>
class ContourFinder {
public:
  /// @brief initialize the finder to work with the given array
  /// @details This is a basic constructor, which stores a reference to the 
  /// working array inside and rewinds the iterator to the first point.
  /// The array can have any number of dimensions, but only two first will be used in the
  /// search (i.e. the contour is a curve, rather than a surface). 
  /// It is also possible to give a central position for the contour search (default is to 
  /// search for a peak). If peak position is given, it should have the same dimensionality as
  /// the array.
  /// @param[in] array an array to work with
  /// @param[in] pred predicate to define the contour
  /// @param[in] peak peak position around which the search is performed (allows to work with
  ///            local optima). The predicate should give false for this point to get a
  ///            sensible result out. Default is IPosition(1,-1), which means to search for
  ///            a maximum and use its position.
  /// @param[in] clip if true, a contour will always be closed by returning edge pixels if
  ///            contour goes beyond the array  
  ContourFinder(const casa::Array<T> &array, const P &pred,
                const casa::IPosition &peak = casa::IPosition(1,-1),
                bool clip = true);
  
  /// @brief default constructor, serves as an end-mark
  /// @details This constructor makes an iterator, which is equivalent to the end method 
  /// for stl containers.
  ContourFinder();
  
  /// @brief Comparison operator
  /// @details It checks whether the iterator reached an end. Only comparison with an end-mark
  /// is allowed.
  /// @param[in] other another iterator (e.g. an end mark)
  /// @return true, if iterators are equal
  bool operator==(const ContourFinder<T,P> &other) const;
  
  /// @brief Comparison operator
  /// @details It checks whether the iterator has more points to iterate over. Only comparison 
  /// with an end-mark is allowed.
  /// @param[in] other another iterator (e.g. an end mark)
  /// @return true, if iterators are not equal
  bool operator!=(const ContourFinder<T,P> &other) const;
  
  /// @brief access operator
  /// @details It returns IPosition for the current point of the contour.
  /// @return const reference to the currnt point of the contour.
  const casa::IPosition& operator*() const;
   
  /// @brief access operator
  /// @details It returns IPosition for the current point of the contour by pointer.
  /// @return const pointer to the currnt point of the contour.
  const casa::IPosition* operator->() const;
  
  /// @brief rewind the iterator
  /// @details This method rewinds the iterator to its initial state. It should not be called
  /// for an end-mark iterator.
  /// @return a reference to itself (to be able to call stl algorithms in a more concise way)
  ContourFinder<T,P>& init();
  
  /// @brief increment operator
  /// @details It makes a step to the next contour point.
  /// @return a reference to itself (to allow chaining)
  ContourFinder<T,P>& operator++();
protected:
  
  /// @brief a helper method to find a contour point searching along the first coordinate
  /// @details This method is called from both init() and operator++() methods.
  /// It modifies itsTestedPosition, which can be left beyond the image on the first coordinate
  /// @return true if a contour point is found, false if the edge is reached.
  bool searchAlongFirstAxis();
    
private:
  /// @brief array to work with
  casa::Array<T> itsArray;
  
  /// @brief predicate
  P itsPredicate;
  
  /// @brief position of the peak 
  /// @details This data member has the same shape as the array
  casa::IPosition itsPeak;
  
  /// @brief true, if contour is clipped at the image edges
  /// @details If false, the contour may be broken if it goes beyond the image dimensions
  bool itsDoClip;
  
  /// @brief increment for each dimension
  /// @details The code first tests the positive increment and then the negative one. 
  casa::IPosition itsIncrements;
  
  /// @brief current tested coordinates
  casa::IPosition itsTestedPosition;
  
  /// @brief true, if this object is an end-mark
  bool itsIsEndMark;
};


} // namespace synthesis

} // namespace askap

#endif // #ifndef CONTOUR_FINDER_H

#include <measurementequation/ContourFinder.tcc>

