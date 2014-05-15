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

#include <casa/Arrays/ArrayMath.h>

#ifndef CONTOUR_FINDER_TCC
#define CONTOUR_FINDER_TCC

namespace askap {

namespace synthesis {

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
template<typename T, typename P> 
ContourFinder<T,P>::ContourFinder(const casa::Array<T> &array, const P &pred,
                const casa::IPosition &peak, bool clip) : itsArray(array),
                itsPredicate(pred), itsPeak(peak), itsDoClip(clip),
                itsIsEndMark(false) 
{
  ASKAPASSERT(array.shape().size());
  ASKAPASSERT(array.shape()[0]>1);
  if (this->itsPeak == casa::IPosition(1,-1)) {
      T minVal, maxVal;
      casa::IPosition minPos, maxPos;
      minMax(minVal,maxVal,minPos,maxPos,this->itsArray);
      this->itsPeak = abs(minVal) > abs(maxVal) ? minPos : maxPos;
  }
  this->init();
}

/// @brief default constructor, serves as an end-mark
/// @details This constructor makes an iterator, which is equivalent to the end method 
/// for stl containers.
template<typename T, typename P> 
ContourFinder<T,P>::ContourFinder() : itsIsEndMark(true) {}

/// @brief Comparison operator
/// @details It checks whether the iterator reached an end. Only comparison with an end-mark
/// is allowed.
/// @param[in] other another iterator (e.g. an end mark)
/// @return true, if iterators are equal
template<typename T, typename P> 
bool ContourFinder<T,P>::operator==(const ContourFinder<T,P> &other) const
{
  if (this->itsIsEndMark) {
      return other.itsEndMark;
  } 
  return !other.itsEndMark;
}

/// @brief Comparison operator
/// @details It checks whether the iterator has more points to iterate over. Only comparison 
/// with an end-mark is allowed.
/// @param[in] other another iterator (e.g. an end mark)
/// @return true, if iterators are not equal
template<typename T, typename P> 
bool ContourFinder<T,P>::operator!=(const ContourFinder<T,P> &other) const
{
  if (this->itsIsEndMark) {
      return !other.itsEndMark;
  } 
  return other.itsEndMark;  
}

/// @brief access operator
/// @details It returns IPosition for the current point of the contour.
/// @return const reference to the currnt point of the contour.
template<typename T, typename P> 
const casa::IPosition& ContourFinder<T,P>::operator*() const
{
  return this->itsTestedPosition;
}
   
/// @brief access operator
/// @details It returns IPosition for the current point of the contour by pointer.
/// @return const pointer to the currnt point of the contour.
template<typename T, typename P> 
const casa::IPosition* ContourFinder<T,P>::operator->() const
{
  return &(this->itsTestedPosition);
}
  
/// @brief rewind the iterator
/// @details This method rewinds the iterator to its initial state. It should not be called
/// for an end-mark iterator.
/// @return a reference to itself (to be able to call stl algorithms in a more concise way)
template<typename T, typename P> 
ContourFinder<T,P>& ContourFinder<T,P>::init()
{
  this->itsIncrements = casa::IPosition(itsArray.shape().size(),1);
  this->itsTestedPosition = itsPeak;
  if (!this->searchAlongFirstAxis()) {
      this->operator++();
  } // if tested position beyond the image and the clipping option is unused
  return *this;
}

/// @brief a helper method to find a contour point searching along the first coordinate
/// @details This method is called from both init() and operator++() methods.
/// It modifies itsTestedPosition, which can be left beyond the image on the first coordinate
/// @return true if a contour point is found, false if the edge is reached.
template<typename T, typename P> 
bool ContourFinder<T,P>::searchAlongFirstAxis()
{
  for (this->itsTestedPosition[0] = this->itsPeak[0]; (this->itsTestedPosition[0] < this->itsArray.shape()[0]) && 
          (this->itsTestedPosition[0]>=0); this->itsTestedPosition[0] += this->itsIncrements[0]) {
       if (this->itsPredicate(this->itsArray(this->itsTestedPosition))) {
           return  true;
       } 
  }
  // reached the end of the image
  if (this->itsDoClip) {
      // result is the edge
      this->itsTestedPosition[0] = this->itsIncrements[0]<0 ? 0 : this->itsArray.shape()[0]-1;
      return true;
  }
  return false;
}

  
/// @brief increment operator
/// @details It makes a step to the next contour point.
/// @return a reference to itself (to allow chaining)
template<typename T, typename P> 
ContourFinder<T,P>& ContourFinder<T,P>::operator++()
{
  this->itsTestedPosition[0] = this->itsPeak[0];
  // loop until all points are found
  // the code would normally exist by a return operator inside the loop
  while (!this->itsIsEndMark) {
     if (this->itsIncrements[0] == 1) {
         // reverse the direction of the search
         this->itsIncrements[0] = -1;
         if (this->searchAlongFirstAxis()) {
             return *this;
         }
     } else {
       this->itsIncrements[0] = 1;
       // need increment of the other coordinate
       if (this->itsArray.ndim() == 1) {
           this->itsIsEndMark = true;
       } else {
           this->itsTestedPosition[1] += this->itsIncrements[1];
           if (this->itsTestedPosition[1] >= this->itsArray.shape()[1]) {
               this->itsTestedPosition[1] = this->itsPeak[1]-1;
               this->itsIncrements[1] = -1;
           } 
           if (this->itsTestedPosition[1]<0) {
               this->itsIsEndMark = true;
           }
       }
       if (!this->itsIsEndMark) {
           if (this->searchAlongFirstAxis()) {
               return *this;
           }       
       }
     }
  }
  return *this;
}


} // namespace synthesis

} // namespace askap

#endif // #ifndef CONTOUR_FINDER_TCC
