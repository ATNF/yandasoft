/// @file
///
/// This file contains utilities to search for a support of the convolution function.
/// They could, in principle, be moved to a higher level (to Base), but left here
/// for now as they are not logically a part of fitting, but would introduce
/// a casacore dependence if moved to askap.
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

#ifndef SUPPORT_SEARCHER_H
#define SUPPORT_SEARCHER_H

#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/IPosition.h>
#include <casa/BasicSL/Complex.h>


namespace askap {

namespace synthesis {

/// @brief A utility class searching for a support
/// @details The methods of these class operate on a 2D matrix and search
/// for a support
class SupportSearcher {
public:
   /// @brief initialize the searcher with some cutoff
   /// @details
   /// @param[in] cutoff the cutoff value. The meaning could be either relative
   /// (with respect to the absolute peak of the image) or absolute, depending
   /// on which method is used.
   explicit SupportSearcher(double cutoff);
   
   /// @brief obtain a cutoff value
   /// @details This method simply returns the current cutoff value
   /// @return cutoff value
   inline double cutoff() const { return itsCutoff;}
   
   /// @brief obtain peak position
   /// @details The method throws exception if no prior support search has
   /// been done.
   /// @return peak position determined during the last search for support
   casa::IPosition peakPos() const;
   
   /// @brief obtain peak value
   /// @details The method throws exception if no prior support search has
   /// been done.
   /// @return peak value determined during the last search for support
   double peakVal() const;
   
   /// @brief obtain the bottom left corner of the support
   /// @details This method returns the bottom left corner of the support. It
   /// throws an exception of no prior search for support has been done.
   /// @return bottom left corner of the support
   casa::IPosition blc() const;
   
   /// @brief obtain the top right corner of the support
   /// @details This method returns the bottom left corner of the support. It
   /// throws an exception of no prior search for support has been done.
   /// @return bottom left corner of the support
   casa::IPosition trc() const;
   
   /// @brief obtain a size of the smallest square support
   /// @details This method essentially returns the largest length across both 
   /// axes (i.e. max(trc-blc)). It throws an exception if no prior search for
   /// support has been done.
   /// @return support size
   casa::uInt support() const;
   
   /// @brief obtain a size of the smallest symmetrical square support
   /// @details This method returns the smallest square support, which is
   /// symmetrical with respect to the centre.
   /// @param[in] shape defines the centre of symmetry (as shape/2)
   casa::uInt symmetricalSupport(const casa::IPosition &shape) const;
   
   /// @brief search assuming the peak is in the centre
   /// @details This search method assumes the peak is in the centre of the
   /// image and has a given value. The search starts at the edges and
   /// terminated as soon as the absolute value higher than 
   /// cutoff*value has been found. Giving value of 1. effectively means 
   /// that the cutoff is an absolute cutoff (default).
   /// @param[in] in input 2D matrix with an image 
   /// @param[in] value assumed peak value
   template<typename T>
   void searchCentered(const casa::Matrix<T> &in, double value = 1.);
   
   /// @brief determine the peak and its position
   /// @details This method fillss only itsPeakPos and itsPeakVal. It is
   /// normally called from one of the search methods, but could be called
   /// separately.
   /// @param[in] in input 2D matrix with an image
   template<typename T>
   void findPeak(const casa::Matrix<T> &in);
   
   /// @brief full search which determines the peak
   /// @details This search method doesn't assume anything about the peak and
   /// searches for its position and peak beforehand. The search starts at the
   /// edges and progresses towards the peak. The edge of the support region
   /// is where the value first time exceeds the cutoff*peakVal, or cutoff*value
   /// if value given as a second parameter is positive
   /// @param[in] in input 2D matrix with an image 
   /// @param[in] value optional peak value, if a positive value is given it will be used
   /// instead of the peak amplitude (although the positon of the peak will still be searched for)
   template<typename T>
   void search(const casa::Matrix<T> &in, const double value = -1.);
protected:
    
   /// @brief do actual support search
   /// @details This method assumes that peak has already been found and
   /// implements the actual search of blc and trc of the support region.
   /// @param[in] in input 2D matrix with an image 
   template<typename T>
   void doSupportSearch(const casa::Matrix<T> &in);

   /// @brief debug method to save the matrix
   /// @details This method is used for debugging only and stores given
   /// complex matrix (i.e. if NaN appears in the CF calculation). It does
   /// nothing for the generic value type.
   /// @param[in] in input 2D matrix with an image
   template<typename T>
   static void debugStoreImage(const casa::Matrix<T> &in) {}
         
private:
   /// @brief relative cutoff level (from the absolute peak)
   double itsCutoff;
   /// @brief Peak position (either assumed or found)
   casa::IPosition itsPeakPos;
   /// @brief peak value (either assumed or found)
   double itsPeakVal;
   /// @brief Bottom left corner of the support
   casa::IPosition itsBLC;
   /// @brief Top right corner of the support
   casa::IPosition itsTRC;
};

} // namespace synthesis

} // namespace askap

// implementation of template functions
#include <gridding/SupportSearcher.tcc>

#endif // #ifndef SUPPORT_SEARCHER_H
