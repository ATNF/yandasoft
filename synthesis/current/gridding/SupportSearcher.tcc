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

#ifndef SUPPORT_SEARCHER_TCC
#define SUPPORT_SEARCHER_TCC

#include <askap/AskapError.h>
#include <profile/AskapProfiler.h>

namespace askap {

namespace synthesis {


/// @brief determine the peak and its position
/// @details This method fills only itsPeakPos and itsPeakVal. It is
/// normally called from one of the search methods, but could be called
/// separately.
/// @param[in] in input 2D matrix with an image
template<typename T>
void SupportSearcher::findPeak(const casa::Matrix<T> &in)
{ 
  ASKAPDEBUGTRACE("SupportSearcher::findPeak");

  itsPeakPos.resize(in.shape().nelements(),casa::False);
  itsPeakPos = 0;
  itsPeakVal = -1;
  #ifdef _OPENMP
  #pragma omp parallel default(shared)
  {
  #pragma omp for
  #endif
  for (int iy=0;iy<int(in.ncolumn());++iy) {
       double tempPeakVal = -1;      
       int tempPeakX = 0, tempPeakY = 0;
       for (int ix=0;ix<int(in.nrow());++ix) {
	    // the following line has been commented out until we find a better work around on the delphinus/minicp bug
	    // See ticket:2307
            //const double curVal = std::abs(in(ix,iy));
	    const double curVal = std::abs(casa::DComplex(in(ix,iy)));
            if(tempPeakVal< curVal) {
               tempPeakX = ix;
               tempPeakY = iy;
               tempPeakVal = curVal;
            }
       }
       #ifdef _OPENMP
       #pragma omp critical
       {
       #endif
       if (itsPeakVal < tempPeakVal) {
           itsPeakPos(0) = tempPeakX;
           itsPeakPos(1) = tempPeakY;
           itsPeakVal = tempPeakVal;
       }
       #ifdef _OPENMP
       }
       #endif
  }
  #ifdef _OPENMP
  }
  #endif
#ifdef ASKAP_DEBUG  
  if (itsPeakVal<0) {
      ASKAPTHROW(CheckError, "An empty matrix has been passed to SupportSearcher::findPeak, please investigate. Shape="<<
                 in.shape());
  }
  if (std::isinf(itsPeakVal) || std::isnan(itsPeakVal)) {
      debugStoreImage(in);
  }
  ASKAPCHECK(!std::isnan(itsPeakVal), "Peak value is not a number, please investigate. itsPeakPos="<<itsPeakPos);
  ASKAPCHECK(!std::isinf(itsPeakVal), "Peak value is infinite, please investigate. itsPeakPos="<<itsPeakPos);
#endif // #ifdef ASKAP_DEBUG

  ASKAPCHECK(itsPeakVal>0.0, "Unable to find peak in the support searcher (in either making a convolution function or fitting the PSF), all values appear to be zero. itsPeakVal=" 
             << itsPeakVal);
}

/// @brief debug method to save the matrix
/// @details This method is used for debugging only and stores given
/// complex matrix (i.e. if NaN appears in the CF calculation). It does
/// nothing for the generic value type.
/// @param[in] in input 2D matrix with an image
template<>
void SupportSearcher::debugStoreImage(const casa::Matrix<casa::Complex> &in);


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
void SupportSearcher::search(const casa::Matrix<T> &in, const double value)
{
  findPeak(in);
  if (value > 0) {
      itsPeakVal = value;
  }
  doSupportSearch(in);
}


/// @brief do actual support search
/// @details This method assumes that peak has already been found and
/// implements the actual search of blc and trc of the support region.
/// @param[in] in input 2D matrix with an image 
template<typename T>
void SupportSearcher::doSupportSearch(const casa::Matrix<T> &in)
{
  ASKAPDEBUGTRACE("SupportSearcher::doSupportSearch");

  ASKAPDEBUGASSERT(in.shape().nelements() == 2);
  ASKAPDEBUGASSERT(itsPeakPos.nelements() == 2);
  itsBLC.resize(2,casa::False);
  itsTRC.resize(2,casa::False);
  itsBLC = -1;
  itsTRC = -1;
  ASKAPCHECK(itsPeakVal>0.0, "A positive peak value of the convolution function is expected inside doSupportSearch, itsPeakVal" << 
                             itsPeakVal);
  ASKAPCHECK(!std::isinf(itsPeakVal), "Peak value is infinite, this shouldn't happen. itsPeakPos="<<itsPeakPos); 
  ASKAPCHECK(itsPeakPos[0]>0 && itsPeakPos[1]>0, "Peak position of the convolution function "<<itsPeakPos<<
             " is too close to the edge, increase maxsupport");
  ASKAPCHECK(itsPeakPos[0] + 1 < int(in.nrow()) && itsPeakPos[1] + 1 < int(in.ncolumn()), 
             "Peak position of the convolution function "<<itsPeakPos<<" is too close to the edge, increase maxsupport");
  
  const double absCutoff = itsCutoff*itsPeakVal;
  #ifdef _OPENMP
  #pragma omp parallel sections
  {
  #pragma omp section
  {
  #endif 
  for (int ix = 0; ix<=itsPeakPos(0); ++ix) {
       if (casa::abs(in(ix, itsPeakPos(1))) > absCutoff) {
           itsBLC(0) = ix;
           break;
       }
  }

  #ifdef _OPENMP
  }
  #pragma omp section
  {
  #endif
  
  for (int iy = 0; iy<=itsPeakPos(1); ++iy) {
       if (casa::abs(in(itsPeakPos(0),iy)) > absCutoff) {
           itsBLC(1) = iy;
           break;
       }
  }

  #ifdef _OPENMP
  }
  #pragma omp section
  {
  #endif

  for (int ix = int(in.nrow())-1; ix>=itsPeakPos(0); --ix) {
       if (casa::abs(in(ix, itsPeakPos(1))) > absCutoff) {
           itsTRC(0) = ix;
           break;
       }
  }
  
  #ifdef _OPENMP
  }
  #pragma omp section
  {
  #endif

  for (int iy = int(in.ncolumn())-1; iy>=itsPeakPos(1); --iy) {
       if (casa::abs(in(itsPeakPos(0),iy)) > absCutoff) {
           itsTRC(1) = iy;
           break;
       }
  }

  #ifdef _OPENMP
  }
  }
  #endif
  
  ASKAPCHECK((itsBLC(0)>=0) && (itsBLC(1)>=0) && (itsTRC(0)>=0) && 
             (itsTRC(1)>=0), "Unable to find the support on one of the coordinates (try decreasing the value of .gridder.cutoff) Effective support is 0. itsBLC="<<itsBLC<<" itsTRC="<<itsTRC<<" itsPeakPos="<<itsPeakPos<<" in.shape()="<<in.shape()<<" absCutoff="<<absCutoff<<" itsPeakVal="<<itsPeakVal);
}


/// @brief search assuming the peak is in the centre
/// @details This search method assumes the peak is in the centre of the
/// image and has a given value. The search starts at the edges and
/// terminated as soon as the absolute value higher than 
/// cutoff*value has been found. Giving value of 1. effectively means 
/// that the cutoff is an absolute cutoff (default).
/// @param[in] in input 2D matrix with an image 
/// @param[in] value assumed peak value
template<typename T>
void SupportSearcher::searchCentered(const casa::Matrix<T> &in, double value)
{
  itsPeakVal = value;
  itsPeakPos.resize(in.shape().nelements(), casa::False);
  ASKAPDEBUGASSERT(itsPeakPos.nelements() == 2);
  itsPeakPos = in.shape();
  itsPeakPos(0)/=2;
  itsPeakPos(1)/=2;
  doSupportSearch(in);  
}


} // namespace synthesis

} // namespace askap


#endif // #ifndef SUPPORT_SEARCHER_TCC


