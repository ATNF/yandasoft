/// @file
/// 
/// @brief utilities related to illumination pattern
/// @details This class is written for experiments with eigenbeams and synthetic beams.
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

#ifndef ILLUMINATION_UTILS_H
#define ILLUMINATION_UTILS_H

#include <boost/shared_ptr.hpp>
#include <gridding/IBasicIllumination.h>
#include <coordinates/Coordinates/CoordinateSystem.h>


#include <string>

namespace askap {

namespace synthutils {

/// @brief utilities related to illumination pattern
/// @details This class is written for experiments with eigenbeams and synthetic beams.
class IlluminationUtils {
public:
   /// @brief constructor
   /// @details 
   /// @param[in] illum illumination pattern to work with   
   /// @param[in] size desired image size
   /// @param[in] cellsize uv-cell size
   /// @param[in] oversample oversampling factor (default 1)
   IlluminationUtils(const boost::shared_ptr<synthesis::IBasicIllumination> &illum,
                     size_t size, double cellsize, size_t oversample);
   
   /// @brief constructor from a parset file 
   /// @details
   /// This version extracts all required parameters from the supplied parset file 
   /// using the same factory, which provides illumination patterns for gridders.
   /// @param[in] parset parset file name 
   IlluminationUtils(const std::string &parset);
   
   /// @brief save the pattern into an image
   /// @details 
   /// @param[in] name file name
   /// @param[in] what type of the image requested, e.g. amplitude (default),
   /// real, imag, phase, complex. Minimum match applies.
   void save(const std::string &name, const std::string &what = "amp");

   /// @brief save the voltage pattern into an image
   /// @details 
   /// @param[in] name file name
   /// @param[in] what type of the image requested, e.g. amplitude (default),
   /// real, imag, phase, complex. Minimum match applies.
   void saveVP(const std::string &name, const std::string &what = "amp");
   
   /// @brief switch to the single element case
   void useSingleElement(); 
   
   /// @brief switch the code to synthetic pattern
   /// @details
   /// @param[in] offsets a matrix with offsets of the elements (number of columns should be 2,
   /// number of rows is the number of elements).
   /// @param[in] weights a vector of complex weights
   void useSyntheticPattern(const casa::Matrix<double> &offsets, 
                            const casa::Vector<casa::Complex> &weights);
protected:
   /// @brief save complex array into an image
   /// @details 
   /// @param[in] name file name
   /// @param[in] coords coordinate system
   /// @param[in] arr array to take the data from
   /// @param[in] what type of the image requested, e.g. amplitude (default),
   /// real, imag, phase, complex. Minimum match applies.
   void saveComplexImage(const std::string &name, const casa::CoordinateSystem &coords,
                         const casa::Array<casa::Complex> &arr,
                         const std::string &what = "amp");
      
private:
   /// @brief illumination pattern corresponding to the single feed
   boost::shared_ptr<synthesis::IBasicIllumination> itsElementIllumination;

   /// @brief illumination pattern to use (may be synthetic)
   boost::shared_ptr<synthesis::IBasicIllumination> itsIllumination;   
   
   /// @brief size of the pattern to work with 
   size_t itsSize;
   
   /// @brief required cell size of the pixellized pattern (wavelengths)
   double itsCellSize;
   
   /// @brief oversampling factor
   size_t itsOverSample;
};

} // namespace synthutils

} // namespace askap

#endif // #ifndef ILLUMINATION_UTILS_H
