/// @file
/// 
/// @brief Composite calibration component  (a sum of two or three others).
/// @details This template acts as a composite effect with the resulting
/// Mueller matrix equal to the sum of input Mueller matrices.
/// @note I currently forsee two ways of dealing with the composite effects,
/// especially sums. First, if the effect is solvable it has to be included
/// in the effect chain and used with CalibrationME. This template is
/// intended for this case. The second way is to have a separate composite
/// equation replacing CalibrationME, which adds some effect to the data.
/// It is more appropriate for simulator, which can add some non-solvable 
/// modifications of the data (e.g. noise). The main benefit of this 
/// second approach is an ability to construct the equations more dynamically.
/// The main drawback is inability to solve for parameters using just the 
/// functionality of wrapped classes.
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

#ifndef SUM_H
#define SUM_H

#include <measurementequation/MEComponent.h>
#include <dataaccess/IConstDataAccessor.h>
#include <measurementequation/BlockCDMOperations.h>

namespace askap {

namespace synthesis {

/// @brief Composite calibration component  (a sum of two or three others).
/// @details This template acts as a composite effect with the resulting
/// Mueller matrix equal to the sum of input Mueller matrices.
/// @note I currently forsee two ways of dealing with the composite effects,
/// especially sums. First, if the effect is solvable it has to be included
/// in the effect chain and used with CalibrationME. This template is
/// intended for this case. The second way is to have a separate composite
/// equation replacing CalibrationME, which adds some effect to the data.
/// It is more appropriate for simulator, which can add some non-solvable 
/// modifications of the data (e.g. noise). The main benefit of this 
/// second approach is an ability to construct the equations more dynamically.
/// The main drawback is inability to solve for parameters using just the 
/// functionality of wrapped classes.
/// @ingroup measurementequation
template<typename Effect1,typename  Effect2,typename  Effect3 = MEComponent<false> >
struct Sum : public MEComponent<Effect1::theirFDPFlag || Effect2::theirFDPFlag || Effect3::theirFDPFlag> {

   /// @brief constructor, store reference to paramters
   /// @param[in] par shared pointer to parameters
   inline explicit Sum(const scimath::Params::ShPtr &par) :  
            itsEffect1(par), itsEffect2(par), itsEffect3(par) {}

   /// @brief main method returning Mueller matrix and derivatives
   /// @details This method has to be overloaded (in the template sense) for
   /// all classes representing various calibration effects. CalibrationME
   /// template will call it when necessary.
   /// @param[in] chunk accessor to work with
   /// @param[in] row row of the chunk to work with
   /// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
   /// this effect
   inline scimath::ComplexDiffMatrix get(const accessors::IConstDataAccessor &chunk, 
                                casa::uInt row) const
   { return BlockCDMOperations<Effect1::theirFDPFlag,Effect2::theirFDPFlag,Effect3::theirFDPFlag>::sum(
              itsEffect1.get(chunk,row),itsEffect2.get(chunk,row),itsEffect3.get(chunk,row)); }
   
private:
   /// @brief buffer for the first effect
   Effect1 itsEffect1;
   /// @brief buffer for the second effect
   Effect2 itsEffect2;
   /// @brief buffer for the third effect
   Effect3 itsEffect3;
};


/// @brief specialization for two items only
template<typename Effect1, typename Effect2>
struct Sum<Effect1, Effect2, MEComponent<false> > : public MEComponent<Effect1::theirFDPFlag || Effect2::theirFDPFlag> {

   /// @brief constructor, store reference to paramters
   /// @param[in] par shared pointer to parameters
   inline explicit Sum(const scimath::Params::ShPtr &par) :  
            itsEffect1(par), itsEffect2(par) {}

   /// @brief main method returning Mueller matrix and derivatives
   /// @details This method has to be overloaded (in the template sense) for
   /// all classes representing various calibration effects. CalibrationME
   /// template will call it when necessary.
   /// @param[in] chunk accessor to work with
   /// @param[in] row row of the chunk to work with
   /// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
   /// this effect
   inline scimath::ComplexDiffMatrix get(const accessors::IConstDataAccessor &chunk, 
                                casa::uInt row) const
   { return BlockCDMOperations<Effect1::theirFDPFlag,Effect2::theirFDPFlag,false>::sum(
              itsEffect1.get(chunk,row),itsEffect2.get(chunk,row)); }
   
private:
   /// @buffer first effect
   Effect1 itsEffect1;
   /// @buffer second effect
   Effect2 itsEffect2;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef SUM_H
