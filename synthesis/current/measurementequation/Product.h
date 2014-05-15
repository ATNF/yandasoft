/// @file
/// 
/// @brief Composite calibration component  (a product of two others).
/// @details This template acts as a composite effect with the resulting
/// Mueller matrix equal to the matrix product of two input Mueller matrices.
/// @note There are plans to extend the interface to several multipliers
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

#ifndef PRODUCT_H
#define PRODUCT_H

#include <measurementequation/MEComponent.h>
#include <dataaccess/IConstDataAccessor.h>
#include <measurementequation/BlockCDMOperations.h>

namespace askap {

namespace synthesis {


/// @brief Composite calibration component  (a product of two others).
/// @details This template acts as a composite effect with the resulting
/// Mueller matrix equal to the matrix product of two input Mueller matrices.
/// @note There are plans to extend the interface to several multipliers
/// @ingroup measurementequation
template<typename Effect1,typename  Effect2,typename  Effect3 = MEComponent<false> >
struct Product  : public MEComponent<Effect1::theirFDPFlag || Effect2::theirFDPFlag || Effect3::theirFDPFlag> {

   /// @brief constructor, store reference to paramters
   /// @param[in] par shared pointer to parameters
   inline explicit Product(const scimath::Params::ShPtr &par) :  
            itsEffect1(par), itsEffect2(par), itsEffect3(par) {}

   /// @brief main method returning Mueller matrix and derivatives
   /// @details This method has to be overloaded (in the template sense) for
   /// all classes representing various calibration effects. CalibrationME
   /// template will call it when necessary. It returns 
   /// @param[in] chunk accessor to work with
   /// @param[in] row row of the chunk to work with
   /// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
   /// this effect
   inline scimath::ComplexDiffMatrix get(const accessors::IConstDataAccessor &chunk, 
                                casa::uInt row) const
   {  return BlockCDMOperations<Effect1::theirFDPFlag,Effect2::theirFDPFlag,Effect3::theirFDPFlag>::product( 
           itsEffect1.get(chunk,row), itsEffect2.get(chunk,row), itsEffect3.get(chunk,row)); }
   
private:
   /// @brief buffer for the first effect
   Effect1 itsEffect1;
   /// @brief buffer for the second effect
   Effect2 itsEffect2;
   /// @brief buffer for the third effect
   Effect3 itsEffect3;
};

/// @brief specialization for two multipliers only
template<typename Effect1, typename Effect2>
struct Product<Effect1, Effect2, MEComponent<false> > : public MEComponent<Effect1::theirFDPFlag || Effect2::theirFDPFlag> {

   /// @brief constructor, store reference to paramters
   /// @param[in] par shared pointer to parameters
   inline explicit Product(const scimath::Params::ShPtr &par) :  
            itsEffect1(par), itsEffect2(par){}

   /// @brief main method returning Mueller matrix and derivatives
   /// @details This method has to be overloaded (in the template sense) for
   /// all classes representing various calibration effects. CalibrationME
   /// template will call it when necessary. It returns 
   /// @param[in] chunk accessor to work with
   /// @param[in] row row of the chunk to work with
   /// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
   /// this effect
   inline scimath::ComplexDiffMatrix get(const accessors::IConstDataAccessor &chunk, 
                                casa::uInt row) const
   {  return BlockCDMOperations<Effect1::theirFDPFlag,Effect2::theirFDPFlag,false>::product( 
           itsEffect1.get(chunk,row), itsEffect2.get(chunk,row)); }
   
private:
   /// @buffer first effect
   Effect1 itsEffect1;
   /// @buffer second effect
   Effect2 itsEffect2;
};



} // namespace synthesis

} // namespace askap

#endif // #ifndef PRODUCT_H
