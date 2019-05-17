/// @file
/// 
/// @brief Helper class working with block ComplexDiffMatrix
/// @details This templated class encapsulates product and summation 
/// of ComplexDiffMatrix. Specialisation allows to bypass unnecessary operations
/// if the effects are simple (i.e. without frequency dependence)
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

#ifndef SYNTHESIS_BLOCK_CDM_OPERATIONS_H
#define SYNTHESIS_BLOCK_CDM_OPERATIONS_H

#include<fitting/ComplexDiffMatrix.h>

namespace askap {

namespace synthesis {

/// @brief helper class to do multiplications
/// @details FDP1..FDP3 describe whether the appropriate ComplexDiffMatrix is a block one
template<bool FDP1, bool FDP2, bool FDP3>
struct BlockCDMOperations {

 /// @brief two-operand product
 /// @param[in] cdm1 first matrix
 /// @param[in] cdm2 second matrix 
 /// @return product matrix
 inline static scimath::ComplexDiffMatrix product(const scimath::ComplexDiffMatrix &cdm1, 
                  const scimath::ComplexDiffMatrix &cdm2);

 /// @brief 3-operand product
 /// @param[in] cdm1 first matrix
 /// @param[in] cdm2 second matrix 
 /// @param[in] cdm2 third matrix 
 /// @return product matrix
 inline static scimath::ComplexDiffMatrix product(const scimath::ComplexDiffMatrix &cdm1, 
                  const scimath::ComplexDiffMatrix &cdm2, const scimath::ComplexDiffMatrix &cdm3);

 /// @brief two-operand sum
 /// @param[in] cdm1 first matrix
 /// @param[in] cdm2 second matrix 
 /// @return sum matrix
 inline static scimath::ComplexDiffMatrix sum(const scimath::ComplexDiffMatrix &cdm1, 
                  const scimath::ComplexDiffMatrix &cdm2);

 /// @brief 3-operand sum
 /// @param[in] cdm1 first matrix
 /// @param[in] cdm2 second matrix 
 /// @param[in] cdm2 third matrix 
 /// @return sum matrix
 inline static scimath::ComplexDiffMatrix sum(const scimath::ComplexDiffMatrix &cdm1, 
                  const scimath::ComplexDiffMatrix &cdm2, const scimath::ComplexDiffMatrix &cdm3);
}; 

/// @brief specialised version of the helper for frequency-independent case
template<> 
struct BlockCDMOperations<false,false,false> {

 /// @brief two-operand product
 /// @param[in] cdm1 first matrix
 /// @param[in] cdm2 second matrix 
 /// @return product matrix
 inline static scimath::ComplexDiffMatrix product(const scimath::ComplexDiffMatrix &cdm1, 
                  const scimath::ComplexDiffMatrix &cdm2);

 /// @brief 3-operand product
 /// @param[in] cdm1 first matrix
 /// @param[in] cdm2 second matrix 
 /// @param[in] cdm2 third matrix 
 /// @return product matrix
 inline static scimath::ComplexDiffMatrix product(const scimath::ComplexDiffMatrix &cdm1, 
                  const scimath::ComplexDiffMatrix &cdm2, const scimath::ComplexDiffMatrix &cdm3);

 /// @brief specialised version of the 2-operand sum for frequency-independent case
 /// @param[in] cdm1 first matrix
 /// @param[in] cdm2 second matrix 
 /// @return sum matrix
 inline static scimath::ComplexDiffMatrix sum(const scimath::ComplexDiffMatrix &cdm1, 
                  const scimath::ComplexDiffMatrix &cdm2);

 /// @brief specialised version of the 3-operand sum for frequency-independent case
 /// @param[in] cdm1 first matrix
 /// @param[in] cdm2 second matrix 
 /// @param[in] cdm2 third matrix 
 /// @return sum matrix
 inline static scimath::ComplexDiffMatrix sum(const scimath::ComplexDiffMatrix &cdm1, 
                  const scimath::ComplexDiffMatrix &cdm2, const scimath::ComplexDiffMatrix &cdm3);

}; 


} // namespace synthesis

} // namespace askap

#include <measurementequation/BlockCDMOperations.tcc>

#endif // #ifndef SYNTHESIS_BLOCK_CDM_OPERATIONS_H

