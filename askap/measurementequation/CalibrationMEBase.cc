/// @file
/// 
/// @brief Base class for generic measurement equation for calibration.
/// @details This is a base class for a template designed to represent any 
/// possible measurement equation we expect to encounter in calibration. 
/// It is a result of evolution of the former GainCalibrationEquation, which 
/// will probably be completely substituted by this template in the future. 
/// See CalibrationME template for more details. This class contains all
/// functionality, which doesn't depend on the template parameter.
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

// casa includes
#include <casacore/casa/complex.h>

// own includes
#include <askap/measurementequation/CalibrationMEBase.h>
#include <askap/profile/AskapProfiler.h>

#include <casacore/casa/Arrays/MatrixMath.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/dataaccess/MemBufferDataAccessor.h>
#include <askap/scimath/fitting/GenericNormalEquations.h>
#include <askap/scimath/fitting/DesignMatrix.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/ComplexDiff.h>
#include <askap/scimath/fitting/ComplexDiffMatrix.h>

// std includes
//#include <algorithm>
//#include <functional>
#include <string>
#include <exception>


using casacore::IPosition;
using askap::scimath::DesignMatrix;
using askap::scimath::ComplexDiffMatrix;


using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

/// @brief Standard constructor using the parameters and the
/// data iterator.
/// @param[in] ip Parameters
/// @param[in] idi data iterator
/// @param[in] ime measurement equation describing perfect visibilities
/// @note In the future, measurement equations will work with accessors
/// only, and, therefore, the dependency on iterator will be removed
CalibrationMEBase::CalibrationMEBase(const askap::scimath::Params& ip,
          const IDataSharedIter& idi, 
          const boost::shared_ptr<IMeasurementEquation const> &ime) :
            scimath::Equation(ip), MultiChunkEquation(idi), askap::scimath::GenericEquation(ip),
            GenericMultiChunkEquation(idi), itsPerfectVisME(ime) {}
  
/// @brief Predict model visibilities for one accessor (chunk).
/// @details This version of the predict method works with
/// a single chunk of data only. It seems that all measurement
/// equations should work with accessors rather than iterators
/// (i.e. the iteration over chunks should be moved to the higher
/// level, outside this class). In the future, I expect that
/// predict() without parameters will be deprecated.
/// @param[in] chunk a read-write accessor to work with
void CalibrationMEBase::predict(IDataAccessor &chunk) const
{ 
  ASKAPTRACE("CalibrationMEBase::predict");
  casacore::Cube<casacore::Complex> &rwVis = chunk.rwVisibility();
  ASKAPDEBUGASSERT(rwVis.nelements());
  ASKAPCHECK(itsPerfectVisME, "Perfect ME should be defined before calling CalibrationMEBase::predict");
  itsPerfectVisME->predict(chunk);
  if (isFrequencyDependent()) {
      for (casacore::uInt row = 0; row < chunk.nRow(); ++row) {
           ComplexDiffMatrix calCDM = buildComplexDiffMatrix(chunk, row);
           //casacore::Matrix<casacore::Complex> thisRow = chunk.visibility().yzPlane(row);
           // because the visibility is in the new order (pol,chan,row), taking the 
           // xyPlane() is the same as yzPlane() in the old order (row,chan,pol). However
           // the order (row,col)  of the Matrix is reversed
           casacore::Matrix<casacore::Complex> thisRow = chunk.visibility().xyPlane(row);
           ASKAPDEBUGASSERT(calCDM.nColumn() == thisRow.nrow() * thisRow.ncolumn());
           ASKAPDEBUGASSERT(calCDM.nRow() == thisRow.nrow());
           //for (casacore::uInt chan=0; chan < thisRow.nrow(); ++chan) {
           for (casacore::uInt chan=0; chan < thisRow.ncolumn(); ++chan) {
                ComplexDiffMatrix thisChanCDM = calCDM.extractBlock(chan * calCDM.nRow(),calCDM.nRow());
                //ComplexDiffMatrix cdm = thisChanCDM * ComplexDiffMatrix(thisRow.row(chan));
                ComplexDiffMatrix cdm = thisChanCDM * ComplexDiffMatrix(thisRow.column(chan));
                for (casacore::uInt pol = 0; pol < chunk.nPol(); ++pol) {
                      //rwVis(row, chan, pol) = cdm[pol].value();
                      rwVis(pol, chan, row) = cdm[pol].value();
                }
           }
      } 
  } else {
     for (casacore::uInt row = 0; row < chunk.nRow(); ++row) {
         //ComplexDiffMatrix cdm = buildComplexDiffMatrix(chunk, row) *
         //     ComplexDiffMatrix(casacore::transpose(chunk.visibility().yzPlane(row)));

         // with the new order (pol,chan,row), the xyPlane(row) is the same as the yzPlane(row)
         // in the older order (row,chan,pol) except that the (row,col) is reversed and hence
         // the statement below does not need the transpose() function call

         ComplexDiffMatrix cdm = buildComplexDiffMatrix(chunk, row) *
              ComplexDiffMatrix(chunk.visibility().xyPlane(row));
         for (casacore::uInt chan = 0; chan < chunk.nChannel(); ++chan) {
             for (casacore::uInt pol = 0; pol < chunk.nPol(); ++pol) {
                // cdm is transposed! because we need a vector for
                // each spectral channel for a proper matrix multiplication
                //rwVis(row, chan, pol) = cdm(pol, chan).value();
                rwVis(pol, chan, row) = cdm(pol, chan).value();
                //std::cout << "rwVis(" << pol << "," << chan << "," << row << ") = "
                //          << rwVis(pol, chan, row) << std::endl;
             }
         }
     }
  }
}

/// @brief form correction matrix
/// @details This is a helper method which forms correction matrix out of
/// ComplexDiffMatrix (to encapsulate the code used both in frequency-dependent
/// and frequency-independent cases).
/// @param[in] cdm a square ComplexDiffMatrix describing the effect
/// @return correction matrix
casacore::Matrix<casacore::Complex> CalibrationMEBase::getCorrectionMatrix(const scimath::ComplexDiffMatrix &cdm)
{
  ASKAPTRACE("CalibrationMEBase::getCorrectionMatrix");
  ASKAPASSERT(cdm.nRow()==cdm.nColumn()); // need to do something about this sometime in the future
  casacore::Matrix<casacore::Complex> effect(cdm.nRow(),cdm.nColumn());

  casacore::Matrix<casacore::Complex> reciprocal;
  casacore::Complex det=0;
  for (casacore::uInt i = 0; i < effect.nrow(); ++i) {
       for (casacore::uInt j = 0; j < effect.ncolumn(); ++j) {
            effect(i, j) = cdm(i, j).value();
       }
  }
  invertSymPosDef(reciprocal, det, effect);
  if (casacore::abs(det)<1e-5) {
      ASKAPTHROW(AskapError, "Unable to apply gains, determinate too close to 0. D="<<casacore::abs(det));           
  }
  return reciprocal;  
}


/// @brief correct model visibilities for one accessor (chunk).
/// @details This method corrects the data in the given accessor
/// (accessed via rwVisibility) for the calibration errors 
/// represented by this measurement equation (i.e. an inversion of
/// the matrix has been performed). 
/// @param[in] chunk a read-write accessor to work with
/// @note Need to think what to do in the inversion is unsuccessful
/// e.g. amend flagging information? This is not yet implemented as
/// existing accessors would throw an exception if flagging info is 
/// changed.
void CalibrationMEBase::correct(IDataAccessor &chunk) const
{
  ASKAPTRACE("CalibrationMEBase::correct");
  casacore::Cube<casacore::Complex> &rwVis = chunk.rwVisibility();
  ASKAPDEBUGASSERT(rwVis.nelements());
  if (isFrequencyDependent()) {
      for (casacore::uInt row = 0; row < chunk.nRow(); ++row) {
           ComplexDiffMatrix thisRowCDM = buildComplexDiffMatrix(chunk, row);
           //casacore::Matrix<casacore::Complex> thisRow = rwVis.yzPlane(row);       
           casacore::Matrix<casacore::Complex> thisRow = rwVis.xyPlane(row);       
           //for (casacore::uInt chan=0; chan < thisRow.nrow(); ++chan) {
           for (casacore::uInt chan=0; chan < thisRow.ncolumn(); ++chan) {
                //ComplexDiffMatrix thisChanCDM = thisRowCDM.extractBlock(chan * thisRow.ncolumn(), thisRow.ncolumn());
                ComplexDiffMatrix thisChanCDM = thisRowCDM.extractBlock(chan * thisRow.nrow(), thisRow.nrow());
                //casacore::Vector<casacore::Complex> visPolVector = thisRow.row(chan);
                casacore::Vector<casacore::Complex> visPolVector = thisRow.column(chan);
                // no need to transpose here because we deal with vectors
                casacore::Matrix<casacore::Complex> correctedPolVector = getCorrectionMatrix(thisRowCDM) * visPolVector;
                visPolVector = correctedPolVector;
           }
      }
  } else {
      for (casacore::uInt row = 0; row < chunk.nRow(); ++row) {
           ComplexDiffMatrix cdm = buildComplexDiffMatrix(chunk, row);
                    
  
           // we need to transpose the inverse matrix to be able to bulk-process all
           // channels which are presented in a cube slice which is nchan x npol matrix 
          casacore::Matrix<casacore::Complex> reciprocal = transpose(getCorrectionMatrix(cdm));
              
          //casacore::Matrix<casacore::Complex> thisRow = rwVis.yzPlane(row);       
          casacore::Matrix<casacore::Complex> thisRow = rwVis.xyPlane(row);       
          //casacore::Matrix<casacore::Complex> temp(thisRow.nrow(),reciprocal.ncolumn(),
          //                           casacore::Complex(0.,0.)); // = thisRow*reciprocal;
          casacore::Matrix<casacore::Complex> temp(thisRow.ncolumn(),reciprocal.nrow(),
                                     casacore::Complex(0.,0.)); // = thisRow*reciprocal;
                                     
          // the code below is just a matrix multiplication doing on-the-fly transpose
          for (casacore::uInt i = 0; i < temp.nrow(); ++i) {
               for (casacore::uInt j = 0; j < temp.ncolumn(); ++j) {
                    //for (casacore::uInt k = 0; k < thisRow.ncolumn(); ++k) {
                    for (casacore::uInt k = 0; k < thisRow.nrow(); ++k) {
                         temp(i,j) += thisRow(i,k)*reciprocal(k,j);
                    }  
               }
          }
       
          thisRow = temp;              
      }
  }
}

/// @brief Calculate the normal equation for one accessor (chunk).
/// @details This version of the method works on a single chunk of
/// data only (one iteration).It seems that all measurement
/// equations should work with accessors rather than iterators
/// (i.e. the iteration over chunks should be moved to the higher
/// level, outside this class). In the future, I expect that
/// the variant of the method without parameters will be deprecated.
/// @param[in] chunk a read-write accessor to work with
/// @param[in] ne Normal equations
void CalibrationMEBase::calcGenericEquations(const IConstDataAccessor &chunk,
                              askap::scimath::GenericNormalEquations& ne) const
{  
  ASKAPTRACE("CalibrationMEBase::calcGenericEquations");
  MemBufferDataAccessor  buffChunk(chunk);
  ASKAPDEBUGASSERT(buffChunk.visibility().nelements());
  ASKAPCHECK(itsPerfectVisME, "Perfect ME should be defined before calling CalibrationMEBase::predict");
  
  itsPerfectVisME->predict(buffChunk);
  const casacore::Cube<casacore::Complex> &measuredVis = chunk.visibility();
  
  if (isFrequencyDependent()) {
      for (casacore::uInt row = 0; row < buffChunk.nRow(); ++row) {
          ComplexDiffMatrix thisRowCDM = buildComplexDiffMatrix(buffChunk, row);
          //casacore::Matrix<casacore::Complex> thisRowPerfectVis = buffChunk.visibility().yzPlane(row);
          casacore::Matrix<casacore::Complex> thisRowPerfectVis = buffChunk.visibility().xyPlane(row);
          //casacore::Matrix<casacore::Complex> thisRowMeasuredVis = measuredVis.yzPlane(row);
          casacore::Matrix<casacore::Complex> thisRowMeasuredVis = measuredVis.xyPlane(row);
          ASKAPDEBUGASSERT(thisRowCDM.nColumn() == thisRowPerfectVis.nrow() * thisRowPerfectVis.ncolumn());
          //for (casacore::uInt chan=0; chan < thisRowPerfectVis.nrow(); ++chan) {
          for (casacore::uInt chan=0; chan < thisRowPerfectVis.ncolumn(); ++chan) {
               //ComplexDiffMatrix thisChanCDM = thisRowCDM.extractBlock(chan * thisRowPerfectVis.ncolumn(),thisRowPerfectVis.ncolumn());
               ComplexDiffMatrix thisChanCDM = thisRowCDM.extractBlock(chan * thisRowPerfectVis.nrow(),thisRowPerfectVis.nrow());
               //casacore::Vector<casacore::Complex> perfectVisPolVector = thisRowPerfectVis.row(chan);
               casacore::Vector<casacore::Complex> perfectVisPolVector = thisRowPerfectVis.column(chan);
               //casacore::Vector<casacore::Complex> measuredVisPolVector = thisRowMeasuredVis.row(chan); 
               casacore::Vector<casacore::Complex> measuredVisPolVector = thisRowMeasuredVis.column(chan); 
               ComplexDiffMatrix cdm = thisChanCDM * ComplexDiffMatrix(perfectVisPolVector);
               
               DesignMatrix designmatrix;
               
               // we can probably add below actual weights taken from the data accessor
               designmatrix.addModel(cdm, measuredVisPolVector, 
                       casacore::Vector<double>(measuredVisPolVector.nelements(),1.));
               /*
               // temporary hack to investigate the impact of cross-pols
               casacore::Vector<double> wtVector(measuredVisPolVector.nelements(),1.);
               ASKAPDEBUGASSERT(wtVector.nelements()>3);
               wtVector[1] = 0.; wtVector[2] = 0.;
               designmatrix.addModel(cdm, measuredVisPolVector, wtVector);
               */
      
               ne.add(designmatrix);               
          }
      }     
  } else {
     // process all frequency channels at once as the effect ComplexDiffMatrix is the same for all of them
     for (casacore::uInt row = 0; row < buffChunk.nRow(); ++row) { 
          //ComplexDiffMatrix cdm = buildComplexDiffMatrix(buffChunk, row) * 
               //ComplexDiffMatrix(casacore::transpose(buffChunk.visibility().yzPlane(row)));
          ComplexDiffMatrix cdm = buildComplexDiffMatrix(buffChunk, row) * ComplexDiffMatrix(buffChunk.visibility().xyPlane(row));
          //casacore::Matrix<casacore::Complex> measuredSlice = transpose(measuredVis.yzPlane(row));
          casacore::Matrix<casacore::Complex> measuredSlice = measuredVis.xyPlane(row);
       
          DesignMatrix designmatrix;
          // we can probably add below actual weights taken from the data accessor
          //designmatrix.addModel(cdm, measuredSlice, 
          //          casacore::Matrix<double>(measuredSlice.nrow(),
          //          measuredSlice.ncolumn(),1.));
          designmatrix.addModel(cdm, measuredSlice, 
                    casacore::Matrix<double>(measuredSlice.nrow(), measuredSlice.ncolumn(),1.));
      
          ne.add(designmatrix);
     }
  }
}                                   

/// @brief determines whether to scale the noise estimate
/// @details This is one of the configuration methods, it controlls
/// whether the noise estimate is scaled aggording to applied calibration
/// factors or not.
/// @param[in] scale if true, the noise will be scaled
void CalibrationMEBase::scaleNoise(bool scale)
{
  ASKAPTHROW(AskapError, "CalibrationMEBase::scaleNoise is not implemented, scale="<<scale);
}
  
/// @brief determines whether to allow data flagging
/// @details This is one of the configuration methods, it controlls
/// the behavior of the correct method in the case when the matrix inversion
/// fails. If data flagging is allowed, corresponding visibilities are flagged
/// otherwise an exception is thrown.
/// @param[in] flag if true, the flagging is allowed
void CalibrationMEBase::allowFlag(bool flag)
{
  ASKAPTHROW(AskapError, "CalibrationMEBase::allowFlag is not implemented, flag="<<flag);
}

/// @brief determines whether beam=0 calibration is used for all beams or not
/// @details It is handy to be able to apply the same solution for all beams. With 
/// this flag set, beam=0 solution will be used for all beams.
/// @param[in] flag if true, beam=0 calibration is applied to all beams
void CalibrationMEBase::beamIndependent(bool flag)
{
  ASKAPTHROW(AskapError, "CalibrationMEBase::beamIndependent is not implemented, flag="<<flag);
}


