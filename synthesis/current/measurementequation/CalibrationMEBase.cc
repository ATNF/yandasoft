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
#include <casa/complex.h>

// own includes
#include <measurementequation/CalibrationMEBase.h>
#include <profile/AskapProfiler.h>

#include <casa/Arrays/MatrixMath.h>
#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <dataaccess/MemBufferDataAccessor.h>
#include <fitting/GenericNormalEquations.h>
#include <fitting/DesignMatrix.h>
#include <fitting/Params.h>
#include <fitting/ComplexDiff.h>
#include <fitting/ComplexDiffMatrix.h>


// std includes
//#include <algorithm>
//#include <functional>
#include <string>
#include <exception>


using casa::IPosition;
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
  casa::Cube<casa::Complex> &rwVis = chunk.rwVisibility();
  ASKAPDEBUGASSERT(rwVis.nelements());
  ASKAPCHECK(itsPerfectVisME, "Perfect ME should be defined before calling CalibrationMEBase::predict");
 
  itsPerfectVisME->predict(chunk);
  if (isFrequencyDependent()) {
      for (casa::uInt row = 0; row < chunk.nRow(); ++row) {
           ComplexDiffMatrix calCDM = buildComplexDiffMatrix(chunk, row);
           casa::Matrix<casa::Complex> thisRow = chunk.visibility().yzPlane(row);
           ASKAPDEBUGASSERT(calCDM.nColumn() == thisRow.nrow() * thisRow.ncolumn());
           ASKAPDEBUGASSERT(calCDM.nRow() == thisRow.ncolumn());
           for (casa::uInt chan=0; chan < thisRow.nrow(); ++chan) {
                ComplexDiffMatrix thisChanCDM = calCDM.extractBlock(chan * calCDM.nRow(),calCDM.nRow());
                ComplexDiffMatrix cdm = thisChanCDM * ComplexDiffMatrix(thisRow.row(chan));
                for (casa::uInt pol = 0; pol < chunk.nPol(); ++pol) {
                     rwVis(row, chan, pol) = cdm[pol].value();
                }
           }
      } 
  } else {
     for (casa::uInt row = 0; row < chunk.nRow(); ++row) {
        ComplexDiffMatrix cdm = buildComplexDiffMatrix(chunk, row) * 
             ComplexDiffMatrix(casa::transpose(chunk.visibility().yzPlane(row)));
       
         for (casa::uInt chan = 0; chan < chunk.nChannel(); ++chan) {
             for (casa::uInt pol = 0; pol < chunk.nPol(); ++pol) {
                // cdm is transposed! because we need a vector for
                // each spectral channel for a proper matrix multiplication
                rwVis(row, chan, pol) = cdm(pol, chan).value();
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
casa::Matrix<casa::Complex> CalibrationMEBase::getCorrectionMatrix(const scimath::ComplexDiffMatrix &cdm)
{
  ASKAPTRACE("CalibrationMEBase::getCorrectionMatrix");
  ASKAPASSERT(cdm.nRow()==cdm.nColumn()); // need to do something about this sometime in the future
  casa::Matrix<casa::Complex> effect(cdm.nRow(),cdm.nColumn());

  casa::Matrix<casa::Complex> reciprocal;
  casa::Complex det=0;
  for (casa::uInt i = 0; i < effect.nrow(); ++i) {
       for (casa::uInt j = 0; j < effect.ncolumn(); ++j) {
            effect(i, j) = cdm(i, j).value();
       }
  }
  invertSymPosDef(reciprocal, det, effect);
  if (casa::abs(det)<1e-5) {
      ASKAPTHROW(AskapError, "Unable to apply gains, determinate too close to 0. D="<<casa::abs(det));           
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
  casa::Cube<casa::Complex> &rwVis = chunk.rwVisibility();
  ASKAPDEBUGASSERT(rwVis.nelements());
  if (isFrequencyDependent()) {
      for (casa::uInt row = 0; row < chunk.nRow(); ++row) {
           ComplexDiffMatrix thisRowCDM = buildComplexDiffMatrix(chunk, row);
           casa::Matrix<casa::Complex> thisRow = rwVis.yzPlane(row);       
           for (casa::uInt chan=0; chan < thisRow.nrow(); ++chan) {
                ComplexDiffMatrix thisChanCDM = thisRowCDM.extractBlock(chan * thisRow.ncolumn(), thisRow.ncolumn());
                casa::Vector<casa::Complex> visPolVector = thisRow.row(chan);
                // no need to transpose here because we deal with vectors
                casa::Matrix<casa::Complex> correctedPolVector = getCorrectionMatrix(thisRowCDM) * visPolVector;
                visPolVector = correctedPolVector;
           }
      }
  } else {
      for (casa::uInt row = 0; row < chunk.nRow(); ++row) {
           ComplexDiffMatrix cdm = buildComplexDiffMatrix(chunk, row);
                    
  
           // we need to transpose the inverse matrix to be able to bulk-process all
           // channels which are presented in a cube slice which is nchan x npol matrix 
          casa::Matrix<casa::Complex> reciprocal = transpose(getCorrectionMatrix(cdm));
              
          casa::Matrix<casa::Complex> thisRow = rwVis.yzPlane(row);       
          casa::Matrix<casa::Complex> temp(thisRow.nrow(),reciprocal.ncolumn(),
                                     casa::Complex(0.,0.)); // = thisRow*reciprocal;
                                     
          // the code below is just a matrix multiplication doing on-the-fly transpose
          for (casa::uInt i = 0; i < temp.nrow(); ++i) {
               for (casa::uInt j = 0; j < temp.ncolumn(); ++j) {
                    for (casa::uInt k = 0; k < thisRow.ncolumn(); ++k) {
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
  const casa::Cube<casa::Complex> &measuredVis = chunk.visibility();
  
  if (isFrequencyDependent()) {
      for (casa::uInt row = 0; row < buffChunk.nRow(); ++row) {
          ComplexDiffMatrix thisRowCDM = buildComplexDiffMatrix(buffChunk, row);
          casa::Matrix<casa::Complex> thisRowPerfectVis = buffChunk.visibility().yzPlane(row);
          casa::Matrix<casa::Complex> thisRowMeasuredVis = measuredVis.yzPlane(row);
          ASKAPDEBUGASSERT(thisRowCDM.nColumn() == thisRowPerfectVis.nrow() * thisRowPerfectVis.ncolumn());
          for (casa::uInt chan=0; chan < thisRowPerfectVis.nrow(); ++chan) {
               ComplexDiffMatrix thisChanCDM = thisRowCDM.extractBlock(chan * thisRowPerfectVis.ncolumn(),thisRowPerfectVis.ncolumn());
               casa::Vector<casa::Complex> perfectVisPolVector = thisRowPerfectVis.row(chan);
               casa::Vector<casa::Complex> measuredVisPolVector = thisRowMeasuredVis.row(chan); 
               ComplexDiffMatrix cdm = thisChanCDM * ComplexDiffMatrix(perfectVisPolVector);
               
               DesignMatrix designmatrix;
               
               // we can probably add below actual weights taken from the data accessor
               designmatrix.addModel(cdm, measuredVisPolVector, 
                       casa::Vector<double>(measuredVisPolVector.nelements(),1.));
               /*
               // temporary hack to investigate the impact of cross-pols
               casa::Vector<double> wtVector(measuredVisPolVector.nelements(),1.);
               ASKAPDEBUGASSERT(wtVector.nelements()>3);
               wtVector[1] = 0.; wtVector[2] = 0.;
               designmatrix.addModel(cdm, measuredVisPolVector, wtVector);
               */
      
               ne.add(designmatrix);               
          }
      }     
  } else {
     // process all frequency channels at once as the effect ComplexDiffMatrix is the same for all of them
     for (casa::uInt row = 0; row < buffChunk.nRow(); ++row) { 
          ComplexDiffMatrix cdm = buildComplexDiffMatrix(buffChunk, row) * 
               ComplexDiffMatrix(casa::transpose(buffChunk.visibility().yzPlane(row)));
          casa::Matrix<casa::Complex> measuredSlice = transpose(measuredVis.yzPlane(row));
       
          DesignMatrix designmatrix;
          // we can probably add below actual weights taken from the data accessor
          designmatrix.addModel(cdm, measuredSlice, 
                    casa::Matrix<double>(measuredSlice.nrow(),
                    measuredSlice.ncolumn(),1.));
      
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
/// @details It is handy to be able to apply the same solution for all beams. 
/// With this flag set, beam=0 solution will be used for all beams.
/// @param[in] flag if true, beam=0 calibration is applied to all beams
void CalibrationMEBase::beamIndependent(bool flag)
{
  ASKAPTHROW(AskapError, "CalibrationMEBase::beamIndependent is not implemented, flag="<<flag);
}

/// @brief normalise gain update amplitudes before application
/// @details It may be useful to only update phases during calibration, for
/// instance during self-calibration with a limited number of antennas.
/// @param[in] flag if true, only phases will be updated when applying
/// calibration solutions
void CalibrationMEBase::phaseOnly(bool flag)
{
  ASKAPTHROW(AskapError, "CalibrationMEBase::phaseOnly is not implemented, flag="<<flag);
}


