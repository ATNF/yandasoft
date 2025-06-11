/// @file
/// @brief utility to extract visibility data into an image (generalisation of fringetest)
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


#include <askap/dataaccess/TableDataSource.h>
// #include <askap/askap_accessors.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".extractdata");

#include <askap/askap/AskapError.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/scimath/utils/PolConverter.h>

#include <askap/dataaccess/TableManager.h>
#include <askap/dataaccess/IDataConverterImpl.h>
#include <askap/dataaccess/ParsetInterface.h>
#include <askap/scimath/fft/FFTWrapper.h>
#include <casacore/images/Images/ImageFITSConverter.h>
#include <casacore/images/Images/PagedImage.h>

// casa
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <askap/scimath/utils/ImageUtils.h>
#include <casacore/casa/Arrays/Cube.h>

#include <Common/ParameterSet.h>
#include <askap/utils/CommandLineParser.h>



// std
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

using namespace askap;
using namespace askap::accessors;

void analyseDelay(const casa::Matrix<casa::Complex> &fringes, const casa::uInt padding, double avgTime, 
                  const accessors::IConstDataAccessor &acc)
{
  ASKAPDEBUGASSERT(acc.nRow() == fringes.ncolumn());
  ASKAPDEBUGASSERT(acc.nChannel() * padding == fringes.nrow());
  for (casa::uInt row = 0; row < acc.nRow(); ++row) {
       
  }
}

casa::Matrix<casa::Complex> flagOutliers(const casa::Matrix<casa::Complex> &in) {
  return in;
  /*
  casa::Matrix<casa::Complex> result(in);
  for (casa::uInt row=0;row<result.nrow(); ++row) {
       for (casa::uInt col=0; col<result.ncolumn(); ++col) {
            if (casa::abs(result(row,col))>1) {
                result(row,col) = 0.;
            }
       }
  }
  return result;
  */
}

casa::Matrix<casa::Complex> replaceFlagsWithZeros(const casa::Matrix<casa::Complex> &in, const casa::Matrix<casa::Bool> &flag)
{
   ASKAPDEBUGASSERT(in.ncolumn()>0);
   ASKAPDEBUGASSERT(in.nrow()>0);
   ASKAPCHECK(in.shape() == flag.shape(), "Flag and visibility matrices should have the same shape");
   casa::Matrix<casa::Complex> result = in.copy();
   for (casa::uInt row=0;row<result.nrow(); ++row) {
        for (casa::uInt col=0; col<result.ncolumn(); ++col) {
             if (flag(row,col)) {
                 result(row,col) = casa::Complex(0.,0.);
             }
        }
   }
   return result;
}

casa::Matrix<casa::Complex> padSecond(const casa::Matrix<casa::Complex> &in, const casa::uInt factor) {
   if (factor == 1) {
       return in;
   }
   ASKAPDEBUGASSERT(factor>0);
   ASKAPDEBUGASSERT(in.ncolumn()>0);
   ASKAPDEBUGASSERT(in.nrow()>0);
   casa::Matrix<casa::Complex> result(in.nrow(), in.ncolumn()*factor,casa::Complex(0.,0.));
   const casa::uInt start = in.ncolumn()*(factor-1)/2;
   result(casa::IPosition(2,0,start), casa::IPosition(2, in.nrow() - 1, start + in.ncolumn() - 1)) = in;
   return result;
}

void process(const IConstDataSource &ds, const LOFAR::ParameterSet &parset) {
  const size_t nAvg = size_t(parset.getUint32("nAvg", 1u));
  const size_t padding = size_t(parset.getUint32("padding", 1u));
  const bool doFFT = parset.getBool("dofft",false);
  
  IDataSelectorPtr sel=ds.createSelector();
  sel<<parset;
  //sel->chooseBaseline(0,1);
  //sel->chooseCrossCorrelations();
  //sel->chooseFeed(1);
  IDataConverterPtr conv=ds.createConverter();  
  conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"MHz");
  conv->setEpochFrame(casa::MEpoch(casa::Quantity(56150.0,"d"),
                      casa::MEpoch::Ref(casa::MEpoch::UTC)),"s");
  conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));                    
  casa::Matrix<casa::Complex> buf;
  double avgTime = 0.;
  size_t counter = 0;
  casa::Cube<casa::Complex> imgBuf;
  const casa::uInt maxSteps = parset.getUint32("maxcycles",2000u);
  const std::string stokesStr = parset.getString("stokes", "XX"); 
  const casa::Vector<casa::Stokes::StokesTypes> stokesVector = scimath::PolConverter::fromString(stokesStr);
  if (stokesVector.nelements() != 1) {
      ASKAPTHROW(AskapError, "Exactly one stokes parameter should be defined, you have "<<stokesStr);
  }
  const bool replaceFlags = parset.getBool("zeroflags",false);
  
  casa::uInt currentStep = 0;
  casa::Vector<casa::uInt> ant1IDs;
  casa::Vector<casa::uInt> ant2IDs;
  casa::Vector<casa::uInt> beam1IDs;
  casa::Vector<casa::Stokes::StokesTypes> stokes;
  casa::uInt polIndex = 0;
    
  for (IConstDataSharedIter it=ds.createConstIterator(sel,conv);it!=it.end();++it) {  
       if (buf.nelements() == 0) {
           buf.resize(it->nRow(),it->frequency().nelements()*padding);
           buf.set(casa::Complex(0.,0.));
           ant1IDs = it->antenna1().copy();
           ant2IDs = it->antenna2().copy();
           beam1IDs = it->feed1().copy();
           stokes = it->stokes();
           polIndex = stokes.nelements();
           for (casa::uInt pol = 0; pol<stokes.nelements();++pol) {
                if (stokes[pol] == stokesVector[0]) {
                    polIndex = pol;
                    break;
                }
           }
           ASKAPCHECK(polIndex < stokes.nelements(), "Requested stokes "<<stokesStr<<" is not found in the dataset");
           for (casa::uInt row = 0; row<it->nRow();++row) {
                std::cout<<"plane "<<row<<" corresponds to "<<ant1IDs[row]<<" - "<<ant2IDs[row]<<" baseline, beam="<<beam1IDs[row]<<std::endl;
           }
           imgBuf.resize(buf.ncolumn(),maxSteps,it->nRow());
           imgBuf.set(casa::Complex(0.,0.));
       } else { 
           ASKAPCHECK(buf.ncolumn() == padding*it->frequency().nelements(), 
                  "Number of channels seem to have been changed, previously "<<buf.ncolumn()<<" now "<<it->frequency().nelements());
           if (imgBuf.nplane() != it->nRow()) {
               std::cerr << "The number of rows in the accessor is "<<it->nRow()<<", previously "<<imgBuf.nplane()<<" - ignoring"<<std::endl;
               continue;
           }
           ASKAPCHECK(imgBuf.nplane() == it->nRow(), "The number of rows in the accessor "<<it->nRow()<<
                      " is different to the maximum number of baselines");
           ASKAPDEBUGASSERT(ant1IDs.nelements() == it->nRow());
           ASKAPDEBUGASSERT(ant2IDs.nelements() == it->nRow());
           for (casa::uInt row = 0; row<it->nRow(); ++row) {
                ASKAPCHECK(ant1IDs[row] == it->antenna1()[row], "Mismatch of antenna 1 index for row "<<row<<
                           " - got "<<it->antenna1()[row]<<" expected "<<ant1IDs[row]);
                ASKAPCHECK(ant2IDs[row] == it->antenna2()[row], "Mismatch of antenna 2 index for row "<<row<<
                           " - got "<<it->antenna2()[row]<<" expected "<<ant2IDs[row]);
                ASKAPCHECK(beam1IDs[row] == it->feed1()[row], "Mismatch of beam 1 index for row "<<row<<
                           " - got "<<it->feed1()[row]<<" expected "<<beam1IDs[row]);
           }
           ASKAPCHECK(it->stokes().nelements() == stokes.nelements(), "Polarisation properties change between different chunks of data, this is not supported");
           for (casa::uInt pol = 0; pol < stokes.nelements(); ++pol) {
                ASKAPCHECK(it->stokes()[pol] == stokes[pol], "Available polarisation products appear to have been changed, expected "<<scimath::PolConverter::toString(stokes)<<
                           " got "<<scimath::PolConverter::toString(it->stokes()))
           }
       }
       ASKAPASSERT(it->nRow() == buf.nrow());
       ASKAPASSERT(it->nChannel()*padding == buf.ncolumn());
       ASKAPASSERT(it->nPol() >= 1);
       ASKAPDEBUGASSERT(polIndex < it->nPol());
       if (replaceFlags) {
           //buf += flagOutliers(padSecond(replaceFlagsWithZeros(it->visibility().xyPlane(polIndex), it->flag().xyPlane(polIndex)),padding));
           buf += flagOutliers(padSecond(replaceFlagsWithZeros(it->visibility().yzPlane(polIndex), it->flag().yzPlane(polIndex)),padding));
       } else {
           //buf += flagOutliers(padSecond(it->visibility().xyPlane(polIndex),padding));
           buf += flagOutliers(padSecond(it->visibility().yzPlane(polIndex),padding));
       }
       avgTime += it->time();
       if (++counter == nAvg) {
           buf /= float(nAvg);
           avgTime /= float(nAvg);
           if (doFFT) {
               for (casa::uInt row = 0; row<buf.nrow(); ++row) {
                    casa::Vector<casa::Complex> curRow = buf.row(row);
                    scimath::fft(curRow, true);
               }
           }
           ASKAPCHECK(currentStep < imgBuf.ncolumn(), "Image buffer is too small (in time axis), increase maxcycles");
           imgBuf.xzPlane(currentStep++) = casa::transpose(buf);
           buf.set(casa::Complex(0.,0.));
           avgTime = 0.;
           counter = 0;
       }
       //cout<<"time: "<<it->time()<<endl;
  }
  if (counter!=0) {
      buf /= float(counter);
      avgTime /= double(counter);
      if (doFFT) {
          for (casa::uInt row = 0; row<buf.nrow(); ++row) {
               casa::Vector<casa::Complex> curRow = buf.row(row);
               scimath::fft(curRow, true);
          }
      }
      ASKAPCHECK(currentStep < imgBuf.ncolumn(), "Image buffer is too small (in time axis)");
      imgBuf.xzPlane(currentStep) = casa::transpose(buf);
  } else if (currentStep > 0) {
      --currentStep;
  }
  //std::cout<<imgBuf.shape()<<" "<<currentStep<<std::endl;
  const bool doDiff = parset.getBool("makediff", false);
  if (doDiff) {
      std::cerr<<"Calculating difference between the data on adjacent integration cycles"<<std::endl;
      ASKAPCHECK(currentStep >= 2, "Need at least two integrations to compute the difference");
      casa::Matrix<casa::Float> wrapCompensation(imgBuf.nrow(), imgBuf.nplane(),0.);
      const float threshold = 3. * casa::C::pi / 2;

      // unwrap phase in time
      for (casa::uInt step = 1; step < currentStep; ++step) {
           casa::Matrix<casa::Complex> curStepMatr = imgBuf.xzPlane(step);
           casa::Matrix<casa::Complex> prevStepMatr = imgBuf.xzPlane(step-1);
           
           for (casa::uInt row = 0; row<prevStepMatr.nrow(); ++row) {
                for (casa::uInt col = 0; col<prevStepMatr.ncolumn(); ++col) {
                     casa::Float prevPhase = real(prevStepMatr(row,col));
                     const casa::Float curPhase = arg(curStepMatr(row,col));
                     const casa::Float prevOrigPhase = prevPhase - wrapCompensation(row,col);
                     const casa::Float diff = curPhase - prevOrigPhase;
                     if (diff >= threshold) {
                         wrapCompensation(row,col) -= 2. * casa::C::pi;
                     } else if (diff <= - threshold) {
                         wrapCompensation(row,col) += 2. * casa::C::pi;
                     }
                     //curStepMatr(row,col) = casa::polar(1.f,curPhase+wrapCompensation(row,col));
                     curStepMatr(row,col) = curPhase+wrapCompensation(row,col);
                }
           }
      }
      
      const float phaseRateUnit = 2. * casa::C::pi / 268435456. / 54e-6; // unit used in fringe rotator
      const float inttime = 4.97664 * nAvg; // integration time
      // calculate the difference
      for (casa::uInt step = 1; step < currentStep; ++step) {
           casa::Matrix<casa::Complex> curStepMatr = imgBuf.xzPlane(step);
           casa::Matrix<casa::Complex> prevStepMatr = imgBuf.xzPlane(step-1);
           for (casa::uInt row = 0; row<prevStepMatr.nrow(); ++row) {
                for (casa::uInt col = 0; col<prevStepMatr.ncolumn(); ++col) {
                     casa::Float prevPhase = real(prevStepMatr(row,col));
                     const casa::Float curPhase = real(curStepMatr(row,col));
                     prevPhase -= curPhase;
                     //prevStepMatr(row,col) = casa::polar(1.f, prevPhase);
                     prevStepMatr(row,col) = prevPhase / phaseRateUnit / inttime;
                }
           }
      }
      

      --currentStep;
  }

  const std::string what2export = parset.getString("datatype","amplitude"); 
  const std::string casaImg = "result.img";
  if (what2export == "amplitude") {
      scimath::saveAsCasaImage(casaImg, casa::amplitude(imgBuf(casa::IPosition(3,0,0,0),
                 casa::IPosition(3,imgBuf.nrow()-1,currentStep,imgBuf.nplane()-1))));
  } else if (what2export == "phase") {
      scimath::saveAsCasaImage(casaImg, casa::phase(imgBuf(casa::IPosition(3,0,0,0),
                 casa::IPosition(3,imgBuf.nrow()-1,currentStep,imgBuf.nplane()-1))));
  } else if (what2export == "real") {
      scimath::saveAsCasaImage(casaImg, casa::real(imgBuf(casa::IPosition(3,0,0,0),
                 casa::IPosition(3,imgBuf.nrow()-1,currentStep,imgBuf.nplane()-1))));
  } else if (what2export == "imag") {
      scimath::saveAsCasaImage(casaImg, casa::imag(imgBuf(casa::IPosition(3,0,0,0),
                 casa::IPosition(3,imgBuf.nrow()-1,currentStep,imgBuf.nplane()-1))));
  } else {
      ASKAPTHROW(AskapError,"Unknown datatype requested: "<<what2export<<", only amplitude and phase are supported");
  }
  {
    std::cerr<<"Written CASA image: "<<casaImg<<std::endl;
    casa::PagedImage<casa::Float> img(casaImg);
    casa::String error;
    if (!casa::ImageFITSConverter::ImageToFITS(error, img, casa::String("result.fits"),64u,casa::False, casa::False, -32,1.0,-1.0,casa::True,casa::False,"extractdata.cc")) {
        std::cerr<<"Error converting CASA image into FITS: "<<error<<std::endl;
    } else {
        std::cerr<<"Successfully written FITS image: result.fits"<<std::endl;
    }
  }
  
  
  // exporting first (or given) row into a dat file
  if ((currentStep>0) || (counter!=0)) {
      const casa::uInt row2export = parset.getUint32("row2export",0u); // choose row to export first row
      ASKAPCHECK((row2export < currentStep) && (row2export < imgBuf.ncolumn()), 
                 "Row "<<row2export<<" selected for exporting into a dat file does not exist");
      std::ofstream os("fringe.dat");
      for (casa::uInt chan=0; chan<imgBuf.nrow(); ++chan) {
           os<<chan<<" ";
           for (casa::uInt baseline_beam = 0; baseline_beam < imgBuf.nrow(); ++baseline_beam) {
                os<<" "<<casa::abs(imgBuf(casa::IPosition(3,chan,row2export,baseline_beam)))<<" "<<casa::arg(imgBuf(casa::IPosition(3,chan,row2export,baseline_beam)))*180./casa::C::pi;
           }
           os<<std::endl;
      }     
  }
  
}


int main(int argc, char **argv) {
  try {
     
     casa::Timer timer;

     timer.mark();
     
     utils::Parser parser; // a command line parser
     
     // command line parameter
     utils::FlaggedParameter<std::string> parsetPar("-c", "");
     // this parameter is optional
     parser.add(parsetPar, utils::Parser::return_default);
     
     utils::GenericParameter<std::string> msNamePar("");
     
     parser.add(msNamePar, utils::Parser::return_default);
     
     parser.process(argc, argv);

     const LOFAR::ParameterSet parset = (parsetPar.defined() ? LOFAR::ParameterSet(parsetPar) : LOFAR::ParameterSet());

     std::string msName = parset.getString("dataset","");
     if (msName == "") {
         msName = msNamePar;
     }
     ASKAPCHECK(msName != "", "Measurement set should be defined");
          
     ASKAPCHECK(parset.isDefined("dataset") == !msNamePar.defined(), 
        "You can only define the measurement set in one place, either in parset or in command line");      
          
     std::cerr<<"Processing measurement set "<<msName<<std::endl;
          
     TableDataSource ds(msName,TableDataSource::MEMORY_BUFFERS);     
     
     std::cerr<<"Initialization: "<<timer.real()<<std::endl;
     timer.mark();     
     
     process(ds, parset);
     
     std::cerr<<"Job: "<<timer.real()<<std::endl;     
  }  catch (const utils::XParser &ex) {
        ASKAPLOG_FATAL_STR(logger, "Command line parser error, wrong arguments " << argv[0]);
        ASKAPLOG_FATAL_STR(logger, "Usage: " << argv[0] << " [-c parsetFile] [msName]");
        return 1;
  } 
  catch(const AskapError &ce) {
     cerr<<"AskapError has been caught. "<<ce.what()<<endl;
     return -1;
  }
  catch(const std::exception &ex) {
     cerr<<"std::exception has been caught. "<<ex.what()<<endl;
     return -1;
  }
  catch(...) {
     cerr<<"An unexpected exception has been caught"<<endl;
     return -1;
  }
  return 0;
}
