/// @file 
/// This is a test file intended to study timing/preformance of preconditioning 
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

#include <iostream>
#include <stdexcept>
#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <casa/Logging/LogIO.h>
#include <askap/Log4cxxLogSink.h>
#include <casa/OS/Timer.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/IPosition.h>
#include <boost/shared_ptr.hpp>

#include <measurementequation/IImagePreconditioner.h>
#include <measurementequation/WienerPreconditioner.h>
#include <measurementequation/RobustPreconditioner.h>
#include <measurementequation/GaussianNoiseME.h>
#include <utils/ImageUtils.h>

#include <askapparallel/AskapParallel.h>

ASKAP_LOGGER(logger, ".tpreconditioner");


using namespace askap;
using namespace askap::synthesis;

/// @brief class for random number generation
/// @details We could have used casa stuff directly
struct RandomGenerator : public GaussianNoiseME 
{
  explicit RandomGenerator(double variance, askap::askapparallel::AskapParallel& comms) :
      GaussianNoiseME(variance, casa::Int(time(0)), 
           casa::Int(comms.rank())) {}
  using GaussianNoiseME::getRandomComplexNumber;
  float operator()() const { return casa::real(getRandomComplexNumber());}
};

/// @brief fill given array with some rather arbitrary values
void fillArray(casa::Array<float> &in, const RandomGenerator &rg)
{
  // first fill array with the noise
  casa::Vector<float> flattened(in.reform(casa::IPosition(1,in.nelements())));
  for (size_t i=0; i<flattened.nelements(); ++i) {
       flattened[i] = rg();
  }
  // add some pattern
  casa::IPosition index(in.shape().nelements(),0);
  ASKAPASSERT(index.nelements()>=2);
  ASKAPASSERT(in.shape()[0]>1);
  ASKAPASSERT(in.shape()[1]>1);

  // central peak
  index[0]=in.shape()[0]/2;
  index[1]=in.shape()[0]/2;
  in(index) = 1.;
  // a ring
  const float radius = 30.;
  for (size_t i=0; i<1000; ++i) {
     const float angle = float(i)/500.*casa::C::pi;
     const float dx = radius*cos(angle);
     const float dy = -radius*sin(angle);
     index[0] = in.shape()[0]/2 + int(dx);
     index[1] = in.shape()[1]/2 + int(dy);
     if ((index[0]>=0) && (index[0]<in.shape()[0]) && (index[1]>=0) && (index[1]<in.shape()[1])) {
         in(index) += 0.1;
     }
  }  
}


int main(int argc, char **argv) {
  try {
     casa::Timer timer;

     timer.mark();
     // Initialize MPI (also succeeds if no MPI available).
     askap::askapparallel::AskapParallel ap(argc, (const char **&)argv);

     // Ensure that CASA log messages are captured
     casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
     casa::LogSink::globalSink(globalSink);

     RandomGenerator rg(0.01, ap);
     // hard coded parameters of the test
     const casa::Int size = 1024;
     const size_t numberOfRuns = 5;
     //
     const casa::IPosition shape(2,size,size);
     
     casa::Array<float> psf(shape);
     fillArray(psf,rg);
     casa::Array<float> img(shape);
     fillArray(img,rg);
          
     std::cerr<<"Image initialization: "<<timer.real()<<std::endl;
     timer.mark();
     
     const float noisepower = 100.;
     WienerPreconditioner wp(noisepower,false);
     //WienerPreconditioner wp(-1);
     //RobustPreconditioner wp(-1);
     
     std::cerr<<"Initialization of preconditioner: "<<timer.real()<<std::endl;            
     timer.mark();     
     
     for (size_t run=0; run<numberOfRuns; ++run) {
          wp.doPreconditioning(psf,img);
     }
     
     std::cerr<<"Preconditioning <"<<numberOfRuns<<" run(s)>: "<<timer.real()<<std::endl;
     
     timer.mark();     
     scimath::saveAsCasaImage("outpsf.casa",psf);
     scimath::saveAsCasaImage("outimg.casa",img);     
     std::cerr<<"Storing results: "<<timer.real()<<std::endl;     
     // just to keep it active
     ap.isParallel();
  }
  catch(const AskapError &ce) {
     std::cerr<<"AskapError has been caught. "<<ce.what()<<std::endl;
     return -1;
  }
  catch(const std::exception &ex) {
     std::cerr<<"std::exception has been caught. "<<ex.what()<<std::endl;
     return -1;
  }
  catch(...) {
     std::cerr<<"An unexpected exception has been caught"<<std::endl;
     return -1;
  }
  return 0;
}
