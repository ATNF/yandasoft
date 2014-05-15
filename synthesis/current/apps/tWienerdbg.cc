/// @file 
/// This is a test file intended to study Wiener filter preconditioning
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
#include <measurementequation/RobustPreconditioner.h>
#include <measurementequation/WienerPreconditioner.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <measurementequation/GaussianNoiseME.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <utils/ImageUtils.h>

#include <askapparallel/AskapParallel.h>

ASKAP_LOGGER(logger, ".tpreconditioner");


using namespace askap;
using namespace askap::synthesis;

int main(int argc, char **argv) {
  try {
     casa::Timer timer;

     timer.mark();
     // Initialize MPI (also succeeds if no MPI available).
     askap::askapparallel::AskapParallel ap(argc, (const char **&)argv);

     // Ensure that CASA log messages are captured
     casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
     casa::LogSink::globalSink(globalSink);


     LOFAR::ParameterSet parset("test.in");
     SynthesisParamsHelper::setUpImageHandler(parset);
     casa::Array<float> psf = SynthesisParamsHelper::imageHandler().read("picmf.psf");
     casa::Array<float> img = SynthesisParamsHelper::imageHandler().read("testimg.img");
     ASKAPASSERT(psf.shape().nonDegenerate() == img.shape().nonDegenerate());
     
          
     std::cerr<<"Image initialization: "<<timer.real()<<std::endl;
     timer.mark();
     
     boost::shared_ptr<WienerPreconditioner> wp = WienerPreconditioner::createPreconditioner(parset);
     ASKAPASSERT(wp);
     
     std::cerr<<"Initialization of preconditioner: "<<timer.real()<<std::endl;            
     timer.mark();     
     
     wp->doPreconditioning(psf,img);
     
     std::cerr<<"Preconditioning: time="<<timer.real()<<std::endl;
     
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
