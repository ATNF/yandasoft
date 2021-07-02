/// @file
/// This is a test file intended to test the beam fitter code
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
#include <casacore/casa/Logging/LogIO.h>
#include <askap/Log4cxxLogSink.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <boost/shared_ptr.hpp>

#include <measurementequation/SynthesisParamsHelper.h>
#include <utils/ImageUtils.h>

#include <askapparallel/AskapParallel.h>

ASKAP_LOGGER(logger, ".tbeamfitter");


using namespace askap;
using namespace askap::synthesis;

int main(int argc, char **argv) {
  try {
     //casa::Timer timer;

     //timer.mark();
     // Initialize MPI (also succeeds if no MPI available).
     askap::askapparallel::AskapParallel ap(argc, (const char **&)argv);

     // Ensure that CASA log messages are captured
     casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
     casa::LogSink::globalSink(globalSink);

     // test.in contains:
     // imagetype                               = fits
     // Images.Names                            = [image.psf.test] # ->image.psf.test.fits
     // cutoff                                  = 0.5

     LOFAR::ParameterSet parset("test.in");

     // set up for reading fits images
     SynthesisParamsHelper::setUpImageHandler(parset);

     // load image into params
     askap::scimath::Params::ShPtr params (new scimath::Params(true));
     SynthesisParamsHelper::loadImages(params,parset.makeSubset("Images."));

     // fit beam
     ASKAPLOG_INFO_STR(logger, "Fitting restoring beam");
     double cutoff = parset.getDouble("cutoff",0.5);
     casa::Vector<casa::Quantum<double> > restoringBeam(3);
     try {
         restoringBeam =
            SynthesisParamsHelper::fitBeam(*params, cutoff,"image.psf.test");
     } catch( const AskapError& x) {
         ASKAPLOG_INFO_STR(logger, "Fit failed, retrying once with 20% higher cutoff");
         restoringBeam =
            SynthesisParamsHelper::fitBeam(*params, cutoff*1.2,"image.psf.test");
     }
     ASKAPLOG_INFO_STR(logger, "Fitted beam: " << restoringBeam[0].getValue("arcsec") <<
        " x "<<restoringBeam[1].getValue("arcsec") <<" arcsec at position angle "<<
        restoringBeam[2].getValue("deg")<<" deg");

     //timer.mark();
     //std::cerr<<"Storing results: "<<timer.real()<<std::endl;

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
