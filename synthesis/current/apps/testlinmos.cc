/// @file
///
/// @brief helper application to support MRO experiments
/// @details In future, linmos would replace this utility.
/// 
///
/// @copyright (c) 2012 CSIRO
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

// Package level header file
#include "askap_synthesis.h"

// System includes
#include <sstream>
#include <boost/shared_ptr.hpp>

// casa includes
#include "casa/OS/Timer.h"
#include <casa/Arrays/Array.h>

// other 3rd party
#include <Common/ParameterSet.h>


// ASKAPsoft includes
#include <askap/AskapError.h>
#include <askap/AskapLogging.h>
#include <askap/Log4cxxLogSink.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <Common/ParameterSet.h>
#include <fitting/Params.h>
#include <askap/AskapUtil.h>
#include <imageaccess/IImageAccess.h>
#include <boost/shared_ptr.hpp>

ASKAP_LOGGER(logger, ".testlinmos");

using namespace askap;
using namespace askap::synthesis;

casa::MVDirection convertDir(const std::string &ra, const std::string &dec) {
  casa::Quantity tmpra,tmpdec;
  casa::Quantity::read(tmpra, ra);
  casa::Quantity::read(tmpdec,dec);
  return casa::MVDirection(tmpra,tmpdec);  
}

// helper method to load beam offsets from the parset file sharing the same format
// as csimulator feed definition. Due to lack of proper phase tracking per beam in the first experiments
// we have to set all beam offsets in the MS to zero. Therefore, this information has to be supplied 
// by other means.
casa::Vector<casa::MVDirection> loadBeamOffsets(const std::string &fname, const casa::MVDirection &centre)
{
   ASKAPLOG_INFO_STR(logger,  "Loading beam offsets from "<<fname);
   LOFAR::ParameterSet parset(fname);
   const std::vector<std::string> beamNames(parset.getStringVector("feeds.names"));
   const casa::uInt nBeams = beamNames.size();
   ASKAPLOG_INFO_STR(logger,  "File contains description for "<<nBeams<<" beams");
   ASKAPCHECK(nBeams > 0, "No beams specified");
   casa::Vector<casa::MVDirection> result(nBeams, centre);
   double spacing = 1.;
   if (parset.isDefined("feeds.spacing")) {
       casa::Quantity qspacing = asQuantity(parset.getString("feeds.spacing"));
       spacing = qspacing.getValue("rad");
       ASKAPLOG_INFO_STR(logger, "Scaling beam offsets by " << qspacing);       
   } 
   for (casa::uInt beam = 0; beam < nBeams; ++beam) {
        const std::string parName = "feeds." + beamNames[beam];
        const std::vector<double> xy(parset.getDoubleVector(parName));
        ASKAPCHECK(xy.size() == 2, "Expect two elements for each offset");
        result[beam].shift(xy[0]*spacing, xy[1]*spacing, casa::True);
   }
   return result;
}

void process() {
  
  // centres for each beam
  const casa::Vector<casa::MVDirection> centres = 
     loadBeamOffsets("beamoffsets.in", convertDir("15h56m58.871","-79.14.04.28"));

     // first attempt at MRO
     //loadBeamOffsets("beamoffsets.in", convertDir("15h39m16.97","-78.42.06.17"));

  /*
  casa::Vector<casa::MVDirection> centres(4);

  centres[0] = convertDir("13:26:51.70","-42.45.38.90");
  centres[1] = convertDir("13:24:08.26","-42.45.38.90");
  centres[2] = convertDir("13:26:52.37","-43.15.38.87");
  centres[3] = convertDir("13:24:07.59","-43.15.38.87");
  */
  /*
  centres[0] = convertDir("15:56:58.87","-79.14.04.28");
  centres[1] = convertDir("16:17:49.28","-77.17.18.49");
  centres[2] = convertDir("16:08:15.09","-78.16.24.53");
  //centres[3] = convertDir("15:55:21.65","-79.40.36.30");
  */
  
  for (casa::uInt beam = 0; beam<centres.nelements(); ++beam) {
       ASKAPLOG_INFO_STR(logger,  "Beam "<<beam + 1<<" pointing centre is "<<printDirection(centres[beam]));
  }
  
  const double cutoff = 2e-1;
  const double fwhm = 1.22*3e8/928e6/12;
  
  accessors::IImageAccess& iacc = SynthesisParamsHelper::imageHandler();
  const casa::IPosition shape = iacc.shape("beam0.img");
  const casa::Vector<casa::Quantum<double> > beamInfo = iacc.beamInfo("beam0.img");
  ASKAPCHECK(beamInfo.nelements()>=3, "beamInfo is supposed to have at least 3 elements");
  const casa::CoordinateSystem cs = iacc.coordSys("beam0.img");
  casa::Vector<casa::Array<float> > pixels(centres.nelements());
  ASKAPASSERT(pixels.nelements()>=1);
  for (size_t beam = 0; beam < pixels.nelements(); ++beam) {
       pixels[beam] = iacc.read("beam" + utility::toString<size_t>(beam)+".img");
       ASKAPASSERT(pixels[beam].shape() == pixels[0].shape());
       ASKAPASSERT(pixels[beam].shape().nonDegenerate().nelements() == 2);
  }
  casa::IPosition curpos(pixels[0].shape());
  for (casa::uInt dim=0; dim<curpos.nelements(); ++dim) {
       curpos[dim] = 0;
  }
  ASKAPASSERT(curpos.nelements()>=2);
  const casa::DirectionCoordinate &dc = cs.directionCoordinate(0);
  casa::Vector<casa::Double> pixel(2,0.);
  casa::MVDirection world;
  for (int x=0; x<pixels[0].shape()[0];++x) {
       for (int y=0; y<pixels[0].shape()[1];++y) {
            pixel[0] = double(x);
            pixel[1] = double(y);
            dc.toWorld(world,pixel);
            const double offsetBeam0 = world.separation(centres[0]);
            const double wt0 = exp(-offsetBeam0*offsetBeam0*4.*log(2.)/fwhm/fwhm);
            double sumsqwt = wt0 * wt0;
            curpos[0] = x;
            curpos[1] = y;
            double resflux = pixels[0](curpos) * wt0;
            for (size_t beam=1; beam<pixels.nelements(); ++beam) {
                 /*
                 if (beam==0) {
                     continue;
                 }
                 */
                 const double offsetThisBeam = world.separation(centres[beam]);
                 const double thisWt = exp(-offsetThisBeam*offsetThisBeam*4.*log(2.)/fwhm/fwhm);
                 sumsqwt += thisWt*thisWt; // * (beam == 3 ? 0.21 : 1.);
                 resflux += pixels[beam](curpos)*thisWt;
            }
            if (sqrt(sumsqwt)<cutoff) {
                resflux = 0.;
            } else {
                resflux /= sqrt(sumsqwt);
            }
            // for weight
            pixels[0](curpos) = sqrt(sumsqwt); //resflux;
            // for image
            //pixels[0](curpos) = resflux;
       }
  }
 
  // write result
  iacc.create("image.result", shape, cs);
  iacc.write("image.result",pixels[0]);
  iacc.setBeamInfo("image.result",beamInfo[0].getValue("rad"), beamInfo[1].getValue("rad"),
                   beamInfo[2].getValue("rad"));
}

// Main function
int main(int argc, const char** argv)
{
    // Now we have to initialize the logger before we use it
    // If a log configuration exists in the current directory then
    // use it, otherwise try to use the programs default one
    std::ifstream config("askap.log_cfg", std::ifstream::in);
    if (config) {
        ASKAPLOG_INIT("askap.log_cfg");
    } else {
        std::ostringstream ss;
        ss << argv[0] << ".log_cfg";
        ASKAPLOG_INIT(ss.str().c_str());
    }

    // Ensure that CASA log messages are captured
    casa::LogSinkInterface* globalSink = new Log4cxxLogSink();
    casa::LogSink::globalSink(globalSink);

    try {
        casa::Timer timer;
        timer.mark();

        SynthesisParamsHelper::setUpImageHandler(LOFAR::ParameterSet());
        
        process();       
        
        ASKAPLOG_INFO_STR(logger,  "Total times - user:   " << timer.user() << " system: " << timer.system()
                << " real:   " << timer.real());
        ///==============================================================================
    } catch (const askap::AskapError& x) {
        ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
        return 1;
    } catch (const std::exception& x) {
        ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
        return 1;
    }

    return 0;
}

        
