/// @file imgscale.cc
///
/// @copyright (c) 2014 CSIRO
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

#include <casa/Arrays/IPosition.h>
#include <images/Images/PagedImage.h>
#include <CommandLineParser.h>
#include <askap/AskapError.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/Coordinate.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Quanta/Quantum.h>
#include <measures/Measures/MDirection.h>
#include <casa/Quanta/MVDirection.h>


#include <stdexcept>
#include <iostream>

using namespace askap;
using namespace std;

void printDirection(ostream &os,const casa::MDirection &dir)  {
    double lngbuf=dir.getValue().getLong("deg").getValue();
    if (lngbuf<0) lngbuf+=360.;
    os<<(dir.getRefString()!="GALACTIC"?casa::MVAngle::Format(casa::MVAngle::TIME):
          casa::MVAngle::Format(casa::MVAngle::ANGLE))<<casa::MVAngle(casa::Quantity(lngbuf,"deg"))<<" "<<
          casa::MVAngle(dir.getValue().getLat("deg"))<<
          " ("<<dir.getRefString()<<")";
}


// Main function
int main(int argc, const char** argv) { 
  try {
     cmdlineparser::Parser parser; // a command line parser
     // command line parameter
	 cmdlineparser::GenericParameter<std::string> imgfile;
	 parser.add(imgfile);

	 // I hope const_cast is temporary here
	 parser.process(argc, const_cast<char**> (argv));
         casa::PagedImage<casa::Float> img(imgfile.getValue());
         casa::IPosition shape = img.shape();
         ASKAPASSERT(shape.nelements()>=2);
         const int xcentre = 269;
         const int ycentre = 247;
         const int boxsz = 20;
         ASKAPASSERT(xcentre+boxsz<shape[0]);
         ASKAPASSERT(ycentre+boxsz<shape[1]);
         ASKAPASSERT(xcentre-boxsz>=0);
         ASKAPASSERT(ycentre-boxsz>=0);
         for (int x=xcentre-boxsz;x<=xcentre+boxsz;++x) {
              for (int y=xcentre-boxsz;y<=xcentre+boxsz;++y) {
                   casa::IPosition pos(shape.nelements(),0);
                   pos[0]=x;
                   pos[1]=shape[1]-y-1;
                   const float val = img(pos);
                   img.putAt(val*10., pos);
              }
         }
  }
  ///==============================================================================
  catch (const cmdlineparser::XParser &ex) {
	 std::cerr << "Usage: " << argv[0] << " imagefile"
			<< std::endl;
  }

  catch (const askap::AskapError& x) {
     std::cerr << "Askap error in " << argv[0] << ": " << x.what()
        << std::endl;
     exit(1);
  } 
  catch (const std::exception& x) {
	 std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what()
			<< std::endl;
	 exit(1);
  }
  exit(0);  
}
