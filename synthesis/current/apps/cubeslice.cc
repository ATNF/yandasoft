#include <casa/Arrays/IPosition.h>
#include <CommandLineParser.h>
#include <askap/AskapError.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/Coordinate.h>
#include <casa/Quanta/MVDirection.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/PagedImage.h>
#include <images/Images/SubImage.h>



#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

using namespace askap;
using namespace std;

// Main function
int main(int argc, const char** argv) { 
  try {
     cmdlineparser::Parser parser; // a command line parser
     // command line parameter
     cmdlineparser::FlaggedParameter<int> nchan("-n",1);
     cmdlineparser::GenericParameter<std::string> imgfile;
     cmdlineparser::GenericParameter<std::string> outfile;
// 	 cmdlineparser::GenericParameter<int> startchan;
     cmdlineparser::FlaggedParameter<int> startchan("-c",-1);
     cmdlineparser::FlaggedParameter<std::string> subsection("-s","[]");
	 
     parser.add(nchan,cmdlineparser::Parser::return_default); // optional
     parser.add(startchan,cmdlineparser::Parser::return_default);
     parser.add(subsection,cmdlineparser::Parser::return_default);
     parser.add(imgfile);
     parser.add(outfile);
     
     
     parser.process(argc, argv);
     casa::PagedImage<casa::Float> img(imgfile.getValue());
     ASKAPCHECK(img.ok(),"Error loading "<<imgfile.getValue());
     ASKAPCHECK(img.shape().nelements()>=3,"Work with at least 3D cubes!");
     
     const casa::IPosition shape = img.shape();
     
     casa::IPosition blc(shape.nelements(),0);
     casa::IPosition trc(shape);

     std::string section = subsection.getValue();

     if(startchan>=0){

       blc[3] = startchan;
       trc[3] = nchan;
       ASKAPCHECK(blc[3]>=0 && blc[3]<shape[3], "Start channel is outside the number of channels or negative, shape: "<<shape);
       ASKAPCHECK(trc[3]+blc[3]<=shape[3], "Subcube extends beyond the original cube, shape:"<<shape);
     
     }
     else if(section!="[]"){
       
       std::vector<std::string> ranges(shape.nelements()); 
       std::string temp;
       std::stringstream ss;
       ss.str(section);
       getline(ss,temp,'[');
       for(size_t i=0;i<shape.nelements()-1;i++){
	 getline(ss,ranges[i],',');
       }
       getline(ss,ranges[shape.nelements()-1],']');

       for(size_t i=0;i<shape.nelements();i++){
	 if(ranges[i] != "*"){
	   // only need to change from default for the case of a:b
	   blc[i] = atoi( ranges[i].substr(0, ranges[i].find(":")).c_str() );
	   trc[i] = atoi( ranges[i].substr(ranges[i].find(":")+1).c_str() ) - blc[i] + 1;
	 }
       }

     }
     else{
       throw cmdlineparser::XParser();
     }

     std::cerr << "blc  = " << blc << ", trc = " << trc << "\n";

     casa::Slicer slc(blc,trc,casa::IPosition(shape.nelements(),1));
     
     casa::SubImage<casa::Float> si = casa::SubImage<casa::Float>(img,slc,casa::AxesSpecifier(casa::True));
     casa::PagedImage<casa::Float> res(si.shape(),si.coordinates(),std::string(outfile.getValue()));
     res.put(si.get());
  }
  ///==============================================================================
  catch (const cmdlineparser::XParser &ex) {
	 std::cerr << "Usage: " << argv[0] << " [-n number_of_chan] [-c start_chan] [-s subsection_string] input_cube output_image"
			<< std::endl
		   << "       One of -c or -s must be used" << std::endl
		   << "       Format for subsection_string is [*,100:200,*,*], where * means entire axis, and range is minimum to maximum pixel value" << std::endl;
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
