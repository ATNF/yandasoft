/// @file cubemerge.cc
///
/// @brief
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

// System includes
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ASKAPsoft includes
#include <casa/Arrays/IPosition.h>
#include <CommandLineParser.h>
#include <askap/AskapError.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/Coordinate.h>
#include <coordinates/Coordinates/LinearCoordinate.h>
#include <casa/Quanta/MVDirection.h>
#include <imageaccess/CasaImageAccess.h>
#include <images/Images/PagedImage.h>

// Using
using namespace askap;
using namespace std;

// Main function
int main(int argc, const char** argv) { 
    try {
        if (argc < 3) {
            throw cmdlineparser::XParser();
        }
        cmdlineparser::Parser parser; // a command line parser
        // command line parameter
        std::vector<cmdlineparser::GenericParameter<std::string> > inputParameters(argc-2);     
        for (std::vector<cmdlineparser::GenericParameter<std::string> >::iterator it = inputParameters.begin();
                it!= inputParameters.end(); ++it) {
            parser.add(*it);
        }
        cmdlineparser::GenericParameter<std::string> outfile;
        parser.add(outfile);

        parser.process(argc, argv);

        std::vector<std::string> inputFiles(inputParameters.size());
        for (size_t i=0;i<inputFiles.size();++i) {
            inputFiles[i] = inputParameters[i].getValue();
            std::cout<<"Input image "<<i<<" is "<<inputFiles[i]<<std::endl;
        }
        std::cout<<"Output will be stored to "<<outfile.getValue()<<std::endl;
        ASKAPCHECK(inputFiles.size()>0, "At least one input image should be defined");

        accessors::CasaImageAccess ia;
        const casa::IPosition shape = ia.shape(inputFiles[0]);
        ASKAPCHECK(shape.nelements()>=2,"Work with at least 2D images!");

        casa::IPosition newShape(shape.nelements()+1);
        for (int i = 0; i<int(shape.nelements()); ++i) {
            newShape[i] = shape[i];
        }
        newShape[shape.nelements()] = int(inputFiles.size());

        casa::CoordinateSystem csys = ia.coordSys(inputFiles[0]);     
        csys.addCoordinate(casa::LinearCoordinate(1));

        casa::PagedImage<float> outimg(casa::TiledShape(newShape), csys, outfile.getValue());
        for (size_t i=0; i<inputFiles.size(); ++i) {
            const casa::Array<float> buf = ia.read(inputFiles[i]);
            ASKAPCHECK(buf.shape().nonDegenerate() == shape.nonDegenerate(), "Image "<<inputFiles[i]<<
                    " has "<<buf.shape()<<" shape which is different from the shape of the first image "<<shape);

#ifdef _OPENMP
            #pragma omp parallel sections
            {
                #pragma omp section
#endif
                {
                    const float peak = casa::abs(casa::max(buf));
                    const float sum = casa::abs(casa::sum(buf));
                    std::cout<<"Image "<<inputFiles[i]<<" has a peak of "<<peak<<", sum of "<<sum << std::endl;
                }

#ifdef _OPENMP
                #pragma omp section
#endif
                { 
                    casa::IPosition where(newShape.nelements(),0);
                    where[shape.nelements()] = int(i);
                    outimg.putSlice(buf, where);
                }
#ifdef _OPENMP
            }
#endif
        }         
    }
    ///==============================================================================
    catch (const cmdlineparser::XParser &ex) {
        std::cerr << "Usage: " << argv[0] << " input_cube1 [input_cube2 ... input_cubeLast] output_image"
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
