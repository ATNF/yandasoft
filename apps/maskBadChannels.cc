/// @file
///
/// Demo program to mask out problematic channels of a cube based on noise statistics
///
/// @copyright (c) 2019 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 3 of the License,
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
/// @author Matthew Whiting <Matthew.Whiting@csiro.au>
///


// Package level header file
#include "askap_synthesis.h"

// System includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdio.h>

// ASKAPsoft includes
#include <askap/Application.h>
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/StatReporter.h>
#include <askapparallel/AskapParallel.h>

#include <imageaccess/ImageAccessFactory.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicMath/Math.h>

ASKAP_LOGGER(logger, ".maskChannels");

using namespace askap;

class MaskChanApp : public askap::Application
{
public:
    virtual int run(int argc, char* argv[])
        {
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;
                LOFAR::ParameterSet subset(config().makeSubset("MaskChannels."));

                std::string statsfile = subset.getString("statsFile","");
                ASKAPLOG_INFO_STR(logger, "Reading stats file " << statsfile);

                std::string image = subset.getString("image","");
                bool editImage = subset.getBool("editImage", false);
                bool editStats = subset.getBool("editStats", false);
                float threshold = subset.getFloat("threshold", 5.);
                
                std::ifstream fin(statsfile.c_str());
                std::string line, name;
                unsigned int size=0;
                
                while (getline(fin, line),
                       !fin.eof()) {
                    if (line[0] != '#') {
                        size++;
                    }
                }
                fin.close();

                ASKAPLOG_INFO_STR(logger, "Will get " << size << " channels of statistics");

                casa::Vector<unsigned int> chan(size);
                casa::Vector<double> freq(size),mean(size),std(size),median(size),madfm(size),onepc(size),minval(size),maxval(size);
                unsigned int nrow=0;

                fin.open(statsfile.c_str());
                std::string title,units;
                getline(fin,title);
                getline(fin,units);
                while (getline(fin, line),
                       !fin.eof()) {
                    if (line[0] != '#') {
                        std::stringstream ss(line);
                        ss >> chan[nrow] >> freq[nrow] >> mean[nrow] >> std[nrow] >> median[nrow] >> madfm[nrow] >> onepc[nrow] >> minval[nrow] >> maxval[nrow];
                        nrow++;
                    }
                }

                ASKAPLOG_INFO_STR(logger, "Successfully read " << nrow << " rows");
                
                double med=casa::median(std);
                double mad=casa::median(abs(std-med));
                double stdThreshold = med + threshold * mad;
                ASKAPLOG_INFO_STR(logger, "Threshold for stdev = " << stdThreshold << " median="<< med << ", madfm="<<mad);
                casa::Vector<double> tmp=(maxval-minval)/madfm;
                std::vector<double> newtmp;
                for(size_t i=0;i<tmp.size();i++){
                    if(! casa::isNaN(tmp[i])){
                        newtmp.push_back(tmp[i]);
                    }
                }
                casa::Vector<double> rangeOverNoise(newtmp);
                med=casa::median(rangeOverNoise);
                mad=casa::median(abs(rangeOverNoise-med));
                double rangeThreshold = med + threshold * mad;
                ASKAPLOG_INFO_STR(logger, "Threshold for (max-min)/median = " << rangeThreshold << " median="<< med << ", madfm="<<mad);

                boost::shared_ptr<accessors::IImageAccess> iacc = accessors::imageAccessFactory(subset);
                casa::IPosition shape = iacc->shape(image);
                ASKAPLOG_INFO_STR(logger, "Shape of input image = " << shape);
                casa::CoordinateSystem coo = iacc->coordSys(image);
                int specAxis = coo.spectralAxisNumber();
                ASKAPLOG_INFO_STR(logger, "Spectral axis = " << specAxis);
                casa::IPosition chanShape=shape;
                ASKAPLOG_INFO_STR(logger, "Shape of a single channel = " << chanShape);
                chanShape[specAxis]=1;
                casa::Vector<float> datavec(chanShape.product(),0.);
                for(size_t i=0;i<datavec.size();i++) {
                    casa::setNaN(datavec[i]);
                }
                casa::Array<float> data(chanShape,datavec.data());
                                        
                unsigned int numBad=0;
                for(unsigned int i=0;i<size;i++){
                    ASKAPLOG_INFO_STR(logger, "Channel " << i << " has values " << std[i] << " and " << rangeOverNoise[i]);
                    // if( (std[i] > stdThreshold) && (rangeOverNoise[i] > rangeThreshold) ) {
                    if( (std[i] > stdThreshold) ) {
                        ASKAPLOG_INFO_STR(logger, "Bad Channel #" << i << ": std = "<<std[i]);
                        numBad++;
                        if (editImage) {
                            casa::IPosition loc(shape.size(), 0);
                            loc[specAxis] = i;
                            ASKAPLOG_INFO_STR(logger, "Writing data of shape " << data.shape() << " to position " << loc);
                            iacc->write(image,data,loc);
                            std::stringstream ss;
                            ss << "Masked channel " << i << " : std.dev="<<std[i]<<", range/median="<<rangeOverNoise[i];
                            //iacc->addHistory(image,ss.str());
                        }
                        if (editStats) {
                            mean[i] = std[i] = median[i] = madfm[i] = onepc[i] = minval[i] = maxval[i] = 0.0;
                        }
                    }
                }
                ASKAPLOG_INFO_STR(logger, "Detected " << numBad << " bad channels to be masked");

                if (editStats) {
                    std::string outputStats = subset.getString("outputStats","");
                    ASKAPLOG_INFO_STR(logger, "Writing new stats file to " << outputStats);
                    // std::ofstream fout(outputStats.c_str());
                    // fout << title << "\n" << units<<"\n";
                    FILE *fout;
                    fout = fopen(outputStats.c_str(),"w");
                    fprintf(fout, "%s\n%s\n",title.c_str(),units.c_str());
                    for(unsigned int i=0;i<size;i++){
                        fprintf(fout, "%8d %15.6f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
                                i,freq[i],mean[i],std[i],median[i],madfm[i],onepc[i],minval[i],maxval[i]);
                        // fout << std::setw(8);
                        // fout << i;
                        // fout << std::fixed << std::right;
                        // fout << std::setw(16) << std::setprecision(6) << freq[i];
                        // fout << std::setw(11) << std::setprecision(3) << mean[i];
                        // fout << std::setw(11) << std::setprecision(3) << std[i];
                        // fout << std::setw(11) << std::setprecision(3) << median[i];
                        // fout << std::setw(11) << std::setprecision(3) << madfm[i];
                        // fout << std::setw(11) << std::setprecision(3) << onepc[i];
                        // fout << std::setw(11) << std::setprecision(3) << minval[i];
                        // fout << std::setw(11) << std::setprecision(3) << maxval[i] << "\n"; 
                    }
//                    fout.close();
                    fclose(fout);
                }
                
                stats.logSummary();
            } catch (const askap::AskapError& x) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
                std::cerr << "Askap error in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            } catch (const std::exception& x) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            }

            return 0;
        }
};




int main(int argc, char *argv[])
{
    MaskChanApp app;
    return app.main(argc, argv);
}

