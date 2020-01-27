/// @file imcontsub.cc
///
/// @brief Image based Continuum subtraction
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
/// @author Mark Wieringa

// Include package level header file
#include <askap/askap_synthesis.h>


// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/Application.h>
#include <askap/StatReporter.h>
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>

// casacore includes
#include <casacore/scimath/Fitting/LinearFitSVD.h>
#include <casacore/scimath/Functionals/Polynomial.h>
#include <casacore/casa/OS/CanonicalConversion.h>

#include <fitsio.h>
// System include
//#include <string>
//#include <iostream>
//#include <stdlib.h>

// robust contsub C++ version
// Need: image accessors that can do parallel write like MPI-IO, reading/writing section of cube in any direction
// E.g., read spectrum for one or more image rows. Current fits image accessor just does one or more 2-d slices of
// contiguous data

// Need:
// 1. Fits image class that can do collective MPI-IO read/write
// 2. median/quartile determination, selection. excluding NaNs
// 3. polynomial fitting & evaluation

ASKAP_LOGGER(logger, ".imcontsub");

using namespace std;
using namespace askap;
//using namespace casacore;
//using namespace askap::synthesis;

class ImContSubApp : public askap::Application
{
public:
    virtual int run(int argc, char* argv[])
    {
        // This class must have scope outside the main try/catch block
        askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

        try {
            StatReporter stats;
            LOFAR::ParameterSet subset(config().makeSubset("imcontsub."));

            ASKAPLOG_INFO_STR(logger, "ASKAP image based continuum subtraction application " << ASKAP_PACKAGE_VERSION);

            const casa::String infile = subset.getString("inputfitscube","");
            casa::String outfile = subset.getString("outputfitscube","");
            if (outfile=="") {
                size_t pos = infile.rfind(".fits");
                outfile = infile.substr(0,pos) + ".contsub.fits";
            }
            float threshold = subset.getFloat("threshold",2.0);
            int order = subset.getInt("order",2);

            // get header and data size, get image dimensions
            casa::Int nx,ny,nz;
            casa::Long headersize;
            decode_header(infile, nx, ny, nz, headersize);

            if (comms.isMaster()) {
                ASKAPLOG_INFO_STR(logger,"In = "<<infile <<", Out = "<<
                                      outfile <<", threshold = "<<threshold << ", order = "<< order);
                ASKAPLOG_INFO_STR(logger,"Input image dimensions: "<< nx <<", "<< ny << ", "<<nz);
                copy_header(infile, outfile, headersize);
            }

            // All wait for header to be written
            //ASKAPLOG_INFO_STR(logger,"Waiting at barrier");
            comms.barrier();

            // Now process the rest of the file in parallel
            MPI_Status status;
            const int myrank = comms.rank();
            const int numprocs = comms.nProcs();
            const int nrow = ny / numprocs;
            ASKAPASSERT(ny == nrow * numprocs);
            const MPI::Offset blocksize = nx * nrow;    // number of floats
            const MPI::Offset bufsize = nx * nrow * nz;  // local number to read

            MPI_Datatype filetype;
            // # blocks, blocklength, stride (distance between blocks)
            MPI_Type_vector(nz, nx*nrow, nx*ny, MPI_FLOAT, &filetype);
            MPI_Type_commit(&filetype);
            const MPI_Offset offset = headersize + myrank * blocksize * sizeof(float);
            MPI_File fh;
            MPI_File_open(MPI_COMM_WORLD, infile.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            MPI_File_set_view(fh, offset, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
            casa::Float* buf = new casa::Float[bufsize];
            casa::Float* bufl = new casa::Float[bufsize];
            // Collective read of the whole cube
            ASKAPLOG_INFO_STR(logger,"Reading the cube");
            MPI_File_read_all(fh, buf, bufsize, MPI_FLOAT, &status);
            // Take care of endianness
            casa::CanonicalConversion::toLocal(bufl, buf, bufsize);
            MPI_File_close(&fh);

            // Process spectrum by spectrum
            ASKAPLOG_INFO_STR(logger,"Process the spectra");
            casa::Vector<casa::Float> vec(casa::IPosition(1,bufsize), bufl, casa::StorageInitPolicy::SHARE);
            for (int i = 0; i< nx*nrow; i++ ) {
                casa::Vector<casa::Float> spec(vec(casa::Slice(i,nz,blocksize)));
                process_spectrum(spec, threshold, order);
            }

            // Write results to output file
            MPI_File_open(MPI_COMM_WORLD, outfile.c_str(), MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
            MPI_File_set_view(fh, offset, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
            casa::CanonicalConversion::fromLocal(buf, bufl, bufsize);
            ASKAPLOG_INFO_STR(logger,"Writing the cube");
            MPI_File_write_all(fh, buf, bufsize, MPI_FLOAT, &status);
            MPI_File_close(&fh);
            delete [] buf;
            delete [] bufl;

            // Add fits padding to make file size multiple of 2880
            if (comms.isMaster()) {
                fits_padding(outfile);
            }

            // All wait for padding to be written
            //ASKAPLOG_INFO_STR(logger,"Waiting at barrier");
            comms.barrier();
            ASKAPLOG_INFO_STR(logger,"Done");
            // Done
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


        void decode_header(const casa::String& infile, casa::Int& nx, casa::Int& ny,
                            casa::Int& nz, casa::Long& headersize)
        {
            fitsfile *infptr, *outfptr;  // FITS file pointers
            int status = 0;  // CFITSIO status value MUST be initialized to zero!

            fits_open_file(&infptr, infile.c_str(), READONLY, &status); // open input image
            if (status) {
                fits_report_error(stderr, status); // print error message
                return;
            }
            LONGLONG headstart, datastart, dataend;
            fits_get_hduaddrll (infptr, &headstart, &datastart, &dataend, &status);
            //ASKAPLOG_INFO_STR(logger,"header starts at: "<<headstart<<" data start: "<<datastart<<" end: "<<dataend);
            headersize = datastart;
            int naxis;
            long naxes[4];
            fits_get_img_dim(infptr, &naxis, &status);  // read dimensions
            //ASKAPLOG_INFO_STR(logger,"Input image #dimensions: " << naxis);
            fits_get_img_size(infptr, 4, naxes, &status);
            nx = naxes[0];
            ny = naxes[1];
            nz = naxes[2];
            if (naxis >3) {
                //ASKAPLOG_INFO_STR(logger,"Input image dimensions: "<< naxes[0] <<", "<< naxes[1] << ", "<<
                //naxes[2] << ", " << naxes[3]);
                if (nz==1) nz = naxes[3];
            } else {
                //ASKAPLOG_INFO_STR(logger,"Input image dimensions: "<< naxes[0] <<", "<< naxes[1] << ", "<<
                //naxes[2]);
            }
            fits_close_file(infptr, &status);
        }

        void copy_header(const casa::String &infile, const casa::String& outfile, casa::Long headersize)
        {
            // create the new output file and copy header
            ASKAPLOG_INFO_STR(logger,"master creates the new output file and copies header");
            //streampos size;
            char * header;
            ifstream file (infile, ios::in|ios::binary);
            if (file.is_open())
            {
                header = new char [headersize];
                file.read (header, headersize);
                file.close();
                ofstream ofile (outfile, ios::out|ios::binary|ios::trunc);
                if (ofile.is_open()) {
                    ofile.write(header,headersize);
                    ofile.close();
                }
                delete[] header;
            }
        }

        void fits_padding(const casa::String& filename)
        {
            std::ifstream infile(filename, std::ios::binary | std::ios::ate);
            const size_t file_size = infile.tellg();
            const size_t padding = 2880 - (file_size % 2880);
            infile.close();
            if (padding != 2880) {
                std::ofstream ofile;
                ofile.open(filename, std::ios::binary | std::ios::app);
                char * buf = new char[padding];
                for (char* bufp = &buf[padding-1]; bufp >= buf; bufp--) *bufp = 0;
                ofile.write(buf,padding);
                ofile.close();
                delete [] buf;
            }
            ASKAPLOG_INFO_STR(logger,"master added "<< padding % 2880 << " bytes of FITS padding to file of size "<<file_size);
        }

        void process_spectrum(casa::Vector<casa::Float>& vec, float threshold, int order) {
            size_t n = vec.nelements();

            // first remove spectral index

            // take every 10th point
            const int inc = casa::min(n/10,10);
            const int n1 = n/inc;
            casa::Vector<casa::Float> y1(n1);
            for (int i=0; i<n1; i++) y1(i) = vec(i*inc);

            // mask out NaNs and extreme outliers
            casa::Float q5 = casa::fractile(y1,0.05f);
            casa::Float q95 = casa::fractile(y1,0.95f);
            casa::Vector<casa::Bool> mask1(n1);
            int ngood = 0;
            for (int i=0; i<n1; i++) {
                mask1(i) = !casa::isNaN(y1(i)) & (y1(i) >= q5) & (y1(i) <= q95);
                if (mask1(i)) ngood++;
            }
            //ASKAPLOG_INFO_STR(logger,"q5="<<q5<<", q95="<<q95<< ", y1="<< y1);
            //ASKAPLOG_INFO_STR(logger,"n1="<<n1<<", ngood="<<ngood);
            casa::Vector<casa::Float> x(ngood), y(ngood);
            for (int i=0, j=0; i<n1; i++) {
                if (mask1(i)) {
                    y(j) = vec(i*inc);
                    x(j++) = i*inc;
                }
            }

            casa::Float xmean = 0;
            casa::Float offset = 0;
            casa::Float slope = 0;

            if (ngood > 0) {
                xmean = mean(x);
                x -= xmean;
                // fit line to spectrum and subtract it
                offset = sum(y) / ngood;
                slope = sum(x * y) / sum(x * x);
            }

            casa::Vector<casa::Float> x_full(n), y_full(n);
            indgen(x_full);
            x_full -= xmean;
            for (int i = 0; i<n; i++) y_full(i) = vec(i) - (offset + slope * x_full(i));
            //y_full = x_full;
            //y_full *= -slope;
            //y_full -= offset;
            //y_full += vec;
            //ASKAPLOG_INFO_STR(logger," offset = "<<offset<<", slope = "<<slope);

            // now select data within threshold of median
            casa::Float q15 = casa::fractile(y_full, 0.15f);
            casa::Float q50 = casa::fractile(y_full, 0.50f);
            casa::Float iqr = (1/1.35)*2*(q50-q15); // just copying python version - iqr should really be 2*(q50-q25)
            casa::Vector<casa::Bool> mask = (y_full >= q50 - threshold * iqr ) & ( y_full <= q50 + threshold * iqr);

            // fit polynomial of given order
            //casa::LinearFitSVD<casa::Float> fitter;

            // using AutoDiff - slow
            // casa::Polynomial<casa::AutoDiff<casa::Float> > func(order+1);
            // for (int i=0; i<order+1; i++) {
            //     func.setCoefficient( i, 1.0);
            // }
            // fitter.setFunction(func);
            // ASKAPLOG_INFO_STR(logger,"Fitting spectrum");
            // casa::Vector<casa::Float> solution = fitter.fit(x_full, vec, &mask);
            // subtract continuum model from data
            // casa::Vector<casa::Float> model(n);
            // ASKAPLOG_INFO_STR(logger,"Subtracting model");
            // bool ok = fitter.residual(model, x_full, casa::True);
            // vec -= model;

            // using coefficients table - segv
            // fitter.set(order+1);
            // casa::Matrix<casa::Float> xx(n, order+1);
            // for (int i=0; i<n; i++) {
            //     xx(i,0) = 1;
            //     for (int j=1; j<order+1; j++) xx(i,j) = xx(i,j-1)*casa::Float(i+1);
            // }
            // casa::Vector<casa::Float> solution = fitter.fit(xx, vec, &mask);
            // casa::Vector<casa::Float> model(n);
            // // ASKAPLOG_INFO_STR(logger,"Subtracting model");
            // bool ok = fitter.residual(model, xx, casa::True);
            // vec -= model;

            // using LSQaips - need invert before solve
            // Note: fitter uses double internally, giving float inputs just results in internal copying
            //       and is slower
            casa::LSQaips fitter(order+1);
            casa::Vector<casa::Double> xx(order+1);
            casa::VectorSTLIterator<casa::Double> it(xx);
            for (int i=0; i<n; i++) {
                if (mask(i)) {
                    xx(0) = 1;
                    for (int j=1; j<order+1; j++) {
                        xx(j) = xx(j-1) * i;
                    }
                    //ASKAPLOG_INFO_STR(logger,"makeNorm "<<i<<" "<<vec(i));
                    fitter.makeNorm(it,1.0,casa::Double(vec(i)));
                    //fitter.makeNorm(it,1.0f,vec(i));
                }
            }
            casa::uInt nr1;
            //cout << "Invert = " << fitter.invert(nr1);
            //cout << ", rank=" << nr1 << endl;
            casa::Vector<casa::Double> solution(order+1);
            fitter.invert(nr1);
            fitter.solve(solution);
            // cout << "Sol : ";
            // for (uInt i=0; i<order+1; i++) {
            //   cout << solution[i] << " ";
            // }
            // cout << endl;

            casa::Vector<casa::Float> model(n);
            for (int i=0; i<n; i++) {
                model(i)=solution(order);
                for (int j=order-1; j>=0; j--) {
                    model(i) = i * model(i) + solution(j);
                }
            }
            vec -= model;
        }

    };

// Main function
int main(int argc, char* argv[])
{
    ImContSubApp app;
    return app.main(argc, argv);
}
