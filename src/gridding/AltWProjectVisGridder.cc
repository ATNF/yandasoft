/// @file AltWProjectVisGridder.cc
///
/// @copyright (c) 2007,2016 CSIRO
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

// Package level header file
#include <askap_synthesis.h>

// System includes
#include <cmath>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <fft/FFTWrapper.h>
#include <profile/AskapProfiler.h>
//ASKAPSoft package includes
#include <gridding/WProjectVisGridder.h>
#include <gridding/SupportSearcher.h>
#include <utils/PaddingUtils.h>
#include <measurementequation/ImageParamsHelper.h>
#include <utils/ImageUtils.h>
#include <fft/FFTWrapper.h>
// Local package includes
#include "AltWProjectVisGridder.h"
#include <boost/lexical_cast.hpp>


ASKAP_LOGGER(logger, ".gridding.AltWProjectVisGridder");

namespace askap {
namespace synthesis {

AltWProjectVisGridder::AltWProjectVisGridder(const double wmax,
                                       const int nwplanes,
                                       const double cutoff,
                                       const int overSample,
                                       const int maxSupport,
                                       const int limitSupport,
                                       const std::string& name,
                                       const float alpha, const bool writeOut) :
        WProjectVisGridder(wmax, nwplanes, cutoff, overSample,maxSupport,limitSupport,name,alpha),
        itsWriteOut(writeOut)

{

}

AltWProjectVisGridder::~AltWProjectVisGridder()
{
}


/// Clone a copy of this Gridder
IVisGridder::ShPtr AltWProjectVisGridder::clone()
{
    return IVisGridder::ShPtr(new AltWProjectVisGridder(*this));
}


/// @brief static method to create gridder
/// @details Each gridder should have a static factory method, which is
/// able to create a particular type of the gridder and initialise it with
/// the parameters taken form the given parset. It is assumed that the
/// method receives a subset of parameters where the gridder name is already
/// taken out.
/// @param[in] parset input parset file
/// @return a shared pointer to the gridder instance
IVisGridder::ShPtr AltWProjectVisGridder::createGridder(const LOFAR::ParameterSet& parset)
{
    const double wmax = parset.getDouble("wmax", 35000.0);
    const int nwplanes = parset.getInt32("nwplanes", 65);
    const double cutoff = parset.getDouble("cutoff", 1e-3);
    const int oversample = parset.getInt32("oversample", 8);
    const int maxSupport = parset.getInt32("maxsupport", 256);
    const int limitSupport = parset.getInt32("limitsupport", 0);
    const string tablename = parset.getString("tablename", "");
    const float alpha=parset.getFloat("alpha", 1.);
    const bool writeOut=parset.getBool("dumpgrid",false);

    ASKAPLOG_INFO_STR(logger, "Gridding using (Alternate) W projection with " << nwplanes << " w-planes");
    boost::shared_ptr<AltWProjectVisGridder> gridder(new AltWProjectVisGridder(wmax, nwplanes,
            cutoff, oversample, maxSupport, limitSupport, tablename, alpha,writeOut));
    gridder->configureGridder(parset);
    gridder->configureWSampling(parset);
    return gridder;
}
void AltWProjectVisGridder::finaliseGrid(casa::Array<double>& out) {
    static int passThrough = 0;
    ASKAPTRACE("AltWProjectVisGridder::finaliseGrid");
    ASKAPLOG_INFO_STR(logger, "Using Alternate Finalise Grid ");
    ASKAPLOG_INFO_STR(logger, "There are " << itsGrid.size() << " grids");

    ASKAPDEBUGASSERT(itsGrid.size() > 0);
    // buffer for result as doubles
    casa::Array<double> dBuffer(itsGrid[0].shape());
    ASKAPDEBUGASSERT(dBuffer.shape().nelements()>=2);
    ASKAPDEBUGASSERT(itsShape == scimath::PaddingUtils::paddedShape(out.shape(),paddingFactor()));

    /// Loop over all grids Fourier transforming and accumulating
    for (unsigned int i=0; i<itsGrid.size(); i++) {
        casa::Array<casa::DComplex> scratch(itsGrid[i].shape());
        casa::convertArray<casa::DComplex,casa::Complex>(scratch, itsGrid[i]);

        if (itsWriteOut == true) {
          // for debugging
          ASKAPLOG_INFO_STR(logger, "Writing out Grids ");
          casa::Array<float> buf(scratch.shape());
          casa::convertArray<float,double>(buf,imag(scratch));
          string name = boost::lexical_cast<std::string>(passThrough) + "." + boost::lexical_cast<std::string>(i) + ".prefft.imag";
          scimath::saveAsCasaImage(name,buf);
          casa::convertArray<float,double>(buf,real(scratch));
          name = boost::lexical_cast<std::string>(passThrough) + "." + boost::lexical_cast<std::string>(i) + ".prefft.real";
          scimath::saveAsCasaImage(name,buf);
          /*
              // adjust values to extract part which gives a real symmetric FT and the remainder
              casa::Matrix<float> bufM(buf.nonDegenerate());
              for (int x=0; x<int(bufM.nrow()); ++x) {
                   for (int y=0; y<int(bufM.ncolumn())/2; ++y) {
                        const float val = 0.5*(bufM(x,y)+bufM(bufM.nrow() - x -1, bufM.ncolumn() - y -1));
                        bufM(x,y) = val;
                        bufM(bufM.nrow() - x -1, bufM.ncolumn() - y -1) = val;
                   }
              }
              scimath::saveAsCasaImage("uvcoverage.sympart",buf);
              casa::Matrix<casa::DComplex> scratchM(scratch.nonDegenerate());
              for (int x=0; x<int(scratchM.nrow()); ++x) {
                   for (int y=0; y<int(scratchM.ncolumn()); ++y) {
                        scratchM(x,y) -= double(bufM(x,y));
                    }
              }
              // as we ignore imaginary part after FT, make scratch hermitian to be fair
              for (int x=0; x<int(scratchM.nrow()); ++x) {
                   for (int y=0; y<int(scratchM.ncolumn())/2; ++y) {
                        const casa::DComplex val = 0.5*(scratchM(x,y)+
                              conj(scratchM(scratchM.nrow() - x -1, scratchM.ncolumn() - y -1)));
                        scratchM(x,y) = val;
                        scratchM(scratchM.nrow() - x -1, scratchM.ncolumn() - y -1) = conj(val);
                   }
              }
              casa::convertArray<float,double>(buf,imag(scratch));
              scimath::saveAsCasaImage("uvcoverage.asympart.imag",buf);
              casa::convertArray<float,double>(buf,real(scratch));
              scimath::saveAsCasaImage("uvcoverage.asympart.real",buf);
              scimath::fft2d(scratch, false);
              casa::convertArray<float,double>(buf,real(scratch));
              scimath::saveAsCasaImage("psf.asympart.real",buf);

              ASKAPCHECK(false, "Debug termination");
          */

        }

        scimath::fft2d(scratch, false);
        if (i==0) {
            toDouble(dBuffer, scratch);
        } else {
            casa::Array<double> work(dBuffer.shape());
            toDouble(work, scratch);
            dBuffer+=work;
        }

    }

    if (itsWriteOut) {
        string name = boost::lexical_cast<std::string>(passThrough) + ".postfft";
        casa::Array<float> buf(dBuffer.shape());
        casa::convertArray<float,double>(buf, dBuffer);
        scimath::saveAsCasaImage(name, buf);
    }

    // Now we can do the convolution correction
    correctConvolution(dBuffer, passThrough);

    if (itsWriteOut) {
        string name = boost::lexical_cast<std::string>(passThrough) + ".post_convolution_correction";
        casa::Array<float> buf(dBuffer.shape());
        casa::convertArray<float,double>(buf, dBuffer);
        scimath::saveAsCasaImage(name, buf);
    }

    dBuffer*=double(dBuffer.shape()(0))*double(dBuffer.shape()(1));
    out = scimath::PaddingUtils::extract(dBuffer,paddingFactor());

    if (itsWriteOut) {
        string name = boost::lexical_cast<std::string>(passThrough) + ".post_padding";
        casa::Array<float> buf(out.shape());
        casa::convertArray<float,double>(buf, out);
        scimath::saveAsCasaImage(name, buf);
    }

    passThrough++;
}



/// @brief assignment operator
/// @details Defined as private, so it can't be called (to enforce usage of the
/// copy constructor
/// @param[in] other input object
/// @return reference to itself
AltWProjectVisGridder& AltWProjectVisGridder::operator=(const AltWProjectVisGridder &)
{
    ASKAPTHROW(AskapError, "This method is not supposed to be called!");
    return *this;
}


} // namespace askap
} // namespace synthesis
