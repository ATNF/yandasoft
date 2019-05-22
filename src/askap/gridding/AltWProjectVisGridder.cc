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
#include <askap/askap_synthesis.h>

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
#include <askap/gridding/WProjectVisGridder.h>
#include <askap/gridding/SupportSearcher.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/measurementequation/ImageParamsHelper.h>
#include <askap/scimath/utils/ImageUtils.h>
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
void AltWProjectVisGridder::finaliseGrid(casacore::Array<double>& out) {
    static int passThrough = 0;
    ASKAPTRACE("AltWProjectVisGridder::finaliseGrid");
    ASKAPLOG_INFO_STR(logger, "Using Alternate Finalise Grid ");
    ASKAPLOG_INFO_STR(logger, "There are " << itsGrid.size() << " grids");

    ASKAPDEBUGASSERT(itsGrid.size() > 0);
    // buffer for result as doubles
    casacore::Array<double> dBuffer(itsGrid[0].shape());
    ASKAPDEBUGASSERT(dBuffer.shape().nelements()>=2);
    ASKAPDEBUGASSERT(itsShape == scimath::PaddingUtils::paddedShape(out.shape(),paddingFactor()));

    /// Loop over all grids Fourier transforming and accumulating
    for (unsigned int i=0; i<itsGrid.size(); i++) {
        casacore::Array<casacore::DComplex> scratch(itsGrid[i].shape());
        casacore::convertArray<casacore::DComplex,casacore::Complex>(scratch, itsGrid[i]);

        if (itsWriteOut == true) {
          // for debugging
          ASKAPLOG_INFO_STR(logger, "Writing out Grids ");
          casacore::Array<float> buf(scratch.shape());
          casacore::convertArray<float,double>(buf,imag(scratch));
          string name = boost::lexical_cast<std::string>(passThrough) + "." + boost::lexical_cast<std::string>(i) + ".prefft.imag";
          scimath::saveAsCasaImage(name,buf);
          casacore::convertArray<float,double>(buf,real(scratch));
          name = boost::lexical_cast<std::string>(passThrough) + "." + boost::lexical_cast<std::string>(i) + ".prefft.real";
          scimath::saveAsCasaImage(name,buf);
          /*
              // adjust values to extract part which gives a real symmetric FT and the remainder
              casacore::Matrix<float> bufM(buf.nonDegenerate());
              for (int x=0; x<int(bufM.nrow()); ++x) {
                   for (int y=0; y<int(bufM.ncolumn())/2; ++y) {
                        const float val = 0.5*(bufM(x,y)+bufM(bufM.nrow() - x -1, bufM.ncolumn() - y -1));
                        bufM(x,y) = val;
                        bufM(bufM.nrow() - x -1, bufM.ncolumn() - y -1) = val;
                   }
              }
              scimath::saveAsCasaImage("uvcoverage.sympart",buf);
              casacore::Matrix<casacore::DComplex> scratchM(scratch.nonDegenerate());
              for (int x=0; x<int(scratchM.nrow()); ++x) {
                   for (int y=0; y<int(scratchM.ncolumn()); ++y) {
                        scratchM(x,y) -= double(bufM(x,y));
                    }
              }
              // as we ignore imaginary part after FT, make scratch hermitian to be fair
              for (int x=0; x<int(scratchM.nrow()); ++x) {
                   for (int y=0; y<int(scratchM.ncolumn())/2; ++y) {
                        const casacore::DComplex val = 0.5*(scratchM(x,y)+
                              conj(scratchM(scratchM.nrow() - x -1, scratchM.ncolumn() - y -1)));
                        scratchM(x,y) = val;
                        scratchM(scratchM.nrow() - x -1, scratchM.ncolumn() - y -1) = conj(val);
                   }
              }
              casacore::convertArray<float,double>(buf,imag(scratch));
              scimath::saveAsCasaImage("uvcoverage.asympart.imag",buf);
              casacore::convertArray<float,double>(buf,real(scratch));
              scimath::saveAsCasaImage("uvcoverage.asympart.real",buf);
              scimath::fft2d(scratch, false);
              casacore::convertArray<float,double>(buf,real(scratch));
              scimath::saveAsCasaImage("psf.asympart.real",buf);

              ASKAPCHECK(false, "Debug termination");
          */

        }

        scimath::fft2d(scratch, false);
        if (i==0) {
            toDouble(dBuffer, scratch);
        } else {
            casacore::Array<double> work(dBuffer.shape());
            toDouble(work, scratch);
            dBuffer+=work;
        }
    }
    // Now we can do the convolution correction
    correctConvolution(dBuffer);
    dBuffer*=double(dBuffer.shape()(0))*double(dBuffer.shape()(1));
    out = scimath::PaddingUtils::extract(dBuffer,paddingFactor());
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
