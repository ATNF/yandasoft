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

#include <askap/gridding/WStackVisGridder.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.wstackgridder");

#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <casacore/casa/BasicSL/Constants.h>
#include <askap/scimath/fft/FFTWrapper.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/profile/AskapProfiler.h>

using namespace askap;

#include <cmath>

namespace askap
{
  namespace synthesis
  {

    WStackVisGridder::WStackVisGridder(const double wmax, const int nwplanes) :
           WDependentGridderBase(wmax,nwplanes) {}

    WStackVisGridder::~WStackVisGridder() {}

    /// @brief copy constructor
    /// @details It is required to decouple internal arrays between
    /// input object and the copy
    /// @param[in] other input object
    WStackVisGridder::WStackVisGridder(const WStackVisGridder &other) :
       IVisGridder(other), WDependentGridderBase(other), itsGMap(other.itsGMap.copy()) {}


    /// Clone a copy of this Gridder
    IVisGridder::ShPtr WStackVisGridder::clone()
    {
      return IVisGridder::ShPtr(new WStackVisGridder(*this));
    }

    /// Initialize the convolution function into the cube. If necessary this
    /// could be optimized by using symmetries.
    void WStackVisGridder::initIndices(const accessors::IConstDataAccessor& acc)
    {
      ASKAPTRACE("WStackVisGridder::initIndices");
      /// We have to calculate the lookup function converting from
      /// row and channel to plane of the w-dependent convolution
      /// function
      const int nSamples = acc.nRow();
      const int nChan = acc.nChannel();
      const int nPol = acc.nPol();

      itsGMap.resize(nSamples, nPol, nChan);
      const casacore::Vector<casacore::RigidVector<double, 3> > &rotatedUVW = acc.rotatedUVW(getTangentPoint());
      for (int i=0; i<nSamples; ++i)
      {
        const double w=(rotatedUVW(i)(2))/(casacore::C::c);
        for (int chan=0; chan<nChan; ++chan)
        {
          for (int pol=0; pol<nPol; pol++)
          {
            const double freq=acc.frequency()[chan];
            /// Calculate the index into the grids
            itsGMap(i, pol, chan) = getWPlane(w*freq);
          }
        }
      }
    }

    void WStackVisGridder::initialiseGrid(const scimath::Axes& axes,
        const casacore::IPosition& shape, const bool dopsf, const bool dopcf)
    {
      ASKAPTRACE("WStackVisGridder::initialiseGrid");
      ASKAPDEBUGASSERT(shape.nelements()>=2);
      itsShape=scimath::PaddingUtils::paddedShape(shape,paddingFactor());

      initialiseCellSize(axes);

      initStokes();
      configureForPSF(dopsf);
      configureForPCF(dopcf);

      initialiseWeightBuilder();

      /// We need one grid for each plane
      itsGrid.resize(nWPlanes());
      for (int i=0; i<nWPlanes(); ++i)
      {
        itsGrid[i].resize(itsShape);
        itsGrid[i].set(0.0);
      }
      if (isPSFGridder())
      {
        // for a proper PSF calculation
		initRepresentativeFieldAndFeed();
      }

      initialiseSumOfWeights();

      initialiseFreqMapping();

      ASKAPLOG_INFO_STR(logger, "Gridding is set up with tangent centre "<<
             printDirection(getTangentPoint())<<" and image centre "<<
             printDirection(getImageCentre()));

    }

    void WStackVisGridder::multiply(casacore::Array<imtypeComplex>& scratch, int i)
    {
      ASKAPDEBUGTRACE("WStackVisGridder::multiply");
      /// These are the actual cell sizes used
      const float cellx=1.0/(float(itsShape(0))*itsUVCellSize(0));
      const float celly=1.0/(float(itsShape(1))*itsUVCellSize(1));

      const int nx=itsShape(0);
      const int ny=itsShape(1);

      const float w=2.0f*casacore::C::pi*getWTerm(i);
      casacore::ArrayIterator<imtypeComplex> it(scratch, 2);
      while (!it.pastEnd())
      {
        casacore::Matrix<imtypeComplex> mat(it.array());

        /// @todo Optimise multiply loop
        for (int iy=0; iy<ny; iy++)
        {
          float y2=float(iy-ny/2)*celly;
          y2*=y2;
          for (int ix=0; ix<nx; ix++)
          {
            if (casacore::abs(mat(ix, iy))>0.0)
            {
              float x2=float(ix-nx/2)*cellx;
              x2*=x2;
              const float r2=x2+y2;
              if (r2<1.0) {
                  const float phase=w*(1.0-sqrt(1.0-r2));
                  mat(ix, iy)*=imtypeComplex(cos(phase), -sin(phase));
              }
            }
          }
        }
        it.next();
      }
    }

    void WStackVisGridder::updatew(casacore::Array<imtypeComplex>& scratch, int i)
    {
      ASKAPDEBUGTRACE("WStackVisGridder::updatew");
      ASKAPDEBUGASSERT(itsUVCellSize(0)>0);

      // estimate the size of the projection in uv (would-be w-proj kernel size)
      //  - need to know if getWTerm returns the full w or smaller due to snapshotting
      const imtype wThetaPix = static_cast<imtype>(fabs(getWTerm(i)) / (itsUVCellSize(0) * itsUVCellSize(0)));
      imtype wKernelPix;
      if (wThetaPix < 1) {
        wKernelPix = 3;
      } else {
        wKernelPix = 6 + 1.14*wThetaPix;
      }
      // PCF gridding should have been done with a boxcar kernel = 1+3i or 1-3i
      wKernelPix /= 3.;

      const int nx=itsShape(0);
      const int ny=itsShape(1);

      const float w=2.0f*casacore::C::pi*getWTerm(i);
      casacore::ArrayIterator<imtypeComplex> it(scratch, 2);
      while (!it.pastEnd())
      {
        casacore::Matrix<imtypeComplex> mat(it.array());

        /// @todo Optimise loop
        for (int iy=0; iy<ny; iy++)
        {
          for (int ix=0; ix<nx; ix++)
          {
            if (casacore::abs(mat(ix, iy))>0.0)
            {
              // update the imaginary part but not the real part
              mat(ix, iy) = imtypeComplex(real(mat(ix, iy)), imag(mat(ix, iy)) * wKernelPix);
            }
          }
        }
        it.next();
      }
    }

    /// This is the default implementation
    void WStackVisGridder::finaliseGrid(casacore::Array<imtype>& out)
    {
      ASKAPTRACE("WStackVisGridder::finaliseGrid");
      if (isPSFGridder()) {
          ASKAPLOG_INFO_STR(logger, "Stacking " << nWPlanes()
                          << " planes of W stack to get final PSF");
      } else if (isPCFGridder()) {
          ASKAPLOG_INFO_STR(logger, "Stacking " << nWPlanes()
                          << " planes of W stack to get final preconditioner fuction");
      } else {
          ASKAPLOG_INFO_STR(logger, "Stacking " << nWPlanes()
                          << " planes of W stack to get final image");
      }
      ASKAPDEBUGASSERT(itsGrid.size()>0);

      if (isPCFGridder()) {

        // buffer for the result as doubles
        casacore::Array<imtypeComplex> cBuffer(itsGrid[0].shape());
        ASKAPDEBUGASSERT(cBuffer.shape().nelements()>=2);
       
        /// Loop over all grids Fourier transforming and accumulating
        bool first=true;
        for (unsigned int i=0; i<itsGrid.size(); i++)
        {
          if (casacore::max(casacore::amplitude(itsGrid[i]))>0.0)
          {
            casacore::Array<imtypeComplex> scratch(itsGrid[i].shape());
            casacore::convertArray<imtypeComplex,casacore::Complex>(scratch,itsGrid[i]);
            // Don't FFT yet. Stack uv grids
            // the w-support terms stored in the imaginary part of itsGrid do not account for w. Update them.
            updatew(scratch, i);

            if (first)  {
              first=false;
              cBuffer = scratch;
            } else {
              cBuffer += scratch;
            }
          }
        }
        // Now FFT
        scimath::fft2d(cBuffer, false);
        casacore::Array<imtype> dBuffer(itsGrid[0].shape());
        dBuffer = real(cBuffer);
        out = scimath::PaddingUtils::extract(dBuffer, paddingFactor());

      } else {

        // buffer for the result as doubles
        casacore::Array<imtype> dBuffer(itsGrid[0].shape());
        ASKAPDEBUGASSERT(dBuffer.shape().nelements()>=2);
       
        /// Loop over all grids Fourier transforming and accumulating
        bool first=true;
        for (unsigned int i=0; i<itsGrid.size(); i++)
        {
          if (casacore::max(casacore::amplitude(itsGrid[i]))>0.0)
          {
            casacore::Array<imtypeComplex> scratch(itsGrid[i].shape());
            casacore::convertArray<imtypeComplex,casacore::Complex>(scratch,itsGrid[i]);
            scimath::fft2d(scratch, false);
            multiply(scratch, i);
       
            if (first)  {
              first=false;
              dBuffer = real(scratch);
            } else {
              dBuffer += real(scratch);
            }
          }
        }
        // Now we can do the convolution correction
        correctConvolution(dBuffer);
        dBuffer *= static_cast<imtype>(dBuffer.shape()(0)*dBuffer.shape()(1));
        out = scimath::PaddingUtils::extract(dBuffer, paddingFactor());

      }

    }

    void WStackVisGridder::initialiseDegrid(const scimath::Axes& axes,
        const casacore::Array<imtype>& in)
    {
      ASKAPTRACE("WStackVisGridder::initialiseDegrid");
      itsShape = scimath::PaddingUtils::paddedShape(in.shape(),paddingFactor());
      configureForPSF(false);
      configureForPCF(false);

      initialiseCellSize(axes);
      initStokes();

      initialiseFreqMapping();

      itsGrid.resize(nWPlanes());
      if (casacore::max(casacore::abs(in))>0.0) {
        itsModelIsEmpty=false;
        ASKAPLOG_INFO_STR(logger, "Filling " << nWPlanes()
                           << " planes of W stack with model");
        casacore::Array<imtype> scratch(itsShape,static_cast<imtype>(0.));
        scimath::PaddingUtils::extract(scratch, paddingFactor()) = in;
        correctConvolution(scratch);
        for (int i=0; i<nWPlanes(); ++i)
        {
          casacore::Array<imtypeComplex> work(itsShape);
          toComplex(work, scratch);
          multiply(work, i);
          /// Need to conjugate to get sense of w correction correct
          work = casacore::conj(work);
          scimath::fft2d(work, true);
          itsGrid[i].resize(itsShape);
          casacore::convertArray<casacore::Complex,imtypeComplex>(itsGrid[i],work);
        }
      } else {
        itsModelIsEmpty=true;
        ASKAPLOG_INFO_STR(logger, "No need to fill W stack: model is empty");
        for (int i=0; i<nWPlanes(); ++i) {
          itsGrid[i].resize(casacore::IPosition(1, 1));
          itsGrid[i].set(casacore::Complex(0.0));
        }
      }
    }

    int WStackVisGridder::gIndex(int row, int pol, int chan)
    {
      const int plane = itsGMap(row, pol, chan);
      if (plane >=0) notifyOfWPlaneUse(plane);
      return plane;
    }

    /// @brief static method to create gridder
    /// @details Each gridder should have a static factory method, which is
    /// able to create a particular type of the gridder and initialise it with
    /// the parameters taken form the given parset. It is assumed that the
    /// method receives a subset of parameters where the gridder name is already
    /// taken out.
    /// @param[in] parset input parset file
    /// @return a shared pointer to the gridder instance
    IVisGridder::ShPtr WStackVisGridder::createGridder(const LOFAR::ParameterSet& parset)
    {
      double wmax=parset.getDouble("wmax", 35000.0);
      int nwplanes=parset.getInt32("nwplanes", 65);
      ASKAPLOG_INFO_STR(logger, "Gridding using W stacking with "<<nwplanes<<" w-planes in the stack");
      boost::shared_ptr<WStackVisGridder> gridder(new WStackVisGridder(wmax, nwplanes));
      gridder->configureWSampling(parset);
      return gridder;
    }

    /// @brief assignment operator
    /// @details It is required as private to avoid being called
    /// @param[in] other input object
    /// @return reference to itself
    WStackVisGridder& WStackVisGridder::operator=(const WStackVisGridder &)
    {
      ASKAPTHROW(AskapError, "This method is not supposed to be called!");
      return *this;
    }

  }
}
