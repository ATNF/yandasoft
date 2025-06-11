/// @file SphFuncVisGridder.cc
///
/// @copyright (c) 2007,2015,2016 CSIRO
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

#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.sphfuncvisgridder");
#include <askap/askap/AskapError.h>

#include <askap/gridding/SphFuncVisGridder.h>
#include <casacore/casa/Arrays/ArrayIter.h>
#include <askap/scimath/fft/FFTWrapper.h>
#include <askap/profile/AskapProfiler.h>

namespace askap
{
  namespace synthesis
  {
    /// @brief Standard two dimensional gridding
    /// @param[in] alpha alpha value for prolate spheroidal function
    /// @param[in] support support size in pixels (spheroidal function with m=2*support will be generated)
    /// @param[in] oversample number of oversampling planes
    SphFuncVisGridder::SphFuncVisGridder(const float alpha,
                                         const int support,
                                         const int oversample) :
        itsSphFunc(casacore::C::pi*support, alpha), itsInterp(true), itsAlpha(static_cast<double>(alpha))
    {
       ASKAPASSERT(support>=3);
       ASKAPASSERT(oversample>=1);
       ASKAPASSERT(alpha>=0. && alpha<=2.);
       itsOverSample = oversample;
       itsSupport = support;
    }

    SphFuncVisGridder::~SphFuncVisGridder()
    {
    }

    /// @brief static method to create gridder
	/// @details Each gridder should have a static factory method, which is
	/// able to create a particular type of the gridder and initialise it with
	/// the parameters taken form the given parset. It is assumed that the
	/// method receives a subset of parameters where the gridder name is already
	/// taken out.
	/// @return a shared pointer to the gridder instance
	IVisGridder::ShPtr SphFuncVisGridder::createGridder(const LOFAR::ParameterSet &parset)
	{
	  const int oversample=parset.getInt32("oversample", 128);
	  const int support=parset.getInt32("support", 3);
	  const float alpha=parset.getFloat("alpha", 1.);
	  ASKAPLOG_INFO_STR(logger, "Setting up spheroidal function gridder with support="<<
	                    support<<" and oversample="<<oversample);

	  return IVisGridder::ShPtr(new SphFuncVisGridder(alpha,support,oversample));
	}

    /// Clone a copy of this Gridder
    IVisGridder::ShPtr SphFuncVisGridder::clone()
    {
      return IVisGridder::ShPtr(new SphFuncVisGridder(*this));
    }

    void SphFuncVisGridder::initIndices(const accessors::IConstDataAccessor&)
    {
    }

    /// Initialize the convolution function into the cube. If necessary this
    /// could be optimized by using symmetries.
    void SphFuncVisGridder::initConvolutionFunction(const accessors::IConstDataAccessor&)
    {
      ASKAPDEBUGTRACE("SphFuncVisGridder::initConvolutionFunction");
      if(itsConvFunc.size() != 0) {
         // a rather poor way of checking that convolution function has already been initialised
         return;
      }
      itsConvFunc.resize(itsOverSample*itsOverSample);

      if (isPCFGridder()) {

        // A simple grid kernel for use in setting up the preconditioner function.
        // Set up as a nearest neighbour gridder (based partly on the Box gridder).
        // Gridding weight is written to the real part, gridding support is written
        // to the imaginary part.

        itsSupport=1;
        const int cSize=2*itsSupport+1;
        const int cCenter=(cSize-1)/2;

        // This gridding function is divided out of the psf & dirty image
        // before preconditioning, so ideally the support should be reduced
        // back the pixel/s it falls in. Settng to 3 to be conservative.
        for (int fracu = 0; fracu < itsOverSample; ++fracu) {
            for (int fracv = 0; fracv < itsOverSample; ++fracv) {
                const int plane=fracu+itsOverSample*fracv;
                ASKAPDEBUGASSERT(plane>=0 && plane<int(itsConvFunc.size()));
                itsConvFunc[plane].resize(cSize, cSize);
                itsConvFunc[plane].set(0.0);
                // are fracu and fracv being correctly used here?
                // I think they should be -ve, since the offset in nux & nuy is +ve.
                const int ix = -float(fracu)/float(itsOverSample);
                const int iy = -float(fracv)/float(itsOverSample);
                itsConvFunc[plane](ix + cCenter, iy + cCenter) =  casacore::Complex(1.0, 3.0);
            }
        }

      } else {

        /// This must be changed for non-MFS

        const int cSize=2*itsSupport + 1; // 7;
        for (int fracv=0; fracv<itsOverSample; ++fracv) {
          for (int fracu=0; fracu<itsOverSample; ++fracu) {
            const int plane=fracu+itsOverSample*fracv;
            ASKAPDEBUGASSERT(plane>=0 && plane<int(itsConvFunc.size()));
            itsConvFunc[plane].resize(cSize, cSize);
            itsConvFunc[plane].set(0.0);
            casacore::Vector<imtypeComplex> bufx(cSize);
            casacore::Vector<imtypeComplex> bufy(cSize);
            for (int ix=0; ix<cSize; ++ix) {
              const double nux=std::abs(double(itsOverSample*(ix-itsSupport)+fracu))/double(itsSupport*itsOverSample);
              bufx(ix)=grdsf(nux)*std::pow(1.0-nux*nux, itsAlpha);
            }
            for (int iy=0; iy<cSize; ++iy) {
              const double nuy=std::abs(double(itsOverSample*(iy-itsSupport)+fracv))/double(itsSupport*itsOverSample);
              bufy(iy)=grdsf(nuy)*std::pow(1.0-nuy*nuy, itsAlpha);
            }
            for (int ix=0; ix<cSize; ++ix) {
              for (int iy=0; iy<cSize; ++iy) {
                itsConvFunc[plane](ix, iy)=bufx(ix)*bufy(iy);
              }
            }
          } // for fracu
        } // for fracv

      }

      // force normalization for all fractional offsets (or planes)
      for (size_t plane = 0; plane<itsConvFunc.size(); ++plane) {
           if (itsConvFunc[plane].nelements() == 0) {
               // this plane of the cache is unused
               continue;
           }
           const double norm = real(sum(casacore::abs(itsConvFunc[plane])));
           ASKAPDEBUGASSERT(norm>0.);
           itsConvFunc[plane]/=casacore::Complex(norm);
      } // for plane

    }

    void SphFuncVisGridder::correctConvolution(casacore::Array<imtype>& grid,
        scimath::SpheroidalFunction& sf, int support, bool interpolate)
    {
        casacore::IPosition shape = grid.shape();
        ASKAPDEBUGASSERT(shape.nelements()>=2);
        ASKAPDEBUGASSERT(shape(0)>1);
        ASKAPDEBUGASSERT(shape(1)>1);

        const casacore::Int xHalfSize = shape(0)/2;
        const casacore::Int yHalfSize = shape(1)/2;
        casacore::Vector<double> ccfx(shape(0));
        casacore::Vector<double> ccfy(shape(1));

        // initialise buffers to enable a filtering of the correction
        // function in Fourier space.
        casacore::Vector<imtypeComplex> bufx(shape(0));
        casacore::Vector<imtypeComplex> bufy(shape(1));

        // note grdsf(1)=0.
        for (int ix=0; ix<shape(0); ++ix)
        {
          const double nux=std::abs(double(ix-xHalfSize))/double(xHalfSize);
          const double val = sf(nux);
          bufx(ix) = imtypeComplex(val,0.0);
        }

        for (int iy=0; iy<shape(1); ++iy)
        {
          const double nuy=std::abs(double(iy-yHalfSize))/double(yHalfSize);
          const double val = sf(nuy);
          bufy(iy) = imtypeComplex(val,0.0);
        }

        if (interpolate) {
          // The spheroidal is undefined and set to zero at nu=1, but that
          // is not the numerical limit. Estimate it from its neighbours.
          interpolateEdgeValues(bufx);
          interpolateEdgeValues(bufy);
        }

        // Fourier filter the spheroidal (crop in Fourier space in line with
        // gridding kernel support size)
        const bool doFiltering = true;
        if (doFiltering) {
           // Some more advanced gridders have support>3 (e.g. w-proj).
           //
           int support = 3;
           const imtypeComplex maxBefore = bufx(shape(0)/2);
           scimath::fft(bufx, true);
           scimath::fft(bufy, true);
           for (int ix=0; ix<shape(0)/2-support; ++ix) {
             bufx(ix) = 0.0;
           }
           for (int ix=shape(0)/2+support+1; ix<shape(0); ++ix) {
             bufx(ix) = 0.0;
           }
           for (int iy=0; iy<shape(1)/2-support; ++iy) {
             bufy(iy) = 0.0;
           }
           for (int iy=shape(1)/2+support+1; iy<shape(1); ++iy) {
             bufy(iy) = 0.0;
           }
           scimath::fft(bufx, false);
           scimath::fft(bufy, false);
           // Normalise after filtering.
           const imtypeComplex normalisation = maxBefore / bufx(shape(0)/2);
           bufx *= normalisation;
           bufy *= normalisation;
        }

        for (int ix=0; ix<shape(0); ++ix) {
          double val = real(bufx(ix));
          ccfx(ix) = casacore::abs(val) > 1e-10 ? 1.0/val : 0.;
        }
        for (int iy=0; iy<shape(1); ++iy) {
          double val = real(bufy(iy));
          ccfy(iy) = casacore::abs(val) > 1e-10 ? 1.0/val : 0.;
        }

        casacore::ArrayIterator<imtype> it(grid, 2);
        while (!it.pastEnd())
        {
            casacore::Matrix<imtype> mat(it.array());
            ASKAPDEBUGASSERT(int(mat.nrow()) <= shape(0));
            ASKAPDEBUGASSERT(int(mat.ncolumn()) <= shape(1));
            for (int iy=0; iy<shape(1); iy++)
            {
                for (int ix=0; ix<shape(0); ix++)
                {
                    mat(ix, iy)*=ccfx(ix)*ccfy(iy);
                }
            }
            it.next();
        }

      }

    void SphFuncVisGridder::correctConvolution(casacore::Array<imtype>& grid)
    {
      ASKAPTRACE("SphFuncVisGridder::correctConvolution");
      if (isPCFGridder()) return;
      ASKAPDEBUGASSERT(itsShape.nelements()>=2);
      ASKAPDEBUGASSERT(itsShape(0)>1);
      ASKAPDEBUGASSERT(itsShape(1)>1);
      ASKAPDEBUGASSERT(itsShape==grid.shape());
      correctConvolution(grid, itsSphFunc, itsSupport, doInterpolation());
    }

    /*
    // find spheroidal function with m = 6, alpha = 1 using the rational
                  // approximations discussed by fred schwab in 'indirect imaging'.
                  // this routine was checked against fred's sphfn routine, and agreed
                  // to about the 7th significant digit.
                  // the gridding function is (1-nu**2)*grdsf(nu) where nu is the distance
                  // to the edge. the grid correction function is just 1/grdsf(nu) where nu
                  // is now the distance to the edge of the image.
                  double SphFuncVisGridder::grdsf1(double nu) const
                  {

                    double top, bot, delnusq, nuend;
                    int k, part;
                    int np, nq;
                    np=4;
                    nq=2;
                    double p[2][5] =
                      { {8.203343e-2, -3.644705e-1, 6.278660e-1, -5.335581e-1, 2.312756e-1},
                        {4.028559e-3, -3.697768e-2, 1.021332e-1, -1.201436e-1, 6.412774e-2}};
                    double q[2][3]=
                    { {1.0000000, 8.212018e-1, 2.078043e-1}, {1.0000000, 9.599102e-1,
                      2.918724e-1}};
                    double value = 0.0;

                    if ((nu>=0.0)&&(nu<0.75))
                    {
                      part = 0;
                      nuend = 0.75;
                    }
                    else if ((nu>=0.75)&&(nu<=1.00))
                    {
                      part = 1;
                      nuend = 1.00;
                    }
                    else
                    {
                      value = 0.0;
                      return value;
                    }

                    top = p[part][0];
                    bot = q[part][0];
                    delnusq = std::pow(nu, 2) - std::pow(nuend, 2);
                    for (k = 1; k<= np; k++)
                    {
                      double factor=std::pow(delnusq, k);
                      top += p[part][k] * factor;
                    }
                    for (k = 1; k<= nq; k++)
                    {
                      double factor=std::pow(delnusq, k);
                      bot += q[part][k] * factor;
                    }
                    if (bot!=0.0)
                    {
                      value = top/bot;
                    }
                    else
                    {
                      value = 0.0;
                    }
                    if (value<0.0)
                    value=0.0;
                    return value;
                  }
                 */

  }
}
