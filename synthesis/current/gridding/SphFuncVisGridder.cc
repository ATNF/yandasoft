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

#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.sphfuncvisgridder");
#include <askap/AskapError.h>

#include <gridding/SphFuncVisGridder.h>
#include <casa/Arrays/ArrayIter.h>
#include <profile/AskapProfiler.h>


namespace askap
{
  namespace synthesis
  {
    /// @brief Standard two dimensional gridding
    /// @param[in] support support size in pixels (spheroidal function with m=2*support will be generated)
    /// @param[in] oversample number of oversampling planes
    SphFuncVisGridder::SphFuncVisGridder(int support, int oversample) : itsSphFunc(casa::C::pi*support, 1.)
    {
       ASKAPASSERT(support>=3);
       ASKAPASSERT(oversample>=1);       
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
	  ASKAPLOG_INFO_STR(logger, "Setting up spheroidal function gridder with support="<<
	                    support<<" and oversample="<<oversample);
	                    
	  return IVisGridder::ShPtr(new SphFuncVisGridder(support,oversample));
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

      const int cSize=2*itsSupport+1; // 7;

      /// This must be changed for non-MFS

      for (int fracv=0; fracv<itsOverSample; ++fracv) {
        for (int fracu=0; fracu<itsOverSample; ++fracu) {
          const int plane=fracu+itsOverSample*fracv;
          ASKAPDEBUGASSERT(plane>=0 && plane<int(itsConvFunc.size()));
          itsConvFunc[plane].resize(cSize, cSize);
          itsConvFunc[plane].set(0.0);
          for (int ix=0; ix<cSize; ++ix) {
            double nux=std::abs(double(itsOverSample*(ix-itsSupport)+fracu))/double(itsSupport*itsOverSample);
            double fx=grdsf(nux)*(1.0-std::pow(nux, 2));
            for (int iy=0; iy<cSize; ++iy) {
              double nuy=std::abs(double(itsOverSample*(iy-itsSupport)+fracv))/double(itsSupport*itsOverSample);
              double fy=grdsf(nuy)*(1.0-std::pow(nuy, 2));
              itsConvFunc[plane](ix, iy)=fx*fy;
            } // for iy
          } // for ix
        } // for fracu
      } // for fracv
      
      // force normalization for all fractional offsets (or planes)
      for (size_t plane = 0; plane<itsConvFunc.size(); ++plane) {
           if (itsConvFunc[plane].nelements() == 0) {
               // this plane of the cache is unused
               continue;
           }
           const double norm = real(sum(casa::abs(itsConvFunc[plane])));
           ASKAPDEBUGASSERT(norm>0.);
           itsConvFunc[plane]/=casa::Complex(norm);
      } // for plane					        
      
    }

    void SphFuncVisGridder::correctConvolution(casa::Array<double>& grid)
    { 
      ASKAPTRACE("SphFuncVisGridder::correctConvolution");
      ASKAPDEBUGASSERT(itsShape.nelements()>=2);      
      const casa::Int xHalfSize = itsShape(0)/2;
      const casa::Int yHalfSize = itsShape(1)/2;
      casa::Vector<double> ccfx(itsShape(0));
      casa::Vector<double> ccfy(itsShape(1));
      ASKAPDEBUGASSERT(itsShape(0)>1);
      ASKAPDEBUGASSERT(itsShape(1)>1);
      
      
      // note grdsf(-1)=0.
      for (int ix=0; ix<itsShape(0); ++ix)
      {
        const double nux=std::abs(double(ix-xHalfSize))/double(xHalfSize);
        const double val = grdsf(nux);
        ccfx(ix) = casa::abs(val) > 1e-10 ? 1.0/val : 0.;             
      }
      for (int iy=0; iy<itsShape(1); ++iy)
      {
        const double nuy=std::abs(double(iy-yHalfSize))/double(yHalfSize);
        const double val = grdsf(nuy);
        ccfy(iy) = casa::abs(val) > 1e-10 ? 1.0/val : 0.;                     
      }

      casa::ArrayIterator<double> it(grid, 2);
      while (!it.pastEnd())
      {
        casa::Matrix<double> mat(it.array());
        ASKAPDEBUGASSERT(int(mat.nrow()) <= itsShape(0));
        ASKAPDEBUGASSERT(int(mat.ncolumn()) <= itsShape(1));        
        for (int ix=0; ix<itsShape(0); ix++)
        {
          for (int iy=0; iy<itsShape(1); iy++)
          {
            mat(ix, iy)*=ccfx(ix)*ccfy(iy);
          }
        }
        it.next();
      }
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
