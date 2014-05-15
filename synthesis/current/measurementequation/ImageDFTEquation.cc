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

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.imagedftequation");

#include <askap/AskapError.h>

#include <dataaccess/SharedIter.h>
#include <fitting/Params.h>
#include <measurementequation/ImageDFTEquation.h>
#include <fitting/GenericNormalEquations.h>
#include <fitting/DesignMatrix.h>
#include <fitting/Axes.h>

#include <scimath/Mathematics/RigidVector.h>
#include <casa/BasicSL/Constants.h>
#include <casa/BasicSL/Complex.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>

#include <stdexcept>

using askap::scimath::Params;
using askap::scimath::Axes;
using askap::scimath::DesignMatrix;

namespace askap
{
  namespace synthesis
  {

    ImageDFTEquation::ImageDFTEquation(const askap::scimath::Params& ip,
      accessors::IDataSharedIter& idi) : scimath::Equation(ip),
                  askap::scimath::GenericEquation(ip), itsIdi(idi) 
      {
        init();
      };
        
    ImageDFTEquation::ImageDFTEquation(accessors::IDataSharedIter& idi) :  
      itsIdi(idi) 
    {
      reference(defaultParameters().clone());
      init();
    }

    ImageDFTEquation::~ImageDFTEquation() 
    {
    }
    
    ImageDFTEquation::ImageDFTEquation(const ImageDFTEquation& other) :
         Equation(other), GenericEquation(other)
    {
      operator=(other);
    }

    ImageDFTEquation& ImageDFTEquation::operator=(const ImageDFTEquation& other)
    {
      if(this!=&other)
      {
        static_cast<askap::scimath::GenericEquation*>(this)->operator=(other);
        itsIdi=other.itsIdi;
      }
      return *this;
    }

    void ImageDFTEquation::init()
    {
    }
    
    askap::scimath::Params ImageDFTEquation::defaultParameters()
    {
      Params ip;
      ip.add("image");
      return ip;
    }
    
    /// Clone this into a shared pointer
    /// @return shared pointer to a copy
    ImageDFTEquation::ShPtr ImageDFTEquation::clone() const
    {
      return ImageDFTEquation::ShPtr(new ImageDFTEquation(*this));
    }

    void ImageDFTEquation::predict() const
    {
      vector<string> completions(parameters().completions("image.i"));
      vector<string>::const_iterator it;

      if(completions.size()==0) {
        ASKAPLOG_WARN_STR(logger, "No parameters appropriate for ImageFFTEquation");
        return;
      }

//      itsIdi.chooseBuffer("model");

      for (itsIdi.init();itsIdi.hasMore();itsIdi.next())
      {

        const casa::Vector<double>& freq=itsIdi->frequency();
        //const double time=itsIdi->time();
        const uint nChan=freq.nelements();
        const uint nRow=itsIdi->nRow();
        casa::Matrix<double> vis(nRow,2*nChan);
        vis.set(0.0);

        for (it=completions.begin();it!=completions.end();it++)
        {

          string imageName("image.i"+(*it));

          const casa::Array<double> imagePixels(parameters().value(imageName));
          const casa::IPosition imageShape(imagePixels.shape());

          Axes axes(parameters().axes(imageName));
          if(!axes.has("RA")||!axes.has("DEC"))
          {
            throw(std::invalid_argument("RA and DEC specification not present for "+imageName));
          }
          double raStart=axes.start("RA");
          double raEnd=axes.end("RA");
          int raCells=imageShape(axes.order("RA"));

          double decStart=axes.start("DEC");
          double decEnd=axes.end("DEC");
          int decCells=imageShape(axes.order("DEC"));

          casa::Matrix<double> noDeriv(0,0);

          this->calcVisDFT(imagePixels, raStart, raEnd, raCells, decStart, decEnd, decCells,
            freq, itsIdi->uvw(), vis, false, noDeriv);

          for (uint row=0;row<nRow;row++)
          {
            for (uint i=0;i<nChan;i++)
            {
              itsIdi->rwVisibility()(row,i,0) += casa::Complex(vis(row,2*i), vis(row,2*i+1));
            }
          }
        }
      }
    };

    void ImageDFTEquation::calcGenericEquations(askap::scimath::GenericNormalEquations& ne) const
    {
// Loop over all completions i.e. all sources
      vector<string> completions(parameters().completions("image.i"));
      vector<string>::iterator it;

      if(completions.size()==0) {
        ASKAPLOG_WARN_STR(logger, "No parameters appropriate for ImageFFTEquation");
        return;
      }

//      itsIdi.chooseOriginal();

      for (itsIdi.init();itsIdi.hasMore();itsIdi.next())
      {

        const casa::Vector<double>& freq=itsIdi->frequency();
        const uint nChan=freq.nelements();
        const uint nRow=itsIdi->nRow();
        //const double time=itsIdi->time();

// Set up arrays to hold the output values
// Row, Two values (complex) per channel, single pol
        casa::Vector<double> residual(2*nRow*nChan);
        casa::Vector<double> weights(2*nRow*nChan);
        weights.set(1.0);
        casa::Matrix<double> vis(nRow,2*nChan);
        vis.set(0.0);

        for (it=completions.begin();it!=completions.end();it++)
        {
          string imageName("image.i"+(*it));
          if(parameters().isFree(imageName)) {

            const casa::Array<double> imagePixels(parameters().value(imageName));
            const casa::IPosition imageShape(imagePixels.shape());
  
            Axes axes(parameters().axes(imageName));
            if(!axes.has("RA")||!axes.has("DEC"))
            {
              throw(std::invalid_argument("RA and DEC specification not present for "+imageName));
            }
            double raStart=axes.start("RA");
            double raEnd=axes.end("RA");
            int raCells=imageShape(axes.order("RA"));
  
            double decStart=axes.start("DEC");
            double decEnd=axes.end("DEC");
            int decCells=imageShape(axes.order("DEC"));
            const uint nPixels=imagePixels.nelements();
  
            DesignMatrix designmatrix; //old parameters: parameters();
            casa::Matrix<double> imageDeriv(2*nRow*nChan,nPixels);
  
            this->calcVisDFT(imagePixels, raStart, raEnd, raCells,
              decStart, decEnd, decCells, freq, itsIdi->uvw(),
              vis, true, imageDeriv);
  
            for (uint row=0;row<itsIdi->nRow();row++)
            {
              for (uint i=0;i<freq.nelements();i++)
              {
                residual(nChan*row+2*i)=real(itsIdi->visibility()(row,i,0))-vis(row,2*i);
                residual(nChan*row+2*i+1)=imag(itsIdi->visibility()(row,i,0))-vis(row,2*i+1);
              }
            }
  
  // Now we can add the design matrix, residual, and weights
            designmatrix.addDerivative(imageName, imageDeriv);
            designmatrix.addResidual(residual, weights);
            ne.add(designmatrix);
          }
        }
      }
    };

    void ImageDFTEquation::calcVisDFT(const casa::Array<double>& imagePixels,
      const double raStart, const double raEnd, const int raCells,
      const double decStart, const double decEnd, const int decCells,
      const casa::Vector<double>& freq,
      const casa::Vector<casa::RigidVector<double, 3> >& uvw,
      casa::Matrix<double>& vis, bool doDeriv, casa::Matrix<double>& imageDeriv)
    {
      double raInc=(raStart-raEnd)/double(raCells);
      double decInc=(decStart-decEnd)/double(decCells);
      const uint nRow=uvw.nelements();
      const uint nChan=freq.nelements();

      vis.set(0.0);

      for (uint row=0;row<nRow;row++)
      {
        uint pixel=0;
        double u=uvw(row)(0);
        double v=uvw(row)(1);
        double w=uvw(row)(2);
        ASKAPDEBUGASSERT(decCells>=0);
        ASKAPDEBUGASSERT(raCells>=0);
        
        for (uint m=0;m<uint(decCells);++m)
        {
          double dec = decStart + m * decInc;
          for (uint l=0;l<uint(raCells);l++)
          {
            double ra = raStart + l * raInc;
            double delay = casa::C::_2pi * (ra * u + dec * v + sqrt(1 - ra * ra - dec * dec) * w)/casa::C::c;
            double flux = imagePixels(casa::IPosition(2, l, m));
            if(doDeriv)
            {
              for (uint i=0;i<nChan;i++)
              {
                double phase = delay * freq(i);
                vis(row,2*i) += flux * cos(phase);
                vis(row,2*i+1) += flux * sin(phase);
                imageDeriv(nChan*row+2*i,pixel) = cos(phase);
                imageDeriv(nChan*row+2*i+1,pixel) = sin(phase);
              }
            }
            else
            {
              for (uint i=0;i<nChan;i++)
              {
                double phase = delay * freq(i);
                vis(row,2*i) += flux * cos(phase);
                vis(row,2*i+1) += flux * sin(phase);
              }
            }
            pixel++;
          }
        }
      }
    }

  }

}
