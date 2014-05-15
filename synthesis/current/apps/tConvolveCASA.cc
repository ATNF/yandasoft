//   This C++ program has been written to demonstrate the
// convolutional resampling algorithm used in radio
// interferometry. It should compile with:
// 	g++ -O2 tConvolve.cc -o tConvolve
// Compiler option -O3 will optimize away the gridding!
//
// The challenge is to minimize the run time - specifically
// the time per grid addition. On a MacBookPro and an Opteron
// this is about 8ns.
//
//   For further details contact Tim.Cornwell@csiro.au
// May 3, 2007
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


#include <iostream>
#include <cmath>
#include <ctime>
#include <complex>
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;
using std::complex;
using std::abs;

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <scimath/Mathematics/RigidVector.h>

// Typedefs for easy testing
// Cost of using double for Coord is low, cost for
// double for Real is also low
typedef casa::Double Coord;
typedef casa::Float Real;
typedef casa::Complex Value;

// Perform standard data independent gridding
//
// u,v,w - components of spatial frequency
// data - values to be gridded
// freq - temporal frequency (inverse wavelengths)
// cellSize - size of one grid cell in wavelengths
// C - convolution function
// support - Total width of convolution function=2*support+1
// overSample - Oversampling factor for the convolution function
// cOffset - offsets into convolution function per data point
// grid - Output grid
//

int generic(const casa::Vector<casa::RigidVector<casa::Double, 3> >& uvw,
casa::Cube<casa::Complex>& data,
const casa::Vector<casa::Double>& freq,
const Coord cellSize,
const casa::Cube<Real>& C,
const int support,
const int overSample,
const casa::Matrix<casa::uInt>& cOffset,
casa::Cube<Value>& grid)
{

  const int gSize = grid.ncolumn();
  std::cout << "Grid size = " << gSize << std::endl;
  const int nSamples = uvw.size();
  const int nChan = freq.size();
  const int nPol = data.shape()(2);

  //int cSize=2*(support+1)*overSample+1;

// Grid
  grid.set(Value(0));

  cout << "+++++ Forward processing +++++" << endl;

  clock_t start,finish;
  double time;

  Real sumwt=0.0;

  start = clock();
// Loop over all samples adding them to the grid
// First scale to the correct pixel location
// Then find the fraction of a pixel to the nearest pixel
// Loop over the entire support, calculating weights from
// the convolution function and adding the scaled
// visibility to the grid.
  for (int i=0;i<nSamples;i++)
  {
    for (int chan=0;chan<nChan;chan++)
    {
      for (int pol=0;pol<nPol;pol++)
      {

        int coff=cOffset(i,chan);

        Coord uScaled=freq[chan]*uvw(i)(0)/cellSize;
        int iu=int(uScaled);
        int fracu=int(overSample*(uScaled-Coord(iu)));
        iu+=gSize/2;

        Coord vScaled=freq[chan]*uvw(i)(0)/cellSize;
        int iv=int(vScaled);
        int fracv=int(overSample*(vScaled-Coord(iv)));
        iv+=gSize/2;

        casa::Complex& d=data(i,chan,pol);
//#define USEPOINTERS 1
#ifdef USEPOINTERS
        casa::Complex *gptr;
        const Real *cptr;
        for (int suppu=-support;suppu<+support;suppu++)
        {
          gptr=&grid(iu+suppu,iv-support,pol);
          cptr=&C(iu+suppu*overSample+fracu,iv-support*overSample+fracv,coff);
          for (int suppv=-support;suppv<+support;suppv++)
          {
//            Real wt=C(iu+suppu*overSample+fracu,iu+suppv*overSample+fracu,coff);
//            grid(iu+suppu,iv+suppv,pol)+=wt*data(i,chan,pol);
            (*gptr)+=(*cptr)*d;
            sumwt+=(*cptr);
            gptr++;
            cptr+=overSample;
          }
        }
#else
        for (int suppu=-support;suppu<+support;suppu++)
        {
          for (int suppv=-support;suppv<+support;suppv++)
          {
            Real wt=C(iu+suppu*overSample+fracu,iv+suppv*overSample+fracv,coff);
            grid(iu+suppu,iv+suppv,pol)+=wt*d;
          }
        }
#endif
      }
    }
  }
  finish = clock();
// Report on timings
  cout << "    Total weight = " << sumwt << endl;
  time = (double(finish)-double(start))/CLOCKS_PER_SEC;
  cout << "    Time " << time << " (s) " << endl;
  cout << "    Time per visibility sample " << 1e6*time/double(nSamples)
    << " (us) " << endl;
  cout << "    Time per visibility spectral sample "
    << 1e6*time/double(nSamples*nChan) << " (us) " << endl;
  cout << "    Time per grid-addition "
    << 1e9*time/(double(nSamples)*double(nChan)*double((2*support)*(2*support+1)))
    << " (ns) " << endl;

  cout << "+++++ Reverse processing +++++" << endl;

// Just run the gridding in reverse
  start = clock();
  for (int i=0;i<nSamples;i++)
  {
    for (int chan=0;chan<nChan;chan++)
    {
      for (int pol=0;pol<nPol;pol++)
      {

        int coff=cOffset(i,chan);

        Coord uScaled=freq[chan]*uvw(i)(0)/cellSize;
        int iu=int(uScaled);
        int fracu=int(overSample*(uScaled-Coord(iu)));
        iu+=gSize/2;

        Coord vScaled=freq[chan]*uvw(i)(1)/cellSize;
        int iv=int(vScaled);
        int fracv=int(overSample*(vScaled-Coord(iv)));
        iv+=gSize/2;

        double sumviswt=0.0;
        casa::Complex& d=data(i,chan,pol);
#ifdef USEPOINTERS
        casa::Complex *gptr;
        const Real *cptr;
        for (int suppu=-support;suppu<+support;suppu++)
        {
          gptr=&grid(iu+suppu,iv-support,pol);
          cptr=&C(iu+suppu*overSample+fracu,iv-support*overSample+fracv,coff);
          for (int suppv=-support;suppv<+support;suppv++)
          {
            d+=(*cptr)*(*gptr);
            sumviswt+=(*cptr);
            gptr++;
            cptr+=overSample;
          }
        }
#else
        for (int suppu=-support;suppu<+support;suppu++)
        {
          for (int suppv=-support;suppv<+support;suppv++)
          {
            Real wt=C(iu+fracu+suppu*overSample,iv+fracv+suppv*overSample,coff);
            d+=wt*grid(iu+suppu,iv+suppv,pol);
            sumviswt+=wt;
          }
        }
#endif
        data(i,chan,pol)=data(i,chan,pol)/casa::Complex(sumviswt);
      }
    }
  }
  finish = clock();
// Report on timings
  time = (double(finish)-double(start))/CLOCKS_PER_SEC;
  cout << "    Time " << time << " (s) " << endl;
  cout << "    Time per visibility sample " << 1e6*time/double(nSamples)
    << " (us) " << endl;
  cout << "    Time per visibility spectral sample "
    << 1e6*time/double(nSamples*nChan) << " (us) " << endl;
  cout << "    Time per grid-addition "
    << 1e9*time/(double(nSamples)*double(nChan)*double((2*support)*(2*support+1)))
    << " (ns) " << endl;

  return 0;
}


int standard(const casa::Vector<casa::RigidVector<casa::Double, 3> >& uvw,
casa::Cube<casa::Complex>& data,
const casa::Vector<casa::Double>& freq,
const Coord cellSize,
casa::Cube<Value>& grid)

{
  cout << "*************************** Standard gridding ***********************" << endl;
  int support=3;                                  // Support for gridding function in pixels
  const int overSample=100;
  cout << "Support = " << support << " pixels" << endl;

// Convolution function
// We take this to be the product of two Gaussian. More often it
// is the product of two prolate spheroidal wave functions
  int cSize=2*(support+1)*overSample+1;

  casa::Cube<Real> C(cSize,cSize,1);

  int cCenter=(cSize-1)/2;

// Keep this symmetrically to streamline index handling later....
  for (int i=0;i<cSize;i++)
  {
    double i2=std::pow(double(i-cCenter)/double(overSample), 2);
    for (int j=0;j<cSize;j++)
    {
      double r2=i2+std::pow(double(j-cCenter)/double(overSample), 2);
      C(i,j,0)=std::exp(-r2);
    }
  }

  const int nSamples = uvw.size();
  const int nChan = freq.size();

  casa::Matrix<uint> cOffset(nSamples, nChan);
  cOffset.set(0);

  return generic(uvw, data, freq, cellSize, C, support, overSample, cOffset, grid);
}


// Perform w projection (data dependent) gridding
//
// u,v,w - components of spatial frequency
// data - values to be gridded
// nSamples - number of visibility samples
// freq - temporal frequency (inverse wavelengths)
// cellSize - size of one grid cell in wavelengths
// gSize - size of grid in pixels (per axis)
// support - Total width of convolution function=2*support+1
// wCellSize - size of one w grid cell in wavelengths
// wSize - Size of lookup table in w
int wprojection(const casa::Vector<casa::RigidVector<casa::Double, 3> >& uvw,
casa::Cube<casa::Complex>& data,
const casa::Vector<casa::Double>& freq,
const Coord cellSize,
const Coord baseline,
const int wSize,
casa::Cube<Value>& grid)
{

  const int nSamples = uvw.size();
  const int nChan = freq.size();

  cout << "************************* W projection gridding *********************" << endl;
  int support=static_cast<int>(1.5*sqrt(abs(baseline)*static_cast<Coord>(cellSize)*freq(0))/cellSize);
  int overSample=8;
  cout << "Support = " << support << " pixels" << endl;
  const Coord wCellSize=2*baseline*freq(0)/wSize;
  cout << "W cellsize = " << wCellSize << " wavelengths" << endl;

// Convolution function. This should be the convolution of the
// w projection kernel (the Fresnel term) with the convolution
// function used in the standard case. The latter is needed to
// suppress aliasing. In practice, we calculate entire function
// by Fourier transformation. Here we take an approximation that
// is good enough.
  int cSize=2*(support+1)*overSample+1;

  int cCenter=(cSize-1)/2;

  casa::Cube<Real> C(cSize, cSize, wSize);

  for (int k=0;k<wSize;k++)
  {
    if(k!=wSize/2)
    {
      double w=double(k-wSize/2);
      double fScale=sqrt(abs(w)*wCellSize*freq[0])/cellSize;
      for (int j=0;j<cSize;j++)
      {
        double j2=std::pow(double(j-cCenter)/double(overSample), 2);
        for (int i=0;i<cSize;i++)
        {
          double r2=j2+std::pow(double(i-cCenter)/double(overSample), 2);
          C(i,j,k)=static_cast<Real>(std::cos(r2/(w*fScale)));
        }
      }
    }
    else
    {
      for (int j=0;j<cSize;j++)
      {
        double j2=std::pow(double(j-cCenter)/double(overSample), 2);
        for (int i=0;i<cSize;i++)
        {
          double r2=j2+std::pow(double(i-cCenter)/double(overSample), 2);
          C(i,j,k)=static_cast<Real>(std::exp(-r2));
        }
      }
    }
  }

  casa::Matrix<uint> cOffset(nSamples, nChan);
  for (int i=0;i<nSamples;i++)
  {
    for (int chan=0;chan<nChan;chan++)
    {

      Coord wScaled=freq(chan)*uvw(1)(2)/wCellSize;
      cOffset(i,chan)=wSize/2+int(wScaled);
    }
  }

  return generic(uvw, data, freq, cellSize, C, support, overSample, cOffset, grid);
}


int main()
{
  const int baseline=2000;                        // Maximum baseline in meters
  const int nSamples=100000;                      // Number of data samples
  const int gSize=512;                            // Size of output grid in pixels
  const Coord cellSize=50;                        // Cellsize of output grid in wavelengths
  const int wSize=64;                             // Number of lookup planes in w projection
  const int nChan=16;                             // Number of spectral channels
  const int nPol=1;                               // Number of polarizations

// Initialize the data to be gridded
  casa::Vector<casa::RigidVector<casa::Double, 3> > uvw(nSamples);
  casa::Cube<casa::Complex> data(nSamples, nChan, nPol);

  for (int i=0;i<nSamples;i++)
  {
    uvw(i)(0)=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;
    uvw(i)(1)=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;
    uvw(i)(2)=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;

    for (int chan=0;chan<nChan;chan++)
    {
      for (int pol=0;pol<nPol;pol++)
      {
        data(i,chan,pol)=casa::Complex(rand()/float(RAND_MAX), rand()/float(RAND_MAX));
      }
    }
  }

// Measure frequency in inverse wavelengths
  casa::Vector<casa::Double> freq(nChan);
  for (int i=0;i<nChan;i++)
  {
    freq(i)=(1.4e9-2.0e5*Coord(i)/Coord(nChan))/2.998e8;
  }

  casa::Cube<Value> grid(gSize, gSize, nPol);

  standard(uvw, data, freq, cellSize, grid);

  wprojection(uvw, data, freq, cellSize, baseline, wSize, grid);

  cout << "Done" << endl;

  return 0;
}
