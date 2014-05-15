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
// June 5, 2007

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

// Typedefs for easy testing
// Cost of using double for Coord is low, cost for
// double for Real is also low
typedef float Coord;
typedef float Real;
typedef std::complex<Real> Value;

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

int generic(const std::vector<Coord>& u,
const std::vector<Coord>& v,
const std::vector<Coord>& w,
std::vector<Value>& data,
const std::vector<Coord>& freq,
const Coord cellSize,
const std::vector<Real>& C,
const int support,
const int overSample,
const std::vector<unsigned int>& cOffset,
std::vector<Value>& grid)
{

  const int gSize = static_cast<int>(std::sqrt(static_cast<float>(grid.size())));
  std::cout << "Grid size = " << gSize << std::endl;
  const int nSamples = u.size();
  const int nChan = freq.size();

  int cSize=2*(support+1)*overSample+1;

  int cCenter=(cSize-1)/2;

// Grid
  grid.assign(grid.size(), Value(0.0));

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

      int find=i*nChan+chan;

      int coff=cOffset[find];

      Coord uScaled=freq[chan]*u[i]/cellSize;
      int iu=int(uScaled);
      int fracu=int(overSample*(uScaled-Coord(iu)));
      iu+=gSize/2;

      Coord vScaled=freq[chan]*v[i]/cellSize;
      int iv=int(vScaled);
      int fracv=int(overSample*(vScaled-Coord(iv)));
      iv+=gSize/2;

      for (int suppv=-support;suppv<+support;suppv++)
      {
        int vind=cSize*(fracv+overSample*suppv+cCenter)+fracu+cCenter+coff;
        int gind=iu+gSize*(iv+suppv);
        for (int suppu=-support;suppu<+support;suppu++)
        {
          Real wt=C[vind+overSample*suppu];
          grid[gind+suppu]+=wt*data[find];
          sumwt+=wt;
        }
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
    << 1e9*time/(double(nSamples)*double(nChan)*
    double((2*support)*(2*support+1)))
    << " (ns) " << endl;

  cout << "+++++ Reverse processing +++++" << endl;

// Just run the gridding in reverse
  start = clock();
  for (int i=0;i<nSamples;i++)
  {
    for (int chan=0;chan<nChan;chan++)
    {

      Real sumviswt=0.0;

      int find=i*nChan+chan;

      int coff=cOffset[find];

      Coord uScaled=freq[chan]*u[i]/cellSize;
      int iu=int(uScaled);
      int fracu=int(overSample*(uScaled-Coord(iu)));
      iu+=gSize/2;

      Coord vScaled=freq[chan]*v[i]/cellSize;
      int iv=int(vScaled);
      int fracv=int(overSample*(vScaled-Coord(iv)));
      iv+=gSize/2;

      for (int suppv=-support;suppv<+support;suppv++)
      {
        int vind=cSize*(fracv+overSample*suppv+cCenter)+fracu+cCenter+coff;
        int gind=iu+gSize*(iv+suppv);
        for (int suppu=-support;suppu<+support;suppu++)
        {
          Real wt=C[vind+overSample*suppu];
          data[find]=data[find]+wt*grid[gind+suppu];
          sumviswt+=wt;
        }
      }
      data[find]=data[find]/sumviswt;
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
    << 1e9*time/(double(nSamples)*double(nChan)*
    double((2*support)*(2*support+1)))
    << " (ns) " << endl;

  return 0;
}


int standard(const std::vector<Coord>& u, const std::vector<Coord>& v,
const std::vector<Coord>& w,
std::vector<Value>& data,
const std::vector<Coord>& freq,
const Coord cellSize,
std::vector<Value>& grid)
{

  cout << "*************************** Standard gridding ***********************"
    << endl;
  int support=3;                                  // Support for gridding function in pixels
  const int overSample=100;
  cout << "Support = " << support << " pixels" << endl;

// Convolution function
// We take this to be the product of two Gaussian. More often it
// is the product of two prolate spheroidal wave functions
  int cSize=2*(support+1)*overSample+1;

  std::vector<Real> C(cSize*cSize);

  int cCenter=(cSize-1)/2;

// Keep this symmetrically to streamline index handling later....
  for (int i=0;i<cSize;i++)
  {
    double i2=std::pow(double(i-cCenter)/double(overSample), 2);
    for (int j=0;j<cSize;j++)
    {
      double r2=i2+std::pow(double(j-cCenter)/double(overSample), 2);
      C[i+cSize*j]=std::exp(-r2);
    }
  }

  std::vector<unsigned int> cOffset;
  cOffset.assign(data.size(),0);

  return generic(u, v, w, data, freq, cellSize, C, support, overSample,
    cOffset, grid);
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
int wprojection(const std::vector<Coord>& u,
const std::vector<Coord>& v,
const std::vector<Coord>& w,
std::vector<Value>& data,
const std::vector<Coord>& freq,
const Coord cellSize,
const Coord baseline,
const int wSize,
std::vector<Value>& grid)
{

  const int nSamples = u.size();
  const int nChan = freq.size();

  cout << "************************* W projection gridding *********************"
    << endl;
  int support=static_cast<int>(1.5*sqrt(abs(baseline)*static_cast<Coord>(cellSize)*freq[0])/cellSize);
  int overSample=8;
  cout << "Support = " << support << " pixels" << endl;
  const Coord wCellSize=2*baseline*freq[0]/wSize;
  cout << "W cellsize = " << wCellSize << " wavelengths" << endl;

// Convolution function. This should be the convolution of the
// w projection kernel (the Fresnel term) with the convolution
// function used in the standard case. The latter is needed to
// suppress aliasing. In practice, we calculate entire function
// by Fourier transformation. Here we take an approximation that
// is good enough.
  int cSize=2*(support+1)*overSample+1;

  int cCenter=(cSize-1)/2;

  std::vector<Real> C(cSize*cSize*wSize);

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
          long int cind=i+cSize*(j+cSize*k);
          C[cind]=static_cast<Real>(std::cos(r2/(w*fScale)));
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
          long int cind=i+cCenter+cSize*(j+cCenter+cSize*k);
          C[cind]=static_cast<Real>(std::exp(-r2));
        }
      }
    }
  }

  std::vector<unsigned int> cOffset(data.size());
  for (int i=0;i<nSamples;i++)
  {
    for (int chan=0;chan<nChan;chan++)
    {

      int find=i*nChan+chan;

      Coord wScaled=freq[chan]*w[i]/wCellSize;
      cOffset[find]=wSize/2+int(wScaled);
    }
  }

  return generic(u, v, w, data, freq, cellSize, C, support, overSample,
    cOffset, grid);
}


int main()
{
  const int baseline=2000;                        // Maximum baseline in meters
  const int nSamples=100000;                      // Number of data samples
  const int gSize=512;                            // Size of output grid in pixels
  const Coord cellSize=50;                        // Cellsize of output grid in wavelengths
  const int wSize=64;                             // Number of lookup planes in w projection
  const int nChan=16;                             // Number of spectral channels

// Initialize the data to be gridded
  std::vector<Coord> u(nSamples);
  std::vector<Coord> v(nSamples);
  std::vector<Coord> w(nSamples);
  std::vector<Value> data(nSamples*nChan);

  for (int i=0;i<nSamples;i++)
  {
    u[i]=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;
    v[i]=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;
    w[i]=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;
    for (int chan=0;chan<nChan;chan++)
    {
      data[i*nChan+chan]=Coord(rand())/Coord(RAND_MAX);
    }
  }

// Measure frequency in inverse wavelengths
  std::vector<Coord> freq(nChan);
  for (int i=0;i<nChan;i++)
  {
    freq[i]=(1.4e9-2.0e5*Coord(i)/Coord(nChan))/2.998e8;
  }

  std::vector<Value> grid(gSize*gSize);

  standard(u, v, w, data, freq, cellSize, grid);

  wprojection(u, v, w, data, freq, cellSize, baseline, wSize, grid);

  cout << "Done" << endl;

  return 0;
}
