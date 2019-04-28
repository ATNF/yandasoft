//   This C++ program has been written to demonstrate the convolutional resampling algorithm used in radio
// interferometry. It should compile with:
//      g++ -O2 -fstrict-aliasing tConvolveBLAS.cc -o tConvolveBLAS
// Enabling BLAS support on OS X:
//      g++ -DUSEBLAS -O2 -fstrict-aliasing -framework vecLib tConvolveBLAS.cc -o tConvolveBLAS
//
// Strict-aliasing tells the compiler that there are no memory locations accessed through aliases.
//
// The challenge is to minimize the run time - specifically the time per grid addition. On a MacBookPro 
// 2GHz Intel Core Duo this is about 6.0ns. 
//
// For further details contact Tim.Cornwell@csiro.au
// November 22, 2007
// - Rewritten from tConvolve to use BLAS, and to be much smarter about not using strides in C

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

#ifdef USEBLAS

#ifdef __APPLE_CC__
#include <vecLib/cblas.h>
#endif

#endif

using std::cout;
using std::endl;
using std::complex;
using std::abs;

// Typedefs for easy testing
// Cost of using double for Coord is low, cost for
// double for Real is also low
typedef double Coord;
typedef float Real;
typedef std::complex<Real> Value;

/////////////////////////////////////////////////////////////////////////////////
// The next two functions are the kernel of the gridding/degridding.
// The data are presented as a vector. Offsets for the convolution function
// and for the grid location are precalculated so that the kernel does
// not need to know anything about world coordinates or the shape of
// the convolution function. The ordering of cOffset and iu, iv is
// random - some presorting might be advantageous.
//
// Perform gridding
//
// data - values to be gridded in a 1D vector
// support - Total width of convolution function=2*support+1
// C - convolution function shape: (2*support+1, 2*support+1, *)
// cOffset - offset into convolution function per data point
// iu, iv - integer locations of grid points
// grid - Output grid: shape (gSize, *)
// gSize - size of one axis of grid

void residKernel(const std::vector<Value>& obsData,
    std::vector<Value>& modelData, const int support,
    const std::vector<Value>& C, const std::vector<unsigned int>& cOffset,
    const std::vector<unsigned int>& iu, const std::vector<unsigned int>& iv,
    const std::vector<Value>& modelGrid, std::vector<Value>& residGrid,
    const int gSize)
{

  int sSize=2*support+1;
  for (unsigned int dind=0; dind<obsData.size(); dind++)
  {
    modelData[dind]=0.0;

    // Nearly all the L2 cache misses originate here in the next
    // two statements
    // The actual grid point from which we offset
    int gind=iu[dind]+gSize*iv[dind]-support;
    // The Convoluton function point from which we offset
    int cind=cOffset[dind];

    for (int suppv=0; suppv<sSize; suppv++)
    {
#ifdef USEBLAS
      Value dot;
      cblas_cdotu_sub(sSize, &modelGrid[gind], 1, &C[cind], 1, &dot);
      modelData[dind]+=dot;
#else
      for (int suppu=0; suppu<sSize; suppu++)
      {
        modelData[dind]+=modelGrid[gind+suppu]*C[cind+suppu];
      }
#endif
      gind+=gSize;
      cind+=sSize;
    }

    Value residData=obsData[dind]-modelData[dind];

    gind=iu[dind]+gSize*iv[dind]-support;
    cind=cOffset[dind];

    for (int suppv=0; suppv<sSize; suppv++)
    {
#ifdef USEBLAS
      cblas_caxpy(sSize, &residData, &C[cind], 1, &residGrid[gind], 1);
#else
      for (int suppu=0; suppu<sSize; suppu++)
      {
        residGrid[gind+suppu]+=residData*C[cind+suppu];
      }
#endif
      gind+=gSize;
      cind+=sSize;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////

// Initialize W project convolution function 
// - This is application specific and should not need any changes.
//
// nSamples - number of visibility samples
// freq - temporal frequency (inverse wavelengths)
// cellSize - size of one grid cell in wavelengths
// gSize - size of grid in pixels (per axis)
// support - Total width of convolution function=2*support+1
// wCellSize - size of one w grid cell in wavelengths
// wSize - Size of lookup table in w
void initC(const int nSamples, const std::vector<Coord>& w,
    const std::vector<Coord>& freq, const Coord cellSize, const Coord baseline,
    const int wSize, const int gSize, int& support, int& overSample,
    Coord& wCellSize, std::vector<Value>& C)
{

  cout << "Initializing W projection convolution function" << endl;
  support=static_cast<int>(1.5*sqrt(abs(baseline) *static_cast<Coord>(cellSize)
      *freq[0])/cellSize);
  overSample=8;
  cout << "Support = " << support << " pixels" << endl;
  wCellSize=2*baseline*freq[0]/wSize;
  cout << "W cellsize = " << wCellSize << " wavelengths" << endl;

  // Convolution function. This should be the convolution of the
  // w projection kernel (the Fresnel term) with the convolution
  // function used in the standard case. The latter is needed to
  // suppress aliasing. In practice, we calculate entire function
  // by Fourier transformation. Here we take an approximation that
  // is good enough.
  int sSize=2*support+1;

  int cCenter=(sSize-1)/2;

  C.resize(sSize*sSize*overSample*overSample*wSize);
  cout << "Size of convolution function = " << sSize*sSize*overSample
      *overSample*wSize*8/(1024*1024) << " MB" << std::endl;
  cout << "Shape of convolution function = [" << sSize << ", " << sSize << ", "
      << overSample << ", " << overSample << ", " << wSize << "]" << std::endl;

  for (int k=0; k<wSize; k++)
  {
    double w=double(k-wSize/2);
    double fScale=sqrt(abs(w)*wCellSize*freq[0])/cellSize;
    for (int osj=0; osj<overSample; osj++)
    {
      for (int osi=0; osi<overSample; osi++)
      {
        for (int j=0; j<sSize; j++)
        {
          double j2=std::pow((double(j-cCenter)+double(osj)/double(overSample)), 2);
          for (int i=0; i<sSize; i++)
          {
            double r2=j2+std::pow((double(i-cCenter)+double(osi)/double(overSample)), 2);
            long int cind=i+sSize*(j+sSize*(osi+overSample*(osj+overSample*k)));
            if (w!=0.0)
            {
              C[cind]=static_cast<Value>(std::cos(r2/(w*fScale)));
            }
            else
            {
              C[cind]=static_cast<Value>(std::exp(-r2));
            }
          }
        }
      }
    }
  }

  // Now normalise the convolution function
  Real sumC=0.0;
  for (int i=0; i<sSize*sSize*overSample*overSample*wSize; i++)
  {
    sumC+=abs(C[i]);
  }

  for (int i=0; i<sSize*sSize*overSample*overSample*wSize; i++)
  {
    C[i]*=Value(wSize*overSample*overSample/sumC);
  }
}
// Initialize Lookup function
// - This is application specific and should not need any changes.
//
// nSamples - number of visibility samples
// freq - temporal frequency (inverse wavelengths)
// cellSize - size of one grid cell in wavelengths
// gSize - size of grid in pixels (per axis)
// support - Total width of convolution function=2*support+1
// wCellSize - size of one w grid cell in wavelengths
// wSize - Size of lookup table in w
void initCOffset(const std::vector<Coord>& u, const std::vector<Coord>& v,
    const std::vector<Coord>& w, const std::vector<Coord>& freq,
    const Coord cellSize, const Coord wCellSize, const Coord baseline,
    const int wSize, const int gSize, const int support, const int overSample,
    std::vector<unsigned int>& cOffset, std::vector<unsigned int>& iu,
    std::vector<unsigned int>& iv)
{

  const int nSamples = u.size();
  const int nChan = freq.size();

  int sSize=2*support+1;

  // Now calculate the offset for each visibility point
  cOffset.resize(nSamples*nChan);
  iu.resize(nSamples*nChan);
  iv.resize(nSamples*nChan);
  for (int i=0; i<nSamples; i++)
  {
    for (int chan=0; chan<nChan; chan++)
    {

      int dind=i*nChan+chan;

      Coord uScaled=freq[chan]*u[i]/cellSize;
      iu[dind]=int(uScaled);
      if (uScaled<Coord(iu[dind]))
      {
        iu[dind]-=1;
      }
      int fracu=int(overSample*(uScaled-Coord(iu[dind])));
      iu[dind]+=gSize/2;

      Coord vScaled=freq[chan]*v[i]/cellSize;
      iv[dind]=int(vScaled);
      if (vScaled<Coord(iv[dind]))
      {
        iv[dind]-=1;
      }
      int fracv=int(overSample*(vScaled-Coord(iv[dind])));
      iv[dind]+=gSize/2;

      // The beginning of the convolution function for this point
      Coord wScaled=freq[chan]*w[i]/wCellSize;
      int woff=wSize/2+int(wScaled);
      cOffset[dind]=sSize*sSize*(fracu+overSample*(fracv+overSample*woff));
    }
  }

}

// Main testing routine
int main()
{
  // Change these if necessary to adjust run time
  const int nSamples=10000; // Number of data samples
  const int wSize=33; // Number of lookup planes in w projection
  const int nChan=16; // Number of spectral channels

  // Don't change any of these numbers unless you know what you are doing!
  const int gSize=512; // Size of output grid in pixels
  const Coord cellSize=40.0; // Cellsize of output grid in wavelengths
  const int baseline=2000; // Maximum baseline in meters

  // Initialize the data to be gridded
  std::vector<Coord> u(nSamples);
  std::vector<Coord> v(nSamples);
  std::vector<Coord> w(nSamples);
  std::vector<Value> data(nSamples*nChan);
  std::vector<Value> outdata(nSamples*nChan);

  for (int i=0; i<nSamples; i++)
  {
    u[i]=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;
    v[i]=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;
    w[i]=baseline*Coord(rand())/Coord(RAND_MAX)-baseline/2;
    for (int chan=0; chan<nChan; chan++)
    {
      data[i*nChan+chan]=1.0;
      outdata[i*nChan+chan]=0.0;
    }
  }

  std::vector<Value> modelGrid(gSize*gSize);
  modelGrid.assign(modelGrid.size(), Value(1.0));

  std::vector<Value> residGrid(gSize*gSize);
  residGrid.assign(residGrid.size(), Value(1.0));

  // Measure frequency in inverse wavelengths
  std::vector<Coord> freq(nChan);
  for (int i=0; i<nChan; i++)
  {
    freq[i]=(1.4e9-2.0e5*Coord(i)/Coord(nChan))/2.998e8;
  }

  // Initialize convolution function and offsets
  std::vector<std::complex<float> > C;
  int support, overSample;
  std::vector<unsigned int> cOffset;
  // Vectors of grid centers
  std::vector<unsigned int> iu;
  std::vector<unsigned int> iv;
  Coord wCellSize;

  initC(nSamples, w, freq, cellSize, baseline, wSize, gSize, support,
      overSample, wCellSize, C);
  initCOffset(u, v, w, freq, cellSize, wCellSize, baseline, wSize, gSize,
      support, overSample, cOffset, iu, iv);
  int sSize=2*support+1;

  // Now we can do the timing
  cout << "+++++ Residual processing +++++" << endl;

  clock_t start, finish;
  double time;

  start = clock();
  residKernel(data, outdata, support, C, cOffset, iu, iv, modelGrid, residGrid,
      gSize);
  finish = clock();

  // Report on timings
  time = (double(finish)-double(start))/CLOCKS_PER_SEC;
  cout << "    Time " << time << " (s) " << endl;
  cout << "    Time per visibility spectral sample " << 1e6*time/double(data.size()) << " (us) " << endl;
  cout << "    Time per gridding   " << 1e9*time/(double(data.size())* double((sSize)*(sSize))) << " (ns) " << endl;

  return 0;
}
