/// @file
/// @brief Helper class to analyse performance for CF generation
/// @details This template is derived from the gridder class
/// under test. The main method forces contunuous recalculation of
/// convolution functions for the selected maximum number of beams and
/// w-planes. It is therefore possible to profile this core operation more
/// accurately. 
///
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>
///

#include <casa/Quanta/MVDirection.h>
#include <casa/Quanta/Quantum.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.testcfgenperformance");

#include <askap/AskapError.h>
#include <gridding/TestCFGenPerformance.h>
#include <gridding/AProjectGridderBase.h>


using namespace askap;
using namespace askap::synthesis;


/// @brief static method to create test gridder
/// @details This method mimics createGridder interface for a-projection gridders
/// @param[in] parset input parset file
/// @return a shared pointer to the instance
boost::shared_ptr<TestCFGenPerformance> TestCFGenPerformance::createGridder(const LOFAR::ParameterSet& parset)
{
   boost::shared_ptr<TestCFGenPerformance> gridder = createAProjectGridder<TestCFGenPerformance>(parset);
   gridder->configureGridder(parset);

   return gridder;
}

/// @brief clone object, e.g. to test performance with a number of copies
boost::shared_ptr<TestCFGenPerformance> TestCFGenPerformance::clone() const
{
   boost::shared_ptr<TestCFGenPerformance> result(new TestCFGenPerformance(*this));
   return result;
}

/// @brief initialise the accessor
void TestCFGenPerformance::init()
{
  ASKAPASSERT(itsNBeams > 0);
  ASKAPLOG_INFO_STR(logger, "Preparing a dummy accessor");
  
  // some hard-coded beam arrangement and pointing. We could've provided more flexibility if it is found useful.  
  const double maxSeparationInRad = 0.087;
  const casa::uInt nAnt = 36;
  //const double maxSeparationInRad = 0.0087;
  //const casa::uInt nAnt = 6;
  const casa::MVDirection dishPointing = getTangentPoint();
  casa::Vector<casa::MVDirection> beamPointings(itsNBeams, dishPointing);
  int nBeamsOnEachSide = int(sqrt(double(itsNBeams)));
  if (nBeamsOnEachSide == 0) {
      nBeamsOnEachSide = itsNBeams;
  }
  for (casa::uInt beam=0; beam < beamPointings.nelements(); ++beam) {
       const double length = double(beam % nBeamsOnEachSide) / double(nBeamsOnEachSide) * maxSeparationInRad;
       const double angle = 2.*casa::C::pi * double(beam / nBeamsOnEachSide) / double(nBeamsOnEachSide);
       beamPointings[beam].shift(length*cos(angle), length*sin(angle));
  } 
  const casa::uInt nSamples = nAnt * (nAnt - 1) / 2 * beamPointings.nelements();
  itsAccessor.itsAntenna1.resize(nSamples);
  itsAccessor.itsAntenna2.resize(nSamples);
  itsAccessor.itsFeed1.resize(nSamples);
  itsAccessor.itsFeed2.resize(nSamples);
  itsAccessor.itsFeed1PA.resize(nSamples);
  itsAccessor.itsFeed1PA.set(0.);
  itsAccessor.itsFeed2PA.resize(nSamples);
  itsAccessor.itsFeed2PA.set(0.);
  itsAccessor.itsPointingDir1.resize(nSamples);
  itsAccessor.itsPointingDir2.resize(nSamples);
  itsAccessor.itsDishPointing1.resize(nSamples);
  itsAccessor.itsDishPointing1.set(dishPointing);
  itsAccessor.itsDishPointing2.resize(nSamples);
  itsAccessor.itsDishPointing2.set(dishPointing);
  itsAccessor.itsVisibility.resize(nSamples,1,1);
  itsAccessor.itsVisibility.set(casa::Complex(0.));
  itsAccessor.itsFlag.resize(nSamples,1,1);
  itsAccessor.itsFlag.set(false);
  itsAccessor.itsNoise.resize(nSamples,1,1);
  itsAccessor.itsNoise.set(casa::Complex(1.,1.));
  itsAccessor.itsTime = 0.;
  itsAccessor.itsFrequency.resize(1);
  itsAccessor.itsFrequency.set(1.4e9);
  itsAccessor.itsUVW.resize(nSamples);
  itsAccessor.itsUVWRotationDelay.resize(nSamples);
  itsAccessor.itsStokes.resize(1);
  itsAccessor.itsStokes.set(casa::Stokes::XX);
  for (casa::uInt ant1 = 0, index = 0; ant1<nAnt; ++ant1) {
       for (casa::uInt ant2 = 0; ant2<ant1; ++ant2) {
            for (casa::uInt beam=0; beam < beamPointings.nelements(); ++beam,++index) {
                 ASKAPDEBUGASSERT(index < nSamples);
                 itsAccessor.itsAntenna1[index] = ant1;
                 itsAccessor.itsAntenna2[index] = ant2;
                 itsAccessor.itsFeed1[index] = beam;
                 itsAccessor.itsFeed2[index] = beam;
                 itsAccessor.itsPointingDir1[index] = beamPointings[beam];
                 itsAccessor.itsPointingDir2[index] = beamPointings[beam];                 
            }
       }
  }
  ASKAPLOG_INFO_STR(logger, "The accessor stub has been initialised for "<<nAnt<<" antennas and "<<beamPointings.nelements()<<" beams");  
}

/// @brief Initialise the parameters and accessor
/// @param axes axes specifications
/// @param shape Shape of output image: u,v,pol,chan
/// @param dopsf Make the psf?
void TestCFGenPerformance::initialiseGrid(const scimath::Axes& axes,
               const casa::IPosition& shape, const bool dopsf)
{
  AWProjectVisGridder::initialiseGrid(axes,shape,dopsf);
  init();
}

/// @brief main method which initialises CFs
/// @details 
/// @param[in] nRuns number of initialisations
void TestCFGenPerformance::run(const int nRuns)
{
  for (int i=0; i<nRuns; ++i) {
       ASKAPLOG_INFO_STR(logger, "CF generation run "<<(i+1));
       initIndices(itsAccessor);
       initConvolutionFunction(itsAccessor);
       // force recalculation
       resetCFCache();
  }
}





