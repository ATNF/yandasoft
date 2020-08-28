/// @file
/// 
/// @brief Calibration effect: antenna gains without cross-pol
/// @details This is a simple effect which can be used in conjunction
/// with the CalibrationME template (as its template argument)
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

#ifndef IONOSPHERIC_TERM_TCC
#define IONOSPHERIC_TERM_TCC

// own includes
#include <askap/scimath/utils/PolConverter.h>
#include <askap/calibaccess/CalParamNameHelper.h>

// casa includes
#include <casacore/measures/Measures/Stokes.h>

// std includes
#include <string>
#include <utility>

namespace askap {

namespace synthesis {

    using utility::toString;

/// @brief main method returning Mueller matrix and derivatives
/// @details This method has to be overloaded (in the template sense) for
/// all classes representing various calibration effects. CalibrationME
/// template will call it when necessary. It returns 
/// @param[in] chunk accessor to work with
/// @param[in] row row of the chunk to work with
/// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
/// this effect
inline scimath::ComplexDiffMatrix IonosphericTerm::get(const accessors::IConstDataAccessor &chunk, 
                                      casacore::uInt row) const
{
   const casacore::uInt nPol = chunk.nPol();
   ASKAPDEBUGASSERT(nPol != 0);   
   const casacore::uInt nChan = chunk.nChannel();
   const casacore::Vector<casacore::Stokes::StokesTypes> stokes = chunk.stokes();   
   ASKAPDEBUGASSERT(stokes.nelements() == nPol);
   ASKAPDEBUGASSERT(!scimath::PolConverter::isStokes(stokes));
   
   const casacore::uInt ant1 = chunk.antenna1()[row];
   const casacore::uInt ant2 = chunk.antenna2()[row];
   
   const casacore::uInt dir1 = chunk.feed1()[row];
   const casacore::uInt dir2 = chunk.feed2()[row];

   // antenna positions for ant1
   const std::string u1name = "antenna.position.x."+toString(ant1);
   const std::string v1name = "antenna.position.y."+toString(ant1);

   // antenna positions for ant2
   const std::string u2name = "antenna.position.x."+toString(ant2);
   const std::string v2name = "antenna.position.y."+toString(ant2);

   // ionospheric parameters: dl/lambda^2 and dm/lambda_2
   const std::string lname = "ionosphere.coeff.l."+toString(dir1);
   const std::string mname = "ionosphere.coeff.m."+toString(dir2);

   const casacore::Complex u = parameters()->scalarValue(u1name) - parameters()->scalarValue(u2name);
   const casacore::Complex v = parameters()->scalarValue(v1name) - parameters()->scalarValue(v2name);
   const scimath::ComplexDiff phase_const = -casacore::C::_2pi * ( u*getParameter(lname) + v*getParameter(mname) );

   const casacore::Vector<casacore::Double>& frequency = parameters()->value("frequency");

   scimath::ComplexDiffMatrix calFactor(nPol, nPol * nChan, 0.);

   for (casacore::uInt chan = 0; chan < nChan; ++chan) {

       const casacore::Float lambda = casacore::C::c / frequency[chan];
       const askap::scimath::ComplexDiff phase = phase_const * lambda*lambda;

       for (casacore::uInt pol=0; pol<nPol; ++pol) {
       
            const casacore::uInt polIndex = scimath::PolConverter::getIndex(stokes[pol]);
            // polIndex is index in the polarisation frame, i.e.
            // XX is 0, XY is 1, YX is 2 and YY is 3
            // we need an index into matrix 
         
            // gains for antenna 1, polarisation X if XX or XY, or Y if YX or YY
            const std::string g1name = accessors::CalParamNameHelper::paramName(ant1, dir1,
                          polIndex / 2 == 0 ? casacore::Stokes::XX : casacore::Stokes::YY);
         
            // gains for antenna 2, polarisation X if XX or YX, or Y if XY or YY
            const std::string g2name = accessors::CalParamNameHelper::paramName(ant2, dir2,
                          polIndex % 2 == 0 ? casacore::Stokes::XX : casacore::Stokes::YY);
         
            /*
            if (row == 0 && pol == 0) {
                std::cout << "IONOCAL IonosphericTerm (npol = "<<nPol<<"): row = "<<row<<
                                 ", ant1 = "<<ant1<<", ant2 = "<<ant2<<
                                 ", dir1 = "<<dir1<<", dir2 = "<<dir2<<
                                 ", exp(-i*2pi * ("<<u<<"*"<<getParameter(lname).value()<<
                                             " +  "<<v<<"*"<<getParameter(mname).value()<<")" << std::endl;
            }
            */

            // we need to think of how to deal with distributed problem on the cluster (i.e. adding
            // some base to the channel number and propagating it through the framework)    
            // exp( j * phase ) ~ 1 + j * phase
            // could also include gains as free ComplexDiff parameters here, but just use their values for now
            calFactor(pol, pol + chan * nPol) = getParameter(g1name).value() * conj(getParameter(g2name).value()) *
                                                ( 1.0 + casacore::Complex(0.0,1.0) * phase );

       }


   }

   return calFactor;
}

inline casacore::Vector<casacore::RigidVector<double,3> >
    IonosphericTerm::getAntennaPositions(const accessors::IConstDataAccessor &acc) const
{
   // do some checks first?

   // rotate to zenith and use acc.rotatedUVW() instead of acc.uvw()?

   // set positions relative to station zero
   casacore::Vector<casacore::RigidVector<double,3> > xyz(1);
   xyz(0)(0) = 0.0;
   xyz(0)(1) = 0.0;
   xyz(0)(2) = 0.0;
   casacore::uInt nAnt = 1;
   casacore::RigidVector<double,3> xyzAve(0.0);
   for (casacore::uInt row=0; row<acc.nRow(); ++row) {
       if (acc.antenna1()[row] > 0) break;
       // don't really need to be this strict, but do so for now
       ASKAPASSERT(acc.antenna1()[row] == 0);
       ASKAPASSERT(acc.antenna2()[row] == row+1);
       const casacore::uInt ant = row+1;
       nAnt++;
       xyz.resize(ant+1,true);
       xyz(ant)(0) = acc.uvw()[row](0);
       xyz(ant)(1) = acc.uvw()[row](1);
       xyz(ant)(2) = acc.uvw()[row](2);
       xyzAve(0) += xyz(ant)(0);
       xyzAve(1) += xyz(ant)(1);
       xyzAve(2) += xyz(ant)(2);
   }

   // subtract average position?
   for (casacore::uInt coord=0; coord<3; ++coord) {
       xyzAve(coord) /= double(nAnt);
   }
   for (casacore::uInt ant=0; ant<nAnt; ++ant) {
       xyz(ant)(0) -= xyzAve(0);
       xyz(ant)(1) -= xyzAve(1);
       xyz(ant)(2) -= xyzAve(2);
   }

   return xyz;
}

} // namespace synthesis

} // namespace askap


#endif // #ifndef IONOSPHERIC_TERM_TCC
