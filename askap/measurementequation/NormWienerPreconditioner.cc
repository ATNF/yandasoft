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

#include <askap/measurementequation/NormWienerPreconditioner.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.normwienerpreconditioner");

#include <askap/askap/AskapError.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/profile/AskapProfiler.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/lattices/Lattices/SubLattice.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#include <casacore/lattices/LEL/LatticeExpr.h>
using namespace casa;

#include <iostream>
#include <cmath>
using std::abs;

namespace askap
{
  namespace synthesis
  {

    NormWienerPreconditioner::NormWienerPreconditioner() :
	    itsRobust(0.0)
    {
    }
    
    NormWienerPreconditioner::NormWienerPreconditioner(const float& noisepower) :
	    itsRobust(noisepower)
    {
    }
        
    IImagePreconditioner::ShPtr NormWienerPreconditioner::clone()
    {
	    return IImagePreconditioner::ShPtr(new NormWienerPreconditioner(*this));
    }
    
    bool NormWienerPreconditioner::doPreconditioning(casacore::Array<float>& psf,
                                                     casacore::Array<float>& dirty,
                                                     casacore::Array<float>& pcf) const
    {
      {
        ASKAPTRACE("NormWienerPreconditioner::doPreconditioning");

	ASKAPLOG_INFO_STR(logger, "Applying Normalised Wiener filter with robustness parameter " << itsRobust);

       float maxPSFBefore=casacore::max(psf);
       ASKAPLOG_INFO_STR(logger, "Peak of PSF before Normalised Wiener filtering = " << maxPSFBefore);
       casacore::ArrayLattice<float> lpsf(psf);
       casacore::ArrayLattice<float> ldirty(dirty);

       const casacore::IPosition shape = lpsf.shape();
       casacore::ArrayLattice<casacore::Complex> scratch(shape);
       //scratch.set(0.);
       scratch.copyData(casacore::LatticeExpr<casacore::Complex>(toComplex(lpsf)));       
       
       LatticeFFT::cfft2d(scratch, True);
       
       // Construct a Wiener filter
       
       casacore::ArrayLattice<casacore::Complex> wienerfilter(shape);
       //wienerfilter.set(0.);
      
      // Normalize relative to the average weight
       const double noisepower(pow(10.0, 2*itsRobust));
       const double np(noisepower*maxPSFBefore);
       wienerfilter.copyData(casacore::LatticeExpr<casacore::Complex>(maxPSFBefore*conj(scratch)/(real(scratch*conj(scratch)) + np*np)));
              
       // Apply the filter to the lpsf
       // (reuse the ft(lpsf) currently held in 'scratch')
       
       // need to rebuild ft(lpsf) with padding, otherwise there is a scaling error
       //scratch.set(0.);
       scratch.copyData(casacore::LatticeExpr<casacore::Complex>(toComplex(lpsf)));       
       LatticeFFT::cfft2d(scratch, True);      
       //
       scratch.copyData(casacore::LatticeExpr<casacore::Complex> (wienerfilter * scratch));
       
       /*
       SynthesisParamsHelper::saveAsCasaImage("dbg.img",casacore::amplitude(scratch.asArray()));       
       //SynthesisParamsHelper::saveAsCasaImage("dbg.img",lpsf.asArray());
       throw AskapError("This is a debug exception");
       */
       
       LatticeFFT::cfft2d(scratch, False);       
       lpsf.copyData(casacore::LatticeExpr<float>(real(scratch)));
       
       float maxPSFAfter=casacore::max(psf);
       ASKAPLOG_INFO_STR(logger, "Peak of PSF after Normalised Wiener filtering  = " << maxPSFAfter); 
       psf*=maxPSFBefore/maxPSFAfter;
       ASKAPLOG_INFO_STR(logger, "Normalized to unit peak");
      
       // Apply the filter to the dirty image
       //scratch.set(0.);
       scratch.copyData(casacore::LatticeExpr<casacore::Complex>(toComplex(ldirty)));       
       
       LatticeFFT::cfft2d(scratch, True);
 
       scratch.copyData(casacore::LatticeExpr<casacore::Complex> (wienerfilter * scratch));
       LatticeFFT::cfft2d(scratch, False);
       
       ldirty.copyData(casacore::LatticeExpr<float>(real(scratch)));
       dirty*=maxPSFBefore/maxPSFAfter;
	  
       return true;
      }

    }

  }
}


