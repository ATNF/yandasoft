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

#include <measurementequation/WienerPreconditioner.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.wienerpreconditioner");

#include <askap/AskapError.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <utils/PaddingUtils.h>
#include <profile/AskapProfiler.h>

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/Vector.h>
#include <lattices/Lattices/SubLattice.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/Lattices/LatticeFFT.h>
#include <lattices/Lattices/LatticeExpr.h>
using namespace casa;

#include <iostream>
#include <cmath>
using std::abs;

namespace askap
{
  namespace synthesis
  {

    WienerPreconditioner::WienerPreconditioner() :
	    itsParameter(0.0), itsDoNormalise(false), itsUseRobustness(false)
    {
    }
    
    /// @brief constructor with explicitly defined noise power
    /// @param[in] noisepower parameter of the
    /// @param[in] normalise if true, PSF is normalised during filter construction
    WienerPreconditioner::WienerPreconditioner(float noisepower, bool normalise) : 
            itsParameter(noisepower), itsDoNormalise(normalise), itsUseRobustness(false) {}

    /// @brief constructor with explicitly defined robustness
    /// @details In this version, the noise power is calculated from
    /// the robustness parameter
    /// @param[in] robustness robustness parameter (roughly matching Briggs' weighting)
    /// @note Normalisation of PSF is always used when noise power is defined via robustness
    WienerPreconditioner::WienerPreconditioner(float robustness) : itsParameter(robustness), 
            itsDoNormalise(true), itsUseRobustness(true)  {}
        
        
    /// @brief copy constructor
    /// @param[in] other an opject to copy from
    WienerPreconditioner::WienerPreconditioner(const WienerPreconditioner &other) :
          IImagePreconditioner(other),
          itsParameter(other.itsParameter), itsDoNormalise(other.itsDoNormalise), 
          itsUseRobustness(other.itsUseRobustness) 
    {
       if (other.itsTaperCache) {
           itsTaperCache.reset(new GaussianTaperCache(*(other.itsTaperCache)));
       }
    }      
        
    IImagePreconditioner::ShPtr WienerPreconditioner::clone()
    {
	    return IImagePreconditioner::ShPtr(new WienerPreconditioner(*this));
    }
    
    bool WienerPreconditioner::doPreconditioning(casa::Array<float>& psf, casa::Array<float>& dirty) const
    {
      ASKAPTRACE("WienerPreconditioner::doPreconditioning");
      if (!itsUseRobustness && (itsParameter < 1e-6)) {
          return false;
      }

      ASKAPCHECK(psf.shape().conform(dirty.shape()), "Dirty image and PSF do not conform - shapes: " << dirty.shape() << psf.shape());

      if (itsUseRobustness) {
          ASKAPLOG_INFO_STR(logger, "Applying Wiener filter with noise power defined via robustness=" << itsParameter);
      } else {
          ASKAPLOG_INFO_STR(logger, "Applying Wiener filter with noise power=" << itsParameter);
      }

      float maxPSFBefore = casa::max(psf);
      ASKAPLOG_INFO_STR(logger, "Peak of PSF before Wiener filtering = " << maxPSFBefore);

      if (itsDoNormalise) {
          ASKAPLOG_INFO_STR(logger, "The PSF will be normalised to 1 before filter construction");
	  psf=psf/maxPSFBefore;
	  //	  dirty=dirty/maxPSFBefore;
	  maxPSFBefore=1.0;
      }
                  
      casa::ArrayLattice<float> lpsf(psf);      
      casa::ArrayLattice<float> ldirty(dirty);

      const casa::IPosition shape = lpsf.shape();

      // Make the scratch array into which we will calculate the Wiener filter
      casa::ArrayLattice<casa::Complex> scratch(shape);
      scratch.copyData(casa::LatticeExpr<casa::Complex>(toComplex(lpsf)));
      if (itsTaperCache) {
          ASKAPLOG_INFO_STR(logger, "Applying Gaussian taper to the Wiener filter in the image domain");
          casa::Array<casa::Complex> taperArray(itsTaperCache->taper(shape));
          casa::ArrayLattice<casa::Complex> taperLattice(taperArray);
          scratch.copyData(casa::LatticeExpr<casa::Complex>(scratch * taperLattice));
      }
      LatticeFFT::cfft2d(scratch, True);
       
      // Make the transfer function
      casa::ArrayLattice<casa::Complex> xfr(shape);
      xfr.copyData(casa::LatticeExpr<casa::Complex>(toComplex(lpsf)));
      LatticeFFT::cfft2d(xfr, True);
       
      // Calculate the Wiener filter
      casa::ArrayLattice<casa::Complex> wienerfilter(shape);
      const float normFactor = itsDoNormalise ? maxPSFBefore : 1.;
      const float noisePower = (itsUseRobustness ? std::pow(10., 4.*itsParameter) : itsParameter)*normFactor*normFactor;
      ASKAPLOG_INFO_STR(logger, "Effective noise power of the Wiener filter = " << noisePower);     
      wienerfilter.copyData(casa::LatticeExpr<casa::Complex>(normFactor*conj(scratch)/(real(scratch*conj(scratch)) + noisePower)));
      
      /*
      // for debugging - to export Wiener filter
      LatticeFFT::cfft2d(wienerfilter, False);       
      lpsf.copyData(casa::LatticeExpr<float>(real(wienerfilter*conj(wienerfilter))));
      SynthesisParamsHelper::saveAsCasaImage("dbg.img",psf);
      throw 1;
      */
      
      // Apply the Wiener filter to the xfr and transform to the filtered PSF
      scratch.copyData(casa::LatticeExpr<casa::Complex> (wienerfilter * xfr));
      LatticeFFT::cfft2d(scratch, False);       
      lpsf.copyData(casa::LatticeExpr<float>(real(scratch)));
      const float maxPSFAfter=casa::max(psf);
      ASKAPLOG_INFO_STR(logger, "Peak of PSF after Wiener filtering  = " << maxPSFAfter); 
      psf *= maxPSFBefore/maxPSFAfter;
      ASKAPLOG_INFO_STR(logger, "Normalized to unit peak");
      
      // Apply the filter to the dirty image
      scratch.copyData(casa::LatticeExpr<casa::Complex>(toComplex(ldirty)));       
      LatticeFFT::cfft2d(scratch, True);
      scratch.copyData(casa::LatticeExpr<casa::Complex> (wienerfilter * scratch));
      LatticeFFT::cfft2d(scratch, False);

      ldirty.copyData(casa::LatticeExpr<float>(real(scratch)));
      dirty *= maxPSFBefore/maxPSFAfter;
	  
      return true;

    }

    /// @brief static factory method to create preconditioner from a parset
    /// @details
    /// @param[in] parset subset of parset file (with preconditioner.Wiener. removed)
    /// @return shared pointer
    boost::shared_ptr<WienerPreconditioner> WienerPreconditioner::createPreconditioner(const LOFAR::ParameterSet &parset) 
    {
      ASKAPCHECK(parset.isDefined("noisepower") != parset.isDefined("robustness"), 
           "Exactly one parameter, either noisepower or robustness parameter must be given. You gave either none or both of them.");

      boost::shared_ptr<WienerPreconditioner> result;
      if (parset.isDefined("noisepower")) {
          const float noisepower = parset.getFloat("noisepower");
          const bool normalise = parset.getBool("normalise",false);
          result.reset(new WienerPreconditioner(noisepower,normalise));
      } else {
      
          ASKAPDEBUGASSERT(parset.isDefined("robustness"));
     
          const float robustness = parset.getFloat("robustness");
          ASKAPCHECK((robustness >= -2.00001) && (robustness <= 2.0001), 
                     "Robustness parameter is supposed to be between -2 and 2, you have = "<<robustness);
          ASKAPCHECK(!parset.isDefined("normalise"), 
                     "Normalise option of the Wiener preconditioner is not compatible with the "
                     "preconditioner definition via robustness (as normalisation of PSF is always done in this case)");
          result.reset(new WienerPreconditioner(robustness));
      }
      ASKAPASSERT(result);
      // configure tapering
      if (parset.isDefined("taper")) {
          const double fwhm = parset.getDouble("taper");
          result->enableTapering(fwhm);
      }
      //
      
      return result;
    }

    /// @brief assignment operator, to ensure it is not called
    WienerPreconditioner& WienerPreconditioner::operator=(const WienerPreconditioner &) 
    {
      ASKAPTHROW(AskapError, "Assignment operator is not supposed to be used");
      return *this;
    }
    
    /// @brief enable Filter tapering 
    /// @details Wiener filter can optionally be tapered in the image domain, so it is not extended over
    /// the whole field of view.
    /// @param[in] fwhm full width at half maximum of the taper given in image cells
    void WienerPreconditioner::enableTapering(double fwhm)
    {
      ASKAPLOG_INFO_STR(logger, "Wiener filter will be tapered by a circular Gaussian with FWHM="<<fwhm<<
                        " pixels in the image plane");
      itsTaperCache.reset(new GaussianTaperCache(fwhm));
    }    

  }
}


