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
ASKAP_LOGGER(logger, ".measurementequation.imagefftequation");

#include <askap/AskapError.h>
//#include <fft/FFTWrapper.h>

#include <dataaccess/SharedIter.h>
#include <dataaccess/MemBufferDataAccessor.h>
#include <fitting/Params.h>
#include <measurementequation/ImageFFTEquation.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <gridding/SphFuncVisGridder.h>
#include <fitting/ImagingNormalEquations.h>
#include <fitting/DesignMatrix.h>
#include <fitting/Axes.h>
#include <profile/AskapProfiler.h>

#include <gridding/SphFuncVisGridder.h>

#include <scimath/Mathematics/RigidVector.h>

#include <casa/BasicSL/Constants.h>
#include <casa/BasicSL/Complex.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>

#include <stdexcept>

using askap::scimath::Params;
using askap::scimath::Axes;
using askap::scimath::ImagingNormalEquations;
using askap::scimath::DesignMatrix;
using namespace askap::accessors;

namespace askap
{
  namespace synthesis
  {

    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params& ip,
        IDataSharedIter& idi) : scimath::Equation(ip),
      askap::scimath::ImagingEquation(ip), itsIdi(idi), itsSphFuncPSFGridder(false)
    {
      itsGridder = IVisGridder::ShPtr(new SphFuncVisGridder());
      init();
    }
    

    ImageFFTEquation::ImageFFTEquation(IDataSharedIter& idi) :
      itsIdi(idi), itsSphFuncPSFGridder(false)
    {
      itsGridder = IVisGridder::ShPtr(new SphFuncVisGridder());
      reference(defaultParameters().clone());
      init();
    }

    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params& ip,
        IDataSharedIter& idi, IVisGridder::ShPtr gridder) : scimath::Equation(ip),
      askap::scimath::ImagingEquation(ip), itsGridder(gridder), itsIdi(idi), itsSphFuncPSFGridder(false)
    {
      init();
    }
    ;

    ImageFFTEquation::ImageFFTEquation(IDataSharedIter& idi,
        IVisGridder::ShPtr gridder) :
      itsGridder(gridder), itsIdi(idi), itsSphFuncPSFGridder(false)
    {
      reference(defaultParameters().clone());
      init();
    }

    ImageFFTEquation::~ImageFFTEquation()
    {
    }
    
    /// @brief define whether the default spheroidal function gridder is used for PSF
    /// @details We have an option to build PSF using the default spheriodal function
    /// gridder, i.e. no w-term and no primary beam is simulated, as an alternative 
    /// to the same user-defined gridder as used for the model. Apart from the speed, 
    /// it probably makes the overall approximation better (i.e. removes some factors 
    /// of spatial dependence).
    /// @param[in] useSphFunc true, if spheroidal function gridder is to be used for PSF
    void ImageFFTEquation::useSphFuncForPSF(bool useSphFunc)
    {
      itsSphFuncPSFGridder = useSphFunc;
      if (itsSphFuncPSFGridder) {
         ASKAPLOG_INFO_STR(logger, "The default spheroidal function gridder will be used to calculate PSF");
      } else {
         ASKAPLOG_INFO_STR(logger, "The PSF will be calculated by the same gridder type as used for the model and residuals");
      }
    }
    

    askap::scimath::Params ImageFFTEquation::defaultParameters()
    {
      Params ip;
      ip.add("image");
      return ip;
    }

    ImageFFTEquation::ImageFFTEquation(const ImageFFTEquation& other) :
          Equation(), ImagingEquation()
    {
      operator=(other);
    }

    ImageFFTEquation& ImageFFTEquation::operator=(const ImageFFTEquation& other)
    {
      if(this!=&other)
      {
        static_cast<askap::scimath::Equation*>(this)->operator=(other);
        itsIdi=other.itsIdi;
        itsGridder = other.itsGridder;
        itsSphFuncPSFGridder = other.itsSphFuncPSFGridder;
        itsVisUpdateObject = other.itsVisUpdateObject;
      }
      return *this;
    }

    void ImageFFTEquation::init()
    {
    }
    
    /// Clone this into a shared pointer
    /// @return shared pointer to a copy
    ImageFFTEquation::ShPtr ImageFFTEquation::clone() const
    {
      return ImageFFTEquation::ShPtr(new ImageFFTEquation(*this));
    }

    /// @brief setup object function to update degridded visibilities
    /// @details For the parallel implementation of the measurement equation we need
    /// inter-rank communication. To avoid introducing cross-dependency of the measurement
    /// equation and the MPI one can use polymorphic object function to sum degridded visibilities 
    /// across all required ranks in the distributed case and do nothing otherwise.
    /// By default, this class doesn't alter degridded visibilities.
    /// @param[in] obj new object function (or an empty shared pointer to turn this option off)
    void ImageFFTEquation::setVisUpdateObject(const boost::shared_ptr<IVisCubeUpdate> &obj)
    {
      itsVisUpdateObject = obj;
    }
    
    /// @brief helper method to verify whether a parameter had been changed 
    /// @details This method checks whether a particular parameter is tracked. If 
    /// yes, its change monitor is used to verify the status since the last call of
    /// the method, otherwise new tracking begins and true is returned (i.e. to 
    /// update all dependent cache).
    /// @param[in] name name of the parameter
    /// @return true if parameter has been updated since the previous call
    bool ImageFFTEquation::notYetDegridded(const std::string &name) const
    {
      std::map<std::string,scimath::ChangeMonitor>::iterator it = itsImageChangeMonitors.find(name);
      bool result = true;
      if (it != itsImageChangeMonitors.end()) {
          result = parameters().isChanged(name,it->second);
          it->second = parameters().monitorChanges(name);
      } else {
          itsImageChangeMonitors[name] = parameters().monitorChanges(name);
      }
      return result;
    }
    

    void ImageFFTEquation::predict() const
    {
      ASKAPTRACE("ImageFFTEquation::predict");
      const vector<string> completions(parameters().completions("image"));

      // To minimize the number of data passes, we keep copies of the gridders in memory, and
      // switch between these. This optimization may not be sufficient in the long run.

      itsIdi.chooseOriginal();
      ASKAPLOG_DEBUG_STR(logger, "Initialising for model degridding");
      for (vector<string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        string imageName("image"+(*it));
        SynthesisParamsHelper::clipImage(parameters(),imageName);

        if(itsModelGridders.count(imageName)==0) {
          itsModelGridders[imageName]=itsGridder->clone();
        }
	    itsModelGridders[imageName]->customiseForContext(*it);

        if (notYetDegridded(imageName)) {
            ASKAPLOG_DEBUG_STR(logger, "Degridding image "<<imageName);
            const Axes axes(parameters().axes(imageName));
            casa::Array<double> imagePixels(parameters().value(imageName).copy());
            const casa::IPosition imageShape(imagePixels.shape());
            itsModelGridders[imageName]->initialiseDegrid(axes, imagePixels);
        }              
      }
      // Loop through degridding the data
      ASKAPLOG_DEBUG_STR(logger, "Starting to degrid model" );
      
      // report every 5000000 degridded rows into log in the debug mode
      #ifdef ASKAP_DEBUG
      unsigned long report_every = 5000000; 
      unsigned long total_rows = 0;
      unsigned long current_rows = 0;
      #endif // #ifdef ASKAP_DEBUG
      
      for (itsIdi.init();itsIdi.hasMore();itsIdi.next())
      {
        itsIdi->rwVisibility().set(0.0);
        /*
        ASKAPDEBUGASSERT(itsIdi->nPol() == 4);
        for (casa::uInt p=1;p<3;++p) {
             itsIdi->rwVisibility().xyPlane(p).set(0.);
        }
        */
        for (vector<string>::const_iterator it=completions.begin();it!=completions.end();it++)
        {
          string imageName("image"+(*it));
          itsModelGridders[imageName]->degrid(*itsIdi);
        }
        #ifdef ASKAP_DEBUG
        const casa::uInt nRow = itsIdi->nRow();
        total_rows += nRow;
        current_rows += nRow;
        if (current_rows > report_every) {
            current_rows = 0;
            ASKAPLOG_DEBUG_STR(logger, "Degridded "<<total_rows<<" rows of data"); 
        }
        #endif // #ifdef ASKAP_DEBUG
      }
      ASKAPLOG_DEBUG_STR(logger, "Finished degridding model" );
    };
    
    /// @brief assign a different iterator
    /// @details This is a temporary method to assign a different iterator.
    /// All this business is a bit ugly, but should go away when all
    /// measurement equations are converted to work with accessors.
    /// @param idi shared pointer to a new iterator
    void ImageFFTEquation::setIterator(IDataSharedIter& idi)
    {
      itsIdi = idi;
    }
    

    // Calculate the residual visibility and image. We transform the model on the fly
    // so that we only have to read (and write) the data once. This uses more memory
    // but cuts down on IO
    void ImageFFTEquation::calcImagingEquations(askap::scimath::ImagingNormalEquations& ne) const
    {
      ASKAPTRACE("ImageFFTEquation::calcImagingEquations");
      
      // We will need to loop over all completions i.e. all sources
      const vector<string> completions(parameters().completions("image"));

      // To minimize the number of data passes, we keep copies of the gridders in memory, and
      // switch between these. This optimization may not be sufficient in the long run.      
      // Set up initial gridders for model and for the residuals. This enables us to 
      // do both at the same time.

      for (vector<string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        const string imageName("image"+(*it));
        SynthesisParamsHelper::clipImage(parameters(),imageName);
        if(itsModelGridders.count(imageName)==0) {
           itsModelGridders[imageName]=itsGridder->clone();
        }
        if(itsResidualGridders.count(imageName)==0) {
          itsResidualGridders[imageName]=itsGridder->clone();
        }
        if(itsPSFGridders.count(imageName)==0) {
          if (itsSphFuncPSFGridder) {
             boost::shared_ptr<SphFuncVisGridder> psfGridder(new SphFuncVisGridder);
             itsPSFGridders[imageName] = psfGridder;
          } else {
             itsPSFGridders[imageName] = itsGridder->clone();
          }
        }
      }
      // Now we initialise appropriately
      ASKAPLOG_DEBUG_STR(logger, "Initialising for model degridding and residual gridding" );
      if (completions.size() == 0) {
          ASKAPLOG_WARN_STR(logger, "Found no free image parameters, this rank will not contribute usefully to normal equations");
      }
      bool somethingHasToBeDegridded = false;
      for (vector<string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        string imageName("image"+(*it));
        const Axes axes(parameters().axes(imageName));
        casa::Array<double> imagePixels(parameters().value(imageName).copy());
        const casa::IPosition imageShape(imagePixels.shape());
        /// First the model
        itsModelGridders[imageName]->customiseForContext(*it);
        itsModelGridders[imageName]->initialiseDegrid(axes, imagePixels);
        if (!itsModelGridders[imageName]->isModelEmpty()) {
            somethingHasToBeDegridded = true;
        }
        /// Now the residual images, dopsf=false
        itsResidualGridders[imageName]->customiseForContext(*it);
        itsResidualGridders[imageName]->initialiseGrid(axes, imageShape, false);
        // and PSF gridders, dopsf=true
        itsPSFGridders[imageName]->customiseForContext(*it);
        itsPSFGridders[imageName]->initialiseGrid(axes, imageShape, true);        
      }
      // synchronise emtpy flag across multiple ranks if necessary
      if (itsVisUpdateObject) {
          itsVisUpdateObject->aggregateFlag(somethingHasToBeDegridded);
      }      
      // Now we loop through all the data
      ASKAPLOG_DEBUG_STR(logger, "Starting degridding model and gridding residuals" );
      size_t counterGrid = 0, counterDegrid = 0;
      for (itsIdi.init();itsIdi.hasMore();itsIdi.next())
      {
        // buffer-accessor, used as a replacement for proper buffers held in the subtable
        // effectively, an array with the same shape as the visibility cube is held by this class
        MemBufferDataAccessor accBuffer(*itsIdi);
         
        // Accumulate model visibility for all models
        accBuffer.rwVisibility().set(0.0);
        if (somethingHasToBeDegridded) {
            for (vector<string>::const_iterator it=completions.begin();it!=completions.end();++it) {
                 const std::string imageName("image"+(*it));
                 const std::map<std::string, IVisGridder::ShPtr>::iterator grdIt = itsModelGridders.find(imageName);
                 ASKAPDEBUGASSERT(grdIt != itsModelGridders.end());
                 const IVisGridder::ShPtr degridder = grdIt->second;
                 ASKAPDEBUGASSERT(degridder);
                 if (!degridder->isModelEmpty()) {
                     degridder->degrid(accBuffer);
                     counterDegrid+=accBuffer.nRow();
                 }
            }
            // optional aggregation of visibilities in the case of distributed model        
            // somethingHasToBeDegridded is supposed to have consistent value across all participating ranks
            if (itsVisUpdateObject) {
                itsVisUpdateObject->update(accBuffer.rwVisibility());
            }
            //            
        }
        accBuffer.rwVisibility() -= itsIdi->visibility();
        accBuffer.rwVisibility() *= float(-1.);

        /// Now we can calculate the residual visibility and image
        size_t tempCounter = 0; 
#ifdef _OPENMP
        #pragma omp parallel default(shared)
        {
           #pragma omp for reduction(+:tempCounter)
#endif
           for (size_t i = 0; i<completions.size(); ++i) {
                const string imageName("image"+completions[i]);
                if (parameters().isFree(imageName)) {
                    #ifdef _OPENMP
                    #pragma omp task
                    #endif
                    itsResidualGridders[imageName]->grid(accBuffer);
                    #ifdef _OPENMP
                    #pragma omp task
                    #endif
                    itsPSFGridders[imageName]->grid(accBuffer);
                    tempCounter += accBuffer.nRow();
                }
           }
#ifdef _OPENMP
        }
#endif
        counterGrid += tempCounter;
      }
      ASKAPLOG_DEBUG_STR(logger, "Finished degridding model and gridding residuals" );
      ASKAPLOG_DEBUG_STR(logger, "Number of accessor rows iterated through is "<<counterGrid<<" (gridding) and "<<
                        counterDegrid<<" (degridding)");

      // We have looped over all the data, so now we have to complete the 
      // transforms and fill in the normal equations with the results from the
      // residual gridders
      ASKAPLOG_DEBUG_STR(logger, "Adding residual image, PSF, and weights image to the normal equations" );
      for (vector<string>::const_iterator it=completions.begin();it!=completions.end();++it)
      {
        const string imageName("image"+(*it));
        const casa::IPosition imageShape(parameters().value(imageName).shape());

        casa::Array<double> imagePSF(imageShape);
        casa::Array<double> imageWeight(imageShape);
        casa::Array<double> imageDeriv(imageShape);

        itsResidualGridders[imageName]->finaliseGrid(imageDeriv);
        
        /*
        // for debugging/research, store grid prior to FFT
        boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(itsPSFGridders[imageName]);
        if (tvg) {
            tvg->storeGrid("uvcoverage"+(*it),0);
        }
        // end debugging code
        */

        itsPSFGridders[imageName]->finaliseGrid(imagePSF);

        itsResidualGridders[imageName]->finaliseWeights(imageWeight);
        /*{ 
          casa::Array<double> imagePSFWeight(imageShape);
          itsPSFGridders[imageName]->finaliseWeights(imagePSFWeight);
          const double maxPSFWeight = casa::max(imagePSFWeight);
          if (maxPSFWeight > 0) {
              const double psfScalingFactor = casa::max(imageWeight)/maxPSFWeight;
              //std::cout<<"psf peak = "<<casa::max(imagePSF)<<" maxPSFWeight = "<<maxPSFWeight<<" factor="<<
              //     psfScalingFactor<<" maxImageWeight="<<casa::max(imageWeight)<<std::endl;
              imagePSF *= psfScalingFactor;
              // now psf has the same peak as the weight image
          } else {
             if (maxPSFWeight < 0) {
                 ASKAPTHROW(AskapError, "PSF weight is supposed to be non-negative, you have "<<maxPSFWeight);
             }
             // do nothing for zero weight, it just means that this part of normal equations is empty. However,
             // we may still be able to have data after summing all parts of the NE
          } 
        }*/
        {
          casa::IPosition reference(4, imageShape(0)/2, imageShape(1)/2, 0, 0);
          casa::IPosition vecShape(1, imagePSF.nelements());
          casa::Vector<double> imagePSFVec(imagePSF.reform(vecShape));
          casa::Vector<double> imageWeightVec(imageWeight.reform(vecShape));
          casa::Vector<double> imageDerivVec(imageDeriv.reform(vecShape));
          ne.addSlice(imageName, imagePSFVec, imageWeightVec, imageDerivVec,
              imageShape, reference);
        }
      }
    }

  }

}
