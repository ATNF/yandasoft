/// @file
///
/// @brief Gridder adapter to do snap-shot imaging 
/// @details We can handle non-coplanarity via snap-shot imaging. For an approximately co-planar 
///     array the effect of w-term at a short time interval is equivalent to a shift. This gridder
///     uses an accessor adapter to monitor changes of the best-fit plane in the u,v,w-spce. If the
///     departure from the previously fitted plane exceeds the tolerance, the image is regridded to
///     a proper coordinate system (taken the shift out). This is an adapter, which can work with
///     any ASKAPsoft gridder. The real gridder, passed as a parameter during construction, does all 
///     the gridding job, so the snap-shot imaging can be combined with w-projection or any other
///     algorithm. The main driver for snap-shot imaging is an attempt to decrease the support size
///     of convolution functions (largely caused by w-projection).  
///
///     See also Ord et al., 2011, PASA (in press); arXiv:1010.1733
///
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

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.snapshotimaginggridderadapter");

#include <gridding/SnapShotImagingGridderAdapter.h>
#include <askap/AskapError.h>
#include <utils/MultiDimArrayPlaneIter.h>
#include <utils/PaddingUtils.h>

#include <casa/OS/Timer.h>
#include <images/Images/ImageRegrid.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <casa/Arrays/Array.h>
#include <scimath/Mathematics/Interpolate2D.h>

//#include <measurementequation/SynthesisParamsHelper.h>

#include <profile/AskapProfiler.h>

using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

#ifdef _OPENMP
/// @brief mutex to deal with lack of thread safety in casa's regrid
boost::mutex SnapShotImagingGridderAdapter::theirMutex;
#endif


/// @brief initialise the adapter
/// @details
/// @param[in] gridder a shared pointer to the gridder to be wrapped by this adapter
/// @param[in] tolerance w-term tolerance in wavelengths (a new fit is performed if the old
/// plane gives w-deviation exceeding this value)
SnapShotImagingGridderAdapter::SnapShotImagingGridderAdapter(const boost::shared_ptr<IVisGridder> &gridder,
                               const double tolerance, const casa::uInt decimate,
                               const casa::Interpolate2D::Method method) :
     itsAccessorAdapter(tolerance), itsDoPSF(false), itsCoeffA(0.), itsCoeffB(0.),
     itsFirstAccessor(true), itsBuffersFinalised(false), itsNumOfImageRegrids(0), itsTimeImageRegrid(0.),
     itsNumOfInitialisations(0), itsLastFitTimeStamp(0.), itsShortestIntervalBetweenFits(3e7),
     itsLongestIntervalBetweenFits(-1.), itsModelIsEmpty(false), itsClippingFactor(0.), itsNoPSFReprojection(true),
     itsDecimationFactor(decimate), itsInterpolationMethod(method)
{
  ASKAPCHECK(gridder, "SnapShotImagingGridderAdapter should only be initialised with a valid gridder");
  itsGridder = gridder->clone();
}

/// @brief copy constructor
/// @details We need this because the gridder doing actual work is held by a shared pointer,
/// which is a non-trivial type
/// @param[in] other an object to copy from
SnapShotImagingGridderAdapter::SnapShotImagingGridderAdapter(const SnapShotImagingGridderAdapter &other) :
    IVisGridder(other), itsAccessorAdapter(other.itsAccessorAdapter.tolerance()), itsDoPSF(other.itsDoPSF),
    itsAxes(other.itsAxes), itsImageBuffer(other.itsImageBuffer.copy()),
    itsWeightsBuffer(other.itsWeightsBuffer.copy()), itsCoeffA(other.itsCoeffA), itsCoeffB(other.itsCoeffB),
    itsFirstAccessor(other.itsFirstAccessor), itsBuffersFinalised(other.itsBuffersFinalised),
    itsNumOfImageRegrids(other.itsNumOfImageRegrids), itsTimeImageRegrid(other.itsTimeImageRegrid),
    itsNumOfInitialisations(other.itsNumOfInitialisations), itsLastFitTimeStamp(other.itsLastFitTimeStamp),
    itsShortestIntervalBetweenFits(other.itsShortestIntervalBetweenFits), 
    itsLongestIntervalBetweenFits(other.itsLongestIntervalBetweenFits), 
    itsTempInImg(), itsTempOutImg(), itsModelIsEmpty(other.itsModelIsEmpty), itsClippingFactor(other.itsClippingFactor),
    itsNoPSFReprojection(other.itsNoPSFReprojection), itsDecimationFactor(other.itsDecimationFactor),
    itsInterpolationMethod(other.itsInterpolationMethod)
{
  ASKAPCHECK(other.itsGridder, 
       "copy constructor of SnapShotImagingGridderAdapter got an object somehow set up with an empty gridder");
  ASKAPCHECK(!other.itsAccessorAdapter.isAssociated(), 
     "An attempt to copy gridder adapter with the accessor adapter associated with some real data accessor. This shouldn't happen.");
  itsGridder = other.itsGridder->clone();  
}

/// @brief destructor just to print some stats
SnapShotImagingGridderAdapter::~SnapShotImagingGridderAdapter()
{
  if (itsNumOfInitialisations>0) {
      ASKAPLOG_INFO_STR(logger, "SnapShotImagingGridderAdapter usage statistics");
      const std::string msg = itsNoPSFReprojection ? "non-PSF " : "";
      ASKAPLOG_INFO_STR(logger, "   The adapter was initialised for "<<msg<<"gridding and degridding "<<
                        itsNumOfInitialisations<<" times");     
      ASKAPLOG_INFO_STR(logger, "   Total time spent doing image plane regridding is "<<
                        itsTimeImageRegrid<<" (s)");
      ASKAPLOG_INFO_STR(logger, "   Number of regridding events is "<<itsNumOfImageRegrids);
      ASKAPLOG_INFO_STR(logger, "   or "<<double(itsNumOfImageRegrids)/double(itsNumOfInitialisations)<<
                        " times per grid/degrid pass");      
      ASKAPLOG_INFO_STR(logger, "   Image clipping factor (clipping during regrids) is "<< itsClippingFactor);
      if (itsNumOfImageRegrids > 0) {
          ASKAPLOG_INFO_STR(logger, "   Average time spent per image plane regridding is "<<
                      itsTimeImageRegrid/double(itsNumOfImageRegrids)<<" (s)");
      } 
      reportAndInitIntervalStats();     
  }
}

/// @brief report current interval stats and initialise them
/// @details We collect and report such statistics like shortest and longest
/// intervals between changes to the best fit plane (and therefore between 
/// image regrids). As the adapter can be reused multiple times, these
/// stats need to be reset every time a new initialisation is done. This 
/// method reports current stats to the log (if there is something to report; the
/// initial values are such that they shouldn't occur in normal operations and can
/// serve as flags) and initialises them for the next pass.
void SnapShotImagingGridderAdapter::reportAndInitIntervalStats() const
{
  // in principle, (itsNumOfInitialisations > 0) condition is redundant
  if ((itsLongestIntervalBetweenFits > 0) && (itsNumOfInitialisations > 0) ) {
      // intervals were initialised, report them
      ASKAPLOG_DEBUG_STR(logger, "Longest observing time interval between image plane regrids is "<<
                         itsLongestIntervalBetweenFits<<" (s)");
      ASKAPLOG_DEBUG_STR(logger, "Shortest observing time interval between image plane regrids is "<<
                         itsShortestIntervalBetweenFits<<" (s)");                         
  } 
  itsLongestIntervalBetweenFits = -1.;
  itsShortestIntervalBetweenFits = 3e7;
}

/// @brief update interval stats for the new fit
/// @details This method updates interval statistics for the new fit
/// @param[in] time current time reported by the accessor triggering fit update
void SnapShotImagingGridderAdapter::updateIntervalStats(const double time) const
{
  const double interval = fabs(time - itsLastFitTimeStamp);
  if ((interval < itsShortestIntervalBetweenFits) || (itsLongestIntervalBetweenFits < 0)) {
       itsShortestIntervalBetweenFits = interval;
  }
  if (interval > itsLongestIntervalBetweenFits) {
      itsLongestIntervalBetweenFits = interval;
  }
  itsLastFitTimeStamp = time;
}



/// @brief clone a copy of this gridder
/// @return shared pointer to the clone
boost::shared_ptr<IVisGridder> SnapShotImagingGridderAdapter::clone()
{
  boost::shared_ptr<SnapShotImagingGridderAdapter> newOne(new SnapShotImagingGridderAdapter(*this));
  return newOne;
}

/// @brief initialise the gridding
/// @details
/// @param[in] axes axes specifications
/// @param[in] shape Shape of output image: cube: u,v,pol,chan
/// @param[in] dopsf Make the psf?
void SnapShotImagingGridderAdapter::initialiseGrid(const scimath::Axes& axes,
                const casa::IPosition& shape, const bool dopsf)
{
  ASKAPTRACE("SnapShotImagingGridderAdapter::initialiseGrid");

  itsDoPSF = dopsf; // other fields are unused for the PSF gridder
  if (dopsf && itsNoPSFReprojection) {
      // do the standard initialisation for the PSF gridder
      ASKAPDEBUGASSERT(itsGridder);
      itsGridder->initialiseGrid(axes,shape,dopsf);
  } else {
      ASKAPDEBUGASSERT(shape.nelements() >= 2);
      reportAndInitIntervalStats();
      ++itsNumOfInitialisations;
      itsAxes = axes;
      // initialise the buffers for final image and weight
      itsImageBuffer.resize(shape);
      itsWeightsBuffer.resize(shape);
      itsImageBuffer.set(0.);
      itsWeightsBuffer.set(0.);
      // the following flag means the gridding will be 
      // initialised when the first accessor is encountered
      itsFirstAccessor = true; 
      // nothing gridded, zero buffers are the correct output
      itsBuffersFinalised = true; 
  }
}

/// @brief grid the visibility data.
/// @param[in] acc const data accessor to work with
void SnapShotImagingGridderAdapter::grid(IConstDataAccessor& acc)
{
  ASKAPTRACE("SnapShotImagingGridderAdapter::grid");

  ASKAPDEBUGASSERT(itsGridder);
  if (isPSFGridder() && itsNoPSFReprojection) {
      itsAccessorAdapter.associate(acc);
      // for PSF gridder we don't do any image-plane regridding in this mode
      itsGridder->grid(itsAccessorAdapter);
      // we don't really need this line
      itsAccessorAdapter.detach();      
  } else {
      itsAccessorAdapter.associate(acc);
      const scimath::ChangeMonitor cm = itsAccessorAdapter.planeChangeMonitor();
      // the call to rotatedUVW method would assess whether the current plane is still
      // fine. The result is cached, so there is no performance penalty.
      itsAccessorAdapter.rotatedUVW(getTangentPoint());
      if ((cm != itsAccessorAdapter.planeChangeMonitor()) || itsFirstAccessor) {
          if (!itsFirstAccessor) {
              // there is nothing to finalise, if this is the first accessor
              finaliseGriddingOfCurrentPlane();
              itsFirstAccessor = true;
              updateIntervalStats(acc.time());
          } else {
              itsLastFitTimeStamp = acc.time();
          }
          // update plane parameters
          itsCoeffA = itsAccessorAdapter.coeffA();
          itsCoeffB = itsAccessorAdapter.coeffB();          
      } 
      if (itsFirstAccessor) {
          scimath::Axes axes = itsAxes;
          // need to patch axes here before passing to initialise grid
          axes.addDirectionAxis(currentPlaneDirectionCoordinate());
          //
          itsGridder->initialiseGrid(axes,itsImageBuffer.shape(),isPSFGridder());
          itsFirstAccessor = false;
      }
      itsGridder->grid(itsAccessorAdapter);
      itsBuffersFinalised = false;
      // we don't really need this line
      itsAccessorAdapter.detach();
  }
}

/// @brief form the final output image
/// @param[in] out output double precision image or PSF
void SnapShotImagingGridderAdapter::finaliseGrid(casa::Array<double>& out)
{
  ASKAPTRACE("SnapShotImagingGridderAdapter::finaliseGrid");

  ASKAPDEBUGASSERT(itsGridder);
  if (isPSFGridder() && itsNoPSFReprojection) {
      itsGridder->finaliseGrid(out);
  } else {
      if (!itsBuffersFinalised) {
          finaliseGriddingOfCurrentPlane();
      }
      out.assign(itsImageBuffer);  
  }
}

/// @brief finalise weights
/// @details Form the sum of the convolution function squared, multiplied by the weights for each
/// different convolution function. This is used in the evaluation of the second derivative.
/// @param[in] out output double precision sum of weights images
void SnapShotImagingGridderAdapter::finaliseWeights(casa::Array<double>& out)
{
  ASKAPTRACE("SnapShotImagingGridderAdapter::finaliseWeights");

  ASKAPDEBUGASSERT(itsGridder);
  if (isPSFGridder() && itsNoPSFReprojection) {
      itsGridder->finaliseWeights(out);
  } else {
      if (!itsBuffersFinalised) {
          finaliseGriddingOfCurrentPlane();
      }
      out.assign(itsWeightsBuffer);  
  }      
}

/// @brief initialise the degridding
/// @param[in] axes axes specifications
/// @param[in] image input image cube: u,v,pol,chan
void SnapShotImagingGridderAdapter::initialiseDegrid(const scimath::Axes& axes,
					const casa::Array<double>& image)
{
  ASKAPTRACE("SnapShotImagingGridderAdapter::initialiseDegrid");

  itsModelIsEmpty = (casa::max(casa::abs(image)) <= 0.);
  if (itsModelIsEmpty) {
      ASKAPLOG_INFO_STR(logger, "No need to degrid: model is empty");
      return;
  }

  ASKAPDEBUGASSERT(itsGridder);
  reportAndInitIntervalStats();
  ++itsNumOfInitialisations;
  itsDoPSF = false;
  itsAxes = axes;
  itsImageBuffer.assign(image);
  // the following flag means the gridding will be 
  // initialised when the first accessor is encountered
  itsFirstAccessor = true; 
}					

/// @brief make context-dependant changes to the gridder behaviour
/// @param[in] context context description
void SnapShotImagingGridderAdapter::customiseForContext(const std::string &context)
{
  ASKAPDEBUGASSERT(itsGridder);
  itsGridder->customiseForContext(context);
}
			
/// @brief set visibility weights
/// @param[in] viswt shared pointer to visibility weights
void SnapShotImagingGridderAdapter::initVisWeights(const IVisWeights::ShPtr &viswt)
{
  ASKAPDEBUGASSERT(itsGridder);
  itsGridder->initVisWeights(viswt);
}

/// @brief degrid the visibility data.
/// @param[in] acc non-const data accessor to work with  
void SnapShotImagingGridderAdapter::degrid(IDataAccessor& acc)
{
  ASKAPTRACE("SnapShotImagingGridderAdapter::degrid");

  if (itsModelIsEmpty) {
      return;
  }
  ASKAPDEBUGASSERT(itsGridder);
  itsAccessorAdapter.associate(acc);
  const scimath::ChangeMonitor cm = itsAccessorAdapter.planeChangeMonitor();
  // the call to rotatedUVW method would assess whether the current plane is still
  // fine. The result is cached, so there is no performance penalty.
  itsAccessorAdapter.rotatedUVW(getTangentPoint());
  if ((cm != itsAccessorAdapter.planeChangeMonitor()) || itsFirstAccessor) {
       if (!itsFirstAccessor) {
           // there is nothing to finalise, if this is the first accessor
           itsGridder->finaliseDegrid();
           itsFirstAccessor = true;
           updateIntervalStats(acc.time());
       } else {
           itsLastFitTimeStamp = acc.time();
       }
       // update plane parameters
       itsCoeffA = itsAccessorAdapter.coeffA();
       itsCoeffB = itsAccessorAdapter.coeffB();          
  } 
  if (itsFirstAccessor) {
      scimath::Axes axes = itsAxes;
      // need to patch axes here before passing to initialise degrid
      axes.addDirectionAxis(currentPlaneDirectionCoordinate());
      //
      casa::Array<double> scratch(itsImageBuffer.shape());
      imageRegrid(itsImageBuffer,scratch, false);
      itsGridder->initialiseDegrid(axes,scratch);
      itsFirstAccessor = false;
  }
  itsGridder->degrid(itsAccessorAdapter);
  // we don't really need this line
  itsAccessorAdapter.detach();  
}

/// @brief finalise degridding
void SnapShotImagingGridderAdapter::finaliseDegrid()
{
  ASKAPTRACE("SnapShotImagingGridderAdapter::finaliseDegrid");

  if (itsModelIsEmpty) {
      return;
  }
  ASKAPDEBUGASSERT(itsGridder);
  ASKAPCHECK(!itsFirstAccessor, 
       "finaliseDegrid is called while the itsFirstAccessor flag is true. This is not supposed to happen");
  itsGridder->finaliseDegrid();
}

/// @brief finalise gridding for the current plane
/// @details We execute the gridder pointed by itsGridder
/// multiple times. Every time the best fitted plane changes
/// we have to finalise gridding with the wrapped gridder,
/// regrid the result into target frame and add it to buffers.
/// The same has to be done for both image and weights. This
/// method encapsulates all these operations.
void SnapShotImagingGridderAdapter::finaliseGriddingOfCurrentPlane()
{
  ASKAPDEBUGTRACE("SnapShotImagingGridderAdapter::finaliseGriddingOfCurrentPlane");
  ASKAPDEBUGASSERT(itsGridder);
  ASKAPCHECK(!itsFirstAccessor, 
       "finaliseGriddingOfCurrentPlane is called while itsFirstAccessor flag is true. This is not supposed to happen");
  const std::string msg = isPSFGridder() ? "PSF" : "dirty image";     
  ASKAPLOG_DEBUG_STR(logger, "Finalising current "<<msg);
  casa::Array<double> scratch(itsImageBuffer.shape());
  itsGridder->finaliseGrid(scratch);
  imageRegrid(scratch, itsImageBuffer, true);
  ASKAPLOG_DEBUG_STR(logger, "Finalising current weights");
  itsGridder->finaliseWeights(scratch);
  imageRegrid(scratch, itsWeightsBuffer, true);  
  itsBuffersFinalised = true;
}

/// @brief direction coordinate corresponding to the current fit plane
/// @details This method forms a direction coordinate corresponding to the
/// current best fit w=Au+Bv from the direction coordinate stored in 
/// itsAxes. This is used to setup image plane regridding and coordinate system
/// of the wrapped gridder during grid/degrid initialisation.
casa::DirectionCoordinate SnapShotImagingGridderAdapter::currentPlaneDirectionCoordinate() const
{
  ASKAPDEBUGASSERT(itsAxes.hasDirection());
  const casa::DirectionCoordinate dc(itsAxes.directionAxis());
  const casa::MDirection::Types directionType = dc.directionType();
  const casa::Vector<casa::Double> refVal = dc.referenceValue();
  ASKAPDEBUGASSERT(refVal.nelements() == 2);
  const casa::Vector<casa::Double> refPix = dc.referencePixel();
  ASKAPDEBUGASSERT(refPix.nelements() == 2);
  const casa::Vector<casa::Double> inc = dc.increment();
  ASKAPDEBUGASSERT(inc.nelements() == 2);
  const casa::Matrix<casa::Double> xform = dc.linearTransform();
  // now patch projection
  casa::Vector<casa::Double> projParams(2);
  projParams[0] = -coeffA();
  projParams[1] = -coeffB();
  const casa::Projection projection(casa::Projection::SIN, projParams);
  //
  return casa::DirectionCoordinate(directionType, projection, refVal[0],refVal[1],
                   inc[0],inc[1],xform,refPix[0],refPix[1]);
}


/// @brief regrid images between frames
/// @details This method does the core regridding procedure. It iterates
/// over 2D planes of the input array, regrids them into the other frame
/// and either adds the result to the appropriate plane of the output array,
/// if the regridding is into the target frame or replaces the result if it is
/// from the target frame. 
/// @param[in] input input array to be regridded
/// @param[out] output output array
/// @param[in] toTarget true, if regridding is from the current frame into the 
/// target frame (for gridding); false if regridding is from the target frame 
/// into the current frame (for degridding)
/// @note The output and input arrays should have the same shape. The iteration
/// over 2D planes is perfromed explicitly to avoid initialising large scratch 
/// buffers. An exception is raised if input and output arrays have different 
/// shapes
void SnapShotImagingGridderAdapter::imageRegrid(const casa::Array<double> &input, 
           casa::Array<double> &output, bool toTarget) const
{
   ASKAPTRACE("SnapShotImagingGridderAdapter::imageRegrid");    
  
   // for stats
   casa::Timer timer;
   timer.mark();
   ++itsNumOfImageRegrids;
   // actual code
   if (toTarget) {
       ASKAPLOG_DEBUG_STR(logger, "Regridding image from the frame corresponding to the fitted plane w = u * "<<
              coeffA()<<" + v * "<<coeffB()<<", into the target frame");
   } else {
       ASKAPLOG_DEBUG_STR(logger, 
           "Regridding image from the input frame into a frame corresponding to the fitted plane w = u * "<<
              coeffA()<<" + v * "<<coeffB());       
   }
   ASKAPCHECK(input.shape() == output.shape(), 
           "The shape of input and output arrays should be identical, input.shape()="<<
              input.shape()<<", output.shape()="<<output.shape());
   ASKAPDEBUGASSERT(input.shape().nelements() >= 2);
   
   // constness is conceptual, we don't do any assignments to the input array
   // the following line doesn't copy the data (reference semantics)
   casa::Array<double> inRef(input);
   // form coordinate systems
   const casa::DirectionCoordinate dcCurrent = currentPlaneDirectionCoordinate();
   const casa::DirectionCoordinate& dcTarget = itsAxes.directionAxis();
   casa::CoordinateSystem csInput;
   casa::CoordinateSystem csOutput;
   if (toTarget) {
      csInput.addCoordinate(dcCurrent);
      csOutput.addCoordinate(dcTarget);
   } else {
      csInput.addCoordinate(dcTarget);
      csOutput.addCoordinate(dcCurrent);
   }
   
   // iterator over planes
   scimath::MultiDimArrayPlaneIter planeIter(input.shape());
   
   // regridder
   casa::ImageRegrid<double> regridder;
   // regridder works with images, so we have to setup temporary 2D images
   // the following may cause an unnecessary copy, there should be a better way
   // of constructing an image out of an array
   const casa::IPosition tempShape = planeIter.planeShape().nonDegenerate();
   if (!itsTempInImg.shape().isEqual(tempShape)) {
       /* 
       // this resizing is temporary replaced with a more convoluted operation
       // as a workaround to avoid a possible casacore bug with TempImage
       itsTempInImg.resize(casa::TiledShape(tempShape));
       itsTempOutImg.resize(casa::TiledShape(tempShape));       
       */
       // +100 forces to use the memory
       const double maxMemoryInMB = double(tempShape.product()*sizeof(double))/1024./1024.+100;
       itsTempInImg = casa::TempImage<double>(casa::TiledShape(tempShape),csInput,maxMemoryInMB);
       itsTempOutImg = casa::TempImage<double>(casa::TiledShape(tempShape),csOutput,maxMemoryInMB);       
   }
   ASKAPDEBUGASSERT(itsTempInImg.shape().isEqual(itsTempOutImg.shape()));
   const bool csSuccess = itsTempInImg.setCoordinateInfo(csInput) && itsTempOutImg.setCoordinateInfo(csOutput);
   ASKAPCHECK(csSuccess, "Error setting either input or output coordinate frame during image plane regridding");
                  
   for (; planeIter.hasMore(); planeIter.next()) {
        itsTempInImg.put(planeIter.getPlane(inRef));
        #ifdef _OPENMP
        { 
          boost::unique_lock<boost::mutex> lock(theirMutex);
        #endif
          regridder.regrid(itsTempOutImg, itsInterpolationMethod,
                  casa::IPosition(2,0,1), itsTempInImg, false, 3);
        #ifdef _OPENMP
        }
        #endif
        // the next line does not do any copying (reference semantics)
        casa::Array<double> outRef(planeIter.getPlane(output).nonDegenerate());
        // create a lattice to benefit from lattice math operators
        casa::ArrayLattice<double> tempOutputLattice(outRef, casa::True);
        if (toTarget) {
            tempOutputLattice += itsTempOutImg;
        } else {
          // just assign the result
          tempOutputLattice.copyData(itsTempOutImg);
        }
        // optional clipping
        imageClip(outRef);
   }
   itsTimeImageRegrid += timer.real();
}

/// @brief obtain the tangent point
/// @details This method extracts the tangent point (reference position) from the
/// coordinate system.
/// @return direction measure corresponding to the tangent point
casa::MVDirection SnapShotImagingGridderAdapter::getTangentPoint() const
{
   // at this stage, just a copy of the method from TableVisGridder. May need some refactoring 
   // in the future
   ASKAPCHECK(itsAxes.hasDirection(),"Direction axis is missing. axes="<<itsAxes);
   const casa::Vector<casa::Double> refVal(itsAxes.directionAxis().referenceValue());
   ASKAPDEBUGASSERT(refVal.nelements() == 2);
   const casa::Quantum<double> refLon(refVal[0], "rad");
   const casa::Quantum<double> refLat(refVal[1], "rad");
   const casa::MVDirection out(refLon, refLat);
   return out;
}

/// @brief set clipping factor
/// @details The image could be optionally clipped during regridding (to avoid edge effects). 
/// This parameter represents the fraction of the image size (on each directional axis) which is
/// zeroed (equally from both sides). It should be a non-negative number less than 1. Set to zero to avoid
/// any clipping (this is the default behavior)
/// @param[in] factor clipping factor
void SnapShotImagingGridderAdapter::setClippingFactor(const float factor) 
{
  ASKAPCHECK((factor >= 0.) && (factor < 1), "Clipping factor should be a non-negative number less than 1, you have "<<factor);
  itsClippingFactor = factor;
}

/// @brief clip image 
/// @details This method clips the image by zeroing the edges according to the
/// assigned clipping factor.
/// @param[in] img array to modify
void SnapShotImagingGridderAdapter::imageClip(casa::Array<double> &img) const
{
  ASKAPDEBUGTRACE("SnapShotImagingGridderAdapter::imageClip");

  ASKAPDEBUGASSERT(itsClippingFactor < 1.);
  const float unpaddingFactor = 1. - itsClippingFactor;
  ASKAPDEBUGASSERT(unpaddingFactor <= 1.);
  const casa::IPosition shapeToPreserve = scimath::PaddingUtils::paddedShape(img.shape(), unpaddingFactor);
  if (shapeToPreserve != img.shape()) {
      ASKAPASSERT(img.shape().nelements() == 2);
      scimath::PaddingUtils::clip(img, shapeToPreserve);
  }
}

/// @brief control whether to do image reprojection for PSF
/// @details By default we bypass image reprojection for the PSF. It can be changed with this configuration method.
/// @param[in] doIt if true, image reprojection will be done for PSF the same way dirty image and weight are processed,
///                 otherwise (the default), the wrapped gridder is used directly without any reprojection
void SnapShotImagingGridderAdapter::setPSFReprojection(const bool doIt)
{
  itsNoPSFReprojection = !doIt;
  if (doIt) {
        ASKAPLOG_INFO_STR(logger, "PSF image will be reprojected the same way as the residual image");      
  } else {
        ASKAPLOG_INFO_STR(logger, "No PSF reprojection will be done");      
  }
}

/// @brief check whether the model is empty
/// @details A simple check allows us to bypass heavy calculations if the input model
/// is empty (all pixels are zero). This makes sense for degridding only.
/// @brief true, if the model is empty
bool SnapShotImagingGridderAdapter::isModelEmpty() const 
{
  return itsModelIsEmpty;
}


