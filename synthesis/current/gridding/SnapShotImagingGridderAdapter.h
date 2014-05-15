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

#ifndef SNAP_SHOT_IMAGING_GRIDDER_ADAPTER_H
#define SNAP_SHOT_IMAGING_GRIDDER_ADAPTER_H

#include <gridding/IVisGridder.h>
#include <boost/shared_ptr.hpp>
#include <dataaccess/BestWPlaneDataAccessor.h>
#include <fitting/Axes.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <scimath/Mathematics/Interpolate2D.h>
#include <measures/Measures/MDirection.h>
#include <images/Images/TempImage.h>

#ifdef _OPENMP
#include <boost/thread/mutex.hpp>
#endif

namespace askap {

namespace synthesis {

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
/// @ingroup gridding
class SnapShotImagingGridderAdapter : virtual public IVisGridder 
{
public:
   /// @brief initialise the adapter
   /// @details
   /// @param[in] gridder a shared pointer to the gridder to be wrapped by this adapter
   /// @param[in] tolerance w-term tolerance in wavelengths (a new fit is performed if the old
   /// plane gives w-deviation exceeding this value)
   /// @param[in] decimate the decimation factor used for coordinate system interpolation
   /// during regridding. A decimation factor of 0 (or 1) implies no decimation, which is the
   /// slowest and most accurate.
   /// @param[in] method the interpolation method to use for regridding.
   SnapShotImagingGridderAdapter(const boost::shared_ptr<IVisGridder> &gridder,
           const double tolerance,
           const casa::uInt decimate = 0,
           const casa::Interpolate2D::Method method = casa::Interpolate2D::CUBIC);


   /// @brief copy constructor
   /// @details We need this because the gridder doing actual work is held by a shared pointer,
   /// which is a non-trivial type
   /// @param[in] other an object to copy from
   SnapShotImagingGridderAdapter(const SnapShotImagingGridderAdapter &other);
   
   /// @brief destructor just to print some stats
   virtual ~SnapShotImagingGridderAdapter();

   /// @brief clone a copy of this gridder
   /// @return shared pointer to the clone
   virtual boost::shared_ptr<IVisGridder> clone();

   /// @brief initialise the gridding
   /// @details
   /// @param[in] axes axes specifications
   /// @param[in] shape Shape of output image: cube: u,v,pol,chan
   /// @param[in] dopsf Make the psf?
   virtual void initialiseGrid(const scimath::Axes& axes,
                const casa::IPosition& shape, const bool dopsf = true);

   /// @brief grid the visibility data.
   /// @param[in] acc const data accessor to work with
   virtual void grid(accessors::IConstDataAccessor& acc);

   /// @brief form the final output image
   /// @param[in] out output double precision image or PSF
   virtual void finaliseGrid(casa::Array<double>& out);

   /// @brief finalise weights
   /// @details Form the sum of the convolution function squared, multiplied by the weights for each
   /// different convolution function. This is used in the evaluation of the second derivative.
   /// @param[in] out output double precision sum of weights images
   virtual void finaliseWeights(casa::Array<double>& out);

   /// @brief initialise the degridding
   /// @param[in] axes axes specifications
   /// @param[in] image input image cube: u,v,pol,chan
   virtual void initialiseDegrid(const scimath::Axes& axes,
					const casa::Array<double>& image);

   /// @brief make context-dependant changes to the gridder behaviour
   /// @param[in] context context description
   virtual void customiseForContext(const std::string &context);
			
   /// @brief set visibility weights
   /// @param[in] viswt shared pointer to visibility weights
   virtual void initVisWeights(const IVisWeights::ShPtr &viswt);

   /// @brief degrid the visibility data.
   /// @param[in] acc non-const data accessor to work with  
   virtual void degrid(accessors::IDataAccessor& acc);

   /// @brief finalise degridding
   virtual void finaliseDegrid();

   /// @brief set clipping factor
   /// @details The image could be optionally clipped during regridding (to avoid edge effects). 
   /// This parameter represents the fraction of the image size (on each directional axis) which is
   /// zeroed (equally from both sides). It should be a non-negative number less than 1. Set to zero to avoid
   /// any clipping (this is the default behavior)
   /// @param[in] factor clipping factor
   void setClippingFactor(const float factor);

   /// @brief control whether to do image reprojection for PSF
   /// @details By default we bypass image reprojection for the PSF. It can be changed with this configuration method.
   /// @param[in] doIt if true, image reprojection will be done for PSF the same way dirty image and weight are processed,
   ///                 otherwise (the default), the wrapped gridder is used directly without any reprojection
   void setPSFReprojection(const bool doIt);
   
   /// @brief check whether the model is empty
   /// @details A simple check allows us to bypass heavy calculations if the input model
   /// is empty (all pixels are zero). This makes sense for degridding only.
   /// @brief true, if the model is empty
   virtual bool isModelEmpty() const; 
   
protected:
   /// @brief check whether this is a psf gridder
   /// @details We pass all calls directly to the wrapped gridder for PSF calculation,
   /// but may do something fancy if the gridder is calculating dirty image
   /// @return true, if this gridder is a psf gridder
   bool inline isPSFGridder() const { return itsDoPSF;}

   /// @brief coefficient A of the fit
   /// @details This is a parameter of the best fitted plane 
   /// w=Au+Bv corresponding to the current job done by the wrapped gridder
   /// @return coefficient A of the fit
   double inline coeffA() const { return itsCoeffA; }

   /// @brief coefficient B of the fit
   /// @details This is a parameter of the best fitted plane 
   /// w=Au+Bv corresponding to the current job done by the wrapped gridder
   /// @return coefficient B of the fit
   double inline coeffB() const { return itsCoeffB; }
   
   /// @brief finalise gridding for the current plane
   /// @details We execute the gridder pointed by itsGridder
   /// multiple times. Every time the best fitted plane changes
   /// we have to finalise gridding with the wrapped gridder,
   /// regrid the result into target frame and add it to buffers.
   /// The same has to be done for both image and weights. This
   /// method encapsulates all these operations.
   void finaliseGriddingOfCurrentPlane();
   
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
   void imageRegrid(const casa::Array<double> &input, casa::Array<double> &output,
                    bool toTarget) const;
   
   
   /// @brief clip image 
   /// @details This method clips the image by zeroing the edges according to the
   /// assigned clipping factor.
   /// @param[in] img array to modify
   void imageClip(casa::Array<double> &img) const;
   
   /// @brief obtain the tangent point
   /// @details This method extracts the tangent point (reference position) from the
   /// coordinate system.
   /// @return direction measure corresponding to the tangent point
   casa::MVDirection getTangentPoint() const;
   
   /// @brief report current interval stats and initialise them
   /// @details We collect and report such statistics like shortest and longest
   /// intervals between changes to the best fit plane (and therefore between 
   /// image regrids). As the adapter can be reused multiple times, these
   /// stats need to be reset every time a new initialisation is done. This 
   /// method reports current stats to the log (if there is something to report; the
   /// initial values are such that they shouldn't occur in normal operations and can
   /// serve as flags) and initialises them for the next pass.
   void reportAndInitIntervalStats() const;
   
   /// @brief update interval stats for the new fit
   /// @details This method updates interval statistics for the new fit
   /// @param[in] time current time reported by the accessor triggering fit update
   void updateIntervalStats(const double time) const;
   
   /// @brief direction coordinate corresponding to the current fit plane
   /// @details This method forms a direction coordinate corresponding to the
   /// current best fit w=Au+Bv from the direction coordinate stored in 
   /// itsAxes. This is used to setup image plane regridding and coordinate system
   /// of the wrapped gridder during grid/degrid initialisation.
   casa::DirectionCoordinate currentPlaneDirectionCoordinate() const;   
   
private:

   /// @brief gridder doing actual job
   boost::shared_ptr<IVisGridder> itsGridder;
   
   /// @brief adapter dealing with plane fitting
   accessors::BestWPlaneDataAccessor itsAccessorAdapter;
   
   /// @brief true if this gridder is computing PSF
   /// @details We pass all calls to the wrapped gridder for PSF calculation
   bool itsDoPSF;
   
   /// @brief target coordinate system visible from outside
   /// @details This gridder can setup the wrapped gridder to work with a 
   /// different coordinate system (corresponding to some fitted plane).
   /// This field stores the target coordinate system for gridding or 
   /// input coordinate system for degridding.
   scimath::Axes itsAxes;
   
   /// @brief buffer for the image
   /// @details It accumulates dirty or model images in the target
   /// coordinate frame
   casa::Array<double> itsImageBuffer;
   
   /// @brief buffer for weights
   /// @details It accumulates weights during gridding
   casa::Array<double> itsWeightsBuffer;
   
   /// @brief coefficient A of the fit
   /// @details This is a parameter of the best fitted plane 
   /// w=Au+Bv corresponding to the current job done by the wrapped gridder
   double itsCoeffA;
   
   /// @brief coefficient B of the fit
   /// @details This is a parameter of the best fitted plane 
   /// w=Au+Bv corresponding to the current job done by the wrapped gridder
   double itsCoeffB;      
   
   /// @brief true, if this is a first accessor per plane
   /// @details This flag essentially tracks when we need to call initialise method
   bool itsFirstAccessor;
   
   /// @brief true if the final image has been finalised
   /// @details finaliseGrid and finaliseWeights can be called in a random order.
   /// We do the same operations to finalise both. This flag shows that the buffers
   /// are ready to be returned to the end user.
   bool itsBuffersFinalised;
   
   /// @brief number of image regridding events
   /// @details This field is used to build stats.
   mutable unsigned long itsNumOfImageRegrids; 
   
   /// @brief total time spent so far doing image regridding
   /// @details This field is used to build stats. Time is in seconds.
   mutable double itsTimeImageRegrid;
   
   /// @brief number of initialisations
   /// @details This field is incremented every time this adapter is initialised
   /// for gridding or degridding. It is used to build stats in the case the object is
   /// reused multiple times. Note, stats are not printed, if this number is zero
   /// (i.e. for unused gridders). We also don't increment it for PSF gridding as 
   /// all calls are passed through to the wrapped gridder.
   unsigned long itsNumOfInitialisations;
   
   /// @brief time corresponding to current solution
   /// @details This is the time reported by the accessor at the time a new
   /// solution for A and B (i.e. a new fit w=Au+Bv) is constructed. It is used to
   /// build stats for the timespan between individual fits
   mutable double itsLastFitTimeStamp; 
   
   /// @brief shortest time interval between fit operations
   /// @note the units are those returned by the accessor (configured at the top level,
   /// currently in seconds).
   mutable double itsShortestIntervalBetweenFits;

   /// @brief shortest time interval between fit operations
   /// @note the units are those returned by the accessor (configured at the top level,
   /// currently in seconds).
   mutable double itsLongestIntervalBetweenFits;
   
   // buffers used in regridding (promoted to class data members to avoid repeated initialisation
   // inside the loop)
   // we don't need to preserve content outside imageRegrid method. However, it may be handy
   // to initialise these images to a proper size inside initialiseGrid or initialiseDegrid
   // (it is not done at the moment, all usage of these data fields is confimed to imageRegrid)
   
   /// @brief temporary input image for regridding
   mutable casa::TempImage<double> itsTempInImg;
   
   /// @brief temporary output image for regridding
   mutable casa::TempImage<double> itsTempOutImg;   

   /// @brief flag that the model is empty for degridding
   /// @details It allows to bypass expensive image regridding
   bool itsModelIsEmpty;
   
   /// @brief clipping factor
   /// @details The image could be optionally clipped during regridding (to avoid edge effects). 
   /// This parameter represents the fraction of the image size (on each directional axis) which is
   /// zeroed (equally from both sides). It should be a non-negative number less than 1.
   float itsClippingFactor;
   
   /// @brief if true, skip image reprojection for PSF
   /// @details By default we bypass image reprojection for the PSF. Although it is an approximation, the
   /// PSF is approximate anyway. However, this field can be set to false to test the effects on deconvolution.
   bool itsNoPSFReprojection;

   /// The decimation factor used for coordinate system interpolation during
   /// regridding. A decimation factor of 0 (or 1) implies no decimation,
   /// which is the slowest and most accurate. 
   casa::uInt itsDecimationFactor; 

   // Interpolation method used for regridding
   casa::Interpolate2D::Method itsInterpolationMethod;

   #ifdef _OPENMP
   /// @brief mutex to deal with lack of thread safety in casa's regrid
   static boost::mutex theirMutex;
   #endif
};
   
} // namespace synthesis

} // namespace askap

#endif // #ifndef SNAP_SHOT_IMAGING_GRIDDER_ADAPTER_H



