///
/// TableVisGridder: Table-based visibility gridder. This is an incomplete
/// class and cannot be used directly. Classes may be derived from this
/// and the unimplemented methods provided. In some cases, it may be 
/// necessary or more efficient to override the provided methods as well.
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///
#ifndef TABLEVISGRIDDER_H_
#define TABLEVISGRIDDER_H_

// own includes
#include <gridding/IVisGridder.h>
#include <gridding/VisGridderWithPadding.h>
#include <dataaccess/IDataAccessor.h>
#include <gridding/FrequencyMapper.h>

// std includes
#include <string>

// casa includes
#include <casa/BasicSL/Complex.h>

#ifdef _OPENMP
// boost includes
#include <boost/thread/mutex.hpp>
#endif

namespace askap
{
  namespace synthesis
  {
      

    /// @brief Incomplete base class for table-based gridding of visibility data.
    ///
    /// TVG supports gridding of visibility data onto one of a number of grids
    /// using one of a number of gridding functions. After gridding the final
    /// grid may be assembled by summing the grids appropriately. The separate
    /// grids may be for separate pointings, or for separate W planes. The summation
    /// process may include Fourier transformation.
    ///
    /// The main work in derived classes is to provide the convolution function
    /// and stacking operations.
    ///
    /// The sensitivity will vary as a function of position within the image.
    /// This is calculated by direct evaluation.
    ///
    /// @ingroup gridding
    class TableVisGridder : public VisGridderWithPadding
    {
  public:

      /// @brief Standard two dimensional gridding using a convolution function
      /// in a table
      TableVisGridder();
      /// @brief Standard two dimensional gridding using a convolution function
      /// in a table
      /// @param[in] overSample Oversampling (currently limited to <=1)
      /// @param[in] support Support of function
      /// @param[in] padding padding factor (default is 1, i.e. no padding)
      /// @param name Name of table to save convolution function into
      TableVisGridder(const int overSample, const int support,
          const float padding = 1., 
          const std::string& name=std::string(""));

      /// @brief copy constructor
      /// @details it is required to decouple arrays between the input object
      /// and the copy.
      /// @param[in] other input object
      TableVisGridder(const TableVisGridder &other);
                  
      virtual ~TableVisGridder();

      /// @brief Save to a table (for debugging)
      /// @param name Name of table
      void save(const std::string& name);
      

      /// @brief Initialise the gridding
      /// @param axes axes specifications
      /// @param shape Shape of output image: u,v,pol,chan
      /// @param dopsf Make the psf?
      virtual void initialiseGrid(const scimath::Axes& axes,
          const casa::IPosition& shape, const bool dopsf=true);

      /// @brief Grid the visibility data.
      /// @param acc const data accessor to work with
      /// @note a non-const adapter is created behind the scene. If no on-the-fly visibility 
      /// correction is performed, this adapter is equivalent to the original const data accessor
      virtual void grid(accessors::IConstDataAccessor& acc);
      
      /// Form the final output image
      /// @param out Output double precision image or PSF
      virtual void finaliseGrid(casa::Array<double>& out);
      
      /// @brief store given grid
      /// @details This is a helper method for debugging, it stores the amplitude of a given
      /// grid into a CASA image (prior to FFT done as part of finaliseGrid)
      /// @param[in] name image name
      /// @param[in] numGrid number of the grid to store
      void storeGrid(const std::string &name, casa::uInt numGrid) const;

      /// @brief Calculate weights image
      /// @details Form the sum of the convolution function squared, 
      /// multiplied by the weights for each different convolution 
      /// function. This is used in the evaluation of the position
      /// dependent sensitivity
      /// @param out Output double precision sum of weights images
      virtual void finaliseWeights(casa::Array<double>& out);

      /// @brief Initialise the degridding
      /// @param axes axes specifications
      /// @param image Input image: cube: u,v,pol,chan
      virtual void initialiseDegrid(const scimath::Axes& axes,
          const casa::Array<double>& image);

      /// @brief Make context-dependant changes to the gridder behaviour
      /// @param context context
      virtual void customiseForContext(const std::string &context);
      
      /// @brief assign weights
      /// @param viswt shared pointer to visibility weights
      virtual void initVisWeights(const IVisWeights::ShPtr &viswt);
      
      /// @brief Degrid the visibility data.
      /// @param[in] acc non-const data accessor to work with  
      virtual void degrid(accessors::IDataAccessor& acc);

      /// @brief Finalise
      virtual void finaliseDegrid();

      /// @brief set or reset flag forcing gridder to use all data for PSF
      /// @details Change itsUseAllDataForPSF
      /// @param[in] useAll new value of the flag
      void inline useAllDataForPSF(const bool useAll) { itsUseAllDataForPSF = useAll;} 
      
      /// @brief set or reset flag forcing gridder to track weights per oversampling plane
      /// @details change itsTrackWeightPerOversamplePlane
      /// @param[in] flag new value of the flag
      void inline trackWeightPerPlane(const bool flag) { itsTrackWeightPerOversamplePlane = flag;}

      /// @brief set the largest angular separation between the pointing centre and the image centre
      /// @details If the threshold is positive, it is interpreted as the largest allowed angular
      /// separation between the beam (feed in the accessor terminology) pointing centre and the
      /// image centre. We need this to allow imaging of a subset of data (i.e. smaller field of view)
      /// and reject all pointings located outside this smaller image. All accessor rows with
      /// pointingDir1 separated from the image centre by more than this threshold will be ignored.
      /// If the threshold is negative (default), no data rejection based on the pointing direction is done.
      /// The gridder is initialised by default with the negative threshold, i.e. all data are used by default.
      /// @param[in] threshold largest allowed angular separation in radians, use negative value to select all data
      void inline maxPointingSeparation(double threshold = -1.) { itsMaxPointingSeparation = threshold; }

      /// @brief set table name to store the CFs to
      /// @details This method makes it possible to enable writing CFs to disk in destructor after the 
      /// gridder is created. The main use case is to allow a better control of this feature in the parallel
      /// environment (we don't want all workers to write CFs)
      /// @param[in] name table name to store the CFs to (or an empty string if CFs are not to be stored)
      void setTableName(const std::string &name) { itsName = name; }
      
      /// @brief check whether the model is empty
      /// @details A simple check allows us to bypass heavy calculations if the input model
      /// is empty (all pixels are zero). This makes sense for degridding only.
      /// @brief true, if the model is empty
      virtual bool isModelEmpty() const; 
      
  protected:
      /// @brief helper method to print CF cache stats in the log
      /// @details This method is largely intended for debugging. It writes down
      /// to the log the support sizes/offsets for all convolution functions in the cache and
      /// summarises the memory taken up by this cache (per gridder).
      void logCFCacheStats() const;
  
  
      /// @brief shape of the grid
      /// @details The could be a number of grids indexed though gIndex (for each row, polarisation and channel). However, all should
      /// have exactly the same shape. 
      /// @return the shape of grid owned by this gridder
      inline const casa::IPosition& shape() const { return itsShape;}


      /// @brief correct visibilities, if necessary
      /// @details This method is intended for on-the-fly correction of visibilities (i.e. 
      /// facet-based correction needed for LOFAR). This method does nothing in this class, but
      /// can be overridden in the derived classes to plug some effect in. The same method is 
      /// used for both gridding and degridding, with the forward parameter used to distinguish
      /// between these two operations. A non-const accessor has to be modified in situ, if a
      /// correction is required. A buffer for read-only visibilities is created on-demand when
      /// rwVisibility method of the accessor is called for the first time.
      /// @param[in] acc non-const accessor with the data to correct, leave it intact if no
      /// correction is required
      /// @param[in] forward true for degridding (image to vis) and false for gridding (vis to image)
      virtual void correctVisibilities(accessors::IDataAccessor &acc, bool forward);
  
      /// @brief initialise sum of weights
      /// @details We keep track the number of times each convolution function is used per
      /// channel and polarisation (sum of weights). This method is made virtual to be able
      /// to do gridder specific initialisation without overriding initialiseGrid.
      /// This method accepts no parameters as itsShape, itsNWPlanes, etc should have already
      /// been initialised by the time this method is called.
      virtual void initialiseSumOfWeights();
      
      /// @brief zero sum of weights
      /// @details This method just sets all values of the current itsSumWeights biffer to zero
      /// without resizing it. It is done like this for a better encapsulation.
      void inline zeroSumOfWeights() { itsSumWeights.set(0.); }
      
      /// @brief resize sum of weights
      /// @details This method is used inside initialiseSumOfWeights and its overrides in 
      /// derived classes. It resizes itsSumWeights to a given number of convolution
      /// functions taking into account channels/polarisations according to itsShape. 
      /// Moving this operation into the separate method allows better data encapsulation
      /// and tracking weights per oversampling plane or per convolution function depending
      /// on the user's choice.
      /// @param[in] numcf number of convolution functions in the cache (before oversampling)
      void resizeSumOfWeights(const int numcf);
      
      /// @brief obtain sum of weights cube
      /// @return const reference to the sum of weights cube
      inline const casa::Cube<double>& sumOfWeights() const { return itsSumWeights;}
      
      /// @brief log unused spectral planes
      /// @details It is handy to write down the channel numbers into log if the sum of weights 
      /// is zero, i.e. if we have no data for this particular channel. This method does it
      /// (it is called from the destructor, unless the gridder object didn't do any gridding).
      void logUnusedSpectralPlanes() const;
      
      /// @brief translate row of the sum of weights cube into convolution function plane
      /// @details If we are tracking weights per oversampling plane, the row of the sum of
      /// weights cube directly corresponds to the plane of the convolution function cache. 
      /// Otherwise, there is a factor of oversampling squared. It is handy to encapsulate
      /// this functionality here, so all derived classes work do not need to make a 
      /// distinction how the weights are tracked.
      /// @param[in] row row of the sum of weights cube
      inline int cfIndexFromSumOfWeightsRow(const int row) const 
          { return itsTrackWeightPerOversamplePlane ? row : itsOverSample*itsOverSample*row; }
      
      /// @brief helper method to initialise frequency mapping
      /// @details Derived gridders may override initialiseGrid and initialiseDegrid. Howerver, 
      /// they still need to be able to initialise frequency axis mapping (between accessor channels
      /// and image cube), which is handled by a private member class. This method initialises the 
      /// mapper using the content of itsShape and itsAxes, which should be set prior to calling this
      /// method.
      void initialiseFreqMapping();
      
      /// @brief helper method to set up cell size
      /// @details Similar action is required to calculate uv-cell size for gridding and degridding.
      /// Moreover, derived gridders may override initialiseGrid and initialiseDegrid and we don't want
      /// to duplicate the code up there. This method calculates uv-cell size for both ra and dec axes
      /// using coordinate information provided. This method also assigns passed axes parameter to itsAxes.
      /// @param[in] axes coordinate system (ra and dec axes are used).
      void initialiseCellSize(const scimath::Axes& axes);
      
      /// @brief gridder configured to calculate PSF?
      /// @details
      /// @return true if this gridder is configured to calculate PSF, false otherwise
      bool inline isPSFGridder() const { return itsDopsf; }
      
      /// @brief configure gridder to calculate PSF or residuals
      /// @details This method is expected to be called from overridden initialiseGrid method
      /// @param[in] dopsf if true, the gridder is configured to calculate PSF, otherwise
      /// a normal residual grid is calculated.
      void inline configureForPSF(bool dopsf) { itsDopsf = dopsf;}
        
      /// @brief obtain the centre of the image
      /// @details This method extracts RA and DEC axes from itsAxes and
      /// forms a direction measure corresponding to the middle of each axis.
      /// @return direction measure corresponding to the image centre
      casa::MVDirection getImageCentre() const;
      
      /// @brief obtain the tangent point
      /// @details For faceting all images should be constructed for the same tangent
      /// point. This method extracts the tangent point (reference position) from the
      /// coordinate system.
      /// @return direction measure corresponding to the tangent point
      casa::MVDirection getTangentPoint() const;
      
      // data members should be made private in the future!

      /// Axes definition for image
      askap::scimath::Axes itsAxes;

      /// Shape of image
      casa::IPosition itsShape;

      /// Cell sizes in wavelengths
      casa::Vector<double> itsUVCellSize;

//temporary comment out
//private:
      /// @brief Sum of weights (axes are index, pol, chan) 
      casa::Cube<double> itsSumWeights;
protected:      

      /// @brief Convolution function
      /// The convolution function is stored as a vector of arrays so that we can
      /// use any of a number of functions. The index is calculated by cIndex.
      std::vector<casa::Matrix<casa::Complex> > itsConvFunc;

      /// @brief Obtain offset for the given convolution function
      /// @details To conserve memory and speed the gridding up, convolution functions stored in the cache
      /// may have an offset (i.e. essentially each CF should be defined on a bigger support and placed at a 
      /// pixel other than the centre of this support). This method returns this offset, which is the
      /// position of the peak of the given CF on a bigger support w.r.t. the centre of this support. 
      /// The value of (0,0) means no offset from the centre (i.e. support is already centred). 
      /// @param[in] cfPlane plane of the convolution function cache to get the offset for
      /// @return a pair with offsets for each axis
      /// @note if there is no offset defined for a given cfPlane (default behavior), this method returns (0,0)
      std::pair<int,int> getConvFuncOffset(int cfPlane) const;
      
      /// @brief initialise convolution function offsets for a given number of planes
      /// @details The vector with offsets is resized and filled with (0,0).
      /// @param[in] nPlanes number of planes in the cache 
      void initConvFuncOffsets(size_t nPlanes);
      
      /// @brief Assign offset to a particular convolution function
      /// @details
      /// @param[in] cfPlane plane of the convolution function cache to assign the offset for
      /// @param[in] x offset in the first coordinate
      /// @param[in] y offset in the second coordinate
      /// @note For this method, cfPlane should be within the range [0..nPlanes-1].
      void setConvFuncOffset(int cfPlane, int x, int y);
      
      /// @brief assign a given offset to the CF plane
      /// @details 

      /// Return the index into the convolution function for a given
      /// row, polarisation, and channel
      /// @param row Row of accessor
      /// @param pol Polarisation
      /// @param chan Channel
      virtual int cIndex(int row, int pol, int chan);

      /// Return the index into the grid for a given
      /// row and channel
      /// @param row Row of accessor
      /// @param pol Polarisation
      /// @param chan Channel
      virtual int gIndex(int row, int pol, int chan);

      /// @brief Initialize the convolution function - this is the key function to override.
      /// @param[in] acc const accessor to work with
      virtual void initConvolutionFunction(const accessors::IConstDataAccessor& acc) = 0;

      /// @brief Initialise the indices
      /// @param[in] acc const accessor to work with
      virtual void initIndices(const accessors::IConstDataAccessor& acc) = 0;

      /// @brief Correct for gridding convolution function
      /// @param image image to be corrected
      virtual void correctConvolution(casa::Array<double>& image) = 0;

      /// @brief Conversion helper function
      /// @details Copies in to out expanding double into complex values and
      /// padding appropriately if necessary (padding is more than 1)
      /// @param[out] out complex output array
      /// @param[in] in double input array
      /// @param[in] padding padding factor
      static void toComplex(casa::Array<casa::DComplex>& out, const casa::Array<double>& in, 
                     const float padding = 1.);

      /// @brief Conversion helper function
      /// @details Copies real part of in into double array and
      /// extracting an inner rectangle if necessary (padding is more than 1)
      /// @param[out] out real output array
      /// @param[in] in complex input array      
      /// @param[in] padding padding factor
      static void toDouble(casa::Array<double>& out, const casa::Array<casa::DComplex>& in,
                    const float padding = 1.);

      /// @brief a helper method to initialize gridding of the PSF
      /// @details The PSF is calculated using the data for a
      /// representative field/feed only. By default, the first encountered
      /// feed/field is chosen. If the same gridder is reused for another
      /// sequence of data points a new representative feed/field have to be
      /// found. This is done by resetting the cache in initialiseGrid. However,
      /// the latter method can be overridden in the derived classes. To avoid
      /// a duplication of the code, this helper method resets the representative
      /// feed/field cache. It is called from initialiseGrid.
      void initRepresentativeFieldAndFeed();
      
      /// @brief set up itsStokes using the information from itsAxes and itsShape
      void initStokes();
      
      /// @brief obtain stokes for each plane of the current grid
      /// @details The output of this method has a meaning only after initialiseGrid or
      /// initialiseDegrid has been called.
      inline const casa::Vector<casa::Stokes::StokesTypes>& getStokes() const {return itsStokes;}

      /// Support of convolution function
      int itsSupport;
      /// Oversampling of convolution function
      int itsOverSample;

      /// Is the model empty? Used to shortcut degridding
      bool itsModelIsEmpty;

      /// The grid is stored as a cube as well so we can index into that as well.
      std::vector<casa::Array<casa::Complex> > itsGrid;
            
  private:

      /// @brief return the table name to store the result to
      /// @details This method could probably be made private as it is not used outside
      /// this class
      /// @return table name to store the CFs to (or an empty string if CFs are not to be stored)
      std::string tableName() const { return itsName; }

      /// Name of table to save to
      std::string itsName;

      /// @brief assignment operator
      /// @details it is required to decouple arrays between the input object
      /// and the copy.
      /// @param[in] other input object
      /// @return reference to itself
      /// @note assignment operator is made private, so wouldn't need to support both copy constructor and
      /// assignment operator. In the case of inadvertent use, compiler should give an error
      TableVisGridder& operator=(const TableVisGridder &other);
  
      /// @brief polarisation frame for the grid
      /// @details Assumed to be the same for all elements of itsGrid vector.
      /// This field is filled in initialiseGrid or initialiseDegrid using the Axes 
      /// object.
      casa::Vector<casa::Stokes::StokesTypes> itsStokes;
  
      /// Number of samples gridded
      double itsSamplesGridded;
      /// Number of samples degridded
      double itsSamplesDegridded;
      /// Number of flagged visibility vectors (all pols.)
      double itsVectorsFlagged;
      /// Number of grid cells gridded
      double itsNumberGridded;
      /// Number of grid cells degridded
      double itsNumberDegridded;
      /// Time for Coordinates
      double itsTimeCoordinates;
      /// Time for convolution functions
      double itsTimeConvFunctions;
      /// Time for gridding
      double itsTimeGridded;
      /// Time for degridding
      double itsTimeDegridded;

      /// @brief is this gridder a PSF gridder?
      bool itsDopsf;
        
      /// Generic grid/degrid - this is the heart of this framework. It should never
      /// need to be overridden
      /// @param[in] acc non-const data accessor to work with.  
      /// @param[in] forward true for the model to visibility transform (degridding),
      /// false for the visibility to dirty image transform (gridding)
      /// @note We have to pass a non-const accessor because this method can either 
      /// write or read. A bit better re-structuring of the code can help to deal with
      /// constness properly.
      void generic(accessors::IDataAccessor& acc, bool forward);

      /// Visibility Weights
      IVisWeights::ShPtr itsVisWeight;

      /// @brief true if no visibilities have been gridded since the last initialize
      /// @details For PSF calculations we need to take just the first feed and 
      /// field (it is an approximation that they all considered the same). To be
      /// able to extend this check over multiple calls of the generic routine this
      /// flag is used. It is set to true in initialise and then reset to false when
      /// the first visibility is gridded. 
      bool itsFirstGriddedVis;

      /// @brief an index of the feed, which provides data for the PSF calculations
      /// @details This data member is initialized when the first visibility is gridded,
      /// only this feed is used to calculate the PSF.
      casa::uInt itsFeedUsedForPSF;

      /// @brief pointing direction, which provides data for the PSF calculations
      /// @details This data member is initialized when the first visibility is gridded,
      /// only this field is used to calculate the PSF
      casa::MVDirection itsPointingUsedForPSF;
      
      /// @brief use all data for PSF calculation
      /// @details By default we use just a representative feed and field to calculate PSF. 
      /// For research purposes we need an option which allows to take all available data into
      /// account. This is a flag showing that itsFeedUsedForPSF and itsPointingUsedForPSF are ignored.
      /// Default value is false.
      bool itsUseAllDataForPSF;     
      
      /// @brief mapping class between image planes and accessor channels
      /// @details Correspondence between planes of the image cube and accessor channels may be
      /// non-trivial. This class takes care of the mapping.
      FrequencyMapper itsFreqMapper;
      
      /// @brief largest angular separation between the pointing centre and the image centre
      /// @details If the value is positive, it is interpreted as the largest allowed angular
      /// separation between the beam (feed in the accessor terminology) pointing centre and the
      /// image centre. It is intended to allow imaging of a subset of data (i.e. smaller field of view)
      /// and reject all pointings located outside this smaller image. All accessor rows with
      /// pointingDir1 separated from the image centre by more than this threshold are ignored.
      /// If the value is negative, no data rejection based on the pointing direction is done.
      /// Values are in radians.
      double itsMaxPointingSeparation;
      
      /// @brief number of rows rejected due to itsMaxPointingSeparation
      /// @details Accumulated to get proper statistics/debugging info.
      long itsRowsRejectedDueToMaxPointingSeparation;
      
      /// @brief offsets of convolution functions
      /// @details To conserve memory and speed the gridding up, convolution functions stored in the cache
      /// may have an offset (i.e. essentially each CF should be defined on a bigger support and placed at a 
      /// pixel other than the centre of this support). This data field defines this offset, which is the
      /// position of the peak of the given CF on a bigger support w.r.t. the centre of this support. 
      /// The value of (0,0) means no offset from the centre (i.e. support is already centred). If the index
      /// of CF is beyond the size of this vector, (0,0) offset is assumed. By default this vector is empty,
      /// which means no offset. 
      std::vector<std::pair<int,int> > itsConvFuncOffsets;
      
      /// @brief true, if itsSumWeights tracks weights per oversampling plane
      bool itsTrackWeightPerOversamplePlane;

      #ifdef _OPENMP
      /// @brief synchronisation mutex
      mutable boost::mutex itsMutex;
      #endif
    };
  }
}
#endif
