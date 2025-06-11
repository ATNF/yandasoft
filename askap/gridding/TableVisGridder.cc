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

#include <askap/gridding/TableVisGridder.h>
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
#include <askap/dataaccess/IDataAccessor.h>
#include <askap/dataaccess/OnDemandBufferDataAccessor.h>
ASKAP_LOGGER(logger, ".gridding.tablevisgridder");

#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/scimath/fft/FFT2DWrapper.h>

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/casa/Arrays/Slice.h>

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/UVWMachine.h>

#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/ParamsCasaTable.h>

#include <askap/gridding/GridKernel.h>

#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/measurementequation/ImageParamsHelper.h>
#include <askap/scimath/utils/ImageUtils.h>

#include <askap/askap/CasaSyncHelper.h>
#include <askap/profile/AskapProfiler.h>

using namespace askap::scimath;
using namespace askap;

#include <ostream>
#include <sstream>
#include <iomanip>

#include <casacore/casa/OS/Timer.h>

namespace askap {
namespace synthesis {

/// @brief a helper method for a deep copy of casa arrays held in
/// stl vectors
/// @param[in] in input array
/// @param[out] out output array (will be resized)
template<typename T>
void deepCopyOfSTDVector(const std::vector<T> &in,
                         std::vector<T> &out)
{
   out.resize(in.size());

   const typename std::vector<T>::const_iterator inEnd = in.end();
   typename std::vector<T>::iterator outIt = out.begin();
   for (typename std::vector<T>::const_iterator inIt = in.begin();
       inIt != inEnd; ++inIt,++outIt) {
       *outIt = inIt->copy();
   }
}

/// @brief required to mediate thread safety issues of the casa cube
utility::CasaSyncHelper syncHelper;

// DDCALTAG -- itsSourceIndex added to all of the constructors
TableVisGridder::TableVisGridder() : itsSumWeights(),
    itsSupport(-1), itsOverSample(-1),
    itsModelIsEmpty(false), itsSamplesGridded(0), itsSamplesDegridded(0),
    itsVectorsFlagged(0), itsVectorsWFlagged(0), itsNumberGridded(0), itsNumberDegridded(0),
    itsTimeCoordinates(0.0), itsTimeConvFunctions(0.0), itsTimeGridded(0.0),
    itsTimeDegridded(0.0), itsDopsf(false), itsDopcf(false),
    itsFirstGriddedVis(true), itsFeedUsedForPSF(0), itsUseAllDataForPSF(false),
    itsMaxPointingSeparation(-1.), itsRowsRejectedDueToMaxPointingSeparation(0),
    itsTrackWeightPerOversamplePlane(false),itsPARotation(false),itsSwapPols(false),
    its2dGrid(),itsVisPols(),itsPolConv(),itsImagePolFrameVis(),itsImagePolFrameNoise(),
    itsPolVector(),itsImageChan(-1),itsGridIndex(-1),itsSourceIndex(0),itsClearGrid(false)
{
}

TableVisGridder::TableVisGridder(const int overSample, const int support,
    const float padding, const std::string& name) : VisGridderWithPadding(padding),
    itsSumWeights(), itsSupport(support), itsOverSample(overSample),
    itsModelIsEmpty(false), itsName(name), itsSamplesGridded(0), itsSamplesDegridded(0),
    itsVectorsFlagged(0), itsVectorsWFlagged(0), itsNumberGridded(0), itsNumberDegridded(0),
    itsTimeCoordinates(0.0), itsTimeConvFunctions(0.0), itsTimeGridded(0.0),
    itsTimeDegridded(0.0), itsDopsf(false), itsDopcf(false),
    itsFirstGriddedVis(true), itsFeedUsedForPSF(0), itsUseAllDataForPSF(false),
    itsMaxPointingSeparation(-1.), itsRowsRejectedDueToMaxPointingSeparation(0),
    itsTrackWeightPerOversamplePlane(false),itsPARotation(false),itsSwapPols(false),
    its2dGrid(),itsVisPols(),itsPolConv(),itsImagePolFrameVis(),itsImagePolFrameNoise(),
    itsPolVector(),itsImageChan(-1),itsGridIndex(-1),itsSourceIndex(0),itsClearGrid(false)
{
   ASKAPCHECK(overSample>0, "Oversampling must be greater than 0");
   ASKAPCHECK(support>=0, "Maximum support must be zero or greater");
   ASKAPCHECK(padding>0, "Padding factor must be greater than 0");
}

    /// @brief copy constructor
    /// @details it is required to decouple arrays between the input object
    /// and the copy.
    /// @param[in] other input object
TableVisGridder::TableVisGridder(const TableVisGridder &other) :
     IVisGridder(other), VisGridderWithPadding(other),
     itsAxes(other.itsAxes), itsShape(other.itsShape),
     itsUVCellSize(other.itsUVCellSize.copy()),
     itsSumWeights(other.itsSumWeights.copy()),
     itsSupport(other.itsSupport), itsOverSample(other.itsOverSample),
     itsModelIsEmpty(other.itsModelIsEmpty),
     itsName(other.itsName), itsSamplesGridded(other.itsSamplesGridded),
     itsSamplesDegridded(other.itsSamplesDegridded),
     itsVectorsFlagged(other.itsVectorsFlagged),
     itsVectorsWFlagged(other.itsVectorsWFlagged),
     itsNumberGridded(other.itsNumberGridded),
     itsNumberDegridded(other.itsNumberDegridded), itsTimeCoordinates(other.itsTimeCoordinates),
     itsTimeConvFunctions(other.itsTimeConvFunctions),
     itsTimeGridded(other.itsTimeGridded),
     itsTimeDegridded(other.itsTimeDegridded),
     itsDopsf(other.itsDopsf),
     itsDopcf(other.itsDopcf),
     // MV: it looks like there is some technical debt / untidy design here. There is clearly an intention to do a real copy for data fields
     // held by reference. However, we usually clone empty gridders (essentially, clone not to copy objects but to create from a template), so
     // there is no duplication of data. With this usage in mind, it makes sense to have reference semantics for the weight accessor and builder
     // fields (i.e. only shared pointers are copied)
     itsUVWeightAccessor(other.itsUVWeightAccessor),
     itsUVWeightBuilder(other.itsUVWeightBuilder),
     //
     itsFirstGriddedVis(other.itsFirstGriddedVis),
     itsFeedUsedForPSF(other.itsFeedUsedForPSF),
     itsPointingUsedForPSF(other.itsPointingUsedForPSF),
     itsUseAllDataForPSF(other.itsUseAllDataForPSF),
     itsFreqMapper(other.itsFreqMapper),
     itsMaxPointingSeparation(other.itsMaxPointingSeparation),
     itsRowsRejectedDueToMaxPointingSeparation(other.itsRowsRejectedDueToMaxPointingSeparation),
     itsConvFuncOffsets(other.itsConvFuncOffsets),
     itsTrackWeightPerOversamplePlane(other.itsTrackWeightPerOversamplePlane),
     itsPARotation(other.itsPARotation), itsPARotAngle(other.itsPARotAngle),
     itsSwapPols(other.itsSwapPols),
     its2dGrid(other.its2dGrid.copy()),itsVisPols(other.itsVisPols.copy()),
     itsPolConv(other.itsPolConv),itsImagePolFrameVis(other.itsImagePolFrameVis.copy()),
     itsImagePolFrameNoise(other.itsImagePolFrameNoise.copy()),itsPolVector(other.itsPolVector.copy()),
     itsImageChan(other.itsImageChan),itsGridIndex(other.itsGridIndex),itsSourceIndex(other.itsSourceIndex),
     itsClearGrid(other.itsClearGrid)
{
   deepCopyOfSTDVector(other.itsConvFunc,itsConvFunc);
   deepCopyOfSTDVector(other.itsGrid, itsGrid);
   if(other.itsVisWeight) {
      itsVisWeight = other.itsVisWeight->clone();
   } else {
      itsVisWeight = other.itsVisWeight;
   }
}

/// @brief assignment operator
/// @details it is required to decouple arrays between the input object
/// and the copy.
/// @param[in] other input object
/// @return reference to itself
TableVisGridder& TableVisGridder::operator=(const TableVisGridder &)
{
  ASKAPTHROW(AskapError, "This method is not supposed to be called!");
  return *this;
}


TableVisGridder::~TableVisGridder() {

    if (itsNumberGridded>0) {
        if (isPSFGridder()) {
            ASKAPLOG_INFO_STR(logger, "TableVisGridder PSF gridding statistics");
            ASKAPLOG_INFO_STR(logger, "   PSF samples gridded       = "
                              << itsSamplesGridded);
            ASKAPLOG_INFO_STR(logger, "   Visibility vectors flagged (psf)     = "
                              << itsVectorsFlagged);
            ASKAPLOG_INFO_STR(logger, "   Visibility vectors with |w|>wmax (psf) = "
                              << itsVectorsWFlagged);
            ASKAPLOG_INFO_STR(logger, "   Total time for PSF gridding   = "
                << itsTimeGridded << " (s)");
            ASKAPLOG_INFO_STR(logger, "   PSF gridding time         = " << 1e6
                 *itsTimeGridded/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   Total time converting for PSF = "
                << itsTimeCoordinates << " (s)");
            ASKAPLOG_INFO_STR(logger, "   Total time building CFs and indices for PSF = "
                << itsTimeConvFunctions << " (s)");
            ASKAPLOG_INFO_STR(logger, "   PSF coord conversion      = "
                  << 1e6 * itsTimeCoordinates/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   PSF CFs and indices      = "
                  << 1e6 * itsTimeConvFunctions/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   " << GridKernel::info());
            ASKAPLOG_INFO_STR(logger, "   Points gridded (psf)        = "
                  << itsNumberGridded);
            ASKAPLOG_INFO_STR(logger, "   Time per point (psf)        = " << 1e9
                  *itsTimeGridded/itsNumberGridded << " (ns)");
            ASKAPLOG_INFO_STR(logger, "   Performance for PSF         = "
                << 8.0 * 1e-9 * itsNumberGridded/itsTimeGridded << " GFlops");
        }
        else if (isPCFGridder()) {
            ASKAPLOG_INFO_STR(logger, "TableVisGridder PreConditioner Function gridding statistics");
            ASKAPLOG_INFO_STR(logger, "   PCF samples gridded       = "
                              << itsSamplesGridded);
            ASKAPLOG_INFO_STR(logger, "   Visibility vectors flagged (pcf)     = "
                              << itsVectorsFlagged);
            ASKAPLOG_INFO_STR(logger, "   Visibility vectors with |w|>wmax (pcf) = "
                              << itsVectorsWFlagged);
            ASKAPLOG_INFO_STR(logger, "   Total time for PCF gridding   = "
                << itsTimeGridded << " (s)");
            ASKAPLOG_INFO_STR(logger, "   PCF gridding time         = " << 1e6
                 *itsTimeGridded/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   Total time converting for PCF = "
                << itsTimeCoordinates << " (s)");
            ASKAPLOG_INFO_STR(logger, "   Total time building CFs and indices for PCF = "
                << itsTimeConvFunctions << " (s)");
            ASKAPLOG_INFO_STR(logger, "   PCF coord conversion      = "
                  << 1e6 * itsTimeCoordinates/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   PCF CFs and indices      = "
                  << 1e6 * itsTimeConvFunctions/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   " << GridKernel::info());
            ASKAPLOG_INFO_STR(logger, "   Points gridded (pcf)        = "
                  << itsNumberGridded);
            ASKAPLOG_INFO_STR(logger, "   Time per point (pcf)        = " << 1e9
                  *itsTimeGridded/itsNumberGridded << " (ns)");
            ASKAPLOG_INFO_STR(logger, "   Performance for PCF         = "
                << 8.0 * 1e-9 * itsNumberGridded/itsTimeGridded << " GFlops");
        }
        else {
            ASKAPLOG_INFO_STR(logger, "TableVisGridder gridding statistics");
            ASKAPLOG_INFO_STR(logger, "   Samples gridded       = "
                         << itsSamplesGridded);
            ASKAPLOG_INFO_STR(logger, "   Visibility vectors flagged       = "
                          << itsVectorsFlagged);
             ASKAPLOG_INFO_STR(logger, "   Visibility vectors with |w|>wmax = "
                           << itsVectorsWFlagged);
            ASKAPLOG_INFO_STR(logger, "   Total time gridding   = "
                         << itsTimeGridded << " (s)");
            ASKAPLOG_INFO_STR(logger, "   Gridding time         = " << 1e6
                  *itsTimeGridded/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   Total time converting = "
                << itsTimeCoordinates << " (s)");
            ASKAPLOG_INFO_STR(logger, "   Total time building CFs and indices = "
                << itsTimeConvFunctions << " (s)");
            ASKAPLOG_INFO_STR(logger, "   Coord conversion      = "
                  << 1e6 * itsTimeCoordinates/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   CFs and indices      = "
                  << 1e6 * itsTimeConvFunctions/itsSamplesGridded << " (us) per sample");
            ASKAPLOG_INFO_STR(logger, "   " << GridKernel::info());
            ASKAPLOG_INFO_STR(logger, "   Points gridded        = "
                << itsNumberGridded);
            ASKAPLOG_INFO_STR(logger, "   Time per point        = " << 1e9
                *itsTimeGridded/itsNumberGridded << " (ns)");
            ASKAPLOG_INFO_STR(logger, "   Performance           = "
                << 8.0 * 1e-9 * itsNumberGridded/itsTimeGridded << " GFlops");
        }
    }
    if (itsNumberDegridded>0) {
        ASKAPLOG_INFO_STR(logger, "TableVisGridder degridding statistics");
        ASKAPLOG_INFO_STR(logger, "   Samples degridded     = "
                << itsSamplesDegridded);
        ASKAPLOG_INFO_STR(logger, "   Total time degridding = "
                << itsTimeDegridded << " (s)");
        ASKAPLOG_INFO_STR(logger, "   Degridding time       = " << 1e6
                *itsTimeDegridded/itsSamplesDegridded << " (us) per sample");
        ASKAPLOG_INFO_STR(logger, "   Total time converting = "
                << itsTimeCoordinates << " (s)");
        ASKAPLOG_INFO_STR(logger, "   Total time building CFs and indices = "
                << itsTimeConvFunctions << " (s)");
        ASKAPLOG_INFO_STR(logger, "   Coord conversion      = "
                  << 1e6 * itsTimeCoordinates/itsSamplesDegridded << " (us) per sample");
        ASKAPLOG_INFO_STR(logger, "   CFs and indices      = "
                  << 1e6 * itsTimeConvFunctions/itsSamplesDegridded << " (us) per sample");
        ASKAPLOG_INFO_STR(logger, "   " << GridKernel::info());
        ASKAPLOG_INFO_STR(logger, "   Points degridded      = "
                << itsNumberDegridded);
        ASKAPLOG_INFO_STR(logger, "   Time per point        = " << 1e9
                *itsTimeDegridded/itsNumberDegridded << " (ns)");
        ASKAPLOG_INFO_STR(logger, "   Performance           = "
                << 8.0 * 1e-9 * itsNumberDegridded/itsTimeDegridded << " GFlops");
    }
    if (itsMaxPointingSeparation > 0.) {
        ASKAPLOG_INFO_STR(logger, "   Samples rejected due to MaxPointingSeparation = "<<
                                   itsRowsRejectedDueToMaxPointingSeparation);
    }
    if(unusedGridder()) {
      ASKAPLOG_WARN_STR(logger, "Unused gridder");
      if (itsRowsRejectedDueToMaxPointingSeparation != 0) {
          ASKAPLOG_WARN_STR(logger, "It looks like all samples were rejected due to MaxPointingSeparation!");
      }
    }
    else {
      ASKAPLOG_INFO_STR(logger, "   Padding factor    = " << paddingFactor());
      if (itsTrackWeightPerOversamplePlane) {
          ASKAPLOG_INFO_STR(logger, "   Weights were tracked per oversampling plane");
      } else {
          ASKAPLOG_INFO_STR(logger, "   First oversampling plane was used to compute weights");
      }
      logUnusedSpectralPlanes();
      logCFCacheStats();
      if(tableName() != "") {
        save(tableName());
      }
    }
}

void TableVisGridder::save(const std::string& name) {
    if (name.find("image:") != 0) {
        askap::scimath::ParamsCasaTable iptable(name, false);
        askap::scimath::Params ip;
        ASKAPLOG_DEBUG_STR(logger, "Saving " << itsConvFunc.size() << " entries in convolution function");
        for (unsigned int i=0; i<itsConvFunc.size(); i++) {

            casacore::Array<float> realC=real(itsConvFunc[i]);
            // ASKAPLOG_DEBUG_STR(logger, "Entry[" <<  i <<  "] has shape " <<  itsConvFunc[i].shape());
            std::ostringstream os;
            os<<"Real.Convolution";
            os.width(5);
            os.fill('0');
            os<<i;
            ip.add(os.str(), realC);
        }
        iptable.setParameters(ip);
    } else {
        if (itsNumberGridded == 0) {
            ASKAPLOG_DEBUG_STR(logger, "Ignore tablename="<<name<<
                               " option as no visibilities were gridded");
            return;
        }
        if (isPSFGridder()) {
            ASKAPLOG_DEBUG_STR(logger, "Ignore tablename="<<name<<
                               " option for the PSF gridder");
             return;
        } else if (isPCFGridder()) {
            ASKAPLOG_DEBUG_STR(logger, "Ignore tablename="<<name<<
                               " option for the Preconditioner Function gridder");
             return;
        }

        const std::string imgName = name.substr(6);
        // number of planes before oversampling
        const unsigned long nPlanes = itsConvFunc.size()/itsOverSample/itsOverSample;
        if (nPlanes > 0) {
            ASKAPLOG_DEBUG_STR(logger, "Saving convolution functions into a cube "<<
                                  imgName<<" with " << nPlanes<<
                                  " planes (first oversampling plane only)");
            ASKAPDEBUGASSERT(itsConvFunc.size()>0);
            int support = -1;
            for (unsigned int plane = 0; plane<nPlanes; ++plane) {
                 const int x = int(itsConvFunc[plane*itsOverSample*itsOverSample].nrow());
                 const int y = int(itsConvFunc[plane*itsOverSample*itsOverSample].ncolumn());
                 if (support < x) {
                     support = x;
                 }
                 if (support < y) {
                     support = y;
                 }
            }
            casacore::Cube<casacore::Float> imgBuffer(support, support, nPlanes);
            imgBuffer.set(0.);
            for (unsigned int plane = 0; plane<nPlanes; ++plane) {
                //unsigned int peakX = 0, peakY = 0;
                casacore::Float peakVal = -1.;
                for (int x = 0; x<int(imgBuffer.nrow()); ++x) {
                     for (int y = 0; y<int(imgBuffer.ncolumn()); ++y) {
                          const casacore::Matrix<casacore::Complex> thisCF = itsConvFunc[plane*itsOverSample*itsOverSample];
                          const int xOff = (support - int(thisCF.nrow()))/2;
                          const int yOff = (support - int(thisCF.ncolumn()))/2;
                          ASKAPDEBUGASSERT((xOff >= 0) && (yOff >= 0));
                          if ( (x - xOff >= int(thisCF.nrow())) || (y - yOff >= int(thisCF.ncolumn()))) {
                               continue;
                          }
                          if ( (x - xOff < 0) || (y - yOff < 0) ) {
                               continue;
                          }
                          imgBuffer(x,y,plane) = casacore::abs(thisCF(x - xOff,y - yOff));
                     }
                }
            }
            scimath::saveAsCasaImage(imgName,imgBuffer);
        }
    }
}

/// @brief helper method to print CF cache stats in the log
/// @details This method is largely intended for debugging. It writes down
/// to the log the support sizes/offsets for all convolution functions in the cache and
/// summarises the memory taken up by this cache (per gridder).
void TableVisGridder::logCFCacheStats() const
{
   // number of planes before oversampling
   ASKAPDEBUGASSERT(itsOverSample>0);
   const unsigned long nPlanes = itsConvFunc.size()/itsOverSample/itsOverSample;
   unsigned long memUsed = 0;
   for (unsigned int plane = 0; plane<nPlanes; ++plane) {
       ASKAPDEBUGASSERT(plane*itsOverSample*itsOverSample < itsConvFunc.size());
       const casacore::IPosition shape = itsConvFunc[plane*itsOverSample*itsOverSample].shape();
       if (shape.nelements() < 2) {
           ASKAPLOG_DEBUG_STR(logger, "CF plane="<<plane<<" (before oversampling) is malformed");
           continue;
       }
       if (shape.product() == 0) {
           // ASKAPLOG_DEBUG_STR(logger, "CF plane="<<plane<<" (before oversampling) is unused");
           memUsed += sizeof(casacore::Matrix<casacore::Complex>);
           continue;
       }
       if ((shape[0] != shape[1]) || (shape[0] % 2 != 1)) {
           ASKAPLOG_DEBUG_STR(logger, "CF plane="<<plane<<
                                      " (before oversampling) has a rectangular support or even size");
           memUsed += sizeof(casacore::Matrix<casacore::Complex>) +
                      sizeof(casacore::Complex)*shape.product()*itsOverSample*itsOverSample;
           continue;
       }
       const int support = (shape[0] - 1) / 2;
       const std::pair<int,int> cfOffset = getConvFuncOffset(plane);
       ASKAPLOG_DEBUG_STR(logger, "CF plane="<<plane<<" (before oversampling): support="<<
                                  support<<", size="<<shape[0]<<
                                  " at offset ("<<cfOffset.first<<","<<cfOffset.second<<")");
       memUsed += sizeof(casacore::Matrix<casacore::Complex>) +
                  sizeof(casacore::Complex)*shape[0]*shape[1]*itsOverSample*itsOverSample;
   }
   if (nPlanes > 0) {
       float effectiveSize = (float(memUsed)-float(sizeof(casacore::Matrix<casacore::Complex>)*nPlanes)) /
                             (sizeof(casacore::Complex)*itsOverSample*itsOverSample*nPlanes);
       ASKAPLOG_DEBUG_STR(logger, "Cache of convolution functions take " <<
                                  float(memUsed)/1024/1024<<" Mb of memory or ");
       ASKAPLOG_DEBUG_STR(logger, float(memUsed)/nPlanes/1024/1024 <<
                                  " Mb of memory per plane (before oversampling)");
       ASKAPDEBUGASSERT(effectiveSize>=0.);
       effectiveSize = sqrt(effectiveSize);
       ASKAPLOG_DEBUG_STR(logger, "Effective CF size (in terms of memory usage) is "<<
                                  long(effectiveSize)<<", effective support="<<
                                  long((effectiveSize-1)/2));
   }
}


/// This is a generic grid/degrid
void TableVisGridder::generic(accessors::IDataAccessor& acc, bool forward) {
   ASKAPDEBUGTRACE("TableVisGridder::generic");
  if (forward&&itsModelIsEmpty)
        return;

   if (forward && isPSFGridder()) {
       ASKAPTHROW(AskapError,
           "Logic error: the gridder is not supposed to be used for degridding in the PSF mode")
   }
   if (forward && isPCFGridder()) {
       ASKAPTHROW(AskapError,
           "Logic error: the gridder is not supposed to be used for degridding in the PCF mode")
   }

   casacore::Timer timer;

   // Time CFs and indices
   timer.mark();

   initIndices(acc);
   initConvolutionFunction(acc);
   if (!forward) {
      ASKAPCHECK(itsSumWeights.nelements()>0, "SumWeights not yet initialised");
   }

   itsTimeConvFunctions += timer.real();

   // Now time the coordinate conversions, etc.
   // some conversion may have already happened during CF calculation
   timer.mark();
   const casacore::MVDirection imageCentre = getImageCentre();
   const casacore::MVDirection tangentPoint = getTangentPoint();

   // its fine to work with the reference in the openmp case because all our current use cases
   // have the same tangent point for all gridders, otherwise we have to move this call
   // inside the section protected by the lock and make a copy of the returned vector
   const casacore::Vector<casacore::RigidVector<double, 3> > &outUVW = acc.rotatedUVW(tangentPoint);

   #ifdef _OPENMP_WORKING
   boost::unique_lock<boost::mutex> lock(itsMutex);
   const casa::Vector<double> delay = acc.uvwRotationDelay(tangentPoint, imageCentre).copy();
   lock.unlock();
   #else
   const casa::Vector<double> &delay = acc.uvwRotationDelay(tangentPoint, imageCentre);
   #endif

   itsTimeCoordinates += timer.real();

   // Now time the gridding
   timer.mark();

   ASKAPCHECK(itsSupport>=0, "Support must be zero or greater");
   ASKAPCHECK(itsUVCellSize.size()==2, "UV cell sizes not yet set");

   const uint nSamples = acc.nRow();
   const uint nChan = acc.nChannel();
   const uint nPol = acc.nPol();
   const casacore::Vector<casacore::Double>& frequencyList = acc.frequency();
   itsFreqMapper.setupMapping(frequencyList);

   // for now setup the converter inside this method, although it would cause a rebuild
   // of the matrices for every accessor. More intelligent caching is possible with a bit
   // more effort (i.e. one has to detect whether polarisation frames change from the
   // previous call). Need to think about parallactic angle dependence.
   // OPENMP case needs work...
   #ifdef _OPENMP_WORKING
   itsPolConv = (forward ? : scimath::PolConverter(getStokes(),syncHelper.copy(acc.stokes()), false)
                             scimath::PolConverter(syncHelper.copy(acc.stokes()), getStokes()));
   //scimath::PolConverter degridPolConv(getStokes(),syncHelper.copy(acc.stokes()), false);
   #else
   if (nPol != itsVisPols.nelements()  || !allEQ(acc.stokes(), itsVisPols)) {
     itsPolConv = (forward ? scimath::PolConverter(getStokes(),acc.stokes(), false) :
                             scimath::PolConverter(acc.stokes(), getStokes()));
     itsVisPols.assign(acc.stokes());
     itsPolVector.resize(nPol);
   }
   #endif


   ASKAPDEBUGASSERT(itsShape.nelements()>=2);
   casacore::IPosition ipStart(4, 0, 0, 0, 0);
   const casacore::IPosition onePlane(2,shape()(0),shape()(1));
   if (its2dGrid.shape()(0)!=onePlane(0) || its2dGrid.shape()(1)!=onePlane(1)) {
       // Initialise the 2d reference
       its2dGrid.resize(onePlane);
   }
   // number of polarisation planes in the grid
   const casacore::uInt nImagePols = (shape().nelements()<=2) ? 1 : shape()[2];
   if (itsImagePolFrameVis.nelements()!=nImagePols) {
        itsImagePolFrameVis.resize(nImagePols);
        itsImagePolFrameNoise.resize(nImagePols);
   }
   // a buffer for the visibility vector in the polarisation frame used for the grid
   //if (itsPolVector.nelements()!=nPol) itsPolVector.resize(nPol);

   // Loop over all samples adding them to the grid
   // First scale to the correct pixel location
   // Then find the fraction of a pixel to the nearest pixel
   // Loop over the entire itsSupport, calculating weights from
   // the convolution function and adding the scaled
   // visibility to the grid.

   ASKAPDEBUGASSERT(casa::uInt(nChan) <= frequencyList.nelements());
   ASKAPDEBUGASSERT(casa::uInt(nSamples) == acc.uvw().nelements());
   const casa::Cube<casa::Bool>& flagCube = acc.flag();
   // We want these for either gridding or degridding and
   // want to avoid calling them in the loop due to virtual function overheads.
   // Use shared_ptr without delete to hold the pointers.
   boost::shared_ptr<casa::Cube<casa::Complex> > visCube;
   boost::shared_ptr<const casa::Cube<casa::Complex> > roVisCube;
   boost::shared_ptr<const casa::Cube<casa::Complex> > roVisNoise;
   if (forward) {
       visCube.reset(&acc.rwVisibility(), utility::NullDeleter());
   } else {
       roVisCube.reset(&acc.visibility(), utility::NullDeleter());
       roVisNoise.reset(&acc.noise(), utility::NullDeleter());
   }

   const uint iDDOffset = itsSourceIndex * nSamples;

   // MV: always create UVWeight object even if traditional weighting is not done / it is not needed for this particular type of gridder.
   // This is the price paid to have a generic code. However, this object is lightweight (effectively only manages a pointer behind the scene +
   // has some basic metadata), so shouldn't be a huge overhead. It can be moved inside the samples loop (although it is not obvious whether
   // this is better).
   UVWeight uvWeight;

   // Use a separate UVWeight object for the optional gridder-based builder (RW - read/write). Doing it this way (as opposed to reusing uvWeight declared above) allows us
   // to both apply some preliminary weight and build accurate one during the same (first) major cycle (if we found this mode interesting). At this stage,
   // it looks like we're unlikely to use gridder-based weight builder long term. So if declaring a separate unused object here is found to be an unacceptable
   // overhead, we can work with the temporary in the inner loop (which won't be invoked anyway if the builder is not associated with the gridder)
   UVWeight uvWeightRW;

   for (uint i=0; i<nSamples; ++i) {
       if (itsMaxPointingSeparation > 0.) {
           // need to reject samples, if too far from the image centre
           const casacore::MVDirection thisPointing  = acc.pointingDir1()(i);
           if (imageCentre.separation(thisPointing) > itsMaxPointingSeparation) {
               ++itsRowsRejectedDueToMaxPointingSeparation;
               continue;
           }
       }

       if (itsFirstGriddedVis && isPSFGridder()) {
           // data members related to representative feed and field are used for
           // reverse problem only (from visibilities to image).
           if (itsUseAllDataForPSF) {
               ASKAPLOG_DEBUG_STR(logger, "All data are used to estimate PSF");
           } else {
               itsFeedUsedForPSF = acc.feed1()(i);
               itsPointingUsedForPSF = acc.dishPointing1()(i);
               ASKAPLOG_DEBUG_STR(logger, "Using the data for feed "<<itsFeedUsedForPSF<<
                  " and field at "<<printDirection(itsPointingUsedForPSF)<<" to estimate the PSF");
           }
           itsFirstGriddedVis = false;
       }
       if (itsPARotation) {
           // make sure we set the swap before calculating the rotation matrix,
           // as we're just swapping rows or columns in the matrix
           itsPolConv.setSwapPols(itsSwapPols);
           // set parallactic angle
           double pa1=acc.feed1PA()(i) + itsPARotAngle;
           double pa2=acc.feed2PA()(i) + itsPARotAngle;
           itsPolConv.setParAngle(pa1,pa2);
       }

       // we need to allow uv-weight application for both PSF and PCF gridders, perhaps such
       // gridders can be used to build weight as well?
       if (!forward && itsUVWeightAccessor) {
           uvWeight = itsUVWeightAccessor->getWeight(acc.feed1()(i), currentFieldIndex(), itsSourceIndex);
           ASKAPDEBUGASSERT(uvWeight.uSize() == shape()(0));
           ASKAPDEBUGASSERT(uvWeight.vSize() == shape()(1));
       }

       if (!forward && !isPSFGridder() && !isPCFGridder() && itsUVWeightBuilder) {
           uvWeightRW = itsUVWeightBuilder->addWeight(acc.feed1()(i), currentFieldIndex(), itsSourceIndex);
           ASKAPDEBUGASSERT(uvWeightRW.uSize() == shape()(0));
           ASKAPDEBUGASSERT(uvWeightRW.vSize() == shape()(1));
       }

       for (uint chan=0; chan<nChan; ++chan) {
           const double reciprocalToWavelength = frequencyList[chan]/casacore::C::c;
           if (chan == 0) {
              // check for ridiculous frequency to pick up a possible error with input file,
              // not essential for processing as such
              ASKAPCHECK((reciprocalToWavelength>0.1) && (reciprocalToWavelength<30000),
                  "Check frequencies in the input file as the order of magnitude is likely to be wrong, "
                  "comment this statement in the code if you're trying something non-standard. Frequency = "<<
                  frequencyList[chan]/1e9<<" GHz");
           }

           /// Scale U,V to integer pixels plus fractional terms
           const double uScaled=reciprocalToWavelength * outUVW(i)(0) / itsUVCellSize(0);
           int iu = askap::nint(uScaled);
           int fracu=askap::nint(itsOverSample*(double(iu)-uScaled));
           if (fracu<0) {
               iu+=1;
               fracu += itsOverSample;
           } else if (fracu>=itsOverSample) {
               iu-=1;
               fracu -= itsOverSample;
           }
           ASKAPCHECK(fracu>-1,
                   "Fractional offset in u is negative, uScaled="<<uScaled<<
                   " iu="<<iu<<" oversample="<<itsOverSample<<" fracu="<<fracu);
           ASKAPCHECK(fracu<itsOverSample,
                   "Fractional offset in u exceeds oversampling, uScaled="<<uScaled<<
                   " iu="<<iu<<" oversample="<<itsOverSample<<" fracu="<<fracu);
           iu+=itsShape(0)/2;

           const double vScaled=reciprocalToWavelength * outUVW(i)(1) / itsUVCellSize(1);
           int iv = askap::nint(vScaled);
           int fracv=askap::nint(itsOverSample*(double(iv)-vScaled));
           if (fracv<0) {
               iv+=1;
               fracv += itsOverSample;
           } else if (fracv>=itsOverSample) {
               iv-=1;
               fracv -= itsOverSample;
           }
           ASKAPCHECK(fracv>-1,
                   "Fractional offset in v is negative, vScaled="<<vScaled<<
                   " iv="<<iv<<" oversample="<<itsOverSample<<" fracv="<<fracv);
           ASKAPCHECK(fracv<itsOverSample,
                   "Fractional offset in v exceeds oversampling, vScaled="<<vScaled<<
                   " iv="<<iv<<" oversample="<<itsOverSample<<" fracv="<<fracv);
           iv+=itsShape(1)/2;

           // Calculate the delay phasor
           const double phase=casacore::C::_2pi * reciprocalToWavelength * delay(i);

           const casacore::Complex phasor(cos(phase), sin(phase));

           bool allPolGood=true;
           for (uint pol=0; pol<nPol; ++pol) {
               //if (flagCube(i, chan, pol)) {
               if (flagCube(pol, chan, i)) {
                   allPolGood=false;
                   break;
               }
           }

           /*
           // temporary for debugging
           if (allPolGood && !itsFreqMapper.isMapped(chan)) {
              ASKAPLOG_DEBUG_STR(logger, "Channel "<<chan<<" is not mapped to the image cube");
           }
           */
           // check if w values are in range
           bool wGood = true;
           // Ensure that we only use unflagged data, incomplete polarisation vectors are
           // ignored
           // @todo Be more careful about matching polarizations

           if (allPolGood && itsFreqMapper.isMapped(chan)) {
               // obtain which channel of the image this accessor channel is mapped to
               const int imageChan = itsFreqMapper(chan);
               ipStart(3) = imageChan;

               // check that imageChan is within the shape of uvWeight grid, also cater for the
               // situation when weighting is not done
               ASKAPDEBUGASSERT(uvWeight.empty() || imageChan < uvWeight.nPlane());
               ASKAPDEBUGASSERT(uvWeightRW.empty() || imageChan < uvWeightRW.nPlane());

               if (!forward) {
                   if (!isPSFGridder() && !isPCFGridder()) {
                       ASKAPDEBUGASSERT(roVisCube!=0);
                       for (uint pol=0; pol<nPol; pol++) itsPolVector(pol) = (*roVisCube)(pol,chan,iDDOffset+i);
                       itsPolConv.convert(itsImagePolFrameVis,itsPolVector);
                   }
                   // we just don't need this quantity for the forward gridder, although there would be no
                   // harm to always compute it
                   ASKAPDEBUGASSERT(roVisNoise!=0);
                   for (uint pol=0; pol<nPol; pol++) itsPolVector(pol) = (*roVisNoise)(pol,chan,i);
                   itsPolConv.noise(itsImagePolFrameNoise,itsPolVector);
               }
               // Now loop over all image polarizations
               for (uint pol=0; pol<nImagePols; ++pol) {
                   // Lookup the portion of grid to be
                   // used for this row, polarisation and channel
                   const int gInd=gIndex(i, pol, chan);
                   if (gInd<0) {
                       // can't use this data - w out of range
                       wGood = false;
                       itsVectorsWFlagged +=1;
                       break;
                   }
                   //ASKAPCHECK(gInd>-1,"Index into image grid is less than zero");
                   ASKAPCHECK(gInd<int(itsGrid.size()), "Index into image grid exceeds number of planes");

                   // Lookup the convolution function to be
                   // used for this row, polarisation and channel
                   // cIndex gives the index for this row, polarization and channel. On top of
                   // that, we need to adjust for the oversampling since each oversampled
                   // plane is kept as a separate matrix.
                   const int beforeOversamplePlaneIndex = cIndex(i,pol,chan);
                   if (beforeOversamplePlaneIndex<0) {
                       // can't use this data - w out of range
                       wGood = false;
                       itsVectorsWFlagged +=1;
                       break;
                   }
                   const int cInd=fracu+itsOverSample*(fracv+itsOverSample*beforeOversamplePlaneIndex);
                   ASKAPCHECK(cInd>-1,"Index into convolution functions is less than zero");
                   ASKAPCHECK(cInd<int(itsConvFunc.size()),
                           "Index into convolution functions exceeds number of planes");
                   casacore::Matrix<casacore::Complex> & convFunc(itsConvFunc[cInd]);

                   // support only square convolution functions at the moment
                   ASKAPDEBUGASSERT(convFunc.nrow() == convFunc.ncolumn());
                   ASKAPCHECK(convFunc.nrow() % 2 == 1,
                           "Expect convolution function with an odd number of pixels for each axis, CF["<<
                           cInd<<"] has shape="<<convFunc.shape());

                   // we now use support size for this given plane in the CF cache; itsSupport is a maximum
                   // support across all CFs (this allows plane-dependent support size)
                   const int support = (int(convFunc.nrow()) - 1) / 2;
                   ASKAPCHECK(support >= 0, "Support must be zero or greater, CF["<<cInd<<"] has shape="<<
                              convFunc.shape()<<" giving a support of "<<support);

                   // This seems to be the quickest way to get a reference to the matrix we want
                   // It assumes itsGrid is contiguous
                   ASKAPDEBUGASSERT(itsGrid[gInd].contiguousStorage());
                   ipStart(2) = pol;
                   // Check if we need to update the grid reference
                   if ( nImagePols>1 || imageChan!=itsImageChan || gInd!=itsGridIndex ) {
                       its2dGrid.takeStorage(onePlane,&itsGrid[gInd](ipStart),casacore::SHARE);
                       itsImageChan = imageChan;
                       itsGridIndex = gInd;
                   }
                   // the following accounts for a possible offset of the convolution function
                   const std::pair<int,int> cfOffset = getConvFuncOffset(beforeOversamplePlaneIndex);
                   const int iuOffset = iu + cfOffset.first;
                   const int ivOffset = iv + cfOffset.second;

                   /*
                   if (isPSFGridder()) {
                       if (i!=0) continue;
                       std::cout<<"sample "<<i<<" iu="<<iu<<" uScaled="<<uScaled<<" fracu="<<fracu<<" "<<
                              "iv="<<iv<<" vScaled="<<vScaled<<" fracv="<<fracv<<
                              " wplane="<<beforeOversamplePlaneIndex<<std::endl;
                   }
                   */

                   /// Need to check if this point lies on the grid (taking into
                   /// account the support)
                   if (((iuOffset-support)>0)&&((ivOffset-support)>0)&&
                       ((iuOffset+support) <itsShape(0))&&((ivOffset+support)<itsShape(1))) {
                       if (forward) {
                           casacore::Complex cVis(0.,0.);
                           GridKernel::degrid(cVis, convFunc, its2dGrid, iuOffset, ivOffset, support);
                           itsSamplesDegridded+=1.0;
                           itsNumberDegridded+=double((2*support+1)*(2*support+1));
                           if (itsVisWeight) {
                               cVis *= itsVisWeight->getWeight(frequencyList[chan]);
                           }
                           itsImagePolFrameVis[pol] = cVis*phasor;
                       } else {
                           const casacore::Complex visComplexNoise = itsImagePolFrameNoise[pol];

                           const float visNoise = casacore::square(casacore::real(visComplexNoise));
                           //const float visNoise = casacore::norm(visComplexNoise);
                           const float visNoiseWt = (visNoise > 0.) ? 1./visNoise : 0.;
                           ASKAPCHECK(visNoiseWt>0., "Weight is supposed to be a positive number; visNoiseWt="<<
                                      visNoiseWt<<" visNoise="<<visNoise<<" visComplexNoise="<<visComplexNoise);

                           // get the full weight, this one can be zero (although the total weight shouldn't be)
                           // For uvWeight, we use iu, iv instead of iuOffset,ivOffset. The latter are just to take care of
                           // zeros in the CF. It may be worth to experiment, however, whether it would make any difference
                           // (or think about the physical meaning). Also note that we effectively assume that weighting is
                           // the same for all polarisations. It is a can of worms if this breaks down. But 1) we assume it
                           // in many places, 2) samples with incomplete polarisations are ignored. So probably ok.
                           const float visWt = uvWeight.empty() ? visNoiseWt : visNoiseWt * uvWeight(iu, iv, imageChan);

                           // row in itsSumWeights to work with
                           const int sumWeightsRow =
                               itsTrackWeightPerOversamplePlane ? cInd : beforeOversamplePlaneIndex;

                           ASKAPCHECK(itsSumWeights.nelements()>0, "Sum of weights not yet initialised");
                           ASKAPDEBUGASSERT(itsSumWeights.shape().nelements() >= 3);
                           ASKAPCHECK(sumWeightsRow < int(itsSumWeights.shape()(0)),
                                      "Index into itsSumWeights of " << sumWeightsRow <<
                                      " is greater than allowed " << int(itsSumWeights.shape()(0)));
                           ASKAPDEBUGASSERT(pol < uint(itsSumWeights.shape()(1)));
                           ASKAPDEBUGASSERT(imageChan < int(itsSumWeights.shape()(2)));

                           if (!isPSFGridder() && !isPCFGridder()) {
                               /// Gridding visibility data onto grid
                               casacore::Complex rVis = phasor*conj(itsImagePolFrameVis[pol])*visWt;
                               if (itsVisWeight) {
                                   rVis *= itsVisWeight->getWeight(frequencyList[chan]);
                               }
                               GridKernel::grid(its2dGrid, convFunc, rVis, iuOffset, ivOffset, support);

                               itsSamplesGridded+=1.0;
                               itsNumberGridded+=double((2*support+1)*(2*support+1));

                               itsSumWeights(sumWeightsRow, pol, imageChan) += visWt; //1.0;

                               if (!uvWeightRW.empty()) {
                                   uvWeightRW(iu,iv,imageChan) += visNoiseWt;
                               }
                           }
                           /// Grid the PSF?
                           if (isPSFGridder() &&
                               (itsUseAllDataForPSF ||
                                ((itsFeedUsedForPSF == acc.feed1()(i)) &&
                                 (itsPointingUsedForPSF.separation(acc.dishPointing1()(i))<1e-6)))) {
                                casacore::Complex uVis(1.,0.);
                                uVis *= visWt;
                                if (itsVisWeight) {
                                    uVis *= itsVisWeight->getWeight(frequencyList[chan]);
                                }

                                GridKernel::grid(its2dGrid, convFunc, uVis, iuOffset, ivOffset, support);
                                itsSamplesGridded+=1.0;
                                itsNumberGridded+=double((2*support+1)*(2*support+1));

                                itsSumWeights(sumWeightsRow, pol, imageChan) += visWt; //1.0;
                           } // end if psf needs to be done
                           /// Grid the preconditioner function?
                           if (isPCFGridder()) {
                                casacore::Complex uVis(1.,0.);
                                uVis *= visWt;
                                // We don't want different preconditioning for different Taylor terms.
                                //if (itsVisWeight) {
                                //    uVis *= itsVisWeight->getWeight(frequencyList[chan]);
                                //}
                                // MV: probably it would be better not to setup the gridder with visWeight rather than duplicate code for PCFGridder - need to understand this part more

                                // storing w information in the imaginary part of the PCF,
                                // so make them add with conjugate symmetry.
                                if ((ivOffset<itsShape(1)/2 && iuOffset>=itsShape(0)/2) ||
                                    (ivOffset<=itsShape(1)/2 && iuOffset<itsShape(0)/2)) {
                                //if (isPCFGridder() && ivOffset<itsShape(1)/2) {
                                  casacore::Matrix<casacore::Complex> conjFunc = conj(convFunc);
                                  GridKernel::grid(its2dGrid, conjFunc, uVis, iuOffset, ivOffset, support);
                                } else {
                                  GridKernel::grid(its2dGrid, convFunc, uVis, iuOffset, ivOffset, support);
                                }

                                itsSamplesGridded+=1.0;
                                itsNumberGridded+=double((2*support+1)*(2*support+1));

                                // these aren't used. Can probably also disable the PSF weights
                                //itsSumWeights(sumWeightsRow, pol, imageChan) += visWt; //1.0;
                           } // end if pcf needs to be done

                       } // end if forward (else case, reverse operation)
                   } // end of on-grid if statement
               } // end of pol loop
               // need to write back the result for degridding
               if (wGood) {
                   if (forward) {
                       ASKAPDEBUGASSERT(visCube!=0)
                       itsPolConv.convert(itsPolVector,itsImagePolFrameVis);
                       for (uint pol=0; pol<nPol; pol++) (*visCube)(pol,chan,iDDOffset+i) += itsPolVector(pol);
                       // visibilities with w out of range are left unchanged during prediction
                       // as long as subsequent imaging uses the same wmax this should work ok
                       // we may want to flag these data to be sure
                   }
               } else {
                   if (!forward) {
                       itsVectorsWFlagged +=1;
                   }
               }
           } else { // if (allPolGood)
               if (!forward) {
                   itsVectorsFlagged+=1;
               }
          }
       } //end of chan loop
   } //end of i loop

   if (forward) {
       itsTimeDegridded+=timer.real();
   } else {
       itsTimeGridded+=timer.real();
   }
}

// DAM TRADITIONAL
void TableVisGridder::setWeights(accessors::IDataAccessor& acc) {
    ASKAPTRACE("TableVisGridder::setweights");

   initIndices(acc);
   initConvolutionFunction(acc);

   // Now time the coordinate conversions, etc.
   // some conversion may have already happened during CF calculation
   const casa::MVDirection imageCentre = getImageCentre();
   const casa::MVDirection tangentPoint = getTangentPoint();

   const casa::Vector<casa::RigidVector<double, 3> > &outUVW = acc.rotatedUVW(tangentPoint);

   const uint nSamples = acc.nRow();
   const uint nChan = acc.nChannel();
   const uint nPol = acc.nPol();
   const casa::Vector<casa::Double>& frequencyList = acc.frequency();
   itsFreqMapper.setupMapping(frequencyList);

   ASKAPDEBUGASSERT(itsShape.nelements()>=2);
   const casa::IPosition onePlane4D(4, itsShape(0), itsShape(1), 1, 1);
   const casa::IPosition onePlane(2, itsShape(0), itsShape(1));

   for (uint i=0; i<nSamples; ++i) {
       for (uint chan=0; chan<nChan; ++chan) {
           //const double reciprocalToWavelength = frequencyList[chan]/casa::C::c;
           /// Scale U,V to integer pixels plus fractional terms
           const double uScaled=frequencyList[chan]*outUVW(i)(0)/(casa::C::c *itsUVCellSize(0));
           int iu = askap::nint(uScaled);
           int fracu=askap::nint(itsOverSample*(double(iu)-uScaled));
           if (fracu<0) {
               iu+=1;
               fracu += itsOverSample;
           } else if (fracu>=itsOverSample) {
               iu-=1;
               fracu -= itsOverSample;
           }
           iu+=itsShape(0)/2;
           const double vScaled=frequencyList[chan]*outUVW(i)(1)/(casa::C::c *itsUVCellSize(1));
           int iv = askap::nint(vScaled);
           int fracv=askap::nint(itsOverSample*(double(iv)-vScaled));
           if (fracv<0) {
               iv+=1;
               fracv += itsOverSample;
           } else if (fracv>=itsOverSample) {
               iv-=1;
               fracv -= itsOverSample;
           }
           iv+=itsShape(1)/2;

           bool allPolGood=true;
           for (uint pol=0; pol<nPol; ++pol) {
               //if (acc.flag()(i, chan, pol))
               if (acc.flag()(pol, chan, i))
                   allPolGood=false;
           }

           if (allPolGood && itsFreqMapper.isMapped(chan)) {
               // obtain which channel of the image this accessor channel is mapped to
               const int imageChan = itsFreqMapper(chan);

               // number of polarisation planes in the grid
               const casa::uInt nImagePols = (shape().nelements()<=2) ? 1 : shape()[2];

               // Now loop over all image polarizations
               //for (uint pol=0; pol<nImagePols; ++pol) {
               {
                   // DAM the weighting below is done for all pols at once. Might be a prob if pols have diff flags.
                   uint pol=0;
                   // Lookup the portion of grid to be
                   // used for this row, polarisation and channel
                   const int gInd=gIndex(i, pol, chan);

                   /// Make a slicer to extract just this plane
                   const casa::IPosition ipStart(4, 0, 0, pol, imageChan);
                   const casa::Slicer slicer(ipStart, onePlane4D);

                   const int beforeOversamplePlaneIndex = cIndex(i,pol,chan);
                   //const int cInd=fracu+itsOverSample*(fracv+itsOverSample*beforeOversamplePlaneIndex);

                   casa::Array<casa::Complex> aGrid(itsGrid[gInd](slicer));
                   casa::Matrix<casa::Complex> grid(aGrid.nonDegenerate());

                   // the following accounts for a possible offset of the convolution function
                   const std::pair<int,int> cfOffset = getConvFuncOffset(beforeOversamplePlaneIndex);
                   const int iuOffset = iu + cfOffset.first;
                   const int ivOffset = iv + cfOffset.second;

                   ASKAPCHECK(iuOffset >= 0 && ivOffset >= 0 && iuOffset < grid.nrow() && ivOffset < grid.ncolumn(),
                    "grid coordinates out of range in setWeight");

                   if ( real(grid(iuOffset, ivOffset)) > 0.0 ) {
                       //casa::Vector<casa::Complex> thisChanNoise = acc.noise().yzPlane(i).row(chan);
                       casa::Vector<casa::Complex> thisChanNoise = acc.noise().xyPlane(i).column(chan);
                       const float rootInvWgt = sqrt(real(grid(iuOffset, ivOffset)));
                       thisChanNoise *= rootInvWgt;
                   }
               }

           }

       }
   }

}

// DAM TRADITIONAL
void TableVisGridder::addConjugates() {
    ASKAPTRACE("TableVisGridder::addConjugates");

    uint i=0;
    uint chan=0;
    uint pol=0;
    const int gInd=gIndex(i, pol, chan);
    const int imageChan = itsFreqMapper(chan);
    const casa::IPosition onePlane4D(4, itsShape(0), itsShape(1), 1, 1);

    /// Make a slicer to extract just this plane
    const casa::IPosition ipStart(4, 0, 0, pol, imageChan);
    const casa::Slicer slicer(ipStart, onePlane4D);

    casa::Matrix<casa::Complex> aGrid(itsGrid[gInd](slicer).nonDegenerate());

// DAM -- it would be faster to do this all in uv, but testing and don't want to deal with off-by-one issues...

ASKAPLOG_INFO_STR(logger, "DAMDAM before conjugates. sum of grid = " << sum(real(aGrid)) );
// Limit number of fft threads to 8 (more is slower for our fft sizes)
    scimath::FFT2DWrapper<casacore::Complex> fft2d(true,8);
    fft2d(aGrid, false);
    aGrid += conj(aGrid);
    fft2d(aGrid, true);
ASKAPLOG_INFO_STR(logger, "DAMDAM after conjugates. sum of grid = " << sum(real(aGrid)) );

}

// DAM TRADITIONAL
void TableVisGridder::setRobustness(const float robustness) {
    ASKAPTRACE("TableVisGridder::setRobustness");

    uint i=0;
    uint chan=0;
    uint pol=0;
    const int gInd=gIndex(i, pol, chan);
    const int imageChan = itsFreqMapper(chan);
    const casa::IPosition onePlane4D(4, itsShape(0), itsShape(1), 1, 1);

    /// Make a slicer to extract just this plane
    const casa::IPosition ipStart(4, 0, 0, pol, imageChan);
    const casa::Slicer slicer(ipStart, onePlane4D);

    casa::Array<casa::Complex> aGrid(itsGrid[gInd](slicer));
    casa::Matrix<casa::Complex> grid(aGrid.nonDegenerate());
    ASKAPLOG_INFO_STR(logger, "DAMDAM before robustness. sum of grid = " << sum(real(aGrid)) );
    ASKAPLOG_INFO_STR(logger, "DAMDAM before robustness. robustness = " << robustness );

    ASKAPLOG_DEBUG_STR(logger, "DAM estimating the average wgt sum");
    casa::Array<double> wgts(aGrid.nonDegenerate().shape());
    casa::convertArray<double,float>(wgts, real(aGrid.nonDegenerate()));
    double aveWgtSum = sum(wgts*wgts) / sum(wgts);
    ASKAPLOG_DEBUG_STR(logger, "DAM average wgt sum estimate: " << aveWgtSum);

    const float noisePower = (1.0/aveWgtSum)*25.0*std::pow(10., -2.0*robustness);
    ASKAPLOG_DEBUG_STR(logger, "DAM noisePower estimate: " << noisePower);

    for (uint iv=0; iv<itsShape(1); ++iv) {
        for (uint iu=0; iu<itsShape(0); ++iu) {
            //grid(iu, iv) = 1.0 / (noisePower*real(grid(iu, iv)) + 1.0);
            grid(iu, iv) = (noisePower*real(grid(iu, iv)) + 1.0);
        }
    }

}

/// @brief correct visibilities, if necessary
/// @details This method is intended for on-the-fly correction of visibilities (i.e.
/// facet-based correction needed for LOFAR). This method does nothing in this class, but
/// can be overridden in the derived classes to plug some effect in. The same method is
/// used for both gridding and degridding, with the forward parameter used to distinguish
/// between these two operations. A non-const accessor has to be modified in situ, if a
/// correction is required. A buffer for read-only visibilities is created on-demand when
/// rwVisibility method of the accessor is called for the first time.
void TableVisGridder::correctVisibilities(accessors::IDataAccessor &, bool) {}

/// @brief Degrid the visibility data.
/// @param[in] acc non-const data accessor to work with
void TableVisGridder::degrid(accessors::IDataAccessor& acc) {
    ASKAPTRACE("TableVisGridder::degrid");
    generic(acc, true);
    correctVisibilities(acc, true);
}

/// @brief Grid the visibility data.
/// @param acc const data accessor to work with
/// @note a non-const adapter is created behind the scene. If no on-the-fly visibility
/// correction is performed, this adapter is equivalent to the original const data accessor
void TableVisGridder::grid(accessors::IConstDataAccessor& acc) {
    ASKAPTRACE("TableVisGridder::grid");
    accessors::OnDemandBufferDataAccessor bufAcc(acc);
    correctVisibilities(bufAcc, false);
    generic(bufAcc, false);
}

/// @brief obtain the centre of the image
/// @details This method extracts RA and DEC axes from itsAxes and
/// forms a direction measure corresponding to the middle of each axis.
/// @return direction measure corresponding to the image centre
casacore::MVDirection TableVisGridder::getImageCentre() const
{
   ASKAPCHECK(itsAxes.hasDirection(),"Direction axis is missing. axes="<<itsAxes);
   casacore::MDirection out;
   casacore::Vector<casacore::Double> centrePixel(2);
   ASKAPDEBUGASSERT(itsShape.nelements()>=2);
   ASKAPDEBUGASSERT(paddingFactor()>0);

   for (size_t dim=0; dim<2; ++dim) {
        centrePixel[dim] = double(itsShape[dim])/2./double(paddingFactor());
   }
   ASKAPCHECK(syncHelper.toWorld(itsAxes.directionAxis(),out, centrePixel),
        "Unable to obtain world coordinates for the centre of the image. Something is wrong with the coordinate system");
   return out.getValue();
}

/// @brief obtain the tangent point
/// @details For faceting all images should be constructed for the same tangent
/// point. This method extracts the tangent point (reference position) from the
/// coordinate system.
/// @return direction measure corresponding to the tangent point
casacore::MVDirection TableVisGridder::getTangentPoint() const
{
   ASKAPCHECK(itsAxes.hasDirection(),"Direction axis is missing. axes="<<itsAxes);
   const casacore::Vector<casacore::Double> refVal(itsAxes.directionAxis().referenceValue());
   ASKAPDEBUGASSERT(refVal.nelements() == 2);
   const casacore::Quantum<double> refLon(refVal[0], "rad");
   const casacore::Quantum<double> refLat(refVal[1], "rad");
   const casacore::MVDirection out(refLon, refLat);
   return out;
}

/// @brief Conversion helper function
/// @details Copies in to out expanding double into complex values and
/// padding appropriately if necessary (itsPaddingFactor is more than 1)
/// @param[out] out complex output array
/// @param[in] in double input array
/// @param[in] padding padding factor
void TableVisGridder::toComplex(casacore::Array<imtypeComplex>& out,
        const casacore::Array<imtype>& in, const float padding) {
    ASKAPDEBUGTRACE("TableVisGridder::toImtypeComplex");

    out.resize(scimath::PaddingUtils::paddedShape(in.shape(),padding));
    out.set(0.);
    casacore::Array<imtypeComplex> subImage = scimath::PaddingUtils::extract(out,padding);
    casacore::convertArray<imtypeComplex, imtype>(subImage, in);
}

/// @brief Conversion helper function
/// @details Copies real part of in into imtype array and
/// extracting an inner rectangle if necessary (itsPaddingFactor is more than 1)
/// @param[out] out real output array
/// @param[in] in complex input array
/// @param[in] padding padding factor
void TableVisGridder::toDouble(casacore::Array<imtype>& out,
        const casacore::Array<imtypeComplex>& in, const float padding) {
  ASKAPDEBUGTRACE("TableVisGridder::toDouble");
  casacore::Array<imtypeComplex> wrapper(in);
  const casacore::Array<imtypeComplex> subImage = scimath::PaddingUtils::extract(wrapper,padding);
  out.resize(subImage.shape());
  out = real(subImage);
}

/// @brief set up itsStokes using the information from itsAxes and itsShape
void TableVisGridder::initStokes()
{
   const int nPol = itsShape.nelements()>=3 ? itsShape[2] : 1;
   if (itsAxes.has("STOKES")) {
       itsStokes = itsAxes.stokesAxis();
   } else {
       itsStokes.resize(1);
       itsStokes[0] = casacore::Stokes::I;
   }
   ASKAPCHECK(int(itsStokes.nelements()) == nPol, "Stokes axis is not consistent with the shape of the grid. There are "<<
              nPol<<" planes in the grid and "<<itsStokes.nelements()<<
              " polarisation descriptors defined by the STOKES axis");
}


void TableVisGridder::initialiseGrid(const scimath::Axes& axes,
        const casacore::IPosition& shape, const bool dopsf, const bool dopcf) {
     ASKAPTRACE("TableVisGridder::initialiseGrid");

     ASKAPDEBUGASSERT(shape.nelements()>=2);
     itsShape = scimath::PaddingUtils::paddedShape(shape,paddingFactor());

     initialiseCellSize(axes);

     initStokes();

     configureForPSF(dopsf);
     configureForPCF(dopcf);

     initialiseWeightBuilder();

     /// We only need one grid
     itsGrid.resize(1);
     itsGrid[0].resize(itsShape);
     itsGrid[0].set(0.0);

     // This ensures we reinitalise the pol converter
     itsVisPols.resize(0);

     if (isPSFGridder()) {
         // for a proper PSF calculation
         initRepresentativeFieldAndFeed();
     }

     initialiseSumOfWeights();
     ASKAPCHECK(itsSumWeights.nelements()>0, "Sum of weights not yet initialised");
     initialiseFreqMapping();
     ASKAPLOG_DEBUG_STR(logger, "Gridding is set up with tangent centre "<<
             printDirection(getTangentPoint())<<" and image centre "<<
             printDirection(getImageCentre()));
}

/// @brief helper method to set up cell size
/// @details Similar action is required to calculate uv-cell size for gridding and degridding.
/// Moreover, derived gridders may override initialiseGrid and initialiseDegrid and we don't want
/// to duplicate the code up there. This method calculates uv-cell size for both ra and dec axes
/// using coordinate information provided. This method also assigns passed axes parameter to itsAxes.
/// @param[in] axes coordinate system (ra and dec axes are used).
void TableVisGridder::initialiseCellSize(const scimath::Axes& axes)
{
    itsAxes=axes;
    ASKAPCHECK(itsAxes.hasDirection(), "Direction axis is missing. itsAxes:"<<itsAxes);
    casacore::Vector<casacore::Double> increments = itsAxes.directionAxis().increment();
    ASKAPCHECK(increments.nelements() == 2, "Expect 2 elements in the increment vector, you have "<<increments);
    itsUVCellSize.resize(2);
    ASKAPDEBUGASSERT(itsShape.nelements()>=2);
    for (size_t dim = 0; dim<2; ++dim) {
         itsUVCellSize[dim] = 1./(increments[dim]*double(itsShape[dim]));
    }
}

/// @brief helper method to set up weight builder
/// @details Similar action is required for a number of gridders to setup weight builder with the
/// grid parameters. It is only expected to be called from initialiseGrid, after the shape is set
/// (and the weight builder assigned).
/// We could've just added this code to initialiseCellSize, but it would be called unnecessary from
/// initialiseDegrid and for PCF/PSF gridders (although, presumably, the builder won't be set in this
/// cases, so no harm). Having a separate method is neater.
void TableVisGridder::initialiseWeightBuilder()
{
  // why to exclude PSF and PCF gridder here? It can be controlled on the user side
  if (itsUVWeightBuilder && !isPSFGridder() && !isPCFGridder()) {
      ASKAPCHECK(itsShape.nelements()>=2, "Shape has not been initialised before initialiseWeightBuilder is called");
      itsUVWeightBuilder->initialise(itsShape[0], itsShape[1], itsShape.nelements() > 3 ? itsShape[3] : 1u);
  }
}


/// @brief initialise sum of weights
/// @details We keep track the number of times each convolution function is used per
/// channel and polarisation (sum of weights). This method is made virtual to be able
/// to do gridder specific initialisation without overriding initialiseGrid.
/// This method accepts no parameters as itsShape, itsNWPlanes, etc should have already
/// been initialised by the time this method is called.
void TableVisGridder::initialiseSumOfWeights()
{
  resizeSumOfWeights(1);
  zeroSumOfWeights();
}

/// @brief resize sum of weights
/// @details This method is used inside initialiseSumOfWeights and its overrides in
/// derived classes. It resizes itsSumWeights to a given number of convolution
/// functions taking into account channels/polarisations according to itsShape.
/// Moving this operation into the separate method allows better data encapsulation
/// and tracking weights per oversampling plane or per convolution function depending
/// on the user's choice.
/// @param[in] numcf number of convolution functions in the cache (before oversampling)
void TableVisGridder::resizeSumOfWeights(const int numcf)
{
  ASKAPDEBUGASSERT(numcf>0);
  const int numRows = itsTrackWeightPerOversamplePlane ? numcf*itsOverSample*itsOverSample : numcf;
  ASKAPDEBUGASSERT(numRows>0);
  itsSumWeights.resize(numRows,itsShape.nelements()>=3 ? itsShape(2) : 1,
                         itsShape.nelements()>=4 ? itsShape(3) : 1);
}

/// @brief a helper method to initialize gridding of the PSF
/// @details The PSF is calculated using the data for a
/// representative field/feed only. By default, the first encountered
/// feed/field is chosen. If the same gridder is reused for another
/// sequence of data points a new representative feed/field have to be
/// found. This is done by resetting the cache in initialiseGrid. However,
/// the latter method can be overridden in the derived classes. To avoid
/// a duplication of the code, this helper method resets the representative
/// feed/field cache. It is called from initialiseGrid.
void TableVisGridder::initRepresentativeFieldAndFeed()
{
  itsFirstGriddedVis = true;

  /*
  // temporary code for debuggig
  std::cout<<"TableVisGridder::initRepresentativeFieldAndFeed"<<std::endl;
  itsFirstGriddedVis = false;
  casacore::Quantity ra(0.,"rad"), dec(0.,"rad");
  casacore::MVAngle::read(ra,"13:32:22.48");
  casacore::MVAngle::read(dec,"-042.16.56.93");
  itsPointingUsedForPSF = casacore::MVDirection(ra,dec);
  ASKAPLOG_DEBUG_STR(logger, "Field override for PSF, will use "<<printDirection(itsPointingUsedForPSF));
  itsFeedUsedForPSF = 0;
  // end of temporary code
  */
}

/// This is the default implementation
void TableVisGridder::finaliseGrid(casacore::Array<imtype>& out) {
    ASKAPTRACE("TableVisGridder::finaliseGrid");
    ASKAPDEBUGASSERT(itsGrid.size() > 0);
    // buffer for result
    casacore::Array<imtype> dBuffer(itsGrid[0].shape());
    ASKAPDEBUGASSERT(dBuffer.shape().nelements()>=2);
    ASKAPDEBUGASSERT(itsShape == scimath::PaddingUtils::paddedShape(out.shape(),paddingFactor()));
    // Limit number of fft threads to 8 (more is slower for our fft sizes)
    scimath::FFT2DWrapper<imtypeComplex> fft2d(true,8);
    /// Loop over all grids Fourier transforming and accumulating
    for (unsigned int i=0; i<itsGrid.size(); i++) {
        #ifdef ASKAP_FLOAT_IMAGE_PARAMS
        casacore::Array<imtypeComplex> scratch;
        if (itsClearGrid) {
            // we'll throw the grid away at the end, so we can use it as scratch space
            scratch.reference(itsGrid[i]);
        } else {
            // we want the pre-fft grid, make a copy
            scratch = itsGrid[i];
        }
        #else
        casacore::Array<imtypeComplex> scratch(itsGrid[i].shape());
        casacore::convertArray<imtypeComplex,casacore::Complex>(scratch, itsGrid[i]);
        #endif

        /*

        // for debugging
        if (isPSFGridder()) {
            #ifdef ASKAP_FLOAT_IMAGE_PARAMS
            scimath::saveAsCasaImage("uvcoverage.imag",imag(scratch));
            casacore::Array<float> buf(real(scratch));
            scimath::saveAsCasaImage("uvcoverage.real",buf);
            // adjust values to extract part which gives a real symmetric FT and the remainder
            casacore::Matrix<float> bufM(buf.nonDegenerate());
            for (int x=0; x<int(bufM.nrow()); ++x) {
                 for (int y=0; y<int(bufM.ncolumn())/2; ++y) {
                      const float val = 0.5*(bufM(x,y)+bufM(bufM.nrow() - x -1, bufM.ncolumn() - y -1));
                      bufM(x,y) = val;
                      bufM(bufM.nrow() - x -1, bufM.ncolumn() - y -1) = val;
                 }
            }
            scimath::saveAsCasaImage("uvcoverage.sympart",buf);
            casacore::Matrix<imtypeComplex> scratchM(scratch.nonDegenerate());
            for (int x=0; x<int(scratchM.nrow()); ++x) {
                 for (int y=0; y<int(scratchM.ncolumn()); ++y) {
                      scratchM(x,y) -= imtype(bufM(x,y));
                  }
            }
            // as we ignore imaginary part after FT, make scratch hermitian to be fair
            for (int x=0; x<int(scratchM.nrow()); ++x) {
                 for (int y=0; y<int(scratchM.ncolumn())/2; ++y) {
                      const imtypeComplex val = 0.5f*(scratchM(x,y)+
                            conj(scratchM(scratchM.nrow() - x -1, scratchM.ncolumn() - y -1)));
                      scratchM(x,y) = val;
                      scratchM(scratchM.nrow() - x -1, scratchM.ncolumn() - y -1) = conj(val);
                 }
            }
            scimath::saveAsCasaImage("uvcoverage.asympart.imag",imag(scratch));
            scimath::saveAsCasaImage("uvcoverage.asympart.real",real(scratch));
            fft2d(scratch, false);
            scimath::saveAsCasaImage("psf.asympart.real",real(scratch));

            #else
            casacore::Array<float> buf(scratch.shape());
            casacore::convertArray<float,double>(buf,imag(scratch));
            scimath::saveAsCasaImage("uvcoverage.imag",buf);
            casacore::convertArray<float,double>(buf,real(scratch));
            scimath::saveAsCasaImage("uvcoverage.real",buf);
            // adjust values to extract part which gives a real symmetric FT and the remainder
            casacore::Matrix<float> bufM(buf.nonDegenerate());
            for (int x=0; x<int(bufM.nrow()); ++x) {
                 for (int y=0; y<int(bufM.ncolumn())/2; ++y) {
                      const float val = 0.5*(bufM(x,y)+bufM(bufM.nrow() - x -1, bufM.ncolumn() - y -1));
                      bufM(x,y) = val;
                      bufM(bufM.nrow() - x -1, bufM.ncolumn() - y -1) = val;
                 }
            }
            scimath::saveAsCasaImage("uvcoverage.sympart",buf);
            casacore::Matrix<imtypeComplex> scratchM(scratch.nonDegenerate());
            for (int x=0; x<int(scratchM.nrow()); ++x) {
                 for (int y=0; y<int(scratchM.ncolumn()); ++y) {
                      scratchM(x,y) -= double(bufM(x,y));
                  }
            }
            // as we ignore imaginary part after FT, make scratch hermitian to be fair
            for (int x=0; x<int(scratchM.nrow()); ++x) {
                 for (int y=0; y<int(scratchM.ncolumn())/2; ++y) {
                      const imtypeComplex val = 0.5*(scratchM(x,y)+
                            conj(scratchM(scratchM.nrow() - x -1, scratchM.ncolumn() - y -1)));
                      scratchM(x,y) = val;
                      scratchM(scratchM.nrow() - x -1, scratchM.ncolumn() - y -1) = conj(val);
                 }
            }
            casacore::convertArray<float,double>(buf,imag(scratch));
            scimath::saveAsCasaImage("uvcoverage.asympart.imag",buf);
            casacore::convertArray<float,double>(buf,real(scratch));
            scimath::saveAsCasaImage("uvcoverage.asympart.real",buf);
            fft2d(scratch, false);
            casacore::convertArray<float,double>(buf,real(scratch));
            scimath::saveAsCasaImage("psf.asympart.real",buf);
            #endif
            ASKAPCHECK(false, "Debug termination");
        }
        // end of debugging code
        */
        fft2d.transformAllHyperPlanes(scratch, false);

        if (i==0) {
            toDouble(dBuffer, scratch);
        } else {
            casacore::Array<imtype> work(dBuffer.shape());
            toDouble(work, scratch);
            dBuffer+=work;
        }
    }
    clearGrid();
    // Now we can do the convolution correction
    correctConvolution(dBuffer);
    dBuffer*=imtype(double(dBuffer.shape()(0))*double(dBuffer.shape()(1)));
    out = scimath::PaddingUtils::extract(dBuffer,paddingFactor());
}

/// @brief store given grid
/// @details This is a helper method for debugging, it stores the amplitude of a given
/// grid into a CASA image (prior to FFT done as part of finaliseGrid)
/// @param[in] name image name
/// @param[in] numGrid number of the grid to store
void TableVisGridder::storeGrid(const std::string &name, casacore::uInt numGrid) const
{
   ASKAPCHECK(numGrid < itsGrid.size(), "Requested grid number of "<<numGrid<<
              " exceeds or equal to "<<itsGrid.size()<<", the number of grids available");
   casacore::Array<float> buffer(itsGrid[numGrid].shape());
   buffer = amplitude(itsGrid[numGrid]);
   scimath::saveAsCasaImage(name, buffer);
}



/// This is the default implementation
void TableVisGridder::finaliseWeights(casacore::Array<imtype>& out) {
   ASKAPTRACE("TableVisGridder::finaliseWeights");
   ASKAPDEBUGASSERT(itsShape.nelements() >= 4);
    ASKAPDEBUGASSERT(itsShape == scimath::PaddingUtils::paddedShape(out.shape(),paddingFactor()));

    int nPol=itsShape(2);
    int nChan=itsShape(3);
    ASKAPDEBUGASSERT(out.shape().nelements() == 4);
    casa::IPosition ipStart(4, 0, 0, 0, 0);
    casa::IPosition ipEnd(4, out.shape()(0)-1, out.shape()(1)-1, 0, 0);

    ASKAPCHECK(itsSumWeights.nelements()>0, "Sum of weights not yet initialised");
    int nZ=itsSumWeights.shape()(0);

    for (int chan=0; chan<nChan; chan++) {
        ipStart(3) = chan;
        ipEnd(3) = chan;
        for (int pol=0; pol<nPol; pol++) {
            double sumwt=0.0;
            for (int iz=0; iz<nZ; iz++) {
              //              float sumConvFunc=real(casacore::sum(casacore::abs(itsConvFunc[iz])));
              //              ASKAPLOG_DEBUG_STR(logger, "Sum of conv func " << sumConvFunc);
                sumwt+=itsSumWeights(iz, pol, chan);
            }
            ipStart(2) = pol;
            ipEnd(2) = pol;
            out(ipStart,ipEnd).set(sumwt);
        }
    }
}

void TableVisGridder::initialiseDegrid(const scimath::Axes& axes,
        const casacore::Array<imtype>& in) {
    ASKAPTRACE("TableVisGridder::initialiseDegrid");
    configureForPSF(false);
    configureForPCF(false);
    itsShape = scimath::PaddingUtils::paddedShape(in.shape(),paddingFactor());

    initialiseCellSize(axes);
    initStokes();

    initialiseFreqMapping();

    /// We only need one grid
    itsGrid.resize(1);
    itsGrid[0].resize(itsShape);

    // Make sure we reinitalise the pol converter
    itsVisPols.resize(0);
    // Limit number of fft threads to 8 (more is slower for our fft sizes)
    scimath::FFT2DWrapper<imtypeComplex> fft2d(true,8);
    if (casacore::max(casacore::abs(in))>0.0) {
        itsModelIsEmpty=false;
        casacore::Array<imtype> scratch(itsShape,static_cast<imtype>(0.));
        scimath::PaddingUtils::extract(scratch, paddingFactor()) = in;
        correctConvolution(scratch);
        #ifdef ASKAP_FLOAT_IMAGE_PARAMS
        toComplex(itsGrid[0], scratch);
        fft2d.transformAllHyperPlanes(itsGrid[0], true);
        #else
        casacore::Array<imtypeComplex> scratch2(itsGrid[0].shape());
        toComplex(scratch2, scratch);
        fft2d.transformAllHyperPlanes(scratch2), true);
        casacore::convertArray<casacore::Complex,imtypeComplex>(itsGrid[0],scratch2);
        #endif
    } else {
        ASKAPLOG_DEBUG_STR(logger, "No need to degrid: model is empty");
        itsModelIsEmpty=true;
        itsGrid[0].set(casacore::Complex(0.0));
    }
}

/// @brief helper method to initialise frequency mapping
/// @details Derived gridders may override initialiseGrid and initialiseDegrid. Howerver,
/// they still need to be able to initialise frequency axis mapping (between accessor channels
/// and image cube), which is handled by a private member class. This method initialises the
/// mapper using the content of itsShape and itsAxes, which should be set prior to calling this
/// method.
void TableVisGridder::initialiseFreqMapping()
{
  if (itsAxes.has("FREQUENCY") && itsShape.nelements()>=4) {
      itsFreqMapper.setupImage(itsAxes, itsShape(3));
  } else {
      ASKAPLOG_DEBUG_STR(logger, "Forced to use single spectral plane gridding (either "
                                "FREQUENCY axis or the number of channels are missing");
      itsFreqMapper.setupSinglePlaneGridding();
  }
}

/// @brief log unused spectral planes
/// @details It is handy to write down the channel numbers into log if the sum of weights
/// is zero, i.e. if we have no data for this particular channel. This method does it
/// (it is called from the destructor, unless the gridder object didn't do any gridding).
void TableVisGridder::logUnusedSpectralPlanes() const
{
   if (itsSumWeights.nelements()>0) {
       std::string planesList;
       for (casacore::uInt plane = 0; plane<itsSumWeights.nplane(); ++plane) {
            const double sumWt = casacore::sum(itsSumWeights.xyPlane(plane));
            if (sumWt<=0.) {
                if (planesList.size()!=0) {
                    planesList+=",";
                }
                planesList += utility::toString<casacore::uInt>(plane);
            }
       }
       if (planesList.size() == 0) {
           planesList = "none";
       }
       ASKAPLOG_DEBUG_STR(logger, "Unused spectral planes: "<<planesList);
   }
}


/// This is the default implementation
void TableVisGridder::finaliseDegrid() {
    /// Nothing to do
    clearGrid();
}

void TableVisGridder::clearGrid() {
    // Free up the grid memory?
    if (itsClearGrid) {
        ASKAPLOG_INFO_STR(logger,"Clearing the grid - free "<<sizeof(casacore::Complex)*its2dGrid.size()/1024/1024<<" MB/plane, #gridplanes "<<itsGrid.size());
        itsGrid.resize(0);
        its2dGrid.resize(0,0);
        itsGridIndex=-1;
    }
}

// This ShPtr should get deep-copied during cloning.
void TableVisGridder::initVisWeights(const IVisWeights::ShPtr &viswt)
{
    itsVisWeight = viswt;
}

// Customize for the specific type of Visibility weight.
// Input string is whatever is after "image" => "image.i.0.xxx" gives ".i.0.xxx "
void TableVisGridder::customiseForContext(const std::string &context)
{

    // RVU : Set up model dependant gridder behaviour
    //       For MFS, gridders for each Taylor term need different VisWeights.
    //  parse the 'context' string, and generate the "order" parameter.
  ASKAPLOG_DEBUG_STR(logger, "Customising gridder for context " << context);
    /*
    char corder[2];
    corder[0] = *(context.data()+3); // read the fourth character to get the order of the Taylor coefficient.
    corder[1] = '\n';
    int order = atoi(corder);
    //    ASKAPLOG_DEBUG_STR(logger, "Customising gridder for context " << context
    //              << " corder " << corder << " order" << order);
    if(order <0 || order >9) order = 0;
    */
    // MV: more general version
    const ImageParamsHelper iph(context);
    const int order = iph.isTaylorTerm() ? iph.order() : 1;

    if(itsVisWeight) {
        itsVisWeight->setParameters(order);
    }
}



/// This is the default implementation
int TableVisGridder::cIndex(int /*row*/, int /*pol*/, int /*chan*/) {
    return 0;
}

/// This is the default implementation
int TableVisGridder::gIndex(int /*row*/, int /*pol*/, int /*chan*/) {
    return 0;
}

/// @brief Obtain offset for the given convolution function
/// @details To conserve memory and speed the gridding up, convolution functions stored in the cache
/// may have an offset (i.e. essentially each CF should be defined on a bigger support and placed at a
/// pixel other than the centre of this support). This method returns this offset, which is the
/// position of the peak of the given CF on a bigger support w.r.t. the centre of this support.
/// The value of (0,0) means no offset from the centre (i.e. support is already centred).
/// @param[in] cfPlane plane of the convolution function cache to get the offset for
/// @return a pair with offsets for each axis
/// @note if there is no offset defined for a given cfPlane (default behavior), this method returns (0,0)
std::pair<int,int> TableVisGridder::getConvFuncOffset(int cfPlane) const
{
  ASKAPDEBUGASSERT(cfPlane>=0);
  if (cfPlane >= int(itsConvFuncOffsets.size())) {
      return std::pair<int,int>(0,0);
  }
  return itsConvFuncOffsets[cfPlane];
}

/// @brief initialise convolution function offsets for a given number of planes
/// @details The vector with offsets is resized and filled with (0,0).
/// @param[in] nPlanes number of planes in the cache
void TableVisGridder::initConvFuncOffsets(size_t nPlanes)
{
  itsConvFuncOffsets.resize(nPlanes);
  for (std::vector<std::pair<int,int> >::iterator it = itsConvFuncOffsets.begin(); it != itsConvFuncOffsets.end(); ++it) {
       it->first = 0;
       it->second = 0;
  }
}

/// @brief Assign offset to a particular convolution function
/// @details
/// @param[in] cfPlane plane of the convolution function cache to assign the offset for
/// @param[in] x offset in the first coordinate
/// @param[in] y offset in the second coordinate
/// @note For this method, cfPlane should be within the range [0..nPlanes-1].
void TableVisGridder::setConvFuncOffset(int cfPlane, int x, int y)
{
  ASKAPDEBUGASSERT(cfPlane>=0);
  ASKAPCHECK(cfPlane < int(itsConvFuncOffsets.size()), "An attempt to set offset for plane (before oversampling) "<<
             cfPlane<<" while the buffer has been initialised to handle "<<itsConvFuncOffsets.size()<<" planes only");
  itsConvFuncOffsets[cfPlane] = std::pair<int,int>(x,y);
}

/// @brief check whether the model is empty
/// @details A simple check allows us to bypass heavy calculations if the input model
/// is empty (all pixels are zero). This makes sense for degridding only.
/// @brief true, if the model is empty
bool TableVisGridder::isModelEmpty() const
{
  return itsModelIsEmpty;
}

}

}
