/// @file
/// @brief Common functionality for all mosaicing gridders
/// @details AProjectGridderBase class encapsulates common operations for all mosaicing 
/// gridders: CF cache support and recalculation statistics, support for the buffer in the uv-space,
/// and the factory of illumination pattrns.
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
ASKAP_LOGGER(logger, ".gridding.aprojectgridderbase");

#include <gridding/AProjectGridderBase.h>
#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <gridding/DiskIllumination.h>
#include <gridding/ATCAIllumination.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <profile/AskapProfiler.h>


using namespace askap;
using namespace askap::synthesis;
using namespace askap::accessors;

/// @brief initialise common part for mosaicing gridders
/// @param[in] maxFeeds Maximum number of feeds allowed
/// @param[in] maxFields Maximum number of fields allowed
/// @param[in] pointingTol Pointing tolerance in radians
/// @param[in] paTol Parallactic angle tolerance in radians
/// @param[in] freqTol Frequency tolerance (relative, threshold for df/f), negative value 
///        means the frequency axis is ignored 
AProjectGridderBase::AProjectGridderBase(const int maxFeeds, const int maxFields, 
                     const double pointingTol, const double paTol, const double freqTol) :
          itsPointingTolerance(pointingTol),  itsParallacticAngleTolerance(paTol),
          itsLastField(-1), itsCurrentField(0),
          itsDone(maxFeeds, maxFields, false), itsPointings(maxFeeds, maxFields, casa::MVDirection()),
          itsNumberOfCFGenerations(0), itsNumberOfIterations(0), 
          itsNumberOfCFGenerationsDueToPA(0), itsCFParallacticAngle(0),
          itsNumberOfCFGenerationsDueToFreq(0), itsFrequencyTolerance(freqTol),
          itsCFInvalidDueToPA(false), itsCFInvalidDueToFreq(false), itsSlopes(2, maxFeeds, maxFields,0.)
{
  ASKAPCHECK(maxFeeds>0, "Maximum number of feeds must be one or more");
  ASKAPCHECK(maxFields>0, "Maximum number of fields must be one or more");
  ASKAPCHECK(pointingTol>0.0, "Pointing tolerance must be greater than 0.0");
}

/// @brief copy constructor
/// @details It is needed because we have a shared pointer as a data member and want to
/// clone the object instead of copying the reference as if it would be by default.
/// @param[in] other input object
AProjectGridderBase::AProjectGridderBase(const AProjectGridderBase &other) : 
    IVisGridder(other),
    itsPointingTolerance(other.itsPointingTolerance),
    itsParallacticAngleTolerance(other.itsParallacticAngleTolerance),
    itsLastField(other.itsLastField), itsCurrentField(other.itsCurrentField),
    itsDone(other.itsDone.copy()), itsPointings(other.itsPointings.copy()), 
    itsNumberOfCFGenerations(other.itsNumberOfCFGenerations),
    itsNumberOfIterations(other.itsNumberOfIterations),
    itsNumberOfCFGenerationsDueToPA(other.itsNumberOfCFGenerationsDueToPA), 
    itsCFParallacticAngle(other.itsCFParallacticAngle),
    itsNumberOfCFGenerationsDueToFreq(other.itsNumberOfCFGenerationsDueToFreq),
    itsFrequencyTolerance(other.itsFrequencyTolerance),
    itsCachedFrequencies(other.itsCachedFrequencies),
    itsCFInvalidDueToPA(other.itsCFInvalidDueToPA),
    itsCFInvalidDueToFreq(other.itsCFInvalidDueToFreq), itsSlopes(other.itsSlopes.copy())
{
  if (other.itsPattern) {
      itsPattern.reset(new UVPattern(*(other.itsPattern)));
  }
}
  
/// @brief destructor
/// @details We print cache usage stats here. No specific destruction is required for any data member
AProjectGridderBase::~AProjectGridderBase()
{
  size_t nUsed = 0;
  for (casa::uInt feed = 0; feed<itsDone.nrow(); ++feed) {
      for (casa::uInt field = 0; field<itsDone.ncolumn(); ++field) {
           if (isCFValid(feed,field)) {
               ++nUsed;
            }
      }
  }  
  if (itsDone.nelements()) {
      ASKAPLOG_INFO_STR(logger, "   AProjectGridderBase: CF cache memory utilisation (last iteration): "<<
              double(nUsed)/double(itsDone.nrow()*itsDone.ncolumn())*100<<"% of maxfeed*maxfield");
  }
  
  if (itsNumberOfIterations != 0) {
      ASKAPLOG_INFO_STR(logger, "   AProjectGridderBase: CFs were rebuilt "<<
             itsNumberOfCFGenerations<<" times for "<<itsNumberOfIterations<<" iterations");
      ASKAPLOG_INFO_STR(logger, "   Last iteration worked with "<<nUsed<<" CFs");        
      if (itsNumberOfCFGenerations != 0) {
          ASKAPLOG_INFO_STR(logger, "   Parallactic angle change caused "<<
                  itsNumberOfCFGenerationsDueToPA<<" of those rebuilds ("<<
                  double(itsNumberOfCFGenerationsDueToPA)/double(itsNumberOfCFGenerations)*100<<
                  " %)");
          ASKAPLOG_INFO_STR(logger, "   Frequency axis change caused "<<
                  itsNumberOfCFGenerationsDueToFreq<<" of those rebuilds ("<<
                  double(itsNumberOfCFGenerationsDueToFreq)/double(itsNumberOfCFGenerations)*100<<
                  " %)");
      }   
      if (nUsed != 0) { 
          // because nUsed is strictly speaking applicable to the last iteration only we need
          // to filter out rediculous values (and warn the user that the result is approximate
          // anyway)
          const double utilisation = (1.-double(itsNumberOfCFGenerations)/
                                  double(itsNumberOfIterations*nUsed));
          if ((utilisation<1.) && (utilisation>0.)) {
              ASKAPLOG_INFO_STR(logger, "   Approximate CF cache utilisation is "<<
                                        utilisation*100.<<" %");
          }
      }
  }
}


/// @brief set up buffer in the uv-space
/// @details To work with illumination patterns we need a buffer. Moving initialisation
/// out of the loop allows to improve the performance. This method is supposed to be called
/// as soon as all necessary parameters are known.
/// @param[in] uSize size in the direction of u-coordinate
/// @param[in] vSize size in the direction of v-coordinate
/// @param[in] uCellSize size of the uv-cell in the direction of 
///            u-coordinate (in wavelengths)
/// @param[in] vCellSize size of the uv-cell in the direction of 
///            v-coordinate (in wavelengths)
/// @param[in] overSample oversampling factor (default is 1)
void AProjectGridderBase::initUVPattern(casa::uInt uSize, casa::uInt vSize, double uCellSize,
                     double vCellSize, casa::uInt overSample)
{
  itsPattern.reset(new UVPattern(uSize,vSize, uCellSize,vCellSize,overSample));
}

/// @brief checks whether the current field has been updated
/// @details See currentField for more detailed description.
/// @param[in] acc input const accessor to analyse
void AProjectGridderBase::indexField(const IConstDataAccessor &acc)
{
  ASKAPDEBUGTRACE("AProjectGridderBase::indexField");
  // Validate cache using first row only
  bool newField = true;
  ASKAPDEBUGASSERT(acc.nRow()>0);

  casa::uInt firstFeed = acc.feed1()(0);
  ASKAPCHECK(firstFeed<itsDone.nrow(), "Too many feeds: increase maxfeeds");
  casa::MVDirection firstPointing = acc.pointingDir1()(0);

  for (int field=itsLastField; field>-1; --field) {
       if (firstPointing.separation(pointing(firstFeed, field))<itsPointingTolerance) {
           itsCurrentField = field;
           newField = false;
           break;
       }
  }
  if (newField) {
      ++itsLastField;
      ASKAPDEBUGASSERT(itsLastField>=0);
      itsCurrentField = itsLastField;
      ASKAPCHECK(itsCurrentField < itsDone.ncolumn(),
              "Too many fields: increase maxfields " << itsDone.ncolumn());
      itsPointings(firstFeed, itsCurrentField) = firstPointing;
      ASKAPLOG_DEBUG_STR(logger, "Found new field " << itsCurrentField<<" at "<<
                printDirection(firstPointing));
  } 
}

/// @brief check whether CF cache is valid
/// @details This methods validates CF cache for one particular iteration. If necessary, 
/// all values in itsDone are set to false. This method also sets some internal flags to
/// update the stats correctly when updateStats is called. 
/// @param[in] acc input const accessor to analyse
/// @param[in] symmetric true, if illumination pattern is symmetric, false otherwise
void AProjectGridderBase::validateCFCache(const IConstDataAccessor &acc, bool symmetric)
{
  ASKAPDEBUGTRACE("AProjectGridderBase::validateCFCache");
  const int nSamples = acc.nRow();
 
  // flags are used to accumulate CF rebuild statistics
  itsCFInvalidDueToPA = false;
    
  if (!symmetric) {
      // need to check parallactic angles here
      const casa::Vector<casa::Float> &feed1PAs = acc.feed1PA();
      ASKAPDEBUGASSERT(feed1PAs.nelements() == casa::uInt(nSamples));
      for (int row = 0; row<nSamples; ++row) {
           if (fabs(feed1PAs[row] - itsCFParallacticAngle) > itsParallacticAngleTolerance) {
               itsCFInvalidDueToPA = true;
               itsCFParallacticAngle = feed1PAs[row];
               itsDone.set(false);
               break;
           }
      }
  }
    
  // the following flag is used to accululate CF rebuild statistics and internal logic
  itsCFInvalidDueToFreq = false;
    
  // don't bother checking if the cache is rebuilt anyway
  if (!itsCFInvalidDueToPA && (itsFrequencyTolerance >= 0.)) {
      const casa::Vector<casa::Double> &freq = acc.frequency();
      if (freq.nelements() != itsCachedFrequencies.nelements()) {
          itsCFInvalidDueToFreq = true;
      } else {
          // we can also write the following using iterators, if necessary
          for (casa::uInt chan = 0; chan<freq.nelements(); ++chan) {
               const casa::Double newFreq = freq[chan];
               ASKAPDEBUGASSERT(newFreq > 0.);
               if ( fabs(itsCachedFrequencies[chan] - newFreq)/newFreq > itsFrequencyTolerance) {
                    itsCFInvalidDueToFreq = true;
                    break;
               }
          }
      } 
      if (itsCFInvalidDueToFreq) {
          itsDone.set(false);
      }
  }
    
  // cache the current frequency axis if the cache is going to be built
  // do nothing if the tolerance is negative 
  if ((itsCFInvalidDueToPA || itsCFInvalidDueToFreq) && (itsFrequencyTolerance >= 0.)) {
      itsCachedFrequencies.assign(acc.frequency().copy());
  } 
}

/// @brief assignment operator (never to be called)
/// @details It is defined as private, so we can't call it and use copy constructor instead.
/// @param[in] other input object
/// @return reference to itself
AProjectGridderBase& AProjectGridderBase::operator=(const AProjectGridderBase &)
{
  ASKAPTHROW(AskapError, "This method is not supposed to be called!");
  return *this;
}


/// @brief update statistics
/// @details This class maintains cache rebuild statistics. It is impossible to update them 
/// directly in validateCFCache because a priori it is not known how many CFs are recalculated
/// following invalidation. It depends on the actual algorithm and the dataset. To keep track
/// of the cache rebuild stats call this method with the exact number of CFs calculated.
/// @param[in] nDone number of convolution functions rebuilt at this iteration
void AProjectGridderBase::updateStats(casa::uInt nDone)
{
  ++itsNumberOfIterations;
  itsNumberOfCFGenerations += nDone;
  if (itsCFInvalidDueToPA) {
      itsNumberOfCFGenerationsDueToPA += nDone;
  }    
  if (itsCFInvalidDueToFreq) {
      itsNumberOfCFGenerationsDueToFreq += nDone;
  }
}

/// @brief a helper factory of illumination patterns
/// @details Illumination model is required for a number of gridders. This
/// method allows to avoid duplication of code and encapsulates all 
/// functionality related to illumination patterns. 
/// @param[in] parset ParameterSet containing description of illumination to use
/// @return shared pointer to illumination interface
boost::shared_ptr<IBasicIllumination> 
AProjectGridderBase::makeIllumination(const LOFAR::ParameterSet &parset)
{
   const std::string illumType = parset.getString("illumination", "disk");
   const double diameter=SynthesisParamsHelper::convertQuantity(parset.getString("diameter"),"m");
   const double blockage=SynthesisParamsHelper::convertQuantity(parset.getString("blockage"),"m");

   if (illumType == "disk") {
   	    ASKAPLOG_INFO_STR(logger,
					"Using disk illumination model, diameter="<<
					diameter<<" metres, blockage="<<blockage<<" metres");
   
       	return boost::shared_ptr<IBasicIllumination>(new DiskIllumination(diameter,blockage));
   } else if (illumType == "ATCA") {
   	    ASKAPLOG_INFO_STR(logger,
					"Using ATCA illumination model, diameter="<<
					diameter<<" metres, blockage="<<blockage<<" metres");

   	    boost::shared_ptr<ATCAIllumination> illum(new ATCAIllumination(diameter,blockage)); 
   	    ASKAPDEBUGASSERT(illum);
   	    if (parset.getBool("illumination.tapering", true)) {
   	        const double maxDefocusingPhase =
   	             SynthesisParamsHelper::convertQuantity(parset.getString("illumination.tapering.defocusing",
   	                           "0rad"),"rad");
	        illum->simulateTapering(maxDefocusingPhase);
	        ASKAPLOG_INFO_STR(logger,"Tapering of the illumination is simulated, maximum defocusing phase = "<<
	                  maxDefocusingPhase/M_PI*180.<<" deg."); 
	    } else {
	        ASKAPLOG_INFO_STR(logger,"Tapering of the illumination is not simulated");
	    }
	    if (parset.getBool("illumination.feedlegs", true)) {
	        const double width = SynthesisParamsHelper::convertQuantity(
	           parset.getString("illumination.feedlegs.width","1.8m"),"m");
	        const double rotation = SynthesisParamsHelper::convertQuantity(
	           parset.getString("illumination.feedlegs.rotation","45deg"),"rad");   
	        const double shadowingFactor = 
	           parset.getDouble("illumination.feedlegs.shadowing",0.75);   
	        illum->simulateFeedLegShadows(width,rotation,shadowingFactor);
	        ASKAPLOG_INFO_STR(logger,"Feed legs are simulated. Width = "<<width<<" metres, rotated at "<<
	           rotation/M_PI*180.<<" deg, shadowing factor (how much attenuation caused) = "<<shadowingFactor);
	        if (parset.getBool("illumination.feedlegs.wedges", true)) {
	            const double defaultWedgeShadowing[2] = {0.6,0.5};
	            std::vector<double> wedgeShadowing = 
	                parset.getDoubleVector("illumination.feedlegs.wedges.shadowing", 
	                std::vector<double>(defaultWedgeShadowing,defaultWedgeShadowing+2));
	            const double angle = SynthesisParamsHelper::convertQuantity(
	                parset.getString("illumination.feedlegs.wedges.angle","15deg"),"rad");    
	            const double startRadius = SynthesisParamsHelper::convertQuantity(
	                parset.getString("illumination.feedlegs.wedges.startradius","3.5m"),"m");
	            ASKAPCHECK(wedgeShadowing.size() && wedgeShadowing.size()<3, 
	                 "illumination.feedlegs.wedges.shadowing can have either 1 or 2 elements only, "
	                 "you have "<<wedgeShadowing.size());
	            if (wedgeShadowing.size() == 1) {
	                wedgeShadowing.push_back(wedgeShadowing[0]);
	            }     
	            ASKAPDEBUGASSERT(wedgeShadowing.size() == 2);    
          	    illum->simulateFeedLegWedges(wedgeShadowing[0],wedgeShadowing[1],angle,startRadius);	            
          	    ASKAPLOG_INFO_STR(logger,"Feed leg wedges are simulated. Shadowing factors are "<<
          	           wedgeShadowing<<", opening angle is "<<angle/M_PI*180.<<" deg, start radius is "<<
          	           startRadius<<" metres");
	        } else {
	            ASKAPLOG_INFO_STR(logger,"Feed leg wedges are not simulated.");
	        }
	    } else {
	       ASKAPLOG_INFO_STR(logger,"Feed legs are not simulated.");
	    }
	    return illum;
   }
   
   ASKAPTHROW(AskapError, "Unknown illumination type "<<illumType);
   return boost::shared_ptr<IBasicIllumination>(); // to keep the compiler happy
}



