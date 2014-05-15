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

#ifndef A_PROJECT_GRIDDER_BASE_H
#define A_PROJECT_GRIDDER_BASE_H

#include <boost/shared_ptr.hpp>
#include <gridding/IBasicIllumination.h>
#include <dataaccess/IConstDataAccessor.h>
#include <gridding/UVPattern.h>
#include <gridding/IVisGridder.h>
#include <Common/ParameterSet.h>
#include <measurementequation/SynthesisParamsHelper.h>

namespace askap {
namespace synthesis {

/// @brief Common functionality for all mosaicing gridders
/// @details AProjectGridderBase class encapsulates common operations for all mosaicing 
/// gridders: CF cache support and recalculation statistics, support for the buffer in the uv-space,
/// and the factory of illumination pattrns.
/// @note This class is still abstract, all real mosaicing gridders should be derived (virtually)
///       from this class
/// @ingroup gridding
struct AProjectGridderBase : virtual public IVisGridder 
{
  /// @brief initialise common part for mosaicing gridders
  /// @param[in] maxFeeds Maximum number of feeds allowed
  /// @param[in] maxFields Maximum number of fields allowed
  /// @param[in] pointingTol Pointing tolerance in radians
  /// @param[in] paTol Parallactic angle tolerance in radians
  /// @param[in] freqTol Frequency tolerance (relative, threshold for df/f), negative value 
  ///        means the frequency axis is ignored 
  explicit AProjectGridderBase(const int maxFeeds=1, const int maxFields=1, 
                               const double pointingTol=0.0001,                              
                               const double paTol=0.01, const double freqTol = 1e-6);
  
  /// @brief copy constructor
  /// @details It is needed because we have a shared pointer as a data member and want to
  /// clone the object instead of copying the reference as if it would be by default.
  /// @param[in] other input object
  AProjectGridderBase(const AProjectGridderBase &other);
  
  /// @brief destructor
  /// @details We print cache usage stats here. No specific destruction is required for any data member
  virtual ~AProjectGridderBase(); 
  
  /// @brief check if given CF is valid
  /// @details
  /// @param[in] feed feed number to query
  /// @param[in] field field number to query
  inline bool isCFValid(int feed, int field) const { return itsDone(feed,field);}
  
  /// @brief pointing for given feed and field
  /// @details
  /// @param[in] feed feed number to query
  /// @param[in] field field number to query
  inline const casa::MVDirection& pointing(int feed, int field) const { return itsPointings(feed,field);}
  
  /// @brief obtain uv-pattern
  /// @return a reference to uv-pattern
  /// @note One has to initialise uv-pattern at least once before calling this method.
  inline UVPattern& uvPattern() const { ASKAPDEBUGASSERT(itsPattern); return *itsPattern;}
 
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
  void initUVPattern(casa::uInt uSize, casa::uInt vSize, double uCellSize,
                     double vCellSize, casa::uInt overSample = 1);
  
  /// @brief obtain current field
  /// @details Although it is not great, we use the fact that only one field (i.e. dish pointing)
  /// can be represented by a single accessor. It is the case in the current implementation, but
  /// is not, strictly speaking, required by the interface. This class encapsulate all 
  /// related functionality to detect the field change. This method returns the field 
  /// corresponding to the accessor passed during the last call to indexField.
  /// @return current field index
  inline casa::uInt currentField() const { return itsCurrentField; }
  
  /// @brief checks whether the current field has been updated
  /// @details See currentField for more detailed description.
  /// @param[in] acc input const accessor to analyse
  void indexField(const accessors::IConstDataAccessor &acc);
  
  /// @brief check whether CF cache is valid
  /// @details This methods validates CF cache for one particular iteration. If necessary, 
  /// all values in itsDone are set to false. This method also sets some internal flags to
  /// update the stats correctly when updateStats is called. 
  /// @param[in] acc input const accessor to analyse
  /// @param[in] symmetric true, if illumination pattern is symmetric, false otherwise
  void validateCFCache(const accessors::IConstDataAccessor &acc, bool symmetric);

  /// @brief toggle the validity flag for a given CF
  /// @details
  /// @param[in] feed feed number to query
  /// @param[in] field field number to query
  inline void makeCFValid(int feed, int field) { itsDone(feed,field) = true;}
  
  /// @brief update statistics
  /// @details This class maintains cache rebuild statistics. It is impossible to update them 
  /// directly in validateCFCache because a priori it is not known how many CFs are recalculated
  /// following invalidation. It depends on the actual algorithm and the dataset. To keep track
  /// of the cache rebuild stats call this method with the exact number of CFs calculated.
  /// @param[in] nDone number of convolution functions rebuilt at this iteration
  void updateStats(casa::uInt nDone);
   
  /// @brief a helper factory of illumination patterns
  /// @details Illumination model is required for a number of gridders. This
  /// method allows to avoid duplication of code and encapsulates all 
  /// functionality related to illumination patterns. 
  /// @param[in] parset ParameterSet containing description of illumination to use
  /// @return shared pointer to illumination interface
  static boost::shared_ptr<IBasicIllumination> 
      makeIllumination(const LOFAR::ParameterSet &parset);
protected:
  /// @brief helper method to reset CF cache
  inline void resetCFCache() { itsDone.set(false); }

  /// @brief actual factory of derived gridders
  /// @details Gridders derived from this class use exactly the same parameters, but doing
  /// the job using slightly different algorithms. This templated method allows to create both
  /// of them without duplication of the code extracting parameters from the parset file.
  /// @note Parset should have the gridder name already removed from all parameter names (i.e.
  /// a proper subset should be formed first, see IVisGridder::createGridder as well)
  /// @param[in] parset parameter set containing description of the gridder to create
  /// @return shared pointer to the gridder
  template<typename GridderType>
  static boost::shared_ptr<GridderType> createAProjectGridder(const LOFAR::ParameterSet &parset);

  /// @brief read-write access to the buffer of slopes
  /// @return non-const reference to the cube of slopes
  casa::Cube<double>& rwSlopes() { return itsSlopes;}
          
private:
  /// @brief assignment operator (never to be called)
  /// @details It is defined as private, so we can't call it and use copy constructor instead.
  /// @param[in] other input object
  /// @return reference to itself
  AProjectGridderBase& operator=(const AProjectGridderBase &other);


  /// Pointing tolerance in radians
  double itsPointingTolerance;
   
  /// @brief parallactic angle tolerance (in radians)
  /// @details If new angle is different from the one used to compute the cache for more
  /// than this value, the cache will be recomputed. Note, a negative value means to 
  /// always recalculate for asymmetric illumination patterns
  double itsParallacticAngleTolerance;
   
  /// Last field processed
  int itsLastField;
  /// Current field processed
  casa::uInt itsCurrentField;
   
  /// @brief flags that CF is valid for given feed and fields
  casa::Matrix<bool> itsDone;
   
  /// Pointing for this feed and field
  casa::Matrix<casa::MVDirection> itsPointings;
   
  /// @brief buffer in the uv-space
  /// @details It is used to compute convolution functions (buffer for illumination pattern)
  boost::shared_ptr<UVPattern> itsPattern;
      
  // stats for CF cache rebuilds
  /// @brief number of iterations when CFs were generated
  /// @details This number is incremented for each accessor which leads to recomputation of the 
  /// CF cache. In the best case (CFs computed once and
  /// reused later) it should be equal to 1. In the worst case (CFs are
  /// recomputed for every iteration) this should be equal to the number iterations 
  /// (to be exact, the number of calls to initConvolutionFunctions).
  casa::uInt itsNumberOfCFGenerations;
      
  /// @brief total number of iterations
  /// @details This number is incremented each time a new accessor is passed to this gridder
  /// (i.e. counts every iteration). 
  casa::uInt itsNumberOfIterations;
   
  /// @brief number of CFs generated due to parallactic angle change
  /// @details This number is incremented each time a CF is recomputed 
  /// due to a change in parallactic angle. 
  casa::uInt itsNumberOfCFGenerationsDueToPA;
      
  /// @brief parallactic angle for which the cache is valid
  /// @details This buffer is only used and filled if the illumination pattern is asymmetric.
  /// We currently don't account for the VLBI case with notably different parallactic angles.
  /// Therefore, only one angle is stored here. 
  casa::Float itsCFParallacticAngle;
      
  /// @brief number of CFs generated due to a change of frequency 
  /// @details This number is incremented each time a CF is recomputed following a change
  /// in frequency setup.
  casa::uInt itsNumberOfCFGenerationsDueToFreq;
   
  /// @brief relative frequency tolerance
  /// @details If abs(df/f) exceeds this value for any spectral channel, the cache of 
  /// the CFs has to be recomputed. A negative value means that checks of frequency axis are
  /// to be bypassed (i.e. frequency is assumed to be always valid). The latter option may be
  /// useful for smaller ATCA datasets to avoid unnecessary overheads 
  double itsFrequencyTolerance;
      
  /// @brief frequency axis corresponding to the cache
  /// @details To be able to determine whether the cache of CFs needs to be recomputed one has
  /// to know which frequencies were used to compute it. This is a copy of the frequency axis
  /// returned by the accessor on the first iteration
  casa::Vector<casa::Double> itsCachedFrequencies;
                   
  /// @brief flag showing that CFs are rebuilt due to parallactic angle
  bool itsCFInvalidDueToPA;
  
  /// @brief flag showing that CFs are rebuilt due to frequency axis change
  bool itsCFInvalidDueToFreq; 

  /// @brief cube of slopes
  casa::Cube<double> itsSlopes;            
  
};

} // namespace synthesis

} // namespace askap

#include <gridding/AProjectGridderBase.tcc>

#endif // #ifndef A_PROJECT_GRIDDER_BASE_H

