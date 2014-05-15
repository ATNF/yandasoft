/// @file
/// @brief Helper class to analyse performance for CF generation
/// @details This template is derived from the gridder class
/// under test. The main method forces contunuous recalculation of
/// convolution functions for the selected maximum number of beams and
/// w-planes. It is therefore possible to profile this core operation more
/// accurately. 
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

#ifndef TEST_CF_GEN_PERFORMANCE 
#define TEST_CF_GEN_PERFORMANCE

#include <gridding/IVisGridder.h>
#include <gridding/AWProjectVisGridder.h>
#include <boost/shared_ptr.hpp>
#include <dataaccess/DataAccessorStub.h>
#include <Common/ParameterSet.h>


namespace askap {

namespace synthesis {

/// @brief Helper class to analyse performance for CF generation
/// @details This template is derived from the gridder class
/// under test. The main method forces contunuous recalculation of
/// convolution functions for the selected maximum number of beams and
/// w-planes. It is therefore possible to profile this core operation more
/// accurately. 
/// @ingroup gridding
class TestCFGenPerformance : public AWProjectVisGridder {
public:
  /// @brief constructor to initialise common part for mosaicing gridders
  /// @details
  /// @param illum  Antenna illumination model
  /// @param wmax Maximum baseline (wavelengths)
  /// @param nwplanes Number of w planes
  /// @param cutoff Cutoff in determining support e.g. 10^-3 of the peak
  /// @param overSample Oversampling (currently limited to <=1)
  /// @param maxSupport Maximum support to be allowed
  /// @param limitSupport Upper limit of support
  /// @param maxFeeds Maximum number of feeds allowed
  /// @param maxFields Maximum number of fields allowed
  /// @param pointingTol Pointing tolerance in radians
  /// @param paTol Parallactic angle tolerance in radians
  /// @param freqTol Frequency tolerance (relative, threshold for df/f), negative value 
  ///        means the frequency axis is ignored       
  /// @param frequencyDependent Frequency dependent gridding?
  /// @param name Name of table to save convolution function into
  TestCFGenPerformance(const boost::shared_ptr<IBasicIllumination const> &illum,
                      const double wmax, const int nwplanes, const double cutoff,
                      const int overSample, const int maxSupport, const int limitSupport,
                      const int maxFeeds=1, const int maxFields=1, const double pointingTol=0.0001,
                      const double paTol=0.01,
                      const double freqTol = 1e-6,          
                      const bool frequencyDependent=true, 
                      const std::string& name=std::string("")) : 
                       AWProjectVisGridder(illum, wmax,nwplanes,cutoff,overSample, maxSupport, limitSupport,
                            maxFeeds, maxFields, pointingTol, paTol, freqTol, frequencyDependent, name),
                       itsNBeams(maxFeeds) {}
                       
  
  /// @brief copy constructor
  /// @details
  /// @param[in] other another instance to copy from
  TestCFGenPerformance(const TestCFGenPerformance &other) : AWProjectVisGridder(other), itsNBeams(other.itsNBeams) {}
  
  /// @brief clone object, e.g. to test performance with a number of copies
  boost::shared_ptr<TestCFGenPerformance> clone() const;

  /// @brief static method to create test gridder
  /// @details This method mimics createGridder interface for a-projection gridders
  /// @param[in] parset input parset file
  /// @return a shared pointer to the instance
  static boost::shared_ptr<TestCFGenPerformance> createGridder(const LOFAR::ParameterSet& parset);

  /// @brief main method which initialises CFs
  /// @details 
  /// @param[in] nRuns number of initialisations
  void run(const int nRuns = 1);
  
  /// @brief Initialise the parameters and accessor
  /// @param axes axes specifications
  /// @param shape Shape of output image: u,v,pol,chan
  /// @param dopsf Make the psf?
  virtual void initialiseGrid(const scimath::Axes& axes,
               const casa::IPosition& shape, const bool dopsf=true);
  
protected:
  /// @brief initialise the accessor
  void init();
    
private:
  /// @brief number of beams
  int itsNBeams;
   
  /// @brief accessor used in tests
  accessors::DataAccessorStub itsAccessor;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef TEST_CF_GEN_PERFORMANCE 

