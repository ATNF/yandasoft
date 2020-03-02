/// @file
///
/// Actual algorithm implementing delay solver tool 
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

#ifndef ASKAP_UTILITIES_DELAY_SOLVER_IMPL_H
#define ASKAP_UTILITIES_DELAY_SOLVER_IMPL_H

// boost
#include <boost/noncopyable.hpp>

// std
#include <utility>

// casa
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/BasicSL/Complex.h>

// own
#include <dataaccess/IConstDataAccessor.h>
#include <utils/DelayEstimator.h>

namespace askap {

namespace utils {

/// @brief class solving for delays
/// @details This class accumulates data averaged down to a given resolution and solves for delays
/// when requested.
/// @ingroup utils
class DelaySolverImpl : public boost::noncopyable {
public:
    /// @brief constructor
    /// @param[in] targetRes target spectral resolution in Hz, data are averaged to match the desired resolution 
    /// note, integral number of channels are averaged.
    /// @param[in] pol polarisation index to use
    /// @param[in] ampCutoff if positive, amplitudes above ampCutoff will be flagged
    /// @param[in] refAnt reference antenna index   
    explicit DelaySolverImpl(double targetRes, casa::Stokes::StokesTypes pol = casa::Stokes::XX, 
                             float ampCutoff = -1., casa::uInt refAnt = 1);
    
    /// @brief set baselines to exclude
    /// @details An empty vector configures the class to take all available baselines into account.
    /// @param[in] bsln vector with pair of indices representing baselines to exclude from the solution
    void excludeBaselines(const casa::Vector<std::pair<casa::uInt, casa::uInt> > &bsln);
     
    /// @brief process one data accessor
    /// @param[in] acc data accessor to process
    void process(const accessors::IConstDataAccessor &acc);
    
    /// @brief solve for antenna-based delays
    /// @details This method estimates delays for all baselines and then solves for
    /// antenna-based delays honouring baselines to be excluded.
    /// @param[in] useFFT if true, FFT-based delay estimator is used. It is less accurate but
    /// is more robust for large delays and less sensitive to flagged data. 
    /// @return a vector with one delay per antenna (antennas are in the index-increasing order).
    casa::Vector<double> solve(bool useFFT) const; 
    
    /// @brief set initial delay approximation
    /// @details The class can optionally remove some a priori known (approximate) delay in 
    /// full resolution data before averaging takes place. This allows to get a coarse delay 
    /// first in full resolution data and then refine it with
    /// averaging. Empty array means zero initial delay for all antennas (i.e. nothing special is done,
    /// this is the default). There should be one value per antenna. An exception is thrown if
    /// antenna ID greater than or equal to the size of the vector is encountered
    /// @param[in] delays in seconds
    void setApproximateDelays(const casa::Vector<double> &delays);
    
    /// @brief initialise the accumulation
    /// @details This method reverts the state of the object to that before the first call to process method.
    /// This allows to repeat delay estimation with a set of approximate delays.
    void init();
    
    /// @brief set target resolution
    /// @details
    /// @param[in] targetRes target spectral resolution in Hz, data are averaged to match the desired resolution 
    /// note, integral number of channels are averaged.
    void setTargetResolution(double targetRes);
    
protected:

    /// @brief helper method to check that all channels/rows are flagged
    /// @param[in] flags matrix with flags
    /// @return true, if all channels and rows are flagged
    static bool checkAllFlagged(const casa::Matrix<bool> &flags);
    
    /// @brief helper method to obtain delay approximation for the given row
    /// @details Zero is returned if itsDelayApproximation is empty. Otherwise,
    /// delay for the given row is extracted (based on metadata contained in the buffers).
    /// Note, delay units are the same as itsDelayApproximation (i.e. seconds).
    /// @param[in] row row of interest
    double delayApproximation(casa::uInt row) const;
        
private:

    /// @brief target spectral resolution in Hz
    /// @details Data are downsampled as close as possible (only integral number of channels are averaged) to this resolution
    double itsTargetRes;
    
    /// @brief polarisation product to use
    casa::Stokes::StokesTypes itsPol;
    
    /// @brief if positive, amplitudes aboves this cutoff are flagged
    float itsAmpCutoff;
    
    /// @brief reference antenna index
    casa::uInt itsRefAnt;

    /// @brief baselines to exclude (empty vector is all are taken into account)
    casa::Vector<std::pair<casa::uInt, casa::uInt> > itsExcludedBaselines;

    /// @brief summed or averaged spectra
    casa::Matrix<casa::Complex> itsSpcBuffer;
    
    /// @brief counts per channel per row
    casa::Matrix<casa::uInt> itsAvgCounts;

    /// @brief number of accessors processed
    casa::uInt itsNAvg;
    
    /// @brief spectral axis (to ensure it doesn't change throughout the dataset)
    casa::Vector<double> itsFreqAxis;
    
    /// @brief first antenna indices (to ensure they're the same for all iterations)
    casa::Vector<casa::uInt> itsAnt1IDs;

    /// @brief second antenna indices (to ensure they're the same for all iterations)
    casa::Vector<casa::uInt> itsAnt2IDs;
    
    /// @brief delay estimator
    scimath::DelayEstimator itsDelayEstimator;
    
    /// @brief number of spectral channels to average
    casa::uInt itsChanToAverage;
    
    /// @brief a priori delay information for all antennas
    /// @details If this vector is not empty, it contains initial approximation of delays
    /// which are taken out of the full resolution data before averaging takes place. This
    /// allows to get a coarse delay first in full resolution data and then refine it with
    /// averaging. Empty array means zero initial delay (i.e. nothing special is done).
    casa::Vector<double> itsDelayApproximation;  
};

} // namespace utils

} // namespace askap

#endif // #ifndef ASKAP_UTILITIES_DELAY_SOLVER_IMPL_H

