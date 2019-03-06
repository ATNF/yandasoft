/// @file
/// 
/// @brief A measurement equation, which generates a gaussian noise.
/// @details It is required for simulations to be able to add noise to
/// the simulated visibilities. To do it via measurement equations, one
/// has to create a composite measurement equation via SumOfTwoMEs class
/// with one of the input measurement equations set to an instance of
/// GaussianNoiseME defined here. If we need various similar classes the
/// approach probably needs to be changed to something similar to
/// CalibrationME template/effect classes.
/// @note The random number generator is a member of this class. In a 
/// parallel environment this would lead to a number of independent
/// generators used and to the same sequence generated in parallel 
/// branches of code. One needs a global solution (with an internode 
/// comminicaton on the cluster for a proper simulation of random numbers).
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

#ifndef GAUSSIAN_NOISE_ME_H
#define GAUSSIAN_NOISE_ME_H

// own includes
#include <measurementequation/IMeasurementEquation.h>
#include <utils/ComplexGaussianNoise.h>

// casa includes
#include <casacore/casa/BasicMath/Random.h>
#include <casacore/casa/BasicSL/Complex.h>

namespace askap {

namespace synthesis {

/// @brief A measurement equation, which generates a gaussian noise.
/// @details It is required for simulations to be able to add noise to
/// the simulated visibilities. To do it via measurement equations, one
/// has to create a composite measurement equation via SumOfTwoMEs class
/// with one of the input measurement equations set to an instance of
/// GaussianNoiseME defined here. If we need various similar classes the
/// approach probably needs to be changed to something similar to
/// CalibrationME template/effect classes.
/// @note The random number generator is a member of this class. In a 
/// parallel environment this would lead to a number of independent
/// generators used and to the same sequence generated in parallel 
/// branches of code. One needs a global solution (with an internode 
/// comminicaton on the cluster for a proper simulation of random numbers).
/// @ingroup measurementequation
struct GaussianNoiseME : public IMeasurementEquation
{
  /// @brief constructor, initializes random distribution required.
  /// @param[in] variance required variance of the noise (same as rms
  /// squared here because the mean is always zero)
  /// @param[in] seed1 a first seed to initialize the random generator
  /// @param[in] seed2 a second seed to initialize the random generator 
  explicit GaussianNoiseME(double variance, casa::Int seed1 = 0, 
                                            casa::Int seed2 = 10);

  /// @brief constructor, initializes random distribution required.
  /// @details The required noise rms is obtained from the accessor. The
  /// object constructed this way can simulate noise with different statistics
  /// for different visibilities. of the noise (same as rms
  /// squared here because the mean is always zero)
  /// @param[in] seed1 a first seed to initialize the random generator
  /// @param[in] seed2 a second seed to initialize the random generator 
  explicit GaussianNoiseME(casa::Int seed1 = 0, casa::Int seed2 = 10);
  
  
  /// @brief Predict model visibilities for one accessor (chunk).
  /// @details This prediction is done for single chunk of data only. 
  /// It seems that all measurement equations should work with accessors 
  /// rather than iterators (i.e. the iteration over chunks should be 
  /// moved to the higher level, outside this class). 
  /// @param[in] chunk a read-write accessor to work with
  virtual void predict(accessors::IDataAccessor &chunk) const;

  /// @brief Calculate the normal equation for one accessor (chunk).
  /// @details This calculation is done for a single chunk of
  /// data only (one iteration).It seems that all measurement
  /// equations should work with accessors rather than iterators
  /// (i.e. the iteration over chunks should be moved to the higher
  /// level, outside this class). 
  /// @param[in] chunk a read-write accessor to work with
  /// @param[in] ne Normal equations
  virtual void calcEquations(const accessors::IConstDataAccessor &chunk,
                          askap::scimath::INormalEquations& ne) const;
protected:

  /// @brief a helper method to obtain a random complex number.
  /// @details It runs the generator twice for real and imaginary part,
  /// composes a complex number and returns it.
  /// @return a random complex number
  inline casa::Complex getRandomComplexNumber() const 
  {
    return itsGen();
  }
  
private:
  /// @brief random number generator
  scimath::ComplexGaussianNoise itsGen;
  
  /// @brief true, if the variance given explicitly is to be used
  /// @details If the noise distribution is the same for all visibilities,
  /// the class can be set up with some value of variance given explicitly
  /// (it is passed to the random number generator). Alternatively, we can
  /// scale the actual random numbers with the noise figure returned by the
  /// accessor, but simulate them with the variance of 1.
  bool itsExplicitVariance;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef GAUSSIAN_NOISE_ME_H
