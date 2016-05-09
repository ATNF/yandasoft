/// @file
///
/// ComponentEquation: Equation for dealing with discrete components such
/// as point sources and Gaussians.
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

#ifndef SYNCOMPONENTEQUATION_H_
#define SYNCOMPONENTEQUATION_H_

// own include
#include <fitting/Params.h>
#include <fitting/DesignMatrix.h>


#include <dataaccess/CachedAccessorField.h>

#include <measurementequation/IParameterizedComponent.h>
#include <measurementequation/IUnpolarizedComponent.h>
#include <measurementequation/GenericMultiChunkEquation.h>
#include <utils/PolConverter.h>

// casa includes
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>

namespace askap
{
  namespace synthesis
  {

    /// @brief Visibility processing for components
    ///
    /// @details This class does predictions and calculates normal equations
    /// for discrete components such as point sources and Gaussians.
    /// Names are flux.{i,q,u,v}, direction.{ra,dec}, shape.{bmaj,bmin,bpa}
    /// etc.
    ///
    /// @ingroup measurementequation
    
    class ComponentEquation : public GenericMultiChunkEquation
    {
      public:

        /// @brief Standard constructor using the parameters and the
        /// data iterator.
        /// @param ip Parameters
        /// @param idi data iterator
        ComponentEquation(const askap::scimath::Params& ip,
          const accessors::IDataSharedIter& idi);

        /// @brief Constructor using default parameters
        /// @param idi data iterator
        ComponentEquation(const accessors::IDataSharedIter& idi);
        
        /// Return the default parameters
        static askap::scimath::Params defaultParameters();
        
        /// @brief Predict model visibilities for one accessor (chunk).
        /// @details This version of the predict method works with
        /// a single chunk of data only. It seems that all measurement
        /// equations should work with accessors rather than iterators
        /// (i.e. the iteration over chunks should be moved to the higher
        /// level, outside this class). In the future, I expect that
        /// predict() without parameters will be deprecated.
        /// @param[in] chunk a read-write accessor to work with
        virtual void predict(accessors::IDataAccessor &chunk) const;

        /// @brief Calculate the normal equation for one accessor (chunk).
        /// @details This version of the method works on a single chunk of
        /// data only (one iteration).It seems that all measurement
        /// equations should work with accessors rather than iterators
        /// (i.e. the iteration over chunks should be moved to the higher
        /// level, outside this class). In the future, I expect that
        /// the variant of the method without parameters will be deprecated.
        /// @param[in] chunk a read-write accessor to work with
        /// @param[in] ne Normal equations
        virtual void calcGenericEquations(const accessors::IConstDataAccessor &chunk,
                              askap::scimath::GenericNormalEquations& ne) const;
        
        using GenericMultiChunkEquation::predict;
        using askap::scimath::GenericEquation::calcEquations;
        using GenericMultiChunkEquation::calcGenericEquations;
       
        /// Clone this into a shared pointer
        /// @return shared pointer to a copy
        virtual ComponentEquation::ShPtr clone() const;
        
      private:
        /// Initialize this object
        virtual void init();
    protected:
        /// a short cut to shared pointer on a parameterized component
        typedef boost::shared_ptr<IParameterizedComponent> IParameterizedComponentPtr;
    
        /// @brief fill the cache of the components
        /// @details This method convertes the parameters into a vector of 
        /// components. It is called on the first access to itsComponents
        void fillComponentCache(std::vector<IParameterizedComponentPtr> &in) const;
        
        /// @brief helper method to return polarisation index in the visibility cube
        /// @details The visibility cube may have various polarisation products and
        /// ways of arranging them. This method extracts an index corresponding to the
        /// given polarisation product. An exception is thrown, if the requested
        /// product is not present.
        /// @param[in] pol polarisation product
        /// @return index (from 0 to nPol()-1)
        casa::uInt polIndex(casa::Stokes::StokesTypes pol) const;
              
        
        /// @brief a helper method to populate a visibility cube
        /// @details This is method computes visibilities for the one given
        /// component and adds them to the cube provided. This is the most
        /// generic method, which iterates over polarisations. An overloaded
        /// version of the method do the same for unpolarised components
        /// (i.e. it doesn't bother to add zeros)
        ///
        /// @param[in] comp component to generate the visibilities for
        /// @param[in] uvw baseline spacings, one triplet for each data row.
        /// @param[in] freq a vector of frequencies (one for each spectral
        ///            channel) 
        /// @param[in] rwVis a non-const reference to the visibility cube to alter
        void addModelToCube(const IParameterizedComponent& comp,
               const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw,
               const casa::Vector<casa::Double>& freq,
               casa::Cube<casa::Complex> &rwVis) const;

        /// @brief a helper method to populate a visibility cube
        /// @details This is method computes visibilities for the one given
        /// component and adds them to the cube provided. This is a second
        /// version of the method. It is intended for unpolarised components
        /// (i.e. it doesn't bother to add zeros)
        ///
        /// @param[in] comp component to generate the visibilities for
        /// @param[in] uvw baseline spacings, one triplet for each data row.
        /// @param[in] freq a vector of frequencies (one frequency for each 
        ///            spectral channel) 
        /// @param[in] rwVis a non-const reference to the visibility cube to alter
        void addModelToCube(const IUnpolarizedComponent& comp,
               const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw,
               const casa::Vector<casa::Double>& freq,
               casa::Cube<casa::Complex> &rwVis) const;
        
        /// @brief a helper method to update design matrix and residuals
        /// @details This method iterates over a given number of polarisation 
        /// products in the visibility cube. It updates the design matrix with
        /// derivatives and subtracts values from the vector of residuals.
        /// The latter is a flattened vector which should have a size of 
        /// 2*nChan*nPol*nRow. Spectral channel is the most frequently varying
        /// index, then follows the polarisation index, and the least frequently
        /// varying index is the row. The number of channels and the number of
        /// rows always corresponds to that of the visibility cube. The number of
        /// polarisations can be less than the number of planes in the cube to
        /// allow processing of incomplete data cubes (or unpolarised components). In contrast to 
        /// 
        /// @param[in] comp component to generate the visibilities for
        /// @param[in] uvw baseline coorindates for each row
        /// @param[in] freq a vector of frequencies (one frequency for each
        ///            spectral channel)
        /// @param[in] dm design matrix to update (to add derivatives to)
        /// @param[in] residual vector of residuals to update 
        void updateDesignMatrixAndResiduals(
                   const IParameterizedComponent& comp,
                   const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw,
                   const casa::Vector<casa::Double>& freq,
                   scimath::DesignMatrix &dm, casa::Vector<casa::Double> &residual) const;
        
        /// @brief read-write access to parameters
        /// @details This method is overridden to invalidate component cache.
        /// @return a non-const reference to Param::ShPtr
        virtual const scimath::Params::ShPtr& rwParameters() const throw();
        
    private:   
        /// @brief vector of components plugged into this component equation
        /// this has nothing to do with data accessor, we just reuse the class
        /// for a cached field
        accessors::CachedAccessorField<std::vector<IParameterizedComponentPtr> > itsComponents;     
        
        /// @brief True if all components are unpolarised
        mutable bool itsAllComponentsUnpolarised;
        
        /// @brief polarisation converter to be used with this component equation
        /// @details Components are defined in the Stokes frame, this class converts them
        /// into the measurement frame
        mutable scimath::PolConverter itsPolConverter;
    };

  }

}

#endif
