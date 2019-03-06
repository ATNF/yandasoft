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

#include <fitting/Params.h>
#include <dataaccess/SharedIter.h>
#include <measurementequation/ComponentEquation.h>
#include <fitting/INormalEquations.h>
#include <fitting/DesignMatrix.h>

#include <casacore/measures/Measures/Stokes.h>
#include <casacore/scimath/Mathematics/RigidVector.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/scimath/Mathematics/AutoDiff.h>
#include <casacore/scimath/Mathematics/AutoDiffMath.h>


#include <measurementequation/UnpolarizedPointSource.h>
#include <measurementequation/UnpolarizedGaussianSource.h>
#include <measurementequation/Calibrator1934.h>
#include <measurementequation/VectorOperations.h>


#include <stdexcept>

using askap::scimath::INormalEquations;
using askap::scimath::DesignMatrix;

namespace askap
{
  namespace synthesis
  {

    ComponentEquation::ComponentEquation(const askap::scimath::Params& ip,
          const accessors::IDataSharedIter& idi) :  scimath::Equation(ip), MultiChunkEquation(idi),  
           askap::scimath::GenericEquation(ip), GenericMultiChunkEquation(idi),
           itsAllComponentsUnpolarised(false)
    {
      init();
    };

    ComponentEquation::ComponentEquation(const accessors::IDataSharedIter& idi) :
           MultiChunkEquation(idi), GenericMultiChunkEquation(idi),
           itsAllComponentsUnpolarised(false) 
    {
      setParameters(defaultParameters());
      init();
    };
    
    void ComponentEquation::init()
    {
    }

askap::scimath::Params ComponentEquation::defaultParameters()
{
// The default parameters serve as a holder for the patterns to match the actual
// parameters. Shell pattern matching rules apply.
      askap::scimath::Params ip;
      ip.add("flux.i");
      ip.add("direction.ra");
      ip.add("direction.dec");
      ip.add("shape.bmaj");
      ip.add("shape.bmin");
      ip.add("shape.bpa");
      return ip;
}

/// @brief fill the cache of the components
/// @details This method converts the parameters into a vector of 
/// components. It is called on the first access to itsComponents
void ComponentEquation::fillComponentCache(
            std::vector<IParameterizedComponentPtr> &in) const
{ 
  const std::vector<std::string> completions(parameters().completions("flux.i"));
  const std::vector<std::string> calCompletions(parameters().completions("calibrator."));
  in.resize(completions.size() + calCompletions.size());
  if (!in.size()) {
     return;
  }
  
  // we will need to change this variable to false in the loop below, when
  // at least one polarised component is implemented.
  itsAllComponentsUnpolarised = true;
  
  // This loop is over all strings that complete the flux.i.* pattern
  // correctly. An exception will be throw if the parameters are not
  // consistent
  std::vector<IParameterizedComponentPtr>::iterator compIt=in.begin();
  for (std::vector<std::string>::const_iterator it=completions.begin();
        it!=completions.end();++it,++compIt)  {
          const std::string &cur = *it;
          const double ra=parameters().scalarValue("direction.ra"+cur);
          const double dec=parameters().scalarValue("direction.dec"+cur);
          const double fluxi=parameters().scalarValue("flux.i"+cur);
          const double bmaj = parameters().has("shape.bmaj"+cur) ? 
                   parameters().scalarValue("shape.bmaj"+cur) : 0.;
          const double bmin = parameters().has("shape.bmin"+cur) ? 
                   parameters().scalarValue("shape.bmin"+cur) : 0.;
          const double bpa = parameters().has("shape.bpa"+cur) ? 
                   parameters().scalarValue("shape.bpa"+cur) : 0.;
          
          if((bmaj>0.0)&&(bmin>0.0)) {
             // this is a gaussian
             compIt->reset(new UnpolarizedGaussianSource(cur,fluxi,ra,dec,bmaj,
                            bmin,bpa));
          } else {
             // this is a point source
             compIt->reset(new UnpolarizedPointSource(cur,fluxi,ra,dec));
          }
  }
  
  // loop over pre-defined calibrators
  for (std::vector<std::string>::const_iterator it=calCompletions.begin();
        it!=calCompletions.end();++it,++compIt)  {
        ASKAPCHECK(*it == "1934-638", "Only 1934-638 is currently supported, you requested "<<*it);
        compIt->reset(new Calibrator1934());
  }
}   

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
void ComponentEquation::addModelToCube(const IParameterizedComponent& comp,
       const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw,
       const casa::Vector<casa::Double>& freq,
       casa::Cube<casa::Complex> &rwVis) const
{
  ASKAPDEBUGASSERT(rwVis.nrow() == uvw.nelements());
  ASKAPDEBUGASSERT(rwVis.ncolumn() == freq.nelements());
                                
  ASKAPDEBUGASSERT(rwVis.nplane() == itsPolConverter.outputPolFrame().nelements());
  
  // flattened buffer for visibilities 
  std::vector<double> vis(2*freq.nelements()); 
 
  
  for (casa::uInt row=0;row<rwVis.nrow();++row) {
       casa::Matrix<casa::Complex> thisRow = rwVis.yzPlane(row);
       for (casa::Vector<casa::Stokes::StokesTypes>::const_iterator polIt = itsPolConverter.inputPolFrame().begin();
            polIt != itsPolConverter.inputPolFrame().end(); ++polIt) {
            // model given input Stokes
            comp.calculate(uvw[row],freq,*polIt,vis);

            // for most typically used transforms some elements will be zeros, use sparseTransform
            // instead of the full matrix to avoid heavy calculations in the loop
            const std::map<casa::Stokes::StokesTypes, casa::Complex> sparseTransform = 
                   itsPolConverter.getSparseTransform(*polIt); 
            for (std::map<casa::Stokes::StokesTypes, casa::Complex>::const_iterator ci=sparseTransform.begin();
                 ci!=sparseTransform.end(); ++ci) {
                 
                 const casa::uInt pol = polIndex(ci->first);
                 ASKAPDEBUGASSERT(pol < thisRow.ncolumn());
                                    
                 /// next command adds model visibilities to the
                 /// appropriate slice of the visibility cube. Conversions
                 /// between complex and two doubles are handled automatically
                 addScaledVector(vis,thisRow.column(pol),ci->second);
            }
       }          
  }
}               

/// @brief a helper method to populate a visibility cube
/// @details This is method computes visibilities for the one given
/// component and adds them to the cube provided. This is a second
/// version of the method. It is intended for unpolarised components
/// (i.e. it doesn't bother to add zeros)
///
/// @param[in] comp component to generate the visibilities for
/// @param[in] uvw baseline spacings, one triplet for each data row.
/// @param[in] freq a vector of frequencies (one for each spectral
///            channel) 
/// @param[in] rwVis a non-const reference to the visibility cube to alter
void ComponentEquation::addModelToCube(const IUnpolarizedComponent& comp,
       const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw,
       const casa::Vector<casa::Double>& freq,
       casa::Cube<casa::Complex> &rwVis) const
{
  ASKAPDEBUGASSERT(rwVis.nrow() == uvw.nelements());
  ASKAPDEBUGASSERT(rwVis.ncolumn() == freq.nelements());
  ASKAPDEBUGASSERT(rwVis.nplane() >= 1);
  
  // flattened buffer for visibilities 
  std::vector<double> vis(2*freq.nelements());
  
  // only Stokes I is of interest for the unpolarised component, use sparse transform
  const std::map<casa::Stokes::StokesTypes, casa::Complex> sparseTransform = 
        itsPolConverter.getSparseTransform(casa::Stokes::I); 

  for (casa::uInt row=0;row<rwVis.nrow();++row) {
       comp.calculate(uvw[row],freq,vis);
       //
       casa::Matrix<casa::Complex> thisRow = rwVis.yzPlane(row);
       
       ASKAPDEBUGASSERT(thisRow.ncolumn() == itsPolConverter.outputPolFrame().nelements());
       
       for (casa::uInt pol = 0; pol < thisRow.ncolumn(); ++pol) {
            const std::map<casa::Stokes::StokesTypes, casa::Complex>::const_iterator ci = 
                 sparseTransform.find(itsPolConverter.outputPolFrame()[pol]);
            if (ci != sparseTransform.end()) {
                        
                /// next command adds model visibilities to the
                /// appropriate slice of the visibility cube. Conversions
                /// between complex and two doubles are handled automatically
                addScaledVector(vis,thisRow.column(pol), ci->second);       
            }
       }
  }           
}

/// @brief helper method to return polarisation index in the visibility cube
/// @details The visibility cube may have various polarisation products and
/// ways of arranging them. This method extracts an index corresponding to the
/// given polarisation product. An exception is thrown, if the requested
/// product is not present.
/// @param[in] pol polarisation product
/// @return index (from 0 to nPol()-1)
casa::uInt ComponentEquation::polIndex(casa::Stokes::StokesTypes pol) const
{
  casa::Vector<casa::Stokes::StokesTypes> visCubeFrame = itsPolConverter.outputPolFrame();
  for (casa::uInt index = 0; index < visCubeFrame.nelements(); ++index) {
       if (visCubeFrame[index] == pol) {
           return index;
       }
  }
  ASKAPTHROW(AskapError, "Unable to find a match for polarisation product "<<pol);
}


/// @brief Predict model visibilities for one accessor (chunk).
/// @details This version of the predict method works with
/// a single chunk of data only. It seems that all measurement
/// equations should work with accessors rather than iterators
/// (i.e. the iteration over chunks should be moved to the higher
/// level, outside this class). In the future, I expect that
/// predict() without parameters will be deprecated.
/// @param chunk a read-write accessor to work with
void ComponentEquation::predict(accessors::IDataAccessor &chunk) const
{
  const std::vector<IParameterizedComponentPtr> &compList = 
         itsComponents.value(*this,&ComponentEquation::fillComponentCache);
  
  const casa::Vector<casa::Double>& freq = chunk.frequency();
  const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw = chunk.uvw();
  casa::Cube<casa::Complex> &rwVis = chunk.rwVisibility();
         
  // reset all visibility cube to 0
  rwVis.set(0.);
  
  // check whether the polarisation converter is valid
  if (!scimath::PolConverter::equal(itsPolConverter.outputPolFrame(),chunk.stokes())) {
      // reset pol converter because either polarisation frame has changed or 
      // this is the first use. The converter will be used inside addModelToCube shortly      
      itsPolConverter = scimath::PolConverter(scimath::PolConverter::canonicStokes(), chunk.stokes(), true);    
  }
         
  // loop over components
  for (std::vector<IParameterizedComponentPtr>::const_iterator compIt = 
       compList.begin(); compIt!=compList.end();++compIt) {
       
       ASKAPDEBUGASSERT(*compIt); 
       // current component
       const IParameterizedComponent& curComp = *(*compIt);
       try {
            const IUnpolarizedComponent &unpolComp = 
              dynamic_cast<const IUnpolarizedComponent&>(curComp);
            addModelToCube(unpolComp,uvw,freq,rwVis);
       }
       catch (const std::bad_cast&) {
            addModelToCube(curComp,uvw,freq,rwVis);
       }
  }
}


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
void ComponentEquation::updateDesignMatrixAndResiduals(
                   const IParameterizedComponent& comp,
                   const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw,
                   const casa::Vector<casa::Double>& freq,
                   scimath::DesignMatrix &dm, casa::Vector<casa::Double> &residual) const
{
  const size_t nParameters = comp.nParameters();
  // the number of polarisations in the visibility cube
  const casa::uInt nPol = itsPolConverter.outputPolFrame().nelements();
  ASKAPCHECK(nPol > 0, "Polarisation converter doesn't seem to be initialised");
  // number of data points  in the flattened vector
  const casa::uInt nData = nPol*uvw.nelements()*freq.nelements()*2;
  ASKAPDEBUGASSERT(nData!=0);
  
  // we only process Stokes I contribution if all components are unpolarised
  const casa::Vector<casa::Stokes::StokesTypes> inputPolFrame = 
             itsAllComponentsUnpolarised ? casa::Vector<casa::Stokes::StokesTypes>(1, casa::Stokes::I) :
             itsPolConverter.inputPolFrame();
  ASKAPDEBUGASSERT(residual.nelements() == nData);
      
  // Define AutoDiffs to buffer the output of a single call to the calculate 
  // method of the component.
  vector<casa::AutoDiff<double> > visDerivBufferSinglePol(2*freq.nelements(),
                          casa::AutoDiff<double>(0.,nParameters));
                          
  vector<vector<casa::AutoDiff<double> > > visDerivBuffer(nPol);
                            
  casa::Array<casa::Double> derivatives(casa::IPosition(2,nData, nParameters));
           
  for (casa::uInt row=0,offset=0; row<uvw.nelements(); ++row) {
       
       const casa::RigidVector<casa::Double, 3> &thisRowUVW = uvw[row];
  
       // initialise buffers for all polarisations
       for (casa::uInt pol = 0; pol<nPol; ++pol) {                
            visDerivBuffer[pol] = vector<casa::AutoDiff<double> >(2*freq.nelements()*nPol,
                          casa::AutoDiff<double>(0.,nParameters));
       }
                          
       for (casa::Vector<casa::Stokes::StokesTypes>::const_iterator polIt = inputPolFrame.begin();
            polIt != inputPolFrame.end(); ++polIt) {
            
            // these are contribitions to derivatives from polarisation *polInt
            comp.calculate(thisRowUVW,freq,*polIt,visDerivBufferSinglePol);

            // for most typically used transforms some elements will be zeros, use sparseTransform
            // instead of the full matrix to avoid heavy calculations in the loop
            const std::map<casa::Stokes::StokesTypes, casa::Complex> sparseTransform = 
                   itsPolConverter.getSparseTransform(*polIt);
                    
            for (std::map<casa::Stokes::StokesTypes, casa::Complex>::const_iterator ci=sparseTransform.begin();
                 ci!=sparseTransform.end(); ++ci) {
                 const casa::uInt index = polIndex(ci->first);
                 ASKAPDEBUGASSERT(index < nPol);
                 // the following is just a complex multiplication in the flattened form. We could probably
                 // done it in the similar way to everything else with expansion of vector operations, but
                 // do it explicitly for now for simplicity
                 const double factorRe = double(real(ci->second));
                 const double factorIm = double(imag(ci->second));
                                  
                 for (casa::uInt elem = 0; elem < 2*freq.nelements(); elem+=2) {
                      visDerivBuffer[index][elem] += visDerivBufferSinglePol[elem] * factorRe -
                            visDerivBufferSinglePol[elem+1] * factorIm; 
                      visDerivBuffer[index][elem + 1] += visDerivBufferSinglePol[elem] * factorIm +
                            visDerivBufferSinglePol[elem+1] * factorRe; 
                 }
            }            
       }    
       // visDerivBuffer is now filled with data for all polarisations
       
       for (casa::uInt pol=0; pol<nPol; ++pol,offset+=2*freq.nelements()) {
            // copy derivatives from the buffer for each parameter
            for (casa::uInt par=0; par<nParameters; ++par) {
                 // copy derivatives for each channel from visDerivBuffer
                 // to the appropriate slice of the derivatives Array
                 // template takes care of the actual types
                 copyDerivativeVector(par,visDerivBuffer[pol],  
                           derivatives(casa::IPosition(2,offset,par), 
                           casa::IPosition(2,offset+2*freq.nelements()-1,par)));                                 
            }
            // subtract contribution from the residuals
            // next command does: residual slice -= visDerivBuffer
            // taking care of all type conversions via templates
            subtractVector(visDerivBuffer[pol],
                  residual(casa::Slice(offset,2*freq.nelements())));
       }
  }
  // Now we can add the design matrix, residual, and weights
  for (casa::uInt par=0; par<nParameters; ++par) {
       dm.addDerivative(comp.parameterName(par), 
              derivatives(casa::IPosition(2,0,par),
                          casa::IPosition(2,nData-1,par)));
  }                
}

/// @brief Calculate the normal equation for one accessor (chunk).
/// @details This version of the method works on a single chunk of
/// data only (one iteration).It seems that all measurement
/// equations should work with accessors rather than iterators
/// (i.e. the iteration over chunks should be moved to the higher
/// level, outside this class). In the future, I expect that
/// the variant of the method without parameters will be deprecated.
/// @param[in] chunk a read-write accessor to work with
/// @param[in] ne Normal equations
void ComponentEquation::calcGenericEquations(const accessors::IConstDataAccessor &chunk,
                   askap::scimath::GenericNormalEquations& ne) const
{
  const std::vector<IParameterizedComponentPtr> &compList = 
         itsComponents.value(*this,&ComponentEquation::fillComponentCache);
  
  const casa::Vector<double>& freq=chunk.frequency();
  ASKAPDEBUGASSERT(freq.nelements()!=0);
  const casa::Vector<casa::RigidVector<casa::Double, 3> > &uvw = chunk.uvw();
  const casa::Cube<casa::Complex> &visCube = chunk.visibility();
                 
  const casa::uInt nPol = chunk.nPol();
  ASKAPDEBUGASSERT(nPol <= chunk.visibility().nplane());
  
  // check whether the polarisation converter is valid
  if (!scimath::PolConverter::equal(itsPolConverter.outputPolFrame(),chunk.stokes())) {
      // reset pol converter because either polarisation frame has changed or 
      // this is the first use. The converter will be used inside updateDesignMatrix shortly      
      itsPolConverter = scimath::PolConverter(scimath::PolConverter::canonicStokes(), chunk.stokes(), true);    
  }
  
      
  // Set up arrays to hold the output values
  // Two values (complex) per row, channel, pol
  const casa::uInt nData=chunk.nRow()*freq.nelements()*2*nPol;
  ASKAPDEBUGASSERT(nData!=0);
  casa::Vector<casa::Double> residual(nData);
      
  // initialize residuals with the observed visibilities          
  for (casa::uInt row=0,offset=0; row < chunk.nRow(); ++row) {
       for (casa::uInt pol=0; pol<nPol; ++pol,offset+=2*freq.nelements()) {
            // the following command copies visibility slice to the appropriate 
            // residual slice converting complex to two doubles automatically
            // via templates 
            copyVector(visCube.xyPlane(pol).row(row),
                  residual(casa::Slice(offset,2*freq.nelements())));      
       } 
  }
                
  DesignMatrix designmatrix; // old parameters: parameters();
  for (std::vector<IParameterizedComponentPtr>::const_iterator compIt = 
            compList.begin(); compIt!=compList.end();++compIt) {
       ASKAPDEBUGASSERT(*compIt); 
       updateDesignMatrixAndResiduals(*(*compIt),uvw,freq,designmatrix,
                           residual);
  }
  casa::Vector<double> weights(nData,1.);
  designmatrix.addResidual(residual, weights);
  ne.add(designmatrix);
}

/// @brief read-write access to parameters
/// @details This method is overridden to invalidate component cache.
/// @return a non-const reference to Param::ShPtr
const scimath::Params::ShPtr& ComponentEquation::rwParameters() const throw()
{ 
  itsComponents.invalidate();
  return scimath::Equation::rwParameters();
}

/// Clone this into a shared pointer
/// @return shared pointer to a copy
ComponentEquation::ShPtr ComponentEquation::clone() const
{
  return ComponentEquation::ShPtr(new ComponentEquation(*this));
}


} // namespace synthesis

} // namespace askap
