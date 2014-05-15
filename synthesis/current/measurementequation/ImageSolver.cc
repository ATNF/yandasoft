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

#include <measurementequation/ImageSolver.h>
#include <utils/PaddingUtils.h>
#include <profile/AskapProfiler.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>

ASKAP_LOGGER(logger, ".measurementequation.imagesolver");

#include <askap/AskapError.h>

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>

#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/Lattices/LatticeFFT.h>
#include <utils/MultiDimArrayPlaneIter.h>

using namespace askap;
using namespace askap::scimath;

#include <iostream>

#include <cmath>
using std::abs;

#include <map>
#include <vector>
#include <string>
#include <stdexcept>

using std::map;
using std::vector;
using std::string;

namespace askap
{
  namespace synthesis
  {

    ImageSolver::ImageSolver() :
      itsZeroWeightCutoffArea(false), itsZeroWeightCutoffMask(true), itsSaveIntermediate(true)
    {
    }

    void ImageSolver::init()
    {
      resetNormalEquations();
    }
    
    // Add a new element to a list of preconditioners
    void ImageSolver::addPreconditioner(IImagePreconditioner::ShPtr pc)
    {
	    // Add a new element to the map of preconditioners
	    int n = itsPreconditioners.size();
	    itsPreconditioners[n+1] = pc;
    }

    /// @brief perform normalization of the dirty image and psf
	/// @details This method divides the PSF and dirty image by the diagonal of the Hessian.
	/// If a non-void shared pointer is specified for the mask parameter, this method assigns
	/// 0. for those elements where truncation of the weights has been performed and 1. 
	/// otherwise. 
	/// @param[in] diag diagonal of the Hessian (i.e. weights), dirty image will be
	///            divided by an appropriate element of the diagonal or by a cutoff
	///            value
	/// @param[in] tolerance cutoff value given as a fraction of the largest diagonal element
	/// @param[in] psf  point spread function, which is normalized to unity
    /// @param[in] psfRefPeak peak value of the reference PSF before normalisation
    ///            negative value means to take max(psf). PSF is normalised to max(psf)/psfRefPeak
	/// @param[in] dirty dirty image which is normalized by truncated weights (diagonal)
	/// @param[out] mask shared pointer to the output mask showing where the truncation has 
	///             been performed.
	/// @return peak of PSF before normalisation (to be used as psfRefPeak, if necessary)
	/// @note although mask is filled in inside this method it should already have a correct 
	/// size before this method is called. Pass a void shared pointer (default) to skip 
	/// mask-related functionality. Hint: use utility::NullDeleter to wrap a shared pointer
	/// over an existing array reference.
    float ImageSolver::doNormalization(const casa::Vector<double>& diag, const float& tolerance, 
                      casa::Array<float>& psf, float psfRefPeak, casa::Array<float>& dirty, 
				      const boost::shared_ptr<casa::Array<float> > &mask) const
    {
        ASKAPTRACE("ImageSolver::doNormalization");

        const double maxDiag(casa::max(diag));
        const double sumDiag(casa::sum(diag));

        ASKAPCHECK(maxDiag>0., "Maximum diagonal element is supposed to be positive, check that at least some data were gridded, maxDiag="
                   <<maxDiag<<" sumDiag="<<sumDiag);
	ASKAPLOG_INFO_STR(logger, "Solid angle = " << sumDiag/maxDiag << " pixels");

        const double cutoff=tolerance*maxDiag;

	    /// The PSF is just an approximation calculated from a subset of the data. So we
	    /// we allowed to normalize the peak to unity.
	
	//	ASKAPLOG_INFO_STR(logger, "Normalizing PSF by maximum of diagonal equal to "<<maxDiag<<
	//			  ", cutoff weight is "<<tolerance*100<<"\% of the largest diagonal element");
	//	ASKAPLOG_INFO_STR(logger, "Peak of PSF before normalization = " << casa::max(psf));                      
	//        psf /= float(maxDiag);
	//        ASKAPLOG_INFO_STR(logger, "Peak of PSF = " << casa::max(psf));
	
        ASKAPLOG_INFO_STR(logger, "Maximum diagonal element " <<maxDiag<<
			  ", cutoff weight is "<<tolerance*100<<"\% of the largest diagonal element");
        const float unnormalisedMaxPSF = casa::max(psf);	  
        if (psfRefPeak<=0.) {
            ASKAPLOG_INFO_STR(logger, "Normalising PSF to unit peak");
        } else {
            ASKAPLOG_INFO_STR(logger, "Normalising PSF to be "<<unnormalisedMaxPSF/psfRefPeak<<
                        " (psfRefPeak = "<<psfRefPeak<<")");
        }
         
	    ASKAPLOG_INFO_STR(logger, "Peak of PSF before normalisation = " << unnormalisedMaxPSF);
        psf /= float(psfRefPeak<=0. ? unnormalisedMaxPSF : psfRefPeak);
        ASKAPLOG_INFO_STR(logger, "Peak of PSF after normalisation = " << casa::max(psf));
	
	    uint nAbove = 0;
        casa::IPosition vecShape(1,diag.nelements());
        casa::Vector<float> dirtyVector(dirty.reform(vecShape));

        if(mask) {
          ASKAPLOG_DEBUG_STR(logger, "Deconvolution will be based on signal to noise ratio");
        }
        else {
          ASKAPLOG_DEBUG_STR(logger, "Deconvolution will be based on flux only");
        }
        
        ASKAPLOG_DEBUG_STR(logger,"Peak of the dirty vector before normalisation "<<casa::max(dirtyVector));

        casa::Vector<float> maskVector(mask ? mask->reform(vecShape) : casa::Vector<float>());
        for (casa::uInt elem=0; elem<diag.nelements(); ++elem) {
             if(diag(elem)>cutoff) {
                dirtyVector(elem)/=diag(elem);
                if (mask) {
                    maskVector(elem)=sqrt(diag(elem)/float(maxDiag));
                }
                nAbove++;
             } else {
                // a number of actions are possible in the weights cutoff area depending
                // on the state of the solver
                if (zeroWeightCutoffArea()) {
                    dirtyVector(elem) = 0.;
                } else {
                    dirtyVector(elem)/=maxDiag;
                }
                if (mask) {
                    if (zeroWeightCutoffMask()) {
                        maskVector(elem)=0.0;
                    } else {
                        maskVector(elem)=sqrt(tolerance);
                    }
                } // if mask required
             }  // if element > cutoff
        } // loop over elements (image pixels)
        ASKAPLOG_INFO_STR(logger, "Normalized dirty image by truncated weights image");
        if (mask) {
            ASKAPLOG_INFO_STR(logger, "Converted truncated weights image to clean mask");
        } // if mask required
        ASKAPLOG_INFO_STR(logger, 100.0*float(nAbove)/float(diag.nelements()) << "% of the pixels were above the cutoff " << cutoff);
        ASKAPLOG_DEBUG_STR(logger,"Peak of the dirty vector after normalisation "<<casa::max(dirtyVector));
        return unnormalisedMaxPSF;
    }

    // Apply all the preconditioners in the order in which they were created.
    bool ImageSolver::doPreconditioning(casa::Array<float>& psf, casa::Array<float>& dirty) const
    {
        ASKAPTRACE("ImageSolver::doPreconditioning");

        casa::Array<float> oldPSF(psf.copy());
	    bool status=false;
	    for(std::map<int, IImagePreconditioner::ShPtr>::const_iterator pciter=itsPreconditioners.begin(); pciter!=itsPreconditioners.end(); pciter++)
	    {
	      status = status | (pciter->second)->doPreconditioning(psf,dirty);
	    }
	    // we could write the result to the file or return it as a parameter (but we need an image name
	    // here to compose a proper parameter name)
	    //	    if (status) {
	    //	        sensitivityLoss(oldPSF, psf);
	    //	    } else {
	    //	        ASKAPLOG_INFO_STR(logger, "No preconditioning has been done, hence sensitivity loss factor is 1.");
	    //	    }
	    return status;
    }

    /// @brief Solve for parameters, updating the values kept internally
    /// The solution is constructed from the normal equations. The parameters named 
    /// image* are interpreted as images and solved for.
    /// @param[in] ip current model (to be updated)
    /// @param[in] quality Solution quality information
    /// @note Solve for update simply by scaling the data vector by the diagonal term of the
    /// normal equations i.e. the residual image
    bool ImageSolver::solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality)
    {
      ASKAPTRACE("ImageSolver::solveNormalEquations");

      ASKAPLOG_INFO_STR(logger, "Calculating principal solution");

      // Solving A^T Q^-1 V = (A^T Q^-1 A) P
      uint nParameters=0;

      // Find all the free parameters beginning with image
      vector<string> names(ip.completions("image"));
      map<string, uint> indices;

      for (vector<string>::const_iterator it=names.begin(); it!=names.end(); it++)
      {
        string name="image"+*it;
        if (ip.isFree(name))
        {
          indices[name]=nParameters;
          nParameters+=ip.value(name).nelements();
        }
      }
      ASKAPCHECK(nParameters>0, "No free parameters in ImageSolver");

      for (map<string, uint>::const_iterator indit=indices.begin(); indit
          !=indices.end(); indit++) {
        // Axes are dof, dof for each parameter
        //casa::IPosition arrShape(itsParams->value(indit->first).shape());
        for (scimath::MultiDimArrayPlaneIter planeIter(ip.value(indit->first).shape());
             planeIter.hasMore(); planeIter.next()) {
        
             ASKAPCHECK(normalEquations().normalMatrixDiagonal().count(indit->first)>0, "Diagonal not present for solution");
             casa::Vector<double>  diag(normalEquations().normalMatrixDiagonal().find(indit->first)->second);
             ASKAPCHECK(normalEquations().dataVector(indit->first).size()>0, "Data vector not present for solution");
             casa::Vector<double> dv = normalEquations().dataVector(indit->first);
	         ASKAPCHECK(normalEquations().normalMatrixSlice().count(indit->first)>0, "PSF Slice not present");
             casa::Vector<double> slice(normalEquations().normalMatrixSlice().find(indit->first)->second);

             if (planeIter.tag() != "") {
                 // it is not a single plane case, there is something to report
                 ASKAPLOG_INFO_STR(logger, "Processing plane "<<planeIter.sequenceNumber()<<
                                           " tagged as "<<planeIter.tag());
             }


	         casa::Array<float> dirtyArray(planeIter.planeShape());
             casa::convertArray<float, double>(dirtyArray, planeIter.getPlane(dv));
             casa::Array<float> psfArray(planeIter.planeShape());
             casa::convertArray<float, double>(psfArray, planeIter.getPlane(slice));
	
             // Do the preconditioning
             const bool wasPreconditioning = doPreconditioning(psfArray,dirtyArray);
             
             // Normalize by the diagonal
             doNormalization(planeIter.getPlaneVector(diag),tol(),psfArray,dirtyArray);
             
             if (wasPreconditioning) {	
                 // Save the new PSF in parameter class to be saved to disk later
               if(saveIntermediate()) {
                 saveArrayIntoParameter(ip,indit->first,planeIter.shape(),"psf.image",
                                        psfArray,planeIter.position());
               }
             } // if - wasPreconditioning

             ASKAPLOG_INFO_STR(logger, "Peak data vector flux (derivative) "<<max(dirtyArray));

             // uncomment the code below to save the residual image
             // This takes up some memory and we have to ship the residual image out inside
             // the parameter class. Therefore, we may not need this functionality in the 
             // production version (or may need to implement it in a different way).

             if(saveIntermediate())
             {
                saveArrayIntoParameter(ip,indit->first,planeIter.shape(),"residual",
                                       dirtyArray,planeIter.position());
             }
             // end of the code storing residual image

	         casa::Vector<double> value(planeIter.getPlaneVector(ip.value(indit->first)));
             const casa::Vector<float> dirtyVector(dirtyArray.reform(value.shape()));
             for (uint elem=0; elem<dv.nelements(); ++elem) {
                  value(elem) += dirtyVector(elem);
             }
        } // iteration over image planes (polarisation, spectral channels)
      } // iteration over free image parameters (e.g. facets)
      quality.setDOF(nParameters);
      quality.setRank(0);
      quality.setCond(0.0);
      quality.setInfo("Scaled residual calculated");
      
      /// Save the PSF and Weight
      saveWeights(ip);
      savePSF(ip);
      return true;
    }
        
    /// @brief helper method to save part of the NE
    /// @details We need to save slice and diagonal of the normal equations as 
    /// PSF and weights image (savePSF and saveWeights methods) with very similar
    /// operations. This method encapsulates the common code to avoid duplication.
    /// It iterates over all parameters with names starting with "image". 
    /// @param[in] ip model (to be updated with the appropriate parameter) 
    /// @param[in] prefix name prefix for stored parameter (image will be replaced with this prefix, i.e.
    /// "psf" or "weights"
    /// @param[in] nePart part of the normal equations to save (map of vectors)    
    void ImageSolver::saveNEPartIntoParameter(askap::scimath::Params& ip, const std::string &prefix,
                      const std::map<std::string, casa::Vector<double> > &nePart) const 
    {
      ASKAPTRACE("ImageSolver::saveNEPartIntoParameter");

      // get list of all image parameters and iterate over it
      const vector<string> names(ip.completions("image"));
      for (vector<string>::const_iterator it=names.begin(); it!=names.end(); ++it) {
           const std::string name = "image" + *it;
           std::map<std::string, casa::Vector<double> >::const_iterator parIt = nePart.find(name);
        
           if (parIt != nePart.end()) {
               std::map<std::string, casa::IPosition>::const_iterator shpIt = normalEquations().shape().find(name);
               ASKAPDEBUGASSERT(shpIt != normalEquations().shape().end());
               const casa::IPosition arrShape(shpIt->second);
               std::string parName = prefix + *it;
               casa::Array<double> temp(parIt->second.reform(arrShape));
               if (ip.has(parName)) {
                   ip.update(parName, temp);
               } else {
                   scimath::Axes axes(ip.axes(name));
                   ip.add(parName, temp, axes);
               }
           }
      }
    }
    
    /// @brief helper method to save a given array
    /// @details This method encapsulates common functionality to store a given array
    /// as a parameter or a part of the parameter. The idea is similar to saveNEPartIntoParameter,
    /// but functionality is slightly different (the key is to allow storage of a part of the parameter)
    /// @param[in] ip model (to be updated with the appropriate parameter) 
    /// @param[in] imgName image parameter name (to take axes from and to for the output name)
    /// @param[in] shape full shape of the output parameter (arr may be a part of the full parameter)
    /// @param[in] prefix name prefix for stored parameter ("image" in imgName will be replaced by prefix)
    /// @param[in] arr array to save
    /// @param[in] pos position to save (arr is a part of the parameter)
    void ImageSolver::saveArrayIntoParameter(askap::scimath::Params& ip, const std::string &imgName, 
              const casa::IPosition &shape, const std::string &prefix, const casa::Array<double> &arr, 
              const casa::IPosition &pos)
    {
      ASKAPDEBUGTRACE("ImageSolver::saveArrayIntoParameter");

      ASKAPDEBUGASSERT(imgName.find("image")==0);
      ASKAPCHECK(imgName.size()>5, 
                 "Image parameter name should have something appended to word image")           
	  const string parName = prefix + imgName.substr(5);
	  if (!ip.has(parName)) {
	      // create an empty parameter with the full shape
          scimath::Axes axes(ip.axes(imgName));
	      ip.add(parName, shape, axes);
	  }
      ip.update(parName, arr, pos);      
    }
    
    /// @brief helper method to save a given array
    /// @details This variant of the method intended for single precision arrays which are expanded to
    /// double precision inside this method and the double precision version of the method is called to
    /// do the work
    /// @param[in] ip model (to be updated with the appropriate parameter) 
    /// @param[in] imgName image parameter name (to take axes from and to for the output name)
    /// @param[in] shape full shape of the output parameter (arr may be a part of the full parameter)
    /// @param[in] prefix name prefix for stored parameter ("image" in imgName will be replaced by prefix)
    /// @param[in] arr array to save
    /// @param[in] pos position to save (arr is a part of the parameter)
    void ImageSolver::saveArrayIntoParameter(askap::scimath::Params& ip, const std::string &imgName, 
              const casa::IPosition &shape, const std::string &prefix, const casa::Array<float> &arr, 
              const casa::IPosition &pos)
    {
      casa::Array<double> tmpArr(arr.shape());
      casa::convertArray<double,float>(tmpArr,arr);
      saveArrayIntoParameter(ip, imgName, shape, prefix, tmpArr, pos);
    }

    /// @note itsPreconditioner is not cloned, only the ShPtr is.
    Solver::ShPtr ImageSolver::clone() const
    {
      return Solver::ShPtr(new ImageSolver(*this));
    }

    /// @return a reference to normal equations object
    /// @note In this class and derived classes the type returned
    /// by this method is narrowed to always provide image-specific 
    /// normal equations objects
    const scimath::ImagingNormalEquations& ImageSolver::normalEquations() const
    {
      try
      {
        return dynamic_cast<const scimath::ImagingNormalEquations&>(
            Solver::normalEquations());
      }
      catch (const std::bad_cast &bc)
      {
        ASKAPTHROW(AskapError, "An attempt to use incompatible normal "
            "equations class with image solver");
      }
    }
    
    /// @brief a helper method to extract the first plane out of the multi-dimensional array
    /// @details This method just uses MultiDimArrayPlaneIter to extract the first plane
    /// out of the array. It accepts a const reference to the array (which is a conseptual const).
    /// @param[in] in const reference to the input array
    /// @return the array with the first plane
    casa::Array<float> ImageSolver::getFirstPlane(const casa::Array<float> &in)
    {
      casa::Array<float> nonConstArray(in);
      return MultiDimArrayPlaneIter::getFirstPlane(nonConstArray);
    }
    
    
    /// @brief estimate sensitivity loss due to preconditioning
    /// @details Preconditioning (i.e. Wiener filter, tapering) makes the synthesized beam look nice, 
    /// but the price paid is a sensitivity loss. This method gives an estimate (accurate calculations require
    /// gridless weights, which we don't have in our current approach). The method just requires the two
    /// PSFs before and after preconditioning.
    /// @param[in] psfOld an array with original psf prior to preconditioning
    /// @param[in] psfNew an array with the psf after preconditioning has been applied
    /// @return sensitivity loss factor (should be grater than or equal to 1)
    double ImageSolver::sensitivityLoss(const casa::Array<float>& psfOld, const casa::Array<float>& psfNew)
    {
       ASKAPTRACE("ImageSolver::sensitivityLoss");
       ASKAPLOG_INFO_STR(logger, "Estimating sensitivity loss due to preconditioning");
       // current code can't handle cases where the noise is not uniform
       // we need to think about a better approach
       // we also assume that input PSFs are normalized to the same peak value (i.e. 1).
       
       // work with the first slice only if the array is multi-dimensional
       if ((psfOld.shape().nonDegenerate().nelements() >= 2) ||
           (psfNew.shape().nonDegenerate().nelements() >= 2)) {
           ASKAPLOG_INFO_STR(logger, "Sensitivity loss estimate will use a single plane of a multi-dimensional PSF image");      
       }
       casa::Array<float> psfOldSlice = getFirstPlane(psfOld);
       casa::Array<float> psfNewSlice = getFirstPlane(psfNew);
       
       casa::IPosition paddedShape = psfOldSlice.shape();
       ASKAPCHECK(paddedShape == psfNewSlice.shape(), 
            "sensitivityLoss: shapes of two PSFs are supposed to be the same, you have "<<paddedShape<<
            " and "<<psfNew.shape());
       ASKAPDEBUGASSERT(paddedShape.nonDegenerate().nelements() >= 2);     
       paddedShape(0) *= 2;     
       paddedShape(1) *= 2;     
       casa::ArrayLattice<casa::Complex> uvOld(paddedShape);
       casa::ArrayLattice<casa::Complex> uvNew(paddedShape);
       
       casa::ArrayLattice<float> lpsfOld(psfOldSlice);
       casa::ArrayLattice<float> lpsfNew(psfNewSlice);
       
       
       scimath::PaddingUtils::inject(uvOld, lpsfOld);
       scimath::PaddingUtils::inject(uvNew, lpsfNew);
       
       // ratio of FTs is an estimate of the gridded imaging weight. We have to use
       // gridded weight because we don't form ungridded one.
       casa::LatticeFFT::cfft2d(uvOld, casa::True);
       casa::LatticeFFT::cfft2d(uvNew, casa::True);
       
       // The following equation is from Dan Briggs' thesis page 41, eq 3.5
       double sumwtNew = 0.;
       double sumwtOld = 0.;
       double sumwt2Old = 0.;
       double sumwt2New = 0.;

       casa::IPosition cursor(paddedShape.nelements(),0);
       for (int nx = 0; nx<paddedShape(0); ++nx) {
            cursor(0) = nx;
            for (int ny = 0; ny<paddedShape(1); ++ny) {
                 cursor(1)=ny;
                 const double wtOld = casa::abs(uvOld(cursor));
                 ASKAPCHECK(!std::isnan(wtOld), "abs(uvOld) is NaN, uvOld(cursor)="<<uvOld(cursor));	    
                                  
                 const double wtNew = casa::abs(uvNew(cursor));
                 ASKAPCHECK(!std::isnan(wtNew), "abs(uvNew) is NaN, uvOld(cursor)="<<uvNew(cursor));	    
                 sumwtOld += wtOld;
                 sumwtNew += wtNew;
                 sumwt2Old += casa::square(wtOld);
                 sumwt2New += casa::square(wtNew);
            }
       }
       ASKAPCHECK(sumwtNew>0, "Sum of weights is zero in ImageSolver::sensitivityLoss");
       const double loss = sqrt(sumwt2New/sumwt2Old)*sumwtOld/sumwtNew;
       ASKAPLOG_INFO_STR(logger, "The estimate of the sensitivity loss is "<<loss);
       return loss;
    }
    
    /// @brief configure basic parameters of the solver
    /// @details This method encapsulates extraction of basic solver parameters from the parset.
    /// @param[in] parset parset's subset (should have solver.Clean or solver.Dirty removed)
    void ImageSolver::configure(const LOFAR::ParameterSet &parset) {
        setTol(parset.getFloat("tolerance", 0.1));
        setVerbose(parset.getBool("verbose", true));
        setSaveIntermediate(parset.getBool("saveintermediate", true));
        zeroWeightCutoffMask(!parset.getBool("weightcutoff.clean",false));        
        const std::string weightCutoff = parset.getString("weightcutoff","truncate");
        if (weightCutoff == "zero") {
            zeroWeightCutoffArea(true);
            ASKAPLOG_INFO_STR(logger, "Solver is configured to zero pixels in the area where weight is below cutoff (tolerance parameter)");
            ASKAPCHECK(zeroWeightCutoffMask() == true, "With weightcutoff="<<weightCutoff<<
                       " only weightcutoff.clean = false makes sense");
        } else if (weightCutoff == "truncate") {
            ASKAPLOG_INFO_STR(logger, "Solver is configured to truncate pixels in the area where weight is below cutoff (tolerance parameter). These pixels will be normalised by the maximum diagonal");
            zeroWeightCutoffArea(false);
        } else {
            ASKAPTHROW(AskapError, "Only 'zero' and 'truncate' are allowed values for weightcutoff parameter, you have "<<weightCutoff);
        }
        if (zeroWeightCutoffMask()) {
            ASKAPLOG_INFO_STR(logger, "In this area, pixels are masked out, and consequently no S/N-based cleaning will be done");
        } else {
            ASKAPLOG_INFO_STR(logger, "S/N-based clean will search optimum of flux * sqrt(tolerance) in this area");
        }
    }

  } // namespace synthesis
} // namespace askap
