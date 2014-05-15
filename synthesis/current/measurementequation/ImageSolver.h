/// @file
///
/// ImageSolver: This solver calculates the dirty image (or equivalent)
/// for all parameters called image*
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
#ifndef SYNIMAGESOLVER_H_
#define SYNIMAGESOLVER_H_

#include <fitting/Solver.h>
#include <fitting/ImagingNormalEquations.h>
#include <fitting/DesignMatrix.h>
#include <fitting/Params.h>

#include <measurementequation/IImagePreconditioner.h>
#include <boost/shared_ptr.hpp>

#include <Common/ParameterSet.h>


#include <map>
using std::map;

namespace askap
{
  namespace synthesis
  {
    /// @brief Base class for solvers of images
    ///
    /// @details This solver takes the normal equations and simply divides
    /// the data vector by the diagonal of the normal matrix. This
    /// is analogous to making the dirty image or a linear mosaic
    /// of dirty images.
    /// @ingroup measurementequation
    class ImageSolver : public askap::scimath::Solver
    {
      public:
        /// @brief typedef of the shared pointer to ImageSolver
	    typedef boost::shared_ptr<ImageSolver> ShPtr;
	
        /// @brief default constructor
        ImageSolver();

        /// @brief Initialize this solver
        virtual void init();

        /// @brief Solve for parameters, updating the values kept internally
        /// The solution is constructed from the normal equations. The parameters named 
        /// image* are interpreted as images and solved for.
        /// @param[in] ip current model (to be updated)
        /// @param[in] quality Solution quality information
        virtual bool solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality);
        
        /// @brief Clone this object
        virtual askap::scimath::Solver::ShPtr clone() const;
        
        /// @return a reference to normal equations object
        /// @note In this class and derived classes the type returned
        /// by this method is narrowed to always provide image-specific 
        /// normal equations objects
        virtual const scimath::ImagingNormalEquations& normalEquations() const;

	/// @brief Setup the preconditioner
	void addPreconditioner(askap::synthesis::IImagePreconditioner::ShPtr pc);

	/// @brief Do the preconditioning
        /// @details 
        /// @param[in] psf point spread function do modify
        /// @param[in] dirty dirty image to modify
	bool doPreconditioning(casa::Array<float>& psf, casa::Array<float>& dirty) const;
   
	/// @brief perform normalization of the dirty image and psf
	/// @details This method divides the PSF and dirty image by the diagonal of the Hessian.
	/// If a non-void shared pointer is specified for the mask parameter, this method assigns
	/// 0. for those elements where truncation of the weights has been performed and 1. 
	/// otherwise. 
	/// @param[in] diag diagonal of the Hessian (i.e. weights), dirty image will be
	///            divided by an appropriate element of the diagonal or by a cutoff
	///            value
	/// @param[in] tolerance cutoff value given as a fraction of the largest diagonal element
	/// @param[in] psf  point spread function, which is normalized 
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
	float doNormalization(const casa::Vector<double>& diag, 
			     const float& tolerance,
			     casa::Array<float>& psf, 
                             float psfRefPeak,
			     casa::Array<float>& dirty,
			     const boost::shared_ptr<casa::Array<float> >& mask = 
			               boost::shared_ptr<casa::Array<float> >()) const;

	/// @brief perform normalization of the dirty image and psf
	/// @details This is an overloaded version of the method. It also
	/// divides the PSF and dirty image by the diagonal of the Hessian.
	/// However, it assumes that the psf should always be noramlised to 1.
	/// If a non-void shared pointer is specified for the mask parameter, this method assigns
	/// 0. for those elements where truncation of the weights has been performed and 1. 
	/// otherwise. 
	/// @param[in] diag diagonal of the Hessian (i.e. weights), dirty image will be
	///            divided by an appropriate element of the diagonal or by a cutoff
	///            value
	/// @param[in] tolerance cutoff value given as a fraction of the largest diagonal element
	/// @param[in] psf  point spread function, which is normalized 
	/// @param[in] dirty dirty image which is normalized by truncated weights (diagonal)
	/// @param[out] mask shared pointer to the output mask showing where the truncation has 
	///             been performed.
    /// @return peak of PSF before normalisation (to be used as psfRefPeak, if necessary)
	/// @note although mask is filled in inside this method it should already have a correct 
	/// size before this method is called. Pass a void shared pointer (default) to skip 
	/// mask-related functionality. Hint: use utility::NullDeleter to wrap a shared pointer
	/// over an existing array reference.
	inline float doNormalization(const casa::Vector<double>& diag, 
           const float& tolerance, casa::Array<float>& psf, casa::Array<float>& dirty,
           const boost::shared_ptr<casa::Array<float> >& mask = 
			               boost::shared_ptr<casa::Array<float> >()) const
      { return doNormalization(diag,tolerance, psf, -1., dirty, mask); }

    /// @brief configure basic parameters of the solver
    /// @details This method encapsulates extraction of basic solver parameters from the parset.
    /// @param[in] parset parset's subset (should have solver.Clean or solver.Dirty removed)
    virtual void configure(const LOFAR::ParameterSet &parset); 

    /// @brief query weight cutoff behavior
    /// @return true if image pixels corresponding to the weight cutoff area are set to zero during
    /// normalisation.
    inline bool zeroWeightCutoffArea() const { return itsZeroWeightCutoffArea; }
    
    /// @brief query weight cutoff clean mask behavior
    /// @return true if mask is set to zero during normalisation for those pixels which are in the 
    /// weight cutoff area (i.e. not to be cleaned for S/N-based clean)
    inline bool zeroWeightCutoffMask() const { return itsZeroWeightCutoffMask; }

    /// @brief set weight cutoff behavior
    /// @param[in] flag true to set image pixels corresponding to the weight cutoff area to zero during
    /// normalisation.
    inline void zeroWeightCutoffArea(bool flag) { itsZeroWeightCutoffArea = flag; }
      
    /// @brief set weight cutoff clean mask behavior
    /// @param[in] flag true to set mask to zero during normalisation for those pixels which are in the 
    /// weight cutoff area (i.e. to ensure that they are not cleaned during S/N-based clean)
    inline void zeroWeightCutoffMask(bool flag) { itsZeroWeightCutoffMask = flag;}      

    /// @brief set save intermediate flag
      /// @details Force intermediate results to parameters for subsequent saving to disk. This
      /// uses memory and thus probably will be used mostly in debugging
    /// @param[in] flag true to force intermediate results (e.g. residuals, mask) to be saved
      inline void setSaveIntermediate(bool flag) { itsSaveIntermediate=flag;}

    /// @brief get save intermediate flag
      /// @details Force intermediate results to parameters for subsequent saving to disk. This
      /// uses memory and thus probably will be used mostly in debugging
      inline bool saveIntermediate() { return itsSaveIntermediate;}

    /// @brief Save the weights as a parameter
    /// @param[in] ip current model (to be updated)        
    inline void saveWeights(askap::scimath::Params& ip) const 
        { saveNEPartIntoParameter(ip,"weights",normalEquations().normalMatrixDiagonal());}

    /// @brief Save the PSFs as a parameter
    /// @param[in] ip current model (to be updated)        
    inline void savePSF(askap::scimath::Params& ip) const
        { saveNEPartIntoParameter(ip,"psf",normalEquations().normalMatrixSlice());}

  protected:
     
    /// @brief estimate sensitivity loss due to preconditioning
    /// @details Preconditioning (i.e. Wiener filter, tapering) makes the synthesized beam look nice, 
    /// but the price paid is a sensitivity loss. This method gives an estimate (accurate calculations require
    /// gridless weights, which we don't have in our current approach). The method just requires the two
    /// PSFs before and after preconditioning.
    /// @param[in] psfOld an array with original psf prior to preconditioning
    /// @param[in] psfNew an array with the psf after preconditioning has been applied
    /// @return sensitivity loss factor (should be grater than or equal to 1)
    static double sensitivityLoss(const casa::Array<float>& psfOld, const casa::Array<float>& psfNew);
    
    /// @brief a helper method to extract the first plane out of the multi-dimensional array
    /// @details This method just uses MultiDimArrayPlaneIter to extract the first plane
    /// out of the array. It accepts a const reference to the array (which is a conseptual const).
    /// @param[in] in const reference to the input array
    /// @return the array with the first plane
    static casa::Array<float> getFirstPlane(const casa::Array<float> &in);
                
    /// @brief helper method to save part of the NE
    /// @details We need to save slice and diagonal of the normal equations as 
    /// PSF and weights image (savePSF and saveWeights methods) with very similar
    /// operations. This method encapsulates the common code to avoid duplication.
    /// It iterates over all parameters with names starting with "image". 
    /// @param[in] ip model (to be updated with the appropriate parameter) 
    /// @param[in] prefix name prefix for stored parameter (image will be replaced with this prefix, i.e.
    /// "psf" or "weights"
    /// @param[in] nePart part of the normal equations to save (map of vectors)
    void saveNEPartIntoParameter(askap::scimath::Params& ip, const std::string &prefix,
                 const std::map<std::string, casa::Vector<double> > &nePart) const;
    
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
    static void saveArrayIntoParameter(askap::scimath::Params& ip, const std::string &imgName, 
              const casa::IPosition &shape, const std::string &prefix, const casa::Array<double> &arr, 
              const casa::IPosition &pos);
    
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
    static void saveArrayIntoParameter(askap::scimath::Params& ip, const std::string &imgName, 
              const casa::IPosition &shape, const std::string &prefix, const casa::Array<float> &arr, 
              const casa::IPosition &pos);
    
                 
private:
	/// Instance of a preconditioner
	// IImagePreconditioner::ShPtr itsPreconditioner;
        //std::map<std::int, boost::shared_ptr<askap::synthesis::IImagePreconditioner> > itsPreconditioners;
        std::map<int, IImagePreconditioner::ShPtr> itsPreconditioners;

    /// @brief controls weight normalisation
    /// @details If true, area outside weight cutoff area is set to zero. Otherwise, the normalisation is
    /// done by dividing to maximum weight. The default is false;
    bool itsZeroWeightCutoffArea;
      
    /// @brief controls mask used for S/N-based clean
    /// @details If true, the mask in the weight cutoff area is set to zero. This ensures that nothing is cleaned
    /// in those areas for the S/N-based clean. Otherwise, the mask is set to sqrt(tolerance), which corresponds to
    /// normalisation done by dividing to maximum weight. The default is true.
    bool itsZeroWeightCutoffMask; 

      /// @brief Control saving of intermediate images
      /// @details If true, then selected images are saved to parameters and thence to disk
      bool itsSaveIntermediate;

    };

  }
}
#endif
