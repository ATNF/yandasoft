/// @file
/// 
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


// ImplCalWeightSolver -  implementation of the algorithm which 
// solves for the best FPA weights for an optimum calibration on a
// given sky model

#include <casa/aips.h>
#include <components/ComponentModels/SkyComponent.h>
#include <components/ComponentModels/ComponentList.h>
#include <measures/Measures/MDirection.h>
#include <casa/Exceptions/Error.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/String.h>
#include <images/Images/PagedImage.h>
#include <casa/BasicSL/Complex.h>


#include <boost/shared_ptr.hpp>
#include <gridding/IBasicIllumination.h>

class ImplCalWeightSolver {
    casa::MDirection pc; // dish pointing centre
    casa::ComponentList cl; // model of sky brightness
    mutable casa::Matrix<casa::Complex> vismatrix; // visibilities for each element (column)
                                           // and measurement (row)
    casa::ImageInterface<casa::Float> *vp_real; // voltage pattern of a single element
    casa::ImageInterface<casa::Float> *vp_imag; // voltage pattern of a single element
       
public:
    //static const double lambda=0.2;   // wavelength in metres
    ImplCalWeightSolver() throw();
    ~ImplCalWeightSolver() throw(casa::AipsError);
    // set up calculation for a given pointing centre and sky model
    void setSky(const casa::MDirection &ipc, 
		const casa::String &clname) throw(casa::AipsError);
		
	/// @brief make synthetic beam
	/// @details This method constructs synthetic primary beam for the given weights.
	/// @param[in] name output image name
	/// @param[in] weights vector of weights
	void makeSyntheticPB(const std::string &name, 
	                     const casa::Vector<casa::Complex> &weights);
		
    // set up the voltage pattern from a disk-based image
    void setVP(const casa::String &namer,const casa::String &namei) 
	      throw(casa::AipsError); 
    // main method
    casa::Matrix<casa::Complex>
         solveWeights(const casa::Matrix<casa::Double> &feed_offsets,
		      const casa::Vector<casa::Double> &uvw) const
	                            throw(casa::AipsError);
   // solve for eigenvectors for the VP matrix. The first vector
   // (column=0) corresponds to the maximum attainable flux,
   // the last one (column=Nfeeds-1) corresponds to the
   // weight set for an optimal rejection of all known sources.
   // pa - parallactic angle to rotate all source offsets (in radians)
   // if skycat!="", a table with this name will be filled with the offsets
   // w.r.t. the dish pointing centre
   casa::Matrix<casa::Complex>
         eigenWeights(const casa::Matrix<casa::Double> &feed_offsets,
	 casa::Double pa = 0., const casa::String &skycat="")
	                        const throw(casa::AipsError);

   // solve for a basis in the space of weights, which is best for calibration
   // in the sense that the flux from known sources is maximized
   // The result matrix contains the basis vectors in its columns
   // ndim is a number of basis vectors required (should be less than or
   // equal to the number of feeds.
   // pa - parallactic angle to rotate all source offsets (in radians)
   // if skycat!="", a table with this name will be filled with the offsets
   // w.r.t. the dish pointing centre
   casa::Matrix<casa::Complex>
         calBasis(const casa::Matrix<casa::Double> &feed_offsets,
	          casa::uInt ndim, casa::Double pa=0.,
		  const casa::String &skycat="")
		                      const throw(casa::AipsError);
				      
protected:
    // calculate visibility matrix for given feed_offsets
    // uvw - a vector with the uvw coordinates (in the units of wavelength)
    // vismatrix will have visibilities for each element (column) and 
    // measurement (row)
    void formVisMatrix(const casa::Matrix<casa::Double> &feed_offsets,
		       const casa::Vector<casa::Double> &uvw) const
	                throw(casa::AipsError);

    // this version fills vismatrix with sum(F_i E_k*E_l^H), where 
    // E is the voltage pattern value at the source position and
    // F_i is the flux of the ith source
    // pa - parallactic angle to rotate all source offsets (in radians)
    // if skycat!="", a table with this name will be filled with the offsets
    // w.r.t. the dish pointing centre
    void formVPMatrix(const casa::Matrix<casa::Double> &feed_offsets,
                casa::Double pa = 0., const casa::String &skycat = "")
	          const  throw(casa::AipsError);
    
    
    // an auxiliary function to extract a Value of the Voltage Pattern
    // at the given offset (in radians). Return True if successful, and
    // False if the requested offset lies outside the model
    casa::Bool getVPValue(casa::Complex &val, casa::Double l, casa::Double m)
                            const throw(casa::AipsError);
private:
    boost::shared_ptr<askap::synthesis::IBasicIllumination> itsIllumination;
};
