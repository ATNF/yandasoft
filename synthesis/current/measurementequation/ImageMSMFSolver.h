/// @file
///
/// ImageMSMFSSolver: This solver does Multi Scale Multi Frequency deconvolution
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
#ifndef SYNIMAGEMSMFSSOLVER_H_
#define SYNIMAGEMSMFSSOLVER_H_

#include <measurementequation/ImageCleaningSolver.h>
#include <askap/AskapUtil.h>


#include <boost/shared_ptr.hpp>
#include <lattices/Lattices/LatticeCleaner.h>
#include <lattices/Lattices/MultiTermLatticeCleaner.h>

#include <map>

namespace askap
{
  namespace synthesis
  {
    /// @brief Multiscale solver for images.
    ///
    /// @details This solver performs multi-scale clean using the 
    /// casa::LatticeCleaner classes
    ///
    /// @ingroup measurementequation
    class ImageMSMFSolver : public ImageCleaningSolver
    {
      public:

        /// @brief default constructor
        ImageMSMFSolver();

        /// @brief Constructor from scales
        /// @param[in] scales Scales to be solved in pixels
        explicit ImageMSMFSolver(const casa::Vector<float>& scales);
          
        /// @brief Initialize this solver
        virtual void init();

	/// @brief Precondition the normal equations prior to solving them
	/// TODO Send in a parset to specify params for diff kinds of preconditioning.
	///virtual void preconditionNormalEquations();

        /// @brief Solve for parameters
        /// The solution is constructed from the normal equations
        /// The solution is constructed from the normal equations. The parameters named 
        /// image* are interpreted as images and solved for.
        /// @param[in] ip current model (to be updated)        
        /// @param[in] quality Solution quality information
        virtual bool solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality);
        
/// @brief Clone this object
        virtual askap::scimath::Solver::ShPtr clone() const;
        
        /// Set the scales
        void setScales(const casa::Vector<float>& scales);
               
      /// @brief switch the speed up on
      /// @param[in] factor speed up factor
      void setSpeedUp(float factor);

      protected:
      
	/// Parse the image parameter name to get out the Stokes information
        inline std::string getStokes(const std::string& paramstring){return paramstring.substr(6,1);};
	/// Parse the image parameter name to get out the Order information
	/// Order == enumeration index (e.x. second-order Taylor coefficient)
        inline int getOrder(const std::string& paramstring){return atoi(paramstring.substr(8,1).data());};
        /// Create an image parameter string, from stokes and order parameters
	inline std::string makeImageString(const std::string& samplestring, const std::string& stokes, const int& order){std::string newstring(samplestring);newstring.replace(6,1,stokes);newstring.replace(8,1,(utility::toString<int>(order)).data());return newstring;};
	

        /// Scales in pixels
        casa::Vector<float> itsScales;

	/// Map of Cleaners - one for each polarisation index
        std::map<std::string, boost::shared_ptr<casa::MultiTermLatticeCleaner<float> > > itsCleaners;

    private:
      /// @brief if true, use speed up factor (default is false)
      bool itsDoSpeedUp;
      
      /// @brief speed up factor
      /// @details This value is used only if itsDoSpeedUp is true
      float itsSpeedUpFactor;
    };

  }
}
#endif
