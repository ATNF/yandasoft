/// @file
///
/// ImageMultiScaleSolver: This solver calculates the dirty image (or equivalent)
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
#ifndef SYNIMAGEMULTISCALESOLVER_H_
#define SYNIMAGEMULTISCALESOLVER_H_

#include <measurementequation/ImageCleaningSolver.h>

#include <lattices/Lattices/LatticeCleaner.h>
#include <utils/FixedSizeCache.h>

#include <map>

namespace askap {
    namespace synthesis {
        /// @brief Multiscale solver for images.
        ///
        /// @details This solver performs multi-scale clean using the
        /// casa::LatticeCleaner classes
        ///
        /// @ingroup measurementequation
        class ImageMultiScaleSolver : public ImageCleaningSolver {
            public:

                /// @brief default constructor
                /// The default scales are 0, 10, 30 pixels
                ImageMultiScaleSolver();

                /// @brief Constructor from scales.
                /// @param scales Scales to be solved in pixels
                ImageMultiScaleSolver(const casa::Vector<float>& scales);

                /// @brief Initialize this solver
                virtual void init();
 
                /// @brief Solve for parameters, updating the values kept internally
                /// The solution is constructed from the normal equations. The parameters named 
                /// image* are interpreted as images and solved for.
                /// @param[in] ip current model (to be updated)
                /// @param[in] quality Solution quality information
                virtual bool solveNormalEquations(askap::scimath::Params& ip, askap::scimath::Quality& quality);

                /// @brief Clone this object
                /// @return a shared pointer to the clone
                virtual askap::scimath::Solver::ShPtr clone() const;

                /// @brief Set the scales
                /// @param[in] scales vector with scales
                void setScales(const casa::Vector<float>& scales);
                
                /// @brief switch the speed up on
                /// @param[in] factor speed up factor
                void setSpeedUp(float factor);

            protected:
                /// @brief Precondition the PSF and the dirty image
                /// @param[in] psf point spread function to precondition (in/out)
                /// @param[in] dirty dirty image to precondition (in/out)
                void preconditionNE(casa::ArrayLattice<float>& psf, casa::ArrayLattice<float>& dirty);

                /// Scales in pixels
                casa::Vector<float> itsScales;

                /// @brief Cache of Cleaners
                scimath::FixedSizeCache<string, casa::LatticeCleaner<float> > itsCleaners;
            
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
