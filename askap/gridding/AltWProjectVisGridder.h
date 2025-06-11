/// @file AltWProjectVisGridder.h
///
/// AltWProjectVisGridder: W projection gridding - but with alternate finaliseGrid

/// @copyright (c) 2007,2016,2018 CSIRO
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
/// @author Stephen Ord <stephen.ord@csiro.au>
#ifndef ASKAP_SYNTHESIS_ALTWPROJECTVISGRIDDER_H_
#define ASKAP_SYNTHESIS_ALTWPROJECTVISGRIDDER_H_

// ASKAPsoft includes
#include <askap/gridding/WDependentGridderBase.h>

// Local package includes
//#include <askap/dataaccess/IConstDataAccessor.h>

namespace askap
{
    namespace synthesis
    {
        /// @brief Visibility gridder using W projection
        /// @details The visibilities are gridded using a convolution
        /// function that implements a Fresnel transform. This corrects
        /// for the w term in the full synthesis measurement equation.
        ///
        /// The convolution function is calculated straightforwardly
        /// by constructing an image of the complex w-dependent
        /// phasor and Fourier transforming. The calculation is
        /// done using a coarse but large grid in image space so
        /// that it is sub-sampled in uv space.
        ///
        /// The scaling is slow in data points, fast in w planes.
        ///
        /// @ingroup gridding
        class AltWProjectVisGridder : public WProjectVisGridder
        {
            public:
                /// @brief Construct a gridder for W projection
                /// @param wmax Maximum baseline (wavelengths)
                /// @param nwplanes Number of w planes
                /// @param cutoff Cutoff in determining support e.g. 10^-3 of the peak
                /// @param overSample Oversampling (currently limited to <=1)
                /// @param maxSupport Maximum support to be allowed
                /// @param limitSupport Upper limit of support
                /// @param name Name of table to save convolution function into
                //
                AltWProjectVisGridder(const double wmax, const int nwplanes,
                        const double cutoff, const int overSample,
                        const int maxSupport, const int limitSupport,
                        const std::string& name=std::string(""),
                        const float alpha=1.,const bool writeOut=false);

                virtual ~AltWProjectVisGridder();


                /// Clone a copy of this Gridder
                virtual IVisGridder::ShPtr clone();

                /// @brief static method to get the name of the gridder
                /// @details We specify parameters per gridder type in the parset file.
                /// This method returns the gridder name which should be used to extract
                /// a subset of parameters for createGridder method.
                static inline std::string gridderName() { return "AltWProject";}

                /// @brief static method to create gridder
                /// @details Each gridder should have a static factory method, which is
                /// able to create a particular type of the gridder and initialise it with
                /// the parameters taken form the given parset. It is assumed that the
                /// method receives a subset of parameters where the gridder name is already
                /// taken out.
                /// @param[in] parset input parset file
                /// @return a shared pointer to the gridder instance
                static IVisGridder::ShPtr createGridder(const LOFAR::ParameterSet& parset);

                // This should override the default implementation
                // I thought I had already done this.

                virtual void finaliseGrid(casacore::Array<imtype>& out);


            private:
                /// @brief assignment operator
                /// @details Defined as private, so it can't be called (to enforce usage of the
                /// copy constructor
                /// @param[in] other input object
                /// @return reference to itself
                AltWProjectVisGridder& operator=(const AltWProjectVisGridder &other);

                /// dump the Grids
                bool itsWriteOut;


        };
    }
}
#endif
