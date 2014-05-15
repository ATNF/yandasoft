/// @file AWProjectVisGridder.h
///
/// AWProjectVisGridder: Grids visibility data using the self-convolution of 
/// the antenna illumination pattern.
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
#ifndef ASKAP_SYNTHESIS_AWPROJECTVISGRIDDER_H_
#define ASKAP_SYNTHESIS_AWPROJECTVISGRIDDER_H_

// ASKAPsoft includes
#include <dataaccess/IConstDataAccessor.h>
#include <boost/shared_ptr.hpp>

// Local package includes
#include <gridding/WProjectVisGridder.h>
#include <gridding/IBasicIllumination.h>
#include <gridding/AProjectGridderBase.h>

namespace askap
{
    namespace synthesis
    {
        /// @brief Gridder that is appropriate for mosaicing. 
        ///
        /// @details The visibilities are gridded using a convolution
        /// function derived from the antenna illumination pattern,
        /// appropriately shifted in position for each feed, and
        /// incorporating the Fresnel term needed to correct for the
        /// w-term in the full measurement equation.
        ///
        /// The scaling is slow in data points, slow in w planes 
        /// (since the calculation of the convolution function
        /// usually dominates).
        ///
        /// @ingroup gridding
        class AWProjectVisGridder : public WProjectVisGridder, virtual protected AProjectGridderBase
        {
            public:

                /// @brief Construct antenna illumination pattern/W term gridder
                /// @param illum  Antenna illumination model
                /// @param wmax Maximum baseline (wavelengths)
                /// @param nwplanes Number of w planes
                /// @param cutoff Cutoff in determining support e.g. 10^-3 of the peak
                /// @param overSample Oversampling (currently limited to <=1)
                /// @param maxSupport Maximum support to be allowed
                /// @param limitSupport Upper limit of support
                /// @param maxFeeds Maximum number of feeds allowed
                /// @param maxFields Maximum number of fields allowed
                /// @param pointingTol Pointing tolerance in radians
                /// @param paTol Parallactic angle tolerance in radians
                /// @param freqTol Frequency tolerance (relative, threshold for df/f), negative value 
                ///        means the frequency axis is ignored       
                /// @param frequencyDependent Frequency dependent gridding?
                /// @param name Name of table to save convolution function into
                AWProjectVisGridder(const boost::shared_ptr<IBasicIllumination const> &illum,
                        const double wmax, const int nwplanes, const double cutoff,
                        const int overSample, const int maxSupport, const int limitSupport,
                        const int maxFeeds=1, const int maxFields=1, const double pointingTol=0.0001,
                        const double paTol=0.01,
                        const double freqTol = 1e-6,          
                        const bool frequencyDependent=true, 
                        const std::string& name=std::string(""));

                /// @brief copy constructor
                /// @details It is required to decouple internal array arrays, otherwise
                /// those arrays are shared between all cloned gridders of this type
                /// @param[in] other input object
                /// @note illumination model is copied as a pointer, so the same model is referenced
                AWProjectVisGridder(const AWProjectVisGridder &other);

                /// Clone a copy of this Gridder
                virtual IVisGridder::ShPtr clone();

                /// Form the sum of the convolution function squared, multiplied by the weights for each
                /// different convolution function. This is used in the evaluation of the second derivative.
                /// @param out Output double precision grid
                virtual void finaliseWeights(casa::Array<double>& out);

                /*
                /// Form the final output image or PSF
                /// @param out Output double precision image or PSF
                virtual void finaliseGrid(casa::Array<double>& out);
                */

                /// @brief Initialise the gridding
                /// @param axes axes specifications
                /// @param shape Shape of output image: u,v,pol,chan
                /// @param dopsf Make the psf?
                virtual void initialiseGrid(const scimath::Axes& axes,
                        const casa::IPosition& shape, const bool dopsf=true);

                /// @brief Initialise the degridding
                /// @param axes axes specifications
                /// @param image Input image: cube: u,v,pol,chan
                virtual void initialiseDegrid(const scimath::Axes& axes,
                        const casa::Array<double>& image);

                /// @brief static method to get the name of the gridder
                /// @details We specify parameters per gridder type in the parset file.
                /// This method returns the gridder name which should be used to extract
                /// a subset of parameters for createGridder method.
                static inline std::string gridderName() { return "AWProject";}				

                /// @brief static method to create gridder
                /// @details Each gridder should have a static factory method, which is
                /// able to create a particular type of the gridder and initialise it with
                /// the parameters taken form the given parset. It is assumed that the 
                /// method receives a subset of parameters where the gridder name is already
                /// taken out. 
                /// @param[in] parset input parset file
                /// @return a shared pointer to the gridder instance					 
                static IVisGridder::ShPtr createGridder(const LOFAR::ParameterSet& parset);

            protected:
                /// @brief initialise sum of weights
                /// @details We keep track the number of times each convolution function is used per
                /// channel and polarisation (sum of weights). This method is made virtual to be able
                /// to do gridder specific initialisation without overriding initialiseGrid.
                /// This method accepts no parameters as itsShape, itsNWPlanes, etc should have already
                /// been initialised by the time this method is called.
                virtual void initialiseSumOfWeights();

                /// @brief Initialise the indices
                /// @param[in] acc const accessor to work with
                virtual void initIndices(const accessors::IConstDataAccessor& acc);

                /// Index into convolution function
                /// @param row Row number
                /// @param pol Polarization id
                /// @param chan Channel number
                virtual int cIndex(int row, int pol, int chan);

                /// Initialize convolution function
                /// @param[in] acc const accessor to work with
                virtual void initConvolutionFunction(const accessors::IConstDataAccessor& acc);

                /// Correct for gridding convolution function
                /// @param image image to be corrected
                virtual void correctConvolution(casa::Array<double>& image);

            private:
                /// @brief assignment operator (not to be called)
                /// @details It is made private, so we can't call it inadvertently
                /// @param[in] other input object
                AWProjectVisGridder& operator=(const AWProjectVisGridder &other);

                /// Reference frequency for illumination pattern. 
                double itsReferenceFrequency;

                /// Antenna illumination model
                boost::shared_ptr<IBasicIllumination const> itsIllumination;

                /// Is the convolution function frequency dependent?
                bool itsFreqDep;

                /// Maximum number of feeds
                int itsMaxFeeds;

                /// Maximum number of fields
                int itsMaxFields;

                /// Cube of slopes
                casa::Cube<double> itsSlopes;
        };

    } // end namespace synthesis
} // end namespace askap
#endif
