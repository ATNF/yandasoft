/// @file DeconvolverBase.h
/// @brief Base class for a deconvolver
/// @details This interface class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
/// @ingroup Deconvolver
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///

#ifndef ASKAP_SYNTHESIS_DECONVOLVERBASE_H
#define ASKAP_SYNTHESIS_DECONVOLVERBASE_H

#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

#include <Common/ParameterSet.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>

#include <askap/deconvolution/DeconvolverState.h>
#include <askap/deconvolution/DeconvolverControl.h>
#include <askap/deconvolution/DeconvolverMonitor.h>

namespace askap {

    namespace synthesis {

        /// @brief Base class for a deconvolver
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a dirty image, psf optionally using a mask and a weights image.
        /// The template argument T is the type, and FT is the transform
        /// e.g. Deconvolver<Double, DComplex>
        /// The interface is by casacore::Array<T>'s holding the various arrays
        /// Usually the arrays are 2-D. However, in the case of e.g. MSMFS the
        /// third axis will be the taylor terms.
        /// @ingroup Deconvolver
        template<class T, class FT> class DeconvolverBase : public boost::noncopyable {

            public:
                typedef boost::shared_ptr<DeconvolverBase<T, FT> > ShPtr;

                virtual ~DeconvolverBase();

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We use a vector
                /// here to allow multiple dirty images and PSFs for the
                /// same model (e.g. as in MFS)
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverBase(casacore::Vector<casacore::Array<T> >& dirty, casacore::Vector<casacore::Array<T> >& psf);

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We keep this
                /// version for compatibility
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverBase(casacore::Array<T>& dirty, casacore::Array<T>& psf);

                /// @brief Get the current dirty image
                /// @detail Get the current dirty image
                casacore::Array<T>& dirty(const casacore::uInt term = 0);

                /// @brief Get the current PSF
                /// @detail Get the current PSF
                casacore::Array<T>& psf(const casacore::uInt term = 0);

                /// @brief Set the initial model
                /// @detail Set the model from which iteration will start
                void setModel(const casacore::Array<T>& model, const casacore::uInt term = 0);

                /// @brief Get the current model
                /// @detail Get the current model
                /// @param[in] term term of interest (i.e. element of the vector)
                /// @return const reference to the array with the model
                const casacore::Array<T>& model(const casacore::uInt term = 0) const;

                /// @brief Non-const version to obtain reference to the current model
                /// @detail Get the current model (for potential update)
                /// @param[in] term term of interest (i.e. element of the vector)
                /// @return reference to the array with the model
                casacore::Array<T>& model(const casacore::uInt term = 0);

                /// @brief Update only the dirty image
                /// @detail Update an existing deconvolver for a changed dirty image
                /// @param[in] newDirty Dirty image (array)
                /// @param[in] term term to update
                virtual void updateDirty(const casacore::Array<T>& newDirty, const casacore::uInt term = 0);

                /// @brief Update only the dirty images
                /// @detail Update an existing deconvolver for a changed dirty images.
                /// @param[in] newDirty Dirty image (vector of arrays)
                virtual void updateDirty(const casacore::Vector<casacore::Array<T> >& newDirty);

                /// @brief Set the weight image
                /// @detail The weights image (actually the sqrt) is used to
                /// aid the deconvolution. The weights image is proportional
                /// to 1/sigma**2
                /// @param[in] weights Weights (array)
                void setWeight(casacore::Array<T> weight, const casacore::uInt term = 0);

                /// @brief Get the weight image
                /// @detail Get the weight
                /// @param[out] weight (array)
                casacore::Array<T> & weight(const casacore::uInt term = 0);

                /// @brief Set the control
                /// @detail The control is used to hold all the information
                /// required to control the algorithm. All decisions
                /// regarding e.g. stopping are delegated to this class.
                /// @param[in] state Shared pointer to the control
                bool setControl(boost::shared_ptr<DeconvolverControl<T> > control);
                boost::shared_ptr<DeconvolverControl<T> > control() const;

                /// @brief Set the monitor
                /// @detail The monitor is used to monitor the algorithm.
                /// All standard monitoring is performed by this class.
                /// @param[in] state Shared pointer to the monitor
                bool setMonitor(boost::shared_ptr<DeconvolverMonitor<T> > monitor);
                boost::shared_ptr<DeconvolverMonitor<T> > monitor() const;


                /// @brief Set the state
                /// @detail The state class is used to communicate to the
                /// monitor class or to other classes.
                /// @param[in] state Shared pointer to the state
                bool setState(boost::shared_ptr<DeconvolverState<T> > state);
                boost::shared_ptr<DeconvolverState<T> > state() const;

                // @brief Perform the deconvolution
                // @detail This is the main deconvolution method.
                virtual bool deconvolve();

                /// @brief Update the residuals
                /// @detail Update the residuals for this model.
                /// This usually requires convolution of the model with
                /// the specified PSF and subtraction from the dirty image.
                virtual void updateResiduals(casacore::Vector<casacore::Array<T> >& model);

                /// @brief Update the residuals: keep for compatibility
                /// @detail Update the residuals for this model.
                /// This usually requires convolution of the model with
                /// the specified PSF and subtraction from the dirty image.
                virtual void updateResiduals(casacore::Array<T>& model);

                /// @brief Restore with specified beam
                /// @detail Restore the model by smoothing
                /// and adding residuals
                virtual bool restore(casacore::Vector<casacore::Array<T> >& restored, casacore::Vector<casacore::Array<T> >& model);

                /// @brief Restore with specified beam
                /// @detail Restore the model by smoothing
                /// and adding residuals
                virtual bool restore(casacore::Vector<casacore::Array<T> >& restored);

                /// @brief Initialize the deconvolution
                /// @detail Initialise e.g. set weights
                virtual void initialise();

                /// @brief Finalise the deconvolution
                /// @detail Finalise any calculations needed at the end
                /// of iteration
                virtual void finalise();

                /// @brief configure basic parameters of the solver
                /// @details This method encapsulates extraction of basic solver parameters from the parset.
                /// @param[in] parset parset
                virtual void configure(const LOFAR::ParameterSet &parset);

                /// @brief obtain the number of terms
                /// @return the number of terms to solve for
                inline casacore::uInt nTerms() const { return itsNumberTerms;}

                /// @brief release working memory not needed between major cycles
                /// @details Deconvolvers can use a lot of memory, try to release as
                /// much as possible without losing essential state
                virtual void releaseMemory();

            private:

                // Number of terms in the expansion > 0
                casacore::uInt itsNumberTerms;

            protected:
                /// @brief obtain peak psf position (for the 0th term if there are many)
                /// @return peak position
                inline casacore::IPosition getPeakPSFPosition() const {return itsPeakPSFPos;}

                // Initialise for both constructors
                void init(casacore::Vector<casacore::Array<T> >& dirty, casacore::Vector<casacore::Array<T> >& psf);

                // Validate the various shapes to ensure consistency
                void validateShapes();

                /// @brief validate PSF, find peak value and position
                /// @details It works with the zero-th term, if there are many
                /// @param[in] slicer optional slicer if only a fraction of the PSF needs to be considered
                /// the default constructed instance of a slicer results in the whole PSF being used.
                /// @note this method updates itsPeakPSFPos and itsPeakPSFVal, this is why it's non-const
                void validatePSF(const casacore::Slicer &slicer = casacore::Slicer());

                casacore::Vector<casacore::Array<T> > itsDirty;

                casacore::Vector<casacore::Array<T> > itsPsf;

                casacore::Vector<casacore::Array<T> > itsModel;

                casacore::Vector<casacore::Array<T> > itsWeight;

                /// The state of the deconvolver
                boost::shared_ptr<DeconvolverState<T> > itsDS;

                /// The control used for the deconvolver
                boost::shared_ptr<DeconvolverControl<T> > itsDC;

                /// @details Find the shape of the PSF to be used, this includes
                /// the effects of psfwidth
                /// @param[in] term term of interest (i.e. element of the vector)
                IPosition findSubPsfShape(const casacore::uInt term = 0) const;

                /// The monitor used for the deconvolver
                boost::shared_ptr<DeconvolverMonitor<T> > itsDM;

            private:
                // Peak and location of peak of PSF(0)
                casacore::IPosition itsPeakPSFPos;
                T itsPeakPSFVal;

                // Beam information, in pixels
                float itsBMaj;
                float itsBMin;
                float itsBPa;

                // Audit the memory in use right now
                void auditAllMemory();

                /// @brief memory occupied by the array set
                /// @details This method iterates over the supplied array set and computes the
                /// total amount of memory used up.
                /// @param[in] vecArray arrays to work with (packed into a casa Vector)
                /// @return amount of memory in bytes
                template<typename Y>
                casacore::uInt auditMemory(casacore::Vector<casacore::Array<Y> >& vecArray);
        };

    } // namespace synthesis

} // namespace askap

#include <askap/deconvolution/DeconvolverBase.tcc>

#endif  // #ifndef I_DECONVOLVERBASE_H
