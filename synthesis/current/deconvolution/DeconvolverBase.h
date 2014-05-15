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
#include <Common/ParameterSet.h>
#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>

#include <deconvolution/DeconvolverState.h>
#include <deconvolution/DeconvolverControl.h>
#include <deconvolution/DeconvolverMonitor.h>

namespace askap {

    namespace synthesis {

        /// @brief Base class for a deconvolver
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a dirty image, psf optionally using a mask and a weights image.
        /// The template argument T is the type, and FT is the transform
        /// e.g. Deconvolver<Double, DComplex>
        /// The interface is by casa::Array<T>'s holding the various arrays
        /// Usually the arrays are 2-D. However, in the case of e.g. MSMFS the
        /// third axis will be the taylor terms.
        /// @ingroup Deconvolver
        template<class T, class FT> class DeconvolverBase {

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
                DeconvolverBase(casa::Vector<casa::Array<T> >& dirty, casa::Vector<casa::Array<T> >& psf);

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We keep this
                /// version for compatibility
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverBase(casa::Array<T>& dirty, casa::Array<T>& psf);

                /// @brief Get the current dirty image
                /// @detail Get the current dirty image
                casa::Array<T>& dirty(const casa::uInt term = 0);

                /// @brief Get the current PSF
                /// @detail Get the current PSF
                casa::Array<T>& psf(const casa::uInt term = 0);

                /// @brief Set the initial model
                /// @detail Set the model from which iteration will start
                void setModel(const casa::Array<T> model, const casa::uInt term = 0);

                /// @brief Get the current model
                /// @detail Get the current model
                /// @param[out] model Model image (array)
                casa::Array<T>& model(const casa::uInt term = 0);

                /// @brief Update only the dirty image
                /// @detail Update an existing deconvolver for a changed dirty image
                /// @param[in] dirty Dirty image (array)
                /// @param[in] term term to update
                virtual void updateDirty(casa::Array<T>& dirty, const casa::uInt term = 0);

                /// @brief Update only the dirty images
                /// @detail Update an existing deconvolver for a changed dirty images.
                /// @param[in] dirty Dirty image (vector of arrays)
                virtual void updateDirty(casa::Vector<casa::Array<T> >& dirty);

                /// @brief Set the weight image
                /// @detail The weights image (actually the sqrt) is used to
                /// aid the deconvolution. The weights image is proportional
                /// to 1/sigma**2
                /// @param[in] weights Weights (array)
                void setWeight(casa::Array<T> weight, const casa::uInt term = 0);

                /// @brief Get the weight image
                /// @detail Get the weight
                /// @param[out] weight (array)
                casa::Array<T> & weight(const casa::uInt term = 0);

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
                virtual void updateResiduals(casa::Vector<casa::Array<T> >& model);

                /// @brief Update the residuals: keep for compatibility
                /// @detail Update the residuals for this model.
                /// This usually requires convolution of the model with
                /// the specified PSF and subtraction from the dirty image.
                virtual void updateResiduals(casa::Array<T>& model);

                /// @brief Restore with specified beam
                /// @detail Restore the model by smoothing
                /// and adding residuals
                virtual bool restore(casa::Vector<casa::Array<T> >& restored, casa::Vector<casa::Array<T> >& model);

                /// @brief Restore with specified beam
                /// @detail Restore the model by smoothing
                /// and adding residuals
                virtual bool restore(casa::Vector<casa::Array<T> >& restored);

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

            protected:

                // Number of terms in the expansion > 0
                casa::uInt itsNumberTerms;

                // Initialise for both constructors
                void init(casa::Vector<casa::Array<T> >& dirty, casa::Vector<casa::Array<T> >& psf);

                // Validate the various shapes to ensure consistency
                void validateShapes();

                casa::Vector<casa::Array<T> > itsDirty;

                casa::Vector<casa::Array<T> > itsPsf;

                casa::Vector<casa::Array<T> > itsModel;

                casa::Vector<casa::Array<T> > itsWeight;

                /// The state of the deconvolver
                boost::shared_ptr<DeconvolverState<T> > itsDS;

                /// The control used for the deconvolver
                boost::shared_ptr<DeconvolverControl<T> > itsDC;

                // Find the shape of the PSF to be used, this includes
                // the effects of psfwidth
                IPosition findSubPsfShape();

                /// The monitor used for the deconvolver
                boost::shared_ptr<DeconvolverMonitor<T> > itsDM;

                // Peak and location of peak of PSF(0)
                casa::IPosition itsPeakPSFPos;
                T itsPeakPSFVal;

                // Beam information, in pixels
                float itsBMaj;
                float itsBMin;
                float itsBPa;

                // Audit the memory in use right now
                void auditAllMemory();
                casa::uInt auditMemory(casa::Vector<casa::Array<T> >& vecArray);
                casa::uInt auditMemory(casa::Vector<casa::Array<FT> >& vecArray);
        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/DeconvolverBase.tcc>

#endif  // #ifndef I_DECONVOLVERBASE_H


