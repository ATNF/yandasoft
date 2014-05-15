/// @file DeconvolverBase.tcc
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

#include <string>

#include <casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <fft/FFTWrapper.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(decbaselogger, ".deconvolution.base");

#include <deconvolution/DeconvolverBase.h>
#include <deconvolution/DeconvolverState.h>

namespace askap {

    namespace synthesis {

        template<class T, class FT>
        DeconvolverBase<T, FT>::~DeconvolverBase()
        {
            auditAllMemory();
        };

        template<class T, class FT>
        DeconvolverBase<T, FT>::DeconvolverBase(Vector<Array<T> >& dirty, Vector<Array<T> >& psf) :
                itsBMaj(0.0), itsBMin(0.0), itsBPa(0.0)
        {
            init(dirty, psf);
        }

        template<class T, class FT>
        DeconvolverBase<T, FT>::DeconvolverBase(Array<T>& dirty, Array<T>& psf) :
                itsBMaj(0.0), itsBMin(0.0), itsBPa(0.0)
        {
            Vector<Array<T> > dirtyVec(1);
            dirtyVec(0) = dirty.nonDegenerate();
            Vector<Array<T> > psfVec(1);
            psfVec(0) = psf.nonDegenerate();
            init(dirtyVec, psfVec);
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::init(Vector<Array<T> >& dirtyVec, Vector<Array<T> >& psfVec)
        {

            ASKAPCHECK(psfVec.nelements() == dirtyVec.nelements(),
                       "Vectors of dirty images and PSF's not same length");

            itsNumberTerms = dirtyVec.nelements();

            itsDirty.resize(itsNumberTerms);
            itsPsf.resize(itsNumberTerms);
            itsModel.resize(itsNumberTerms);
            itsWeight.resize(itsNumberTerms);

            ASKAPLOG_INFO_STR(decbaselogger, "There are " << itsNumberTerms << " dirty images");

            for (uInt term = 0; term < itsNumberTerms; term++) {

                ASKAPASSERT(dirtyVec(term).nonDegenerate().shape().nelements() == 2);
                ASKAPASSERT(psfVec(term).nonDegenerate().shape().nelements() == 2);

                this->itsDirty(term) = dirtyVec(term).nonDegenerate().copy();
                this->itsPsf(term) = psfVec(term).nonDegenerate().copy();

                ASKAPASSERT(this->itsPsf(term).shape().conform(this->itsDirty(term).shape()));

                ASKAPLOG_INFO_STR(decbaselogger, "Dirty image(" << term << ") has shape: "
                                      << this->dirty(term).shape());

                this->model(term).resize(this->dirty(term).shape());
                this->model(term).set(T(0.0));
            }

            casa::IPosition minPos;
            casa::IPosition maxPos;
            T minVal, maxVal;
            ASKAPLOG_INFO_STR(decbaselogger, "Validating PSF");
            casa::minMax(minVal, maxVal, minPos, maxPos, this->psf(0));

            const Int nx(this->psf(0).shape()(0));
            const Int ny(this->psf(0).shape()(1));

            ASKAPLOG_INFO_STR(decbaselogger, "Maximum of PSF(0) = " << maxVal << " at " << maxPos);
            ASKAPLOG_INFO_STR(decbaselogger, "Minimum of PSF(0) = " << minVal << " at " << minPos);

            if ((maxPos(0) != nx / 2) || (maxPos(1) != ny / 2)) {
                ASKAPTHROW(AskapError, "Peak of PSF(0) is at " << maxPos << ": not at centre pixel: [" << nx / 2 << "," << ny / 2 << "]");
            }

            this->itsPeakPSFVal = maxVal;
            this->itsPeakPSFPos = maxPos;

            itsDS = boost::shared_ptr<DeconvolverState<T> >(new DeconvolverState<T>());
            ASKAPASSERT(itsDS);
            itsDC = boost::shared_ptr<DeconvolverControl<T> >(new DeconvolverControl<T>());
            ASKAPASSERT(itsDC);
            itsDM = boost::shared_ptr<DeconvolverMonitor<T> >(new DeconvolverMonitor<T>());
            ASKAPASSERT(itsDM);

            this->validateShapes();

            auditAllMemory();

        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::configure(const LOFAR::ParameterSet& parset)
        {
            this->itsDC->configure(parset);
            this->itsDM->configure(parset);

            // Get the beam information
            const casa::Vector<float> beam = parset.getFloatVector("beam");
            ASKAPCHECK(beam.size() == 3, "Need three elements for beam. You have " << beam);
            ASKAPLOG_INFO_STR(logger, "Restore solver will convolve with the 2D gaussian: " << beam(0) <<
                              " x " << beam(1) << " pixels at position angle " << beam(2) << " degrees");
            itsBMaj = beam(0);
            itsBMin = beam(1);
            itsBPa = beam(2);
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::setModel(const Array<T> model, const uInt term)
        {
            ASKAPCHECK(term < itsNumberTerms, "Term " << term << " greater than allowed " << itsNumberTerms);
            this->itsModel(term) = model.nonDegenerate().copy();
            this->validateShapes();
        }

        template<class T, class FT>
        Array<T> & DeconvolverBase<T, FT>::model(const uInt term)
        {
            ASKAPCHECK(term < itsNumberTerms, "Term " << term << " greater than allowed " << itsNumberTerms);
            return itsModel(term);
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::updateDirty(Array<T>& dirty, const uInt term)
        {
            ASKAPCHECK(term < itsNumberTerms, "Term " << term << " greater than allowed " << itsNumberTerms);
            if (!dirty.shape().nonDegenerate().conform(this->dirty(term).shape())) {
                throw(AskapError("Updated dirty image has different shape"));
            }
            this->itsDirty(term) = dirty.nonDegenerate().copy();
            this->validateShapes();
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::updateDirty(Vector<Array<T> >& dirtyVec)
        {
            if (dirtyVec.nelements() != this->itsDirty.nelements()) {
                throw(AskapError("Updated dirty image has different shape"));
            }
            this->itsDirty.resize(dirtyVec.nelements());
            for (uInt term = 0; term < dirtyVec.nelements(); term++) {
                if (!dirtyVec(term).nonDegenerate().shape().conform(this->itsDirty(term).nonDegenerate().shape())) {
                    throw(AskapError("Updated dirty image has different shape from original"));
                }
                this->itsDirty(term) = dirtyVec(term).nonDegenerate().copy();
            }
            this->validateShapes();
        }

        template<class T, class FT>
        bool DeconvolverBase<T, FT>::deconvolve()
        {
            throw(AskapError("Called base class deconvolver"));
        }

        template<class T, class FT>
        Array<T> & DeconvolverBase<T, FT>::dirty(const uInt term)
        {
            ASKAPCHECK(term < itsNumberTerms, "Term " << term << " greater than allowed " << itsNumberTerms);
            return itsDirty(term);
        }

        template<class T, class FT>
        Array<T> & DeconvolverBase<T, FT>::psf(const uInt term)
        {
            ASKAPCHECK(term < itsNumberTerms, "Term " << term << " greater than allowed " << itsNumberTerms);
            return itsPsf(term);
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::setWeight(Array<T> weight, const uInt term)
        {
            ASKAPCHECK(term < itsNumberTerms, "Term " << term << " greater than allowed " << itsNumberTerms);
            this->itsWeight(term) = weight.nonDegenerate().copy();
        }

        template<class T, class FT>
        Array<T> & DeconvolverBase<T, FT>::weight(const uInt term)
        {
            ASKAPCHECK(term < itsNumberTerms, "Term " << term << " greater than allowed " << itsNumberTerms);
            return itsWeight(term);
        }

        template<class T, class FT>
        boost::shared_ptr<DeconvolverControl<T> > DeconvolverBase<T, FT>::control() const
        {
            ASKAPASSERT(itsDC);
            return itsDC;
        }

        template<class T, class FT>
        bool DeconvolverBase<T, FT>::setControl(boost::shared_ptr<DeconvolverControl<T> > DC)
        {
            itsDC = DC;
            ASKAPASSERT(itsDC);
            return True;
        }

        template<class T, class FT>
        boost::shared_ptr<DeconvolverMonitor<T> > DeconvolverBase<T, FT>::monitor() const
        {
            ASKAPASSERT(itsDM);
            return itsDM;
        }

        template<class T, class FT>
        bool DeconvolverBase<T, FT>::setMonitor(boost::shared_ptr<DeconvolverMonitor<T> > DM)
        {
            itsDM = DM;
            ASKAPASSERT(itsDM);
            return True;
        }

        template<class T, class FT>
        boost::shared_ptr<DeconvolverState<T> > DeconvolverBase<T, FT>::state() const
        {
            ASKAPASSERT(itsDS);
            return itsDS;
        }

        template<class T, class FT>
        bool DeconvolverBase<T, FT>::setState(boost::shared_ptr<DeconvolverState<T> > DS)
        {
            itsDS = DS;
            ASKAPASSERT(itsDS);
            return True;
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::validateShapes()
        {
            casa::IPosition minPos;
            casa::IPosition maxPos;
            T minVal, maxVal;
            casa::minMax(minVal, maxVal, minPos, maxPos, this->psf(0));
            const Int nx(this->psf(0).shape()(0));
            const Int ny(this->psf(0).shape()(1));
            if ((maxPos(0) != nx / 2) || (maxPos(1) != ny / 2)) {
                ASKAPTHROW(AskapError, "Peak of PSF(0) is at " << maxPos << ": not at centre pixel: [" << nx / 2 << "," << ny / 2 << "]");
            }

            for (uInt term = 0; term < itsNumberTerms; term++) {
                if (!(this->dirty(term).shape().size())) {
                    ASKAPTHROW(AskapError, "Dirty image has zero size");
                }
                if (!(this->psf(term).shape().size())) {
                    ASKAPTHROW(AskapError, "PSF image has zero size");
                }
                if (!(this->psf(term).shape()[0] == this->dirty().shape()[0]) || !(this->psf(term).shape()[1] == this->dirty().shape()[1])) {
                    ASKAPTHROW(AskapError, "PSF has different shape from dirty image");
                }

                // The model and dirty image shapes only need to agree on the
                // first two axes
                if (!(this->model(term).shape().size())) {
                    ASKAPTHROW(AskapError, "Model has zero size");
                }
                if (!(this->model(term).shape()[0] == this->dirty().shape()[0]) || !(this->model(term).shape()[1] == this->dirty().shape()[1])) {
                    ASKAPTHROW(AskapError, "Model has different shape from dirty image");
                }
            }
            if (!this->itsPeakPSFPos.size()) {
                ASKAPTHROW(AskapError, "Position of PSF peak not defined " << this->itsPeakPSFPos);
            }
            if (this->itsPeakPSFVal == 0.0) {
                ASKAPTHROW(AskapError, "PSF peak is zero");
            }
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::initialise()
        {
            ASKAPLOG_INFO_STR(decbaselogger, "Initialising weight images");

            // Always check shapes on initialise
            this->validateShapes();

        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::finalise()
        {
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::updateResiduals(Vector<Array<T> >& model)
        {
            ASKAPCHECK(model.shape() == itsNumberTerms, "Number of terms in model " << model.shape()
                           << " not same as number of terms specified "
                           << itsNumberTerms);

            for (uInt term = 0; term < itsNumberTerms; term++) {
                Array<FT> xfr;
                xfr.resize(psf(term).shape());
                casa::setReal(xfr, psf(term));
                scimath::fft2d(xfr, true);
                Array<FT> work;
                // Find residuals for current model model
                work.resize(model(term).shape());
                work.set(FT(0.0));
                casa::setReal(work, model(term));
                scimath::fft2d(work, true);
                work = xfr * work;
                scimath::fft2d(work, false);
                this->dirty(term) = this->dirty(term) - real(work);
            }
        }

        template<class T, class FT>
        bool DeconvolverBase<T, FT>::restore(Vector<Array<T> >& restored)
        {
            return restore(restored, itsModel);
        }

        template<class T, class FT>
        bool DeconvolverBase<T, FT>::restore(Vector<Array<T> >& restored, Vector<Array<T> >& model)
        {
            if (itsBMaj <= 0.0) {
                ASKAPLOG_INFO_STR(decbaselogger, "Beam not specified - no restoration");
                return false;
            }

            ASKAPCHECK(model.shape() == itsNumberTerms, "Number of terms in model " << model.shape()
                           << " not same as number of terms specified "
                           << itsNumberTerms);

            ASKAPCHECK(restored.shape() == itsNumberTerms, "Number of terms in restored image " << restored.shape()
                           << " not same as number of terms specified "
                           << itsNumberTerms);

            const int nx(model(0).shape()(0));
            const int ny(model(0).shape()(1));
            const IPosition centre(2, nx / 2, ny / 2);
            Matrix<FT> gaussian(nx, ny);
            gaussian.set(0.0);

            const float scalex(sqrt(4.0*log(2.0)) / itsBMaj);
            const float scaley(sqrt(4.0*log(2.0)) / itsBMin);
            const float theta(itsBPa*C::pi / 180.0);
            const float cs(cos(theta));
            const float sn(sin(theta));
            for (int y = 0; y < nx; y++) {
                const float dy = (y - ny / 2) * scaley;
                for (int x = 0; x < nx; x++) {
                    const float dx = (x - nx / 2) * scalex;
                    const float rsq = pow((dx * cs + dy * sn), 2) + pow((-dx * sn + dy * cs), 2);
                    if (rsq < 20.0) {
                        gaussian(x, y) = exp(-rsq);
                    }
                }
            }
            const float volume(sum(real(gaussian)));

            scimath::fft2d(gaussian, true);

            ASKAPLOG_INFO_STR(logger, "Volume of PSF = " << volume << " pixels");

            for (uInt term = 0; term < itsNumberTerms; term++) {
                Array<FT> vis(model(term).shape());
                vis.set(FT(0.0));
                casa::setReal(vis, model(term));
                scimath::fft2d(vis, true);
                vis = vis * gaussian;
                scimath::fft2d(vis, false);
                restored(term).resize(model(term).shape());
                restored(term) = this->dirty(term) + real(vis);
            }
            return true;
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::updateResiduals(Array<T> & model)
        {
            Vector<Array<T> > modelVec(1);
            modelVec(0) = model;
            updateResiduals(modelVec);
        }

        template<class T, class FT>
        uInt DeconvolverBase<T, FT>::auditMemory(Vector<casa::Array<T> >& vecArray)
        {
            uInt memory = 0;
            for (uInt term = 0; term < vecArray.nelements(); term++) {
                memory += sizeof(T) * vecArray(term).nelements();
            }
            return memory;
        }

        template<class T, class FT>
        uInt DeconvolverBase<T, FT>::auditMemory(Vector<casa::Array<FT> >& vecArray)
        {
            uInt memory = 0;
            for (uInt term = 0; term < vecArray.nelements(); term++) {
                memory += sizeof(FT) * vecArray(term).nelements();
            }
            return memory;
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::auditAllMemory()
        {
            ASKAPLOG_DEBUG_STR(decbaselogger, "Dirty images  " << auditMemory(itsDirty));
            ASKAPLOG_DEBUG_STR(decbaselogger, "PSFs          " << auditMemory(itsPsf));
            ASKAPLOG_DEBUG_STR(decbaselogger, "Models        " << auditMemory(itsModel));
            ASKAPLOG_DEBUG_STR(decbaselogger, "Weight images " << auditMemory(itsWeight));
        }

        template<class T, class FT>
        IPosition DeconvolverBase<T, FT>::findSubPsfShape()
        {
            IPosition subPsfShape(2, this->model().shape()(0), this->model().shape()(1));
            // Only use the specified psfWidth if it makes sense
            if (this->control()->psfWidth() > 0) {
                uInt psfWidth = this->control()->psfWidth();
                if ((psfWidth < uInt(this->model().shape()(0))) && (psfWidth < uInt(this->model().shape()(1)))) {
                    subPsfShape(0) = psfWidth;
                    subPsfShape(1) = psfWidth;
                }
            }
            return subPsfShape;
        }
    } // namespace synthesis

} // namespace askap


