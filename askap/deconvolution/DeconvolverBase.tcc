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

#include <casacore/casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(decbaselogger, ".deconvolution.base");

#include <askap/deconvolution/DeconvolverBase.h>
#include <askap/deconvolution/DeconvolverState.h>

namespace askap {

    namespace synthesis {

        template<class T, class FT>
        DeconvolverBase<T, FT>::~DeconvolverBase()
        {
            auditAllMemory();
        };

        template<class T, class FT>
        DeconvolverBase<T, FT>::DeconvolverBase(Vector<Array<T> >& dirty, Vector<Array<T> >& psf) :
                itsPeakPSFVal(0.), itsBMaj(0.0), itsBMin(0.0), itsBPa(0.0)
        {
            init(dirty, psf);
        }

        template<class T, class FT>
        DeconvolverBase<T, FT>::DeconvolverBase(Array<T>& dirty, Array<T>& psf) :
                itsPeakPSFVal(0.), itsBMaj(0.0), itsBMin(0.0), itsBPa(0.0)
        {
            Vector<Array<T> > dirtyVec(1);
            dirtyVec(0) = dirty.nonDegenerate();
            Vector<Array<T> > psfVec(1);
            psfVec(0) = psf.nonDegenerate();
            init(dirtyVec, psfVec);
        }

        /// @brief validate PSF, find peak value and position
        /// @details It works with the zero-th term, if there are many
        /// @param[in] slicer optional slicer if only a fraction of the PSF needs to be considered
        /// the default constructed instance of a slicer results in the whole PSF being used.
        /// @note this method updates itsPeakPSFPos and itsPeakPSFVal, this is why it's non-const
        template<typename T, typename FT>
        void DeconvolverBase<T, FT>::validatePSF(const casacore::Slicer &slicer)
        {
            ASKAPLOG_INFO_STR(decbaselogger, "Validating PSF");
            casacore::Array<T> psfArr = slicer == casacore::Slicer() ? psf(0) : psf(0).nonDegenerate()(slicer);

            casacore::IPosition minPos;
            casacore::IPosition maxPos;
            T minVal, maxVal;
            casacore::minMax(minVal, maxVal, minPos, maxPos, psfArr);

            ASKAPASSERT(psfArr.shape().nelements() >= 2);
            const Int nx(psfArr.shape()(0));
            const Int ny(psfArr.shape()(1));

            ASKAPLOG_INFO_STR(decbaselogger, "Maximum of PSF(0) = " << maxVal << " at " << maxPos);
            ASKAPLOG_INFO_STR(decbaselogger, "Minimum of PSF(0) = " << minVal << " at " << minPos);

            ASKAPDEBUGASSERT(maxPos.nelements() >= 2);
            if ((maxPos(0) != nx / 2) || (maxPos(1) != ny / 2)) {
                ASKAPTHROW(AskapError, "Peak of PSF(0) is at " << maxPos << ": not at centre pixel: [" << nx / 2 << "," << ny / 2 << "]");
            }

            itsPeakPSFVal = maxVal;
            itsPeakPSFPos = maxPos;
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::init(Vector<Array<T> >& dirtyVec, Vector<Array<T> >& psfVec)
        {

            ASKAPCHECK(psfVec.nelements() == dirtyVec.nelements(),
                       "Vectors of dirty images and PSF's not same length");

            itsNumberTerms = dirtyVec.nelements();

            itsDirty.resize(nTerms());
            itsPsf.resize(nTerms());
            itsModel.resize(nTerms());
            itsWeight.resize(nTerms());

            ASKAPLOG_INFO_STR(decbaselogger, "There are " << nTerms() << " dirty images");

            for (uInt term = 0; term < nTerms(); ++term) {

                ASKAPASSERT(dirtyVec(term).nonDegenerate().shape().nelements() == 2);
                ASKAPASSERT(psfVec(term).nonDegenerate().shape().nelements() == 2);

                itsDirty(term).reference(dirtyVec(term).nonDegenerate());
                itsPsf(term).reference(psfVec(term).nonDegenerate());

                ASKAPASSERT(itsPsf(term).shape().conform(itsDirty(term).shape()));

                ASKAPLOG_INFO_STR(decbaselogger, "Dirty image(" << term << ") has shape: "
                                      << dirty(term).shape());

                model(term).resize(dirty(term).shape());
                model(term).set(T(0.0));
            }

            validatePSF();

            itsDS = boost::shared_ptr<DeconvolverState<T> >(new DeconvolverState<T>());
            ASKAPASSERT(itsDS);
            itsDC = boost::shared_ptr<DeconvolverControl<T> >(new DeconvolverControl<T>());
            ASKAPASSERT(itsDC);
            itsDM = boost::shared_ptr<DeconvolverMonitor<T> >(new DeconvolverMonitor<T>());
            ASKAPASSERT(itsDM);

            validateShapes();

            auditAllMemory();

        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::configure(const LOFAR::ParameterSet& parset)
        {
            itsDC->configure(parset);
            itsDM->configure(parset);

            // Get the beam information
            if (parset.isDefined("beam")) {
                const std::vector<float> beam = parset.getFloatVector("beam");
                ASKAPCHECK(beam.size() == 3u, "Need three elements for beam. You have " << beam);
                ASKAPLOG_INFO_STR(decbaselogger, "Restore solver will convolve with the 2D gaussian: " << beam[0] <<
                                  " x " << beam[1] << " pixels at position angle " << beam[2] << " degrees");
                itsBMaj = beam[0];
                itsBMin = beam[1];
                itsBPa = beam[2];
            }
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::setModel(const Array<T>& model, const uInt term)
        {
            ASKAPCHECK(term < nTerms(), "Term " << term << " greater than allowed " << nTerms());
            itsModel(term) = model.nonDegenerate();
            validateShapes();
        }

        template<class T, class FT>
        const Array<T> & DeconvolverBase<T, FT>::model(const uInt term) const
        {
            ASKAPCHECK(term < nTerms(), "Term " << term << " greater than allowed " << nTerms());
            return itsModel(term);
        }

        template<class T, class FT>
        Array<T> & DeconvolverBase<T, FT>::model(const uInt term)
        {
            ASKAPCHECK(term < nTerms(), "Term " << term << " greater than allowed " << nTerms());
            return itsModel(term);
        }


        template<class T, class FT>
        void DeconvolverBase<T, FT>::updateDirty(const Array<T>& newDirty, const uInt term)
        {
            ASKAPLOG_DEBUG_STR(decbaselogger,"updateDirty(dirty,term)");
            ASKAPCHECK(term < nTerms(), "Term " << term << " greater than allowed " << nTerms());
            ASKAPCHECK(dirty(term).nelements()==0 || newDirty.shape().nonDegenerate()==dirty(term).shape(),
                "Updated dirty image has different shape from original");
            itsDirty(term) = newDirty.nonDegenerate();
            validateShapes();
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::updateDirty(const Vector<Array<T> >& dirtyVec)
        {
            ASKAPLOG_DEBUG_STR(decbaselogger,"updateDirty(dirtyVec)");
            ASKAPCHECK(dirtyVec.nelements() == itsDirty.nelements(),
                "Number of taylor terms in updated dirty image differs from original");
            for (uInt term = 0; term < dirtyVec.nelements(); term++) {
                ASKAPCHECK(dirty(term).nelements()==0 || dirtyVec(term).shape().nonDegenerate()==dirty(term).shape(),
                    "Updated dirty image has different shape from original");
                itsDirty(term) = dirtyVec(term).nonDegenerate();
            }
            validateShapes();
        }

        template<class T, class FT>
        bool DeconvolverBase<T, FT>::deconvolve()
        {
            throw(AskapError("Called base class deconvolver"));
        }

        template<class T, class FT>
        Array<T> & DeconvolverBase<T, FT>::dirty(const uInt term)
        {
            ASKAPCHECK(term < nTerms(), "Term " << term << " greater than allowed " << nTerms());
            return itsDirty(term);
        }

        template<class T, class FT>
        Array<T> & DeconvolverBase<T, FT>::psf(const uInt term)
        {
            ASKAPCHECK(term < nTerms(), "Term " << term << " greater than allowed " << nTerms());
            return itsPsf(term);
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::setWeight(Array<T> weight, const uInt term)
        {
            ASKAPCHECK(term < nTerms(), "Term " << term << " greater than allowed " << nTerms());
            itsWeight(term).reference(weight.nonDegenerate());
        }

        template<class T, class FT>
        Array<T> & DeconvolverBase<T, FT>::weight(const uInt term)
        {
            ASKAPCHECK(term < nTerms(), "Term " << term << " greater than allowed " << nTerms());
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
            for (uInt term = 0; term < nTerms(); ++term) {
                 ASKAPCHECK(dirty(term).shape().size() > 0, "Dirty image has zero size for term="<<term);
                 ASKAPCHECK(psf(term).shape().size() > 0, "PSF image has zero size for term="<<term);
                 ASKAPCHECK(psf(term).shape() == dirty().shape(), "PSF has different shape from dirty image for term="<<term);

                 // The model and dirty image shapes only need to agree on the first two axes
                 ASKAPCHECK(model(term).shape().size() > 0, "Model has zero size for term="<<term);
                 ASKAPCHECK(model(term).shape().getFirst(2) == dirty().shape().getFirst(2), "Model has different shape from dirty image for term="<<term);
            }
            ASKAPCHECK(itsPeakPSFPos.size() > 0, "Position of PSF peak not defined " << itsPeakPSFPos);
            ASKAPCHECK(itsPeakPSFVal > 0., "PSF peak is supposed to be positive");
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::initialise()
        {
            ASKAPLOG_INFO_STR(decbaselogger, "Initialising weight images");

            // Always check shapes on initialise
            validateShapes();
            // MV: after I removed another psf peak search from validateShapes some unit tests started to fail
            // because an exception is not thrown. It looks like this method can be another entry point (for no good reason),
            // running psf validation here explicitly. There is some technical debt here, perhaps more thoughts are needed on
            // how to design interfaces of these classes
            validatePSF();
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::finalise()
        {
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::updateResiduals(Vector<Array<T> >& model)
        {
            ASKAPCHECK(model.nelements() == nTerms(), "Number of terms in model " << model.nelements()
                           << " not same as number of terms specified "
                           << nTerms());
            ASKAPASSERT(nTerms() > 0);
            // MV: there is some technical debt here in handling arrays/matrices/degenerate dimenions
            // I don't think this code would've worked correctly anyway if we had extra dimensions
            // (support of which has been removed from the FFT wrapper a while back, but is only apparent now
            // when we need a more rigid interface with the new wrapper). So I am changing Array to Matrix and
            // it will throw exception if we ever encounter extra dimensions
            //
            // uninitialized construction doesn't work for complex type any more in casacore-3.4
            // default constructor of complex will be called for each element
            Matrix<FT> xfr(psf(0).shape());
            Matrix<FT> work(model(0).shape());

            // Don't reorder harmonics as with the original wrapper (hence, pass false to the wrapper), it seems possible to
            // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
            // Limit number of fft threads to 8 (more is slower for our fft sizes)
            scimath::FFT2DWrapper<FT> fft2d(true,8); // MHW changed to true - work image was shifted by half size, circularly in x and y

            for (uInt term = 0; term < nTerms(); ++term) {
                 const Array<T>& thisTermPSF = psf(term);
                 xfr.resize(thisTermPSF.shape());
                 xfr.set(0.);
                 casacore::setReal(xfr, thisTermPSF);
                 fft2d(xfr, true);
                 // Find residuals for current model model
                 const Array<T>& thisTermModel = model(term);
                 work.resize(thisTermModel.shape());
                 work.set(0.);
                 casacore::setReal(work, thisTermModel);
                 fft2d(work, true);
                 work *= xfr;
                 fft2d(work, false);
                 dirty(term) -= real(work);
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

            ASKAPCHECK(model.nelements() == nTerms(), "Number of terms in model " << model.nelements()
                           << " not same as number of terms specified "
                           << nTerms());

            ASKAPCHECK(restored.nelements() == nTerms(), "Number of terms in restored image " << restored.nelements()
                           << " not same as number of terms specified "
                           << nTerms());

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

            // Don't reorder harmonics as with the original wrapper (hence, pass false to the wrapper), it seems possible to
            // skip it here as we use FFT to do convolutions and don't care about particular harmonic placement in the Fourier space
            // Limit number of fft threads to 8 (more is slower for our fft sizes)
            scimath::FFT2DWrapper<FT> fft2d(true,8); //MHW changed to true - components were restored in corners instead of centre

            fft2d(gaussian, true);

            ASKAPLOG_INFO_STR(decbaselogger, "Volume of PSF = " << volume << " pixels");

            for (uInt term = 0; term < nTerms(); ++term) {
                Matrix<FT> vis(model(term).shape());
                vis.set(FT(0.0));
                casacore::setReal(vis, model(term));
                fft2d(vis, true);
                vis *= gaussian;
                fft2d(vis, false);
                restored(term).resize(model(term).shape());
                restored(term) = dirty(term) + real(vis);
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

        /// @brief memory occupied by the array set
        /// @details This method iterates over the supplied array set and computes the
        /// total amount of memory used up.
        /// @param[in] vecArray arrays to work with (packed into a casa Vector)
        /// @return amount of memory in bytes
        template<class T, class FT>
        template<typename Y>
        uInt DeconvolverBase<T, FT>::auditMemory(Vector<casacore::Array<Y> >& vecArray)
        {
            uInt memory = 0;
            for (uInt term = 0; term < vecArray.nelements(); term++) {
                memory += sizeof(Y) * vecArray(term).nelements();
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
        IPosition DeconvolverBase<T, FT>::findSubPsfShape(const casacore::uInt term) const
        {
            const casacore::IPosition modelShape = model(term).shape();
            ASKAPDEBUGASSERT(modelShape.nelements() >= 2);
            IPosition subPsfShape(2, modelShape(0), modelShape(1));
            // Only use the specified psfWidth if it makes sense
            // Also make sure it is even, otherwise the fft of the PSF is incorrect and the clean diverges

            if (control()->psfWidth() > 1) {
                const uInt psfWidth = static_cast<casacore::uInt>(control()->psfWidth()/2)*2;
                if ((psfWidth < static_cast<casacore::uInt>(modelShape(0))) && (psfWidth < static_cast<casacore::uInt>(modelShape(1)))) {
                    subPsfShape(0) = psfWidth;
                    subPsfShape(1) = psfWidth;
                }
            }
            return subPsfShape;
        }

        template<class T, class FT>
        void DeconvolverBase<T, FT>::releaseMemory()
        {
            uInt memory = 0;
            for (uInt term = 0; term < nTerms(); term++) {
                memory += sizeof(T) * itsDirty(term).nelements();
                itsDirty(term).resize();
            }
            ASKAPLOG_DEBUG_STR(decbaselogger,"DeconvolverBase released "<<memory/1024/1024<<" MB");
        }

    } // namespace synthesis

} // namespace askap
