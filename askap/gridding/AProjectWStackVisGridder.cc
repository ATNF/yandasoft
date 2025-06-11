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

#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/casa/Quanta/MVDirection.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/casa/Quanta/MVTime.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>

ASKAP_LOGGER(logger, ".gridding.aprojectwstackgridder");

#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>
#include <askap/profile/AskapProfiler.h>

#include <askap/gridding/AProjectWStackVisGridder.h>
#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/gridding/IBasicIllumination.h>

#include <askap/scimath/utils/PaddingUtils.h>

namespace askap {
namespace synthesis {

AProjectWStackVisGridder::AProjectWStackVisGridder(const boost::shared_ptr<IBasicIllumination const> &illum,
        const double wmax, const int nwplanes, const double,
        const int overSample, const int maxSupport, const int limitSupport,
        const int maxFeeds, const int maxFields,
        const double pointingTol, const double paTol,
        const double freqTol, const bool frequencyDependent,
        const bool spheroidalTaper, const double spheroidalWeightsCutoff,
        const std::string& name) :
    AProjectGridderBase(maxFeeds,maxFields,pointingTol, paTol, freqTol),
    WStackVisGridder(wmax, nwplanes),
    itsReferenceFrequency(0.0),
    itsIllumination(illum),
    itsMaxFeeds(maxFeeds), itsMaxFields(maxFields),
    itsFreqDep(frequencyDependent),
    itsSpheroidalTaper(spheroidalTaper),
    itsSpheroidalWeightsCutoff(spheroidalWeightsCutoff)
{
    ASKAPCHECK(overSample>0, "Oversampling must be greater than 0");
    ASKAPCHECK(maxSupport>0, "Maximum support must be greater than 0")
    ASKAPDEBUGASSERT(itsIllumination);
    itsSupport=0;
    itsOverSample=overSample;
    itsMaxSupport=maxSupport;
    itsLimitSupport=limitSupport;
    setTableName(name);
}

/// @brief copy constructor
/// @details It is required to decouple internal arrays between input object
/// and this copy.
/// @param[in] other input object
/// @note illumination pattern is copied as a shared pointer, hence referencing
/// the same model

AProjectWStackVisGridder::AProjectWStackVisGridder(const AProjectWStackVisGridder &other) :
    IVisGridder(other),
    AProjectGridderBase(other), WStackVisGridder(other),
    itsReferenceFrequency(other.itsReferenceFrequency),
    itsIllumination(other.itsIllumination), itsMaxFeeds(other.itsMaxFeeds),
    itsMaxFields(other.itsMaxFields),
    itsFreqDep(other.itsFreqDep),
    itsSpheroidalTaper(other.itsSpheroidalTaper),
    itsSpheroidalWeightsCutoff(other.itsSpheroidalWeightsCutoff),
    itsMaxSupport(other.itsMaxSupport),
    itsLimitSupport(other.itsLimitSupport),
    itsCMap(other.itsCMap.copy())
{
}


/// Clone a copy of this Gridder
IVisGridder::ShPtr AProjectWStackVisGridder::clone() {
    return IVisGridder::ShPtr(new AProjectWStackVisGridder(*this));
}

/// @brief initialise sum of weights
/// @details We keep track the number of times each convolution function is used per
/// channel and polarisation (sum of weights). This method is made virtual to be able
/// to do gridder specific initialisation without overriding initialiseGrid.
/// This method accepts no parameters as itsShape, itsNWPlanes, etc should have already
/// been initialised by the time this method is called.
void AProjectWStackVisGridder::initialiseSumOfWeights()
{
    // this method is hopefully just a temporary stub until we figure out a better way of
    // managing a cache of convolution functions. It skips initialisation if itsSupport is
    // not zero, which means that some initialisation has been done before.
    // Note, it is not a very good way of doing things!
    if (itsSupport == 0) {
        WStackVisGridder::initialiseSumOfWeights();
    }
    // Reset the weights
    zeroSumOfWeights();
}

/// Initialize the indices into the cube.
void AProjectWStackVisGridder::initIndices(const accessors::IConstDataAccessor& acc) {
    ASKAPTRACE("AProjectWStackVisGridder::initIndices");

    // this calculates current field id
    indexField(acc);

    const int nSamples = acc.nRow();
    const int nChan = acc.nChannel();
    const int nPol = acc.nPol();

    itsCMap.resize(nSamples, nPol, nChan);
    itsCMap.set(0);

    itsGMap.resize(nSamples, nPol, nChan);
    itsGMap.set(0);

    const casacore::Vector<casacore::RigidVector<double, 3> > &rotatedUVW = acc.rotatedUVW(getTangentPoint());

    for (int i=0; i<nSamples; ++i) {
        const int feed=acc.feed1()(i);
        ASKAPCHECK(feed<itsMaxFeeds, "Exceeded specified maximum number of feeds");
        ASKAPCHECK(feed>-1, "Illegal negative feed number");

        const double w=(rotatedUVW(i)(2))/(casacore::C::c);

        for (int chan=0; chan<nChan; ++chan) {
            const double freq=acc.frequency()[chan];
            for (int pol=0; pol<nPol; pol++) {
                int index = -1;
                /// Order is (chan, feed)
                if(itsFreqDep) {
                    index = chan+nChan*(feed+itsMaxFeeds*currentField());
                    ASKAPCHECK(index<itsMaxFields*itsMaxFeeds*nChan, "CMap index too large");
                } else {
                    index = (feed+itsMaxFeeds*currentField());
                    ASKAPCHECK(index < itsMaxFields*itsMaxFeeds, "CMap index too large");
                }

                ASKAPCHECK(index>-1, "CMap index less than zero");

                itsCMap(i, pol, chan) = index;

                /// Calculate the index into the grids
                itsGMap(i, pol, chan) = getWPlane(w*freq);
            }
        }
    }
}
/// @brief Initialise the gridding
/// @param axes axes specifications
/// @param shape Shape of output image: u,v,pol,chan
/// @param dopsf Make the psf?
void AProjectWStackVisGridder::initialiseGrid(const scimath::Axes& axes,
                                              const casacore::IPosition& shape,
                                              const bool dopsf,
                                              const bool dopcf)
{
    ASKAPTRACE("AProjectWStackVisGridder::initialiseGrid");
    WStackVisGridder::initialiseGrid(axes,shape,dopsf,dopcf);

    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casacore::uInt nx=std::min(itsMaxSupport, int(itsShape(0)));
    const casacore::uInt ny=std::min(itsMaxSupport, int(itsShape(1)));

    ASKAPLOG_INFO_STR(logger, "Shape for calculating gridding convolution function = "
            << nx << " by " << ny << " pixels");

    // this is just a buffer in the uv-space
    initUVPattern(nx,ny, itsUVCellSize(0),itsUVCellSize(1),itsOverSample);
}

/// @brief Initialise the degridding
/// @param axes axes specifications
/// @param image Input image: cube: u,v,pol,chan
void AProjectWStackVisGridder::initialiseDegrid(const scimath::Axes& axes,
        const casacore::Array<imtype>& image)
{
    ASKAPTRACE("AProjectWStackVisGridder::initialiseDegrid");
    WStackVisGridder::initialiseDegrid(axes,image);
    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casacore::uInt nx=std::min(itsMaxSupport, int(itsShape(0)));
    const casacore::uInt ny=std::min(itsMaxSupport, int(itsShape(1)));

    ASKAPLOG_INFO_STR(logger, "Shape for calculating degridding convolution function = "
            << nx << " by " << ny << " pixels");

    // this is just a buffer in the uv-space
    initUVPattern(nx,ny, itsUVCellSize(0),itsUVCellSize(1),itsOverSample);
}

/// Initialize the convolution function into the cube. If necessary this
/// could be optimized by using symmetries.
/// @todo Make initConvolutionFunction more robust
void AProjectWStackVisGridder::initConvolutionFunction(const accessors::IConstDataAccessor& acc) {

    ASKAPTRACE("AProjectWStackVisGridder::initConvolutionFunction");
    ASKAPDEBUGASSERT(itsIllumination);
    // just to avoid a repeated call to a virtual function from inside the loop
    const bool hasSymmetricIllumination = itsIllumination->isSymmetric();
    // check whether the output pattern is image based. If so, inverse FFT is not needed
    const bool imageBasedPattern = itsIllumination->isImageBased();
    const int nSamples = acc.nRow();

    validateCFCache(acc, hasSymmetricIllumination);

    casacore::MVDirection out = getImageCentre();

    /// We have to calculate the lookup function converting from
    /// row and channel to plane of the w-dependent convolution
    /// function
    const int nChan = itsFreqDep ? acc.nChannel() : 1;

    if(itsSupport==0) {
        ASKAPLOG_INFO_STR(logger, "Resizing convolution function to "
                << itsOverSample << "*" << itsOverSample << "*" << itsMaxFeeds << "*" << itsMaxFields << "*" << nChan << " entries");
        itsConvFunc.resize(itsOverSample*itsOverSample*itsMaxFeeds*itsMaxFields*nChan);

        ASKAPLOG_INFO_STR(logger, "Resizing sum of weights to " << itsMaxFeeds << "*" << itsMaxFields << "*" << nChan << " entries");
        resizeSumOfWeights(itsMaxFeeds*itsMaxFields*nChan);
        zeroSumOfWeights();
    }

    UVPattern &pattern = uvPattern();
    const casacore::uInt nx = pattern.uSize();
    const casacore::uInt ny = pattern.vSize();

    const casacore::uInt qnx = nx / itsOverSample;
    const casacore::uInt qny = ny / itsOverSample;

    casacore::Vector<double> ccfx;
    casacore::Vector<double> ccfy;
    if (itsSpheroidalTaper) {
        // Include spheroidal for anti-aliasing and general kernel robustness
        ccfx.resize(qnx);
        ccfy.resize(qny);
        for (casacore::uInt qix = 0; qix < qnx; ++qix) {
            const double nux = std::abs(double(qix) - double(qnx / 2)) / double(qnx / 2);
            ccfx(qix) = grdsf(nux);
        }
        for (casacore::uInt qiy = 0; qiy < qny; ++qiy) {
            const double nuy = std::abs(double(qiy) - double(qny / 2)) / double(qny / 2);
            ccfy(qiy) = grdsf(nuy);
        }
        if (doInterpolation()) {
          // The spheroidal is undefined and set to zero at nu=1, but that
          // is not the numerical limit. Estimate it from its neighbours.
          interpolateEdgeValues(ccfx);
          interpolateEdgeValues(ccfy);
        }
    }


    casacore::Matrix<imtypeComplex> cplane(pattern.pattern().shape());
    // Limit number of fft threads to 8 (more is slower for our fft sizes)
    scimath::FFT2DWrapper<imtypeComplex> fft2d(true,8);
    int nDone=0;
    for (int row=0; row<nSamples; ++row) {
        const int feed=acc.feed1()(row);

        if (!isCFValid(feed, currentField())) {
            makeCFValid(feed, currentField());
            nDone++;
            casacore::MVDirection offset(acc.pointingDir1()(row).getAngle());
            const double parallacticAngle = hasSymmetricIllumination ? 0. : acc.feed1PA()(row);

            for (int chan=0; chan<nChan; chan++) {
                /// Extract illumination pattern for this channel
                itsIllumination->getPattern(acc.frequency()[chan], pattern,
                                            out, offset, parallacticAngle, isPSFGridder() || isPCFGridder());
                /// Now convolve the disk with itself using an FFT
                if( !imageBasedPattern ) {
                    fft2d(pattern.pattern(), false);
                }

                double peak=0.0;
                cplane.set(0.);
                for (casacore::uInt iy=0; iy<ny; ++iy) {
                    const int qiy = iy + int(qny)/2 - int(ny)/2;
                    for (casacore::uInt ix=0; ix<nx; ++ix) {
                        const int qix = ix + int(qnx)/2 - int(nx)/2;
                        if (!itsSpheroidalTaper || ((qix>=0) && (qix<qnx) && (qiy>=0) && (qiy<qny))) {
                            const imtype taper = itsSpheroidalTaper ? ccfx(qix) * ccfy(qiy) : 1.0;
                            cplane(ix, iy) = pattern(ix, iy) * conj(pattern(ix,iy)) * taper;
                            if (casacore::abs(cplane(ix,iy))>peak) {
                                peak=casacore::abs(cplane(ix,iy));
                            }
                        }
                    }
                }

                if(peak>0.0) {
                    cplane *= imtypeComplex(1.0/peak);
                }
                // The maximum will be 1.0
                ASKAPLOG_DEBUG_STR(logger, "Max of FT of convolution function = " << casacore::max(cplane));
                fft2d(cplane, true);
                // Now correct for normalization of FFT
                cplane *= imtypeComplex(1.0/(double(nx)*double(ny)));
                ASKAPLOG_DEBUG_STR(logger, "Sum of convolution function before support extraction and decimation = " << casacore::sum(cplane));

                if (itsSupport==0) {
                    // we probably need a proper support search here
                    // it can be encapsulated in a method of the UVPattern class
                    itsSupport = pattern.maxSupport() / (2 * itsOverSample);
                    if (itsSupport < 3) {
                        itsSupport = 3;
                    }

                    ASKAPCHECK(itsSupport>0,
                            "Unable to determine support of convolution function");
                    ASKAPCHECK(itsSupport*itsOverSample<int(nx)/2,
                            "Overflowing convolution function - increase maxSupport or decrease overSample, "<<
                            "Current support size = " << itsSupport << " oversampling factor=" << itsOverSample <<
                            " image size nx=" << nx)
                    if (itsLimitSupport > 0  &&  itsSupport > itsLimitSupport) {
                        ASKAPLOG_INFO_STR(logger, "Convolution function support = "
                                << itsSupport << " pixels exceeds upper support limit; "
                                << "set to limit = " << itsLimitSupport << " pixels");
                        itsSupport = itsLimitSupport;
                    }

                    const int cSize=2*itsSupport+1;
                    // just for logging
                    const double cell = std::abs(pattern.uCellSize())*(casacore::C::c/acc.frequency()[chan]);
                    ASKAPLOG_DEBUG_STR(logger, "Convolution function support = "
                            << itsSupport << " pixels, size = " << cSize
                            << " pixels");
                    ASKAPLOG_DEBUG_STR(logger, "Maximum extent = "<< itsSupport
                            *cell << " (m) sampled at "<< cell
                            << " (m)");
                    ASKAPLOG_DEBUG_STR(logger, "Number of planes in convolution function = "
                            << itsConvFunc.size()<<" or "<<itsConvFunc.size()/itsOverSample/itsOverSample<<
                            " before oversampling with factor "<<itsOverSample);
                } // if itsSupport uninitialized
                int zIndex=chan+nChan*(feed+itsMaxFeeds*currentField());

                // Since we are decimating, we need to rescale by the
                // decimation factor
                const imtype rescale=imtype(itsOverSample*itsOverSample);
                const int cSize=2*itsSupport+1;
                for (int fracu=0; fracu<itsOverSample; fracu++) {
                    for (int fracv=0; fracv<itsOverSample; fracv++) {
                        int plane=fracu+itsOverSample*(fracv+itsOverSample*zIndex);
                        ASKAPDEBUGASSERT(plane>=0 && plane<int(itsConvFunc.size()));
                        itsConvFunc[plane].resize(cSize, cSize);
                        itsConvFunc[plane].set(0.0);
                        // Now cut out the inner part of the convolution function and
                        // insert it into the cache
                        for (int iy=-itsSupport; iy<itsSupport; iy++) {
                            for (int ix=-itsSupport; ix<itsSupport; ix++) {
                                itsConvFunc[plane](ix+itsSupport, iy+itsSupport)
                                    = rescale * cplane(itsOverSample*ix+fracu+nx/2,
                                            itsOverSample*iy+fracv+ny/2);
                            } // for ix
                        } // for iy
                        //
                        //ASKAPLOG_DEBUG_STR(logger, "convolution function for channel "<<chan<<
                        //   " plane="<<plane<<" has an integral of "<<sum(itsConvFunc[plane]));
                        //
                    } // for fracv
                } // for fracu
            } // for chan
        } // if !isDone
    } // for row

    ASKAPCHECK(itsSupport>0, "Support not calculated correctly");
    updateStats(nDone);
}

// To finalize the transform of the weights, we use the following steps:
// 1. For each plane of the convolution function, transform to image plane
// and multiply by conjugate to get abs value squared.
// 2. Sum all planes weighted by the weight for that convolution function.
void AProjectWStackVisGridder::finaliseWeights(casacore::Array<imtype>& out) {

    ASKAPTRACE("AProjectWStackVisGridder::finaliseWeights");
    ASKAPLOG_DEBUG_STR(logger, "Calculating sum of weights image");
    ASKAPDEBUGASSERT(itsShape.nelements()>=3);

    const int nx=itsShape(0);
    const int ny=itsShape(1);
    const int nPol=itsShape(2);
    const int nChan=itsShape(3);

    const int nZ = sumOfWeights().nrow();

    /// We must pad the convolution function to full size, reverse transform
    /// square, and sum multiplied by the corresponding weight
    const int cnx=std::min(itsMaxSupport, nx);
    const int cny=std::min(itsMaxSupport, ny);
    const int ccenx = cnx/2;
    const int cceny = cny/2;

    casacore::Vector<double> ccfx;
    casacore::Vector<double> ccfy;

    if (itsSpheroidalTaper) {
        // the spheroidal ccf needs to be removed from the weights before they are stored
        // it will also be removed from the image in correctConvolution
        ccfx.resize(cnx);
        ccfy.resize(cny);
        for (int ix=0; ix<cnx; ++ix) {
             const double nux = std::abs(double(ix)-double(ccenx))/double(ccenx);
             const double val = grdsf(nux);
             //ccfx(ix) = val; //casacore::abs(val) > 1e-10 ? 1./val : 0.;
             ccfx(ix) = val > itsSpheroidalWeightsCutoff ? val : 0.;
        }
        for (int iy=0; iy<cny; ++iy) {
             const double nuy = std::abs(double(iy)-double(cceny))/double(cceny);
             const double val = grdsf(nuy);
             //ccfy(iy) = val; //casacore::abs(val) > 1e-10 ? 1./val : 0.;
             ccfy(iy) = val > itsSpheroidalWeightsCutoff ? val : 0.;
        }
        // this isn't really needed unless itsSpheroidalWeightsCutoff==0. But shouldn't hurt
        if (doInterpolation()) {
          // The spheroidal is undefined and set to zero at nu=1, but that
          // is not the numerical limit. Estimate it from its neighbours.
          interpolateEdgeValues(ccfx);
          interpolateEdgeValues(ccfy);
        }
    }


    /// This is the output array before sinc padding
    casacore::Array<imtype> cOut(casacore::IPosition(4, cnx, cny, nPol, nChan));
    cOut.set(0.0);
    // Limit number of fft threads to 8 (more is slower for our fft sizes)
    scimath::FFT2DWrapper<imtypeComplex> fft2d(true,8);

    // for debugging
    double totSumWt = 0.;

    for (int iz=0; iz<nZ; ++iz) {
        const int plane = cfIndexFromSumOfWeightsRow(iz);

        bool hasData=false;
        for (int chan=0; chan<nChan; ++chan) {
            for (int pol=0; pol<nPol; ++pol) {
                const double wt = sumOfWeights()(iz, pol, chan);
                ASKAPCHECK(!std::isnan(wt), "sumOfWeights returns NaN for row="<<iz<<
                           " pol="<<pol<<" chan="<<chan);
                if (wt > 0.0) {
                    hasData=true;
                    totSumWt += wt;
                    //break;
                }
            }
        }
        if(hasData) {

            // Now fill the inner part of the uv plane with the convolution function
            // and transform to obtain the image. The uv sampling is fixed here
            // so the total field of view is itsOverSample times larger than the
            // original field of view.
            /// Work space
            casacore::Matrix<imtypeComplex> thisPlane(cnx, cny);
            thisPlane.set(0.0);
            const std::pair<int,int> cfOffset = getConvFuncOffset(iz);

            ASKAPDEBUGASSERT((int(itsConvFunc[plane].nrow()) - 1) / 2 == itsSupport);
            ASKAPDEBUGASSERT(itsConvFunc[plane].nrow() % 2 == 1);
            ASKAPDEBUGASSERT(itsConvFunc[plane].nrow() == itsConvFunc[plane].ncolumn());

            for (int iy=-itsSupport; iy<+itsSupport; ++iy) {
                for (int ix=-itsSupport; ix<+itsSupport; ++ix) {
                	 const int xPos = ix + ccenx + cfOffset.first;
	                 const int yPos = iy + cceny + cfOffset.second;
	                 if ((xPos<0) || (yPos<0) || (xPos>=int(thisPlane.nrow())) || (yPos>=int(thisPlane.ncolumn()))) {
	                     continue;
	                 }
                     thisPlane(xPos, yPos)=itsConvFunc[plane](ix+itsSupport, iy+itsSupport);
                }
            }

            //	  	  ASKAPLOG_DEBUG_STR(logger, "Convolution function["<< iz << "] peak = "<< peak);
            fft2d(thisPlane, false);
            thisPlane*=imtypeComplex(nx*ny);
            const double peak=real(casacore::max(casacore::abs(thisPlane)));
            // ASKAPLOG_DEBUG_STR(logger, "Transform of convolution function["<< iz << "] peak = "<< peak);

            if(peak>0.0) {
                thisPlane*=imtypeComplex(1.0/peak);
            }

            // Now we need to cut out only the part inside the field of view
            for (int chan=0; chan<nChan; ++chan) {
                for (int pol=0; pol<nPol; ++pol) {
                    casacore::IPosition ip(4, 0, 0, pol, chan);
                    const double wt = sumOfWeights()(iz, pol, chan);
                    ASKAPCHECK(!std::isnan(wt), "sumWeights returns NaN for row="<<iz<<
                               " pol="<<pol<<" chan="<<chan);

                    if (itsSpheroidalTaper) {
                       // if itsConvFunc has an additional taper, remove it before updating the weights.
                       for (int ix = 0; ix < cnx; ix++) {
                           ip(0) = ix;
                           for (int iy = 0; iy < cny; iy++) {
                               ip(1) = iy;
                               imtype taper = ccfx(ix) * ccfy(iy);
                               const imtypeComplex val = taper > 0 ? thisPlane(ix, iy) / taper : 0.;
                               cOut(ip) += double(wt) * casacore::real(val * conj(val));
                           }
                       }
                    } else {

                        for (int ix=0; ix<cnx; ++ix) {
                            ip(0)=ix;
                            for (int iy=0; iy<cny; ++iy) {
                                ip(1)=iy;
                                const imtypeComplex val = thisPlane(ix, iy);
                                cOut(ip) += double(wt) * casacore::real(val * conj(val));
                            }
                        }
                    }
                }
            }
        } // if has data
    } // loop over convolution functions
    scimath::PaddingUtils::fftPad(cOut, out, paddingFactor());
    ASKAPLOG_DEBUG_STR(logger,
            "Finished finalising the weights, the sum over all convolution functions is "<<totSumWt);
}

int AProjectWStackVisGridder::cIndex(int row, int pol, int chan) {
    return itsCMap(row, pol, chan);
}

void AProjectWStackVisGridder::correctConvolution(casacore::Array<imtype>& image) {
    // if there are any kernel factors that are not A or W terms, remove them now
    if (itsSpheroidalTaper) {
        SphFuncVisGridder::correctConvolution(image);
    }
}

/// @brief static method to create gridder
/// @details Each gridder should have a static factory method, which is
/// able to create a particular type of the gridder and initialise it with
/// the parameters taken form the given parset. It is assumed that the
/// method receives a subset of parameters where the gridder name is already
/// taken out.
/// @param[in] parset input parset file
/// @return a shared pointer to the gridder instance
IVisGridder::ShPtr AProjectWStackVisGridder::createGridder(const LOFAR::ParameterSet& parset)
{
  return createAProjectGridder<AProjectWStackVisGridder>(parset);
}

/// @brief assignment operator (not to be called)
/// @details It is defined as private, so we can't call it inadvertently
/// @param[in] other input object
/// @return reference to itself
AProjectWStackVisGridder& AProjectWStackVisGridder::operator=(const AProjectWStackVisGridder &)
{
  ASKAPTHROW(AskapError, "This method is not supposed to be called!");
  return *this;
}

}
}
