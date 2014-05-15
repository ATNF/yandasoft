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

#include <casa/BasicSL/Complex.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>
#include <measures/Measures/MDirection.h>
#include <casa/Quanta/MVDirection.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Quanta/MVTime.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>

ASKAP_LOGGER(logger, ".gridding.aprojectwstackgridder");

#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <profile/AskapProfiler.h>

#include <gridding/AProjectWStackVisGridder.h>
#include <fft/FFTWrapper.h>
#include <gridding/IBasicIllumination.h>

#include <utils/PaddingUtils.h>

namespace askap {
namespace synthesis {

AProjectWStackVisGridder::AProjectWStackVisGridder(const boost::shared_ptr<IBasicIllumination const> &illum,
        const double wmax, const int nwplanes, const double,
        const int overSample, const int maxSupport, const int limitSupport, 
        const int maxFeeds,
        const int maxFields, 
        const double pointingTol, const double paTol,
        const double freqTol,       
        const bool frequencyDependent, const std::string& name) :
    AProjectGridderBase(maxFeeds,maxFields,pointingTol, paTol, freqTol),
    WStackVisGridder(wmax, nwplanes), 
    itsReferenceFrequency(0.0),
    itsIllumination(illum),
    itsMaxFeeds(maxFeeds), itsMaxFields(maxFields),
    itsFreqDep(frequencyDependent)
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
    itsFreqDep(other.itsFreqDep), itsMaxSupport(other.itsMaxSupport),
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

    const casa::Vector<casa::RigidVector<double, 3> > &rotatedUVW = acc.rotatedUVW(getTangentPoint());

    for (int i=0; i<nSamples; ++i) {
        const int feed=acc.feed1()(i);
        ASKAPCHECK(feed<itsMaxFeeds, "Exceeded specified maximum number of feeds");
        ASKAPCHECK(feed>-1, "Illegal negative feed number");

        const double w=(rotatedUVW(i)(2))/(casa::C::c);

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
void AProjectWStackVisGridder::initialiseGrid(const scimath::Axes& axes,  const casa::IPosition& shape, const bool dopsf)
{
    ASKAPTRACE("AProjectWStackVisGridder::initialiseGrid");
    WStackVisGridder::initialiseGrid(axes,shape,dopsf);

    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casa::uInt nx=std::min(itsMaxSupport, int(itsShape(0)));
    const casa::uInt ny=std::min(itsMaxSupport, int(itsShape(1)));

    ASKAPLOG_INFO_STR(logger, "Shape for calculating gridding convolution function = "
            << nx << " by " << ny << " pixels");

    // this is just a buffer in the uv-space
    initUVPattern(nx,ny, itsUVCellSize(0),itsUVCellSize(1),itsOverSample);
}

/// @brief Initialise the degridding
/// @param axes axes specifications
/// @param image Input image: cube: u,v,pol,chan
void AProjectWStackVisGridder::initialiseDegrid(const scimath::Axes& axes,
        const casa::Array<double>& image)
{
    ASKAPTRACE("AProjectWStackVisGridder::initialiseDegrid");
    WStackVisGridder::initialiseDegrid(axes,image);      
    /// Limit the size of the convolution function since
    /// we don't need it finely sampled in image space. This
    /// will reduce the time taken to calculate it.
    const casa::uInt nx=std::min(itsMaxSupport, int(itsShape(0)));
    const casa::uInt ny=std::min(itsMaxSupport, int(itsShape(1)));

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
    const int nSamples = acc.nRow();

    validateCFCache(acc, hasSymmetricIllumination);
        
    casa::MVDirection out = getImageCentre();

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
    const casa::uInt nx = pattern.uSize();
    const casa::uInt ny = pattern.vSize();


    int nDone=0;
    for (int row=0; row<nSamples; ++row) {
        const int feed=acc.feed1()(row);

        if (!isCFValid(feed, currentField())) {
            makeCFValid(feed, currentField());
            nDone++;
            casa::MVDirection offset(acc.pointingDir1()(row).getAngle());
            rwSlopes()(0, feed, currentField()) = isPSFGridder() ? 0. : sin(offset.getLong()
                    -out.getLong()) *cos(offset.getLat());
            rwSlopes()(1, feed, currentField())= isPSFGridder() ? 0. : sin(offset.getLat())
                *cos(out.getLat()) - cos(offset.getLat())*sin(out.getLat())
                *cos(offset.getLong()-out.getLong());

            const double parallacticAngle = hasSymmetricIllumination ? 0. : acc.feed1PA()(row);

            for (int chan=0; chan<nChan; chan++) {
                /// Extract illumination pattern for this channel
                itsIllumination->getPattern(acc.frequency()[chan], pattern,
                        rwSlopes()(0, feed, currentField()),
                        rwSlopes()(1, feed, currentField()), 
                        parallacticAngle);

                /// Now convolve the disk with itself using an FFT
                scimath::fft2d(pattern.pattern(), false);

                double peak=0.0;
                for (casa::uInt ix=0; ix<nx; ++ix) {
                    for (casa::uInt iy=0; iy<ny; ++iy) {
                        pattern(ix, iy)=pattern(ix,iy)*conj(pattern(ix,iy));
                        if (casa::abs(pattern(ix,iy))>peak) {
                            peak=casa::abs(pattern(ix,iy));
                        }
                    }
                }
                if(peak>0.0) {
                    pattern.pattern()*=casa::DComplex(1.0/peak);
                }
                // The maximum will be 1.0
                ASKAPLOG_DEBUG_STR(logger, "Max of FT of convolution function = " << casa::max(pattern.pattern()));
                scimath::fft2d(pattern.pattern(), true);	
                // Now correct for normalization of FFT
                pattern.pattern()*=casa::DComplex(1.0/(double(nx)*double(ny)));
                ASKAPLOG_DEBUG_STR(logger, "Sum of convolution function before support extraction and decimation = " << casa::sum(pattern.pattern()));

                if (itsSupport==0) {
                    // we probably need a proper support search here
                    // it can be encapsulated in a method of the UVPattern class
                    itsSupport = pattern.maxSupport();
                    ASKAPCHECK(itsSupport>0,
                            "Unable to determine support of convolution function");
                    ASKAPCHECK(itsSupport*itsOverSample<int(nx)/2,
                            "Overflowing convolution function - increase maxSupport or decrease overSample");
                    if (itsLimitSupport > 0  &&  itsSupport > itsLimitSupport) {
                        ASKAPLOG_INFO_STR(logger, "Convolution function support = "
                                << itsSupport << " pixels exceeds upper support limit; "
                                << "set to limit = " << itsLimitSupport << " pixels");
                        itsSupport = itsLimitSupport;
                    }

                    const int cSize=2*itsSupport+1;
                    // just for logging
                    const double cell = std::abs(pattern.uCellSize())*(casa::C::c/acc.frequency()[chan]);
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
                const double rescale=double(itsOverSample*itsOverSample);
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
                                    = rescale * pattern(itsOverSample*ix+fracu+nx/2,
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
void AProjectWStackVisGridder::finaliseWeights(casa::Array<double>& out) {

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

    /// This is the output array before sinc padding
    casa::Array<double> cOut(casa::IPosition(4, cnx, cny, nPol, nChan));
    cOut.set(0.0);

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
            casa::Matrix<casa::DComplex> thisPlane(cnx, cny);
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
            scimath::fft2d(thisPlane, false);
            thisPlane*=casa::DComplex(nx*ny);
            const double peak=real(casa::max(casa::abs(thisPlane)));
            // ASKAPLOG_DEBUG_STR(logger, "Transform of convolution function["<< iz << "] peak = "<< peak);

            if(peak>0.0) {
                thisPlane*=casa::DComplex(1.0/peak);
            }

            // Now we need to cut out only the part inside the field of view
            for (int chan=0; chan<nChan; ++chan) {
                for (int pol=0; pol<nPol; ++pol) {
                    casa::IPosition ip(4, 0, 0, pol, chan);
                    const double wt = sumOfWeights()(iz, pol, chan);
                    ASKAPCHECK(!std::isnan(wt), "sumWeights returns NaN for row="<<iz<<
                               " pol="<<pol<<" chan="<<chan);
                    for (int ix=0; ix<cnx; ++ix) {
                        ip(0)=ix;
                        for (int iy=0; iy<cny; ++iy) {
                            ip(1)=iy;
                            cOut(ip)+=wt*casa::real(thisPlane(ix, iy)*conj(thisPlane(ix, iy)));
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

void AProjectWStackVisGridder::correctConvolution(casa::Array<double>& /*grid*/) {
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
