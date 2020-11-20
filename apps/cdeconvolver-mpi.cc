/// @file cdeconvolver.cc
///
/// @brief Image deconvolution program
///
/// Control parameters are passed in from a LOFAR ParameterSet file.
///
/// @copyright (c) 2007, 2020 CSIRO
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

// Package level header file
#include <askap/askap_synthesis.h>

// System includes
#include <stdexcept>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".cdeconvolver");
#include <askap/askap/AskapError.h>
#include <askap/askap/Application.h>
#include <askap/askap/StatReporter.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/WienerPreconditioner.h>
#include <askap/deconvolution/DeconvolverBase.h>
#include <askap/deconvolution/DeconvolverFactory.h>
#include <askap/deconvolution/DeconvolverHelpers.h>
#include <askap/deconvolution/DeconvolverBasisFunction.h>
#include <askap/deconvolution/DeconvolverMultiTermBasisFunction.h>
#include <askap/deconvolution/DeconvolverHogbom.h>
#include <askap/distributedimager/CubeBuilder.h>
#include <askap/askapparallel/AskapParallel.h>
#include <askap/scimath/utils/MultiDimArrayPlaneIter.h>
#include <askap/scimath/fft/FFTWrapper.h>
#include <askap/scimath/utils/SpheroidalFunction.h>
#include <askap/gridding/VisGridderFactory.h>
#include <askap/gridding/SphFuncVisGridder.h>

#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/ImageSummary.h>

#include <Blob/BlobArray.h>
#include <Blob/BlobSTL.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobIBufVector.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobOBufVector.h>

using namespace askap;
using namespace askap::synthesis;



class CdeconvolverApp : public askap::Application
{
    public:
    
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsPsfCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsResidualCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsModelCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsRestoredCube;
    
    
        static void correctConvolution(casacore::Array<casacore::Double>&, int, int, bool);

        const bool itsInterp = true;
        static void interpolateEdgeValues(casacore::Vector<casacore::DComplex> &func);
    
        void getPSF(casacore::Array<casacore::Float> &psfArray,casacore::Array<casacore::Complex> &buffer );
    
        void doTheWork(const LOFAR::ParameterSet, casacore::Array<casacore::Complex> &buffer, casacore::Array<casacore::Float> &psfArray,casacore::Array<casacore::Complex> &pcfArray, casacore::Array<casacore::Float> &outpsf, casacore::Array<casacore::Float> &dirty, casacore::Array<casacore::Float> &model, casacore::Array<casacore::Float> &restored);

        std::pair<int, int> get_channel_allocation(askap::askapparallel::AskapParallel &comms, int nchannels)
        {
            auto rank = comms.rank();
            auto nranks = comms.nProcs();
            auto div = nchannels / nranks;
            auto rem = nchannels % nranks;

            if (rem > 0) {
                ASKAPLOG_WARN_STR(logger,"Unbalanced allocation: num of ranks:" << nranks << " not a factor of number of channels: "<< nchannels);
            }
            // Simple round-robin: the first `rem` ranks receive an extra item
            // when rem > 0. That means that:
            //  if rank < rem:  first_chan = (div + 1) * rank
            //                  num_chans = div + 1
            //  if rank >= rem: first_chan = (div + 1) * rem + (div * (rank - rem))
            //                  num_chans = div
            // and that reduces to what's below
            auto first_chan = rank * div + (rank < rem ? div : rem);
            auto num_chans = div + (rank < rem);
            return std::make_pair(first_chan, num_chans);
        }

        virtual int run(int argc, char* argv[])
        {
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));
            try {
                return _run(argc, argv, comms);
            } catch (const std::exception &e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected error: " << e.what());
                comms.abort();
                return 1;
            }
        }

        int _run(int argc, char *argv[], askap::askapparallel::AskapParallel &comms)
        {
            StatReporter stats;
           

            const LOFAR::ParameterSet subset(config().makeSubset("Cdeconvolver."));
           
            ASKAPLOG_INFO_STR(logger, "ASKAP image (MPI) deconvolver " << ASKAP_PACKAGE_VERSION);
            
            // ASKAPCHECK(comms.nProcs() == 1,"Currently only SERIAL mode supported");
                
            
            // Need some metadata for the output cube constructions
            
            // Lets get the grid,pcf and psf cube names from the parset
            
            const std::string gridCubeName = subset.getString("grid");
            const std::string pcfCubeName = subset.getString("pcf");
            const std::string psfGridCubeName = subset.getString("psfgrid");
            
            // ok lets set up some output cubes
            const std::string outPsfCubeName = subset.getString("psf","psf");
            const std::string outResidCubeName = subset.getString("residual","residual");
            const std::string outModelCubeName = subset.getString("model","model");
            const std::string outRestoredCubeName = subset.getString("restored","restored");
            
            
           
            
            // WorkArrays
            casacore::Array<casacore::Float> psfArray;
            casacore::Array<casacore::Complex> pcfArray;
            casacore::Array<casacore::Complex> psfGridBuffer;
            casacore::Array<casacore::Complex> buffer;
            
            
            
            // Lets load in a cube
            casacore::PagedImage<casacore::Complex> grid(gridCubeName);
            casacore::PagedImage<casacore::Complex> pcf(pcfCubeName);
            casacore::PagedImage<casacore::Complex> psfgrid(psfGridCubeName);
            
            const casacore::IPosition shape = grid.shape();
            casacore::IPosition blc(shape.nelements(),0);
            casacore::IPosition trc(shape);
            int nchanCube = trc[3];
            
            if (comms.isMaster()) { // only the master makes the output
            
                
                Int pixelAxis,worldAxis,coordinate;
                CoordinateUtil::findSpectralAxis(pixelAxis,worldAxis,coordinate,grid.coordinates());
                const SpectralCoordinate &sc = grid.coordinates().spectralCoordinate(coordinate);
                casacore::Double baseFreq, nextFreq, freqInc;
                casacore::Double pixelVal=0;
                
                sc.toWorld(baseFreq,pixelVal);
                sc.toWorld(nextFreq,pixelVal+1);
                
                freqInc = nextFreq-baseFreq;
                
                Quantity f0(baseFreq, "Hz");
                Quantity cdelt(freqInc, "Hz");
                
                ASKAPLOG_INFO_STR(logger,"Base Freq " << f0);
                ASKAPLOG_INFO_STR(logger,"Freq inc (CDELT) " << cdelt);
                
                // create the output cubes
                itsPsfCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, cdelt,outPsfCubeName));
                itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, cdelt,outResidCubeName));
                itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, cdelt,outModelCubeName));
                itsRestoredCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, cdelt,outRestoredCubeName));
            }
            else {
                // this should work fine as the cubes will exist by the time
                // they are needed.
                itsPsfCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,outPsfCubeName));
                itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,outResidCubeName));
                itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,outModelCubeName));
                itsRestoredCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,outRestoredCubeName));
            }

            // What fraction of the full problem does a rank have
            int firstChannel, numChannelsLocal;
            std::tie(firstChannel, numChannelsLocal) = get_channel_allocation(comms, nchanCube);
            ASKAPLOG_INFO_STR(logger,"Rank " << comms.rank() << " - RankAllocation starts at " << firstChannel << " and is " << numChannelsLocal << " in size");
            bool firstPassForMaster = true;

            for (int channel = firstChannel; channel < firstChannel + numChannelsLocal; channel++) {

                //FIXME: this is just looping over each channel of the allocation
                
                ASKAPLOG_INFO_STR(logger,"Input image shape " << shape);
                ASKAPLOG_INFO_STR(logger,"Processing Channel " << channel);
                
                casacore::IPosition inblc(shape.nelements(),0); // input bottom left corner of this allocation
                casacore::IPosition intrc(shape); // get the top right
                
                inblc[3] = channel;
                intrc[0] = intrc[0]-1;
                intrc[1] = intrc[1]-1;
                intrc[2] = intrc[2]-1;
                intrc[3] = channel;
                
                const casacore::Slicer slicer(inblc, intrc, casacore::Slicer::endIsLast);
                ASKAPLOG_INFO_STR(logger,"Slicer is " << slicer);
                
                
                
                psfgrid.getSlice(psfGridBuffer,slicer);
                pcf.getSlice(pcfArray,slicer);
                grid.doGetSlice(buffer, slicer);
     
                // do the work
                scimath::MultiDimArrayPlaneIter planeIter(buffer.shape());
                
                for ( ; planeIter.hasMore(); planeIter.next()) {
                    /// FIXME: this is supposed to loop over the polarisations as well as channels
                    /// FIXME: but i have not sorted out the output indexes for this to work
                    
                    
                    casacore::IPosition curpos = planeIter.position();
                    ASKAPLOG_INFO_STR(logger, "Processing from position: " << curpos);
                    // the inputs
                    casacore::Array<casacore::Complex> thisBuffer = planeIter.getPlane(buffer, curpos);
                    casacore::Array<casacore::Complex> thisPCFBuffer = planeIter.getPlane(pcfArray, curpos);
                    casacore::Array<casacore::Complex> thisPSFGridBuffer = planeIter.getPlane(psfGridBuffer, curpos);
                    casacore::Array<casacore::Float> thisPSFBuffer(psfGridBuffer.shape());
                    
                    // the outputs
                    casacore::Array<casacore::Float> psfout;
                    casacore::Array<casacore::Float> dirty;
                    casacore::Array<casacore::Float> model;
                    casacore::Array<casacore::Float> restored;
                    
                    getPSF(thisPSFBuffer,thisPSFGridBuffer);

                    doTheWork(subset, thisBuffer, thisPSFBuffer, thisPCFBuffer, psfout, dirty, model, restored);
                    
                    if (comms.isMaster() && firstPassForMaster == true) {

                      ASKAPLOG_INFO_STR(logger, "Ensuring serial access to cubes");
                      firstPassForMaster = false; 

                    }
                    else { // this is essentially a serializer - it is required for CASA image types
                    // but not FITS
                      int buf;
                      int from = comms.rank() - 1;
                      if (from < 0) {
                            from = comms.nProcs();
                      }
                      ASKAPLOG_INFO_STR(logger, "Waiting for trigger from rank " << from);
                      comms.receive((void *) &buf,sizeof(int),from);
                    }
                    // write out the slice
                    // FIXME: THis is the issue with npol I need to use some position
                    
                    
                    
                    
                    itsPsfCube->writeSlice(psfout, channel);
                    itsResidualCube->writeSlice(dirty, channel);
                    itsModelCube->writeSlice(model, channel);
                    itsRestoredCube->writeSlice(restored, channel);

                    if (comms.rank() < comms.nProcs()-1) { // last rank doesnot use this method
                      int buf;
                      int to = comms.rank()+1;
                      comms.send((void *) &buf,sizeof(int),to);
                    }
                    else {
                      int buf;
                      int to = 0;
                      comms.send((void *) &buf,sizeof(int),to);
                    }
                    
                }
                
                
            }

            stats.logSummary();
            comms.barrier();
            return 0;
        }
};
            
void CdeconvolverApp::getPSF(casacore::Array<casacore::Float> &psfArray,casacore::Array<casacore::Complex> &buffer ) {
   
    casacore::Array<double> dBuffer(buffer.shape());

    //askap::scimath::fft2d(buffer,false);
    //casacore::convertArray<casacore::Float, casacore::Float>(psfArray,real(buffer));

    casacore::Array<casacore::DComplex> scratch(buffer.shape());
    casacore::convertArray<casacore::DComplex,casacore::Complex>(scratch, buffer);
    askap::scimath::fft2d(scratch, false);
    casacore::convertArray<casacore::Float, casacore::Double>(psfArray,real(scratch));
    
    //psfArray = psfArray * psfArray.nelements();
    ASKAPLOG_INFO_STR(logger,"Max PSF array:" << max(psfArray));
    
}

void CdeconvolverApp::doTheWork(const LOFAR::ParameterSet subset, casacore::Array<casacore::Complex> &buffer, casacore::Array<casacore::Float> &psfArray,casacore::Array<casacore::Complex> &pcfArray, casacore::Array<casacore::Float> &outpsf, casacore::Array<casacore::Float> &dirty, casacore::Array<casacore::Float> &model, casacore::Array<casacore::Float> &restored) {
    
    ASKAPLOG_INFO_STR(logger,"Array Shape: " << buffer.shape());
    //double norm = buffer.nelements()/max(psfArray);
    double norm = buffer.nelements();
                
    //askap::scimath::fft2d(buffer,false);
    casacore::Array<casacore::DComplex> scratch(buffer.shape());
    casacore::convertArray<casacore::DComplex,casacore::Complex>(scratch, buffer);
    askap::scimath::fft2d(scratch, false);
    casacore::convertArray<casacore::Complex, casacore::DComplex>(buffer,scratch);

    //buffer = buffer * norm;
    //ASKAPLOG_INFO_STR(logger,"Normalisation Factor is " << norm);
    
    // Convolution correction probably should pull the support from the PARSET
    // alpha is usually 1
    // support is ususally 3
    int alpha = 1;
    int support = 3;
    casacore::Array<casacore::Double> dBuffer(buffer.shape());
    casacore::Array<casacore::Float> fBuffer(buffer.shape());

    casacore::convertArray<casacore::Double, casacore::Float> (dBuffer,real(buffer));
    correctConvolution(dBuffer,support,alpha,true);
    dBuffer = dBuffer * double(norm);
    casacore::convertArray<casacore::Float, casacore::Double>(fBuffer,dBuffer);

    casacore::convertArray<casacore::Double, casacore::Float> (dBuffer,psfArray);
    correctConvolution(dBuffer,support,alpha,true);
    dBuffer = dBuffer * double(norm);
    casacore::convertArray<casacore::Float, casacore::Double>(psfArray,dBuffer);
    
    // *** Preconditioning ***
    if (subset.isDefined("preconditioner.Names")) {
        ASKAPLOG_INFO_STR(logger,"Preparing for preconditioning");

        // Preconditioning Assuming they only want Wiener.
        boost::shared_ptr<WienerPreconditioner>
            wp = WienerPreconditioner::createPreconditioner(subset.makeSubset("preconditioner.Wiener."));
        ASKAPASSERT(wp);
        
        // The preconditioner assumes that the PCF is accumulated in the image domain so FFT it.
        // I'm not normalising it as I dont think I need to.
        
        //askap::scimath::fft2d(pcfArray,false);
        casacore::convertArray<casacore::DComplex,casacore::Complex>(scratch, pcfArray);
        askap::scimath::fft2d(scratch, false);
        casacore::convertArray<casacore::Complex, casacore::DComplex>(pcfArray,scratch);

        casacore::Array<casacore::Float> pcfReal = real(pcfArray) * norm;

        wp->doPreconditioning(psfArray,fBuffer,pcfReal);
        //casacore::Array<casacore::Float> psfArrayUnfiltered = psfArray;
        //wp->doPreconditioning(psfArrayUnfiltered,psfArray,pcfReal);
        //wp->doPreconditioning(psfArrayUnfiltered,fBuffer,pcfReal);

    } else {
        ASKAPLOG_INFO_STR(logger,"No preconditioning");
    }

    // *** Normalisation ***
    const float unnormalisedMaxPSF = max(psfArray);
    if (unnormalisedMaxPSF<=0.) {
        ASKAPTHROW(AskapError, "PSF Error. Peak is " << unnormalisedMaxPSF);
    }

    ASKAPLOG_INFO_STR(logger, "Normalising PSF");
    ASKAPLOG_INFO_STR(logger, "Peak of PSF before normalisation = " << unnormalisedMaxPSF);
    psfArray /= unnormalisedMaxPSF;
    ASKAPLOG_INFO_STR(logger, "Peak of PSF after normalisation = " << max(psfArray));
    ASKAPLOG_INFO_STR(logger, "Normalising Dirty image");
    fBuffer /= unnormalisedMaxPSF;

    // *** Deconvolution ***
    // We have a normalised corrected preconditioned image
    DeconvolverBase<Float, Complex>::ShPtr deconvolver;
    string algorithm = subset.getString("solver.Clean.algorithm", "Basisfunction");
    
    if (algorithm == "Basisfunction") {
        ASKAPLOG_INFO_STR(logger, "Constructing Basisfunction Clean solver");
        deconvolver.reset(new DeconvolverBasisFunction<Float, Complex>(fBuffer, psfArray));
        ASKAPASSERT(deconvolver);
    } else if (algorithm == "MultiTermBasisfunction") {
        ASKAPLOG_INFO_STR(logger, "Constructing MultiTermBasisfunction Clean solver");
        deconvolver.reset(new DeconvolverMultiTermBasisFunction<Float, Complex>(fBuffer, psfArray));
        ASKAPASSERT(deconvolver);
    } else if (algorithm == "Hogbom") {
        ASKAPLOG_INFO_STR(logger, "Constructing Hogbom Clean deconvolver");
        deconvolver.reset(new DeconvolverHogbom<Float, Complex>(fBuffer, psfArray));
        ASKAPASSERT(deconvolver);
    } else {
        ASKAPTHROW(AskapError, "Unknown Clean algorithm " << algorithm);
    }
    
    deconvolver->configure(subset.makeSubset("solver.Clean."));
    
    deconvolver->deconvolve();
    
    dirty = deconvolver->dirty().copy();
    model = deconvolver->model().copy();

    outpsf = psfArray;
    
    casacore::Vector< casacore::Array<casacore::Float> > restored_vec(1); // some times there are more that one restore
    if(deconvolver->restore(restored_vec))
        restored = restored_vec(0).copy();
    
}

void CdeconvolverApp::correctConvolution(casacore::Array<double>& grid, int support=3, int alpha = 1, bool itsInterp = true)
{
    casacore::IPosition itsShape = grid.shape();
    ASKAPDEBUGASSERT(itsShape.nelements()>=2);
    const casacore::Int xHalfSize = itsShape(0)/2;
    const casacore::Int yHalfSize = itsShape(1)/2;
    casacore::Vector<double> ccfx(itsShape(0));
    casacore::Vector<double> ccfy(itsShape(1));
    ASKAPDEBUGASSERT(itsShape(0)>1);
    ASKAPDEBUGASSERT(itsShape(1)>1);

    scimath::SpheroidalFunction itsSphFunc(casacore::C::pi*support, alpha);
  
    // initialise buffers to enable a filtering of the correction
    // function in Fourier space.
    casacore::Vector<casacore::DComplex> bufx(itsShape(0));
    casacore::Vector<casacore::DComplex> bufy(itsShape(1));

    // note grdsf(1)=0.
    for (int ix=0; ix<itsShape(0); ++ix)
    {
        const double nux=std::abs(double(ix-xHalfSize))/double(xHalfSize);
        const double val = itsSphFunc(nux);
        bufx(ix) = casacore::DComplex(val,0.0);
    }

    for (int iy=0; iy<itsShape(1); ++iy)
    {
        const double nuy=std::abs(double(iy-yHalfSize))/double(yHalfSize);
        const double val = itsSphFunc(nuy);
        bufy(iy) = casacore::DComplex(val,0.0);
    }
    if (itsInterp) {
    // The spheroidal is undefined and set to zero at nu=1, but that
    // is not the numerical limit. Estimate it from its neighbours.
        interpolateEdgeValues(bufx);
        interpolateEdgeValues(bufy);
    }
    // Fourier filter the spheroidal (crop in Fourier space in line with
    // gridding kernel support size)
    const bool doFiltering = true;
    if (doFiltering) {
     // Some more advanced gridders have support>3 (e.g. w-proj).
     //
        int support = 3;
        const casacore::DComplex maxBefore = bufx(itsShape(0)/2);
        scimath::fft(bufx, true);
        scimath::fft(bufy, true);
        for (int ix=0; ix<itsShape(0)/2-support; ++ix) {
            bufx(ix) = 0.0;
        }
        for (int ix=itsShape(0)/2+support+1; ix<itsShape(0); ++ix) {
            bufx(ix) = 0.0;
        }
        for (int iy=0; iy<itsShape(1)/2-support; ++iy) {
            bufy(iy) = 0.0;
        }
        for (int iy=itsShape(1)/2+support+1; iy<itsShape(1); ++iy) {
            bufy(iy) = 0.0;
        }
        scimath::fft(bufx, false);
        scimath::fft(bufy, false);
        // Normalise after filtering.
        const casacore::DComplex normalisation = maxBefore / bufx(itsShape(0)/2);
        bufx *= normalisation;
        bufy *= normalisation;
    }

    for (int ix=0; ix<itsShape(0); ++ix) {
        double val = real(bufx(ix));
        ccfx(ix) = casacore::abs(val) > 1e-10 ? 1.0/val : 0.;
    }
    for (int iy=0; iy<itsShape(1); ++iy) {
        double val = real(bufy(iy));
        ccfy(iy) = casacore::abs(val) > 1e-10 ? 1.0/val : 0.;
    }

    casacore::ArrayIterator<double> it(grid, 2);
    while (!it.pastEnd())
    {
        casacore::Matrix<double> mat(it.array());
        ASKAPDEBUGASSERT(int(mat.nrow()) <= itsShape(0));
        ASKAPDEBUGASSERT(int(mat.ncolumn()) <= itsShape(1));
        for (int ix=0; ix<itsShape(0); ix++)
        {
            for (int iy=0; iy<itsShape(1); iy++)
            {
                mat(ix, iy)*=ccfx(ix)*ccfy(iy);
            }
        }
        it.next();
    }

    
}

/// @brief estimate the spheroidal function at nu=1
/// @param[in] func vector to be interpolated
void CdeconvolverApp::interpolateEdgeValues(casacore::Vector<casacore::DComplex> &func)
{
  ASKAPDEBUGASSERT(itsInterp);

  int length = func.shape()[0];
  ASKAPASSERT(length>3);

  //func(0) = func(1) + (func(1)-func(2)) + (func(1)-2.0*func(2)+func(3));
  func(0) = 3. * (func(1) - func(2)) + func(3);
  if (length%2==1) { // nu=1 for the last element as well
    func(length-1) = func(0);
  }

}

// Main function
int main(int argc, char* argv[])
{
    CdeconvolverApp app;
    return app.main(argc, argv);
}
