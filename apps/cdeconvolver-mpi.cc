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

// ASKAPsoft includes
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".cdeconvolver");
#include <askap/AskapError.h>
#include <askap/Application.h>
#include <askap/StatReporter.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/deconvolution/DeconvolverBase.h>
#include <askap/deconvolution/DeconvolverFactory.h>
#include <askap/deconvolution/DeconvolverHelpers.h>
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



using namespace askap;
using namespace askap::synthesis;


class CdeconvolverApp : public askap::Application
{
    public:
    
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsModelCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsResidualCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsRestoredCube;
    
        boost::shared_ptr<IVisGridder> itsGridder;
    
        static void correctConvolution(casacore::Array<casacore::Double>&, int, int, bool);
    
       
        
        virtual int run(int argc, char* argv[])
        {
            StatReporter stats;

            const LOFAR::ParameterSet subset(config().makeSubset("Cdeconvolver."));
            
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));
            
            ASKAPLOG_INFO_STR(logger, "ASKAP image (MPI) deconvolver " << ASKAP_PACKAGE_VERSION);
            // Need some metadata for the output cube constructions
            
            // Lets get the grid,pcf and psf cube names from the parset
            
            const std::string gridCubeName = subset.getString("grid");
            const std::string pcfCubeName = subset.getString("pcf");
            const std::string psfCubeName = subset.getString("psf");
            
            // ok lets set up some output cubes
            const std::string outModelCubeName = subset.getString("model","model");
            const std::string outResidCubeName = subset.getString("residual","residual");
            const std::string outRestoredCubeName = subset.getString("restored","restored");
            
           
            
            // Lets load in a cube
            casacore::PagedImage<casacore::Complex> grid(gridCubeName);
       
            const casa::IPosition shape = grid.shape();
            casa::IPosition blc(shape.nelements(),0);
            casa::IPosition trc(shape);
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
                itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, cdelt,outModelCubeName));
                itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, cdelt,outResidCubeName));
                itsRestoredCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, cdelt,outRestoredCubeName));
                
            }
            
              
            
            
            // lets calculate the allocations ...

           
            
            vector<IPosition> inShapeVec;
            vector<CoordinateSystem> inCoordSysVec;
            // What fraction of the full problem does a rank have
            // THe units here are channels - not polarisations.
            int myFullAllocationSize = 0;
            int myFullAllocationStart = 0;
            int myFullAllocationStop = 0;
            // Where a rank is in its allocation
            int myAllocationSize = 1;
            int myAllocationStart = 0;
            int myAllocationStop = myAllocationStart + myAllocationSize;
            
            if (nchanCube % comms.nProcs() != 0) {
                    ASKAPLOG_WARN_STR(logger,"Unbalanced allocation: num of ranks:" << comms.nProcs() << " not a factor of number of channels: "<< nchanCube);
            }
            
            if (comms.nProcs() >= nchanCube) {
                myFullAllocationSize = 1;
            }
            else {
                myFullAllocationSize = nchanCube/comms.nProcs();
            }
                
                    
            myFullAllocationStart = comms.rank()*myFullAllocationSize;
            myFullAllocationStop = myFullAllocationStart + myFullAllocationSize;

            // unless last rank
            if (comms.rank() == comms.nProcs()-1) {
                myFullAllocationSize = trc[3] - myFullAllocationStart; // we are using End is Last
            }
            
            ASKAPLOG_INFO_STR(logger,"RankAllocation starts at " << myFullAllocationStart << " and is " << myFullAllocationSize << " in size");
            
            int support = 3;
            int alpha = 1;
              // this is private to an inherited class so have to make a new one
                    
            for (myAllocationStart = myFullAllocationStart; myAllocationStart < myFullAllocationStop; myAllocationStart = myAllocationStart +  myAllocationSize) {
                    
                ASKAPLOG_INFO_STR(logger,"Input image shape " << shape);
                ASKAPLOG_INFO_STR(logger,"Processing Channel " << myAllocationStart);
                    
                casa::IPosition inblc(shape.nelements(),0); // input bottom left corner of this allocation
                casa::IPosition intrc(shape); // get the top right
                myAllocationStop = myAllocationStart + myAllocationSize;
                
                inblc[3] = myAllocationStart;
                intrc[0] = intrc[0]-1;
                intrc[1] = intrc[1]-1;
                intrc[2] = intrc[2]-1;
                intrc[3] = myAllocationStart + myAllocationSize-1;
                    
                const casacore::Slicer slicer(inblc, intrc, casacore::Slicer::endIsLast);
                ASKAPLOG_INFO_STR(logger,"Slicer is " << slicer);
                    
                casacore::Array<casacore::Complex> buffer;
                grid.doGetSlice(buffer, slicer);
                
                askap::scimath::fft2d(buffer,false);
                casacore::Array<casacore::Double> dBuffer(buffer.shape());
                casacore::Array<casacore::Float> fBuffer(buffer.shape());
                casacore::convertArray<casacore::Double, casacore::Float> (dBuffer,real(buffer));
                // template this
                correctConvolution(dBuffer,support,alpha,true);
                casacore::convertArray<casacore::Float, casacore::Double>(fBuffer,dBuffer);
                itsModelCube->writeSlice(fBuffer,myAllocationStart);
                
            
                
            }

            // output Model cube
            // this->itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, freqinc, "model"));
            // output Residual cube
            // this->itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset, nchanCube, f0, freqinc,"residual"));
            // boost::shared_ptr<DeconvolverBase<casacore::Float, casacore::Complex> > deconvolver::make(subset));
            
            // deconvolver->deconvolve();

            // Now write the model and residual to disk using the names specified in the 
            // parset. We simply copy the dirty image and then write the array into 
            // the resulting image. 
            // DeconvolverHelpers::putArrayToImage(deconvolver->model(), "model", "dirty", subset);
            // DeconvolverHelpers::putArrayToImage(deconvolver->dirty(), "residual", "dirty", subset);

            // Vector<Array<float> > restored(1);
            // if(deconvolver->restore(restored)) {
            //    DeconvolverHelpers::putArrayToImage(restored(0), "restored", "dirty", subset);
            // }

            stats.logSummary();
            return 0;
        }
};
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
/*
    if (itsInterp) {
    // The spheroidal is undefined and set to zero at nu=1, but that
    // is not the numerical limit. Estimate it from its neighbours.
        interpolateEdgeValues(bufx);
        interpolateEdgeValues(bufy);
    }
*/
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
// Main function
int main(int argc, char* argv[])
{
    CdeconvolverApp app;
    return app.main(argc, argv);
}
