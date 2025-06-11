/// @file cdeconvolver-mpi.cc
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
#include <utility>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/Application.h>
#include <askap/StatReporter.h>
#include <askap/utils/StatsAndMask.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/measurementequation/WienerPreconditioner.h>
#include <askap/measurementequation/GaussianTaperPreconditioner.h>
#include <askap/deconvolution/DeconvolverFactory.h>
#include <askap/deconvolution/DeconvolverBasisFunction.h>
#include <askap/deconvolution/DeconvolverMultiTermBasisFunction.h>
#include <askap/deconvolution/DeconvolverHogbom.h>
#include <askap/distributedimager/CubeBuilder.h>
#include <askap/imageaccess/BeamLogger.h>
#include <askap/imageaccess/WeightsLog.h>
#include <askap/askapparallel/AskapParallel.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>
#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/scimath/utils/SpheroidalFunction.h>
#include <askap/gridding/SphFuncVisGridder.h>

#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/images/Images/ImageOpener.h>


ASKAP_LOGGER(logger, ".cdeconvolver");

using namespace askap;
using namespace askap::synthesis;



class CdeconvolverApp : public askap::Application
{
    public:

        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsPsfCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsResidualCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsModelCube;
        boost::shared_ptr<askap::cp::CubeBuilder<casacore::Float> > itsRestoredCube;
        BeamList itsBeamList;
        /// @brief The list of weights information. Each element of the map is a float,
        /// referenced by the channel number.
        std::map<unsigned int, float> itsWeightsList;
        int itsBeamReferenceChannel;
        LOFAR::ParameterSet itsParset;

        void getRealFFT(casacore::Matrix<casacore::Float> &fArray, casacore::Matrix<casacore::Complex> &cArray);

        // Precondition and deconvolve the inputs to produce the outputs, note inputs are modified (used as scratch)
        void doTheWork(const LOFAR::ParameterSet subset,
                       casacore::Array<casacore::Float> &dirtyIn,
                       casacore::Array<casacore::Float> &psfIn,
                       casacore::Array<casacore::Float> &pcfIn,
                       casacore::Array<casacore::Float> &outpsf,
                       casacore::Array<casacore::Float> &dirty,
                       casacore::Array<casacore::Float> &model,
                       casacore::Array<casacore::Float> &restored,
                       bool writeRestored,
                       float oversampling,
                       casacore::Vector<casacore::Double>& fov,
                       casacore::Vector<casacore::Quantum<double>> & beam);

        void initialiseBeamList(const unsigned int numChannels);
        void writeBeamInfo(askap::askapparallel::AskapParallel &comms);
        void initialiseWeightsList(const unsigned int numChannels);
        void writeWeightsInfo(askap::askapparallel::AskapParallel &comms);

        std::pair<int, int> get_channel_allocation(askap::askapparallel::AskapParallel &comms, int nchannels)
        {
            auto rank = comms.rank();
            auto nranks = comms.nProcs();
            auto div = nchannels / nranks;
            auto rem = nchannels % nranks;

            if (rem > 0) {
                ASKAPLOG_WARN_STR(logger,"Unbalanced allocation: num of ranks:" << nranks <<
                                         " not a factor of number of channels: "<< nchannels);
            }
            // Simple round-robin: the first `rem` ranks receive an extra item
            // when rem > 0. That means that:
            //  if rank < rem:  first_chan = (div + 1) * rank
            //                  num_chans = div + 1
            //  if rank >= rem: first_chan = (div + 1) * rem + (div * (rank - rem))
            //                  num_chans = div
            // and that reduces to what's below
            auto first_chan = rank * div + (rank < rem ? rank : rem);
            auto num_chans = div + (rank < rem);
            return std::make_pair(first_chan, num_chans);
        }

        int run(int argc, char* argv[]) final
        {
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));
            try {
                bool result = _run(argc, argv, comms);
                // Need to destroy the cubes, to force flushing of data
                itsPsfCube.reset();
                itsResidualCube.reset();
                itsModelCube.reset();
                itsRestoredCube.reset();
                return result;
            } catch (const std::exception &e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected error: " << e.what());
                comms.abort();
                return 1;
            }
        }

        int _run(int argc, char *argv[], askap::askapparallel::AskapParallel &comms)
        {
            StatReporter stats;

            LOFAR::ParameterSet subset(config().makeSubset("Cdeconvolver."));
            itsParset = subset;

            ASKAPLOG_INFO_STR(logger, "ASKAP image (MPI) deconvolver " << ASKAP_PACKAGE_VERSION);

            // Need some metadata for the output cube constructions

            // Lets get the grid,pcf and psf cube names from the parset

            const std::vector<std::string> visGridCubeNames = subset.getStringVector("visgrid",{},true);
            uInt nCubes = visGridCubeNames.size();
            ASKAPCHECK(nCubes > 0,"No input cube provided");
            std::vector<std::string> psfGridCubeNames(nCubes);
            std::vector<std::string> pcfGridCubeNames(nCubes);
            for (int i = 0; i < nCubes; i++) {
                psfGridCubeNames[i] = visGridCubeNames[i];
                psfGridCubeNames[i].replace(0,3,"psf");
                pcfGridCubeNames[i] = visGridCubeNames[i];
                pcfGridCubeNames[i].replace(0,3,"pcf");
            }

            const string imageType = subset.getString("imagetype","fits");
            // Now set it in the parset to avoid change of default at lower levels
            subset.replace("imagetype",imageType);
            // CASA images need to be written one process at a time, for fits we have a choice
            // options "serial", "parallel"
            const bool serialWrite = subset.getString("imageaccess.write","serial")=="serial" ||
                     imageType=="casa";

            // what outputs do we want?
            const bool writeResidual = subset.getBool("write.residualimage",false);
            const bool writePsf = subset.getBool("write.psfimage",false);
            const bool writeModel = subset.getBool("write.modelimage",false);
            const bool restore = subset.getBool("restore",true);
            const bool writeRestored = subset.getBool("write.restoredimage",restore);
            ASKAPCHECK(writeResidual||writePsf||writeModel||writeRestored,"Need to request at least one output image");

            // get the calcstats flag from the parset. if it is true, then this task also calculates the image statistics
            const bool calcstats = subset.getBool("calcstats", false);
            // file to store the statistics (optional)
            const std::string outputStats = subset.getString("outputStats","");

            // are we applying custom weights per grid?
            const std::vector<float> weights(subset.getFloatVector("weights",std::vector<float>(nCubes,1.0f)));
            ASKAPCHECK(weights.size() == nCubes,"Number of weights values must match number of input grids");

            // Check if we're loading real/imag fits cubes
            bool combineRealImag = false;
            bool imagePlaneInput = false;
            if (imageType == "fits") {
                casacore::File test(visGridCubeNames[0]+".fits");
                if (!test.exists()) {
                    casacore::File realPart(visGridCubeNames[0]+".real.fits");
                    casacore::File imagPart(visGridCubeNames[0]+".imag.fits");
                    combineRealImag = realPart.isRegular() && imagPart.isRegular();
                    ASKAPLOG_INFO_STR(logger,"Assuming real+imag uv-grid FITS input with RA/Dec coordinates");
                } else {
                    imagePlaneInput = true;
                    ASKAPLOG_INFO_STR(logger,"Assuming (dirty, psf) image FITS input");
                }
            } else if (imageType == "adios") {
                combineRealImag = true;
                ASKAPLOG_INFO_STR(logger, "Assuming CASA ADIOS2 real+imag uv-grid input");
            } else {
                ASKAPLOG_INFO_STR(logger,"Trying to work out casa image data type for "<<visGridCubeNames[0]);
                shared_ptr<casacore::LatticeBase> lattp(casacore::ImageOpener::openImage(visGridCubeNames[0]));
                imagePlaneInput = (lattp->dataType() == casacore::DataType::TpFloat);
                if (imagePlaneInput) {
                    ASKAPLOG_INFO_STR(logger,"Assuming casa image input");
                } else {
                    ASKAPLOG_INFO_STR(logger,"Assuming casa complex uv-grid input");
                }
            }

            // ok lets set up some output cubes

            // WorkArrays
            casacore::Array<casacore::Float> psfImage;
            casacore::Array<casacore::Float> pcfImage;
            casacore::Array<casacore::Float> dirtyImage;

            casacore::Array<casacore::Complex> pcfGrid;
            casacore::Array<casacore::Complex> psfGrid;
            casacore::Array<casacore::Complex> visGrid;


            boost::shared_ptr<accessors::IImageAccess<casacore::Float> > iaccF;
            if (imageType == "adios") {
                iaccF = imageAccessFactory(subset, comms);
            } else {
                iaccF = imageAccessFactory(subset);
            }
            boost::shared_ptr<accessors::IImageAccess<casacore::Complex> > iaccC;
            // Do we need a complex image accessor?
            if (imageType == "casa" && !imagePlaneInput) {
                iaccC.reset(new CasaImageAccess<casacore::Complex>());
            }


            // Lets load in a cube

            // First get the shape and coordinates
            // The input grids may vary in size - for uv grids we can deal with this by finding the largest
            // one. Image grids will have to match.
            casacore::IPosition shape;
            casacore::Vector<casacore::IPosition> shapes(nCubes);
            casacore::CoordinateSystem coordSys;
            // gaussian taper wants to know field of view
            casacore::Vector<casacore::Double> fov(2,0.0);
            uInt iMax = 0;
            for (uInt i = 0; i < nCubes; i++) {
                if (combineRealImag) {
                    // in this case we only support image plane coordinate system
                    shapes(i) = iaccF->shape(visGridCubeNames[i]+".real");
                } else if (imagePlaneInput) {
                    shapes(i) = iaccF->shape(visGridCubeNames[i]);
                } else{
                    // probably UV coordinate system - use parset
                    shapes(i) = iaccC->shape(visGridCubeNames[i]);
                }
                if (shape.nelements()==0) {
                    shape = shapes(i);
                    iMax = i;
                } else if (imagePlaneInput) {
                    ASKAPCHECK(shape == shapes(i),"All input images should have the same shape");
                } else if (shape[0] < shapes(i)[0]) {
                    shape = shapes(i);
                    iMax = i;
                }
            }
            if (combineRealImag) {
                coordSys = iaccF->coordSys(visGridCubeNames[iMax]+".real");
                fov(0) = shape[0] * abs(coordSys.increment()(0));
                fov(1) = shape[1] * abs(coordSys.increment()(1));
            } else if (imagePlaneInput) {
                coordSys = iaccF->coordSys(visGridCubeNames[iMax]);
                fov(0) = shape[0] * abs(coordSys.increment()(0));
                fov(1) = shape[1] * abs(coordSys.increment()(1));
            } else {
                coordSys = iaccC->coordSys(visGridCubeNames[iMax]);
                ASKAPLOG_INFO_STR(logger,"Setting Images.shape=["+std::to_string(shape[0])+","+std::to_string(shape[1])+"]");
                subset.replace("Images.shape","["+std::to_string(shape[0])+","+std::to_string(shape[1])+"]");
            }

            if (fov(0)>0) {
                ASKAPLOG_INFO_STR(logger,"Image Field of View "<<fov(0)<<" radians"<<" or "<<fov(0)*57.3<<" degrees");
            }

            // sort out oversampling
            float oversampling = subset.getFloat("oversampling",1.0);
            ASKAPCHECK(oversampling>=1,"oversampling parameter must be >=1");
            if (oversampling > 1.0) {
                ASKAPLOG_INFO_STR(logger,"Will oversample grid by factor "<<oversampling);
                subset.add("Images.extraoversampling", utility::toString(oversampling));
            }

            casacore::IPosition blc(shape.nelements(),0);
            casacore::IPosition trc(shape);
            int nChanCube = trc[3];

            ASKAPCHECK(trc[2]==1,"Cannot handle more than 1 polarisation in the cube yet");

            // Define reference channel for giving restoring beam
            std::string reference = itsParset.getString("restore.beamReference", "mid");
            if (reference == "mid") {
              itsBeamReferenceChannel = nChanCube / 2;
            } else if (reference == "first") {
              itsBeamReferenceChannel = 0;
            } else if (reference == "last") {
              itsBeamReferenceChannel = nChanCube - 1;
            } else { // interpret reference as a 0-based channel nuumber
              unsigned int num = std::stoi(reference);
              if (num < nChanCube) {
                itsBeamReferenceChannel = num;
              } else {
                ASKAPLOG_WARN_STR(logger, "beamReference value (" << reference
                << ") not valid. Using middle value of " << nChanCube / 2);
                itsBeamReferenceChannel = nChanCube / 2;
              }
            }

            // create/open the output cubes
            if (comms.isMaster() && imageType != "adios") { // only the master makes the output
                if (imageType == "casa" && !imagePlaneInput) {
                    // we have uv-grids with UV coordinates attached
                    // using the parset to get the image parameters
                    Int pixelAxis,worldAxis,coordinate;
                    CoordinateUtil::findSpectralAxis(pixelAxis,worldAxis,coordinate,coordSys);
                    const SpectralCoordinate &sc = coordSys.spectralCoordinate(coordinate);
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
                    if (writePsf) {
                        itsPsfCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,
                        nChanCube, f0, cdelt, "psf.image"));
                    }
                    if (writeResidual) {
                        itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,
                        nChanCube, f0, cdelt, "residual"));
                    }
                    if (writeModel) {
                        itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,
                        nChanCube, f0, cdelt, "image"));
                    }
                    if (writeRestored) {
                        itsRestoredCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,
                        nChanCube, f0, cdelt, "restored"));
                    }
                } else {
                    // FITS case: use the coordinates of the input images
                    if (writePsf) {
                        itsPsfCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,
                        shape, coordSys, "psf.image"));
                    }
                    if (writeResidual) {
                        itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,
                        shape, coordSys, "residual"));
                    }
                    if (writeModel) {
                        itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,
                        shape, coordSys, "image"));
                    }
                    if (writeRestored) {
                        itsRestoredCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,
                        shape, coordSys, "restored"));
                    }
                }
                initialiseBeamList(nChanCube);
                initialiseWeightsList(nChanCube);
            } else if (imageType == "adios") {
                // we have uv-grids with UV coordinates attached
                // using the parset to get the image parameters
                Int pixelAxis,worldAxis,coordinate;
                CoordinateUtil::findSpectralAxis(pixelAxis,worldAxis,coordinate,coordSys);
                const SpectralCoordinate &sc = coordSys.spectralCoordinate(coordinate);
                casacore::Double baseFreq, nextFreq, freqInc;
                casacore::Double pixelVal=0;

                sc.toWorld(baseFreq,pixelVal);
                sc.toWorld(nextFreq,pixelVal+1);

                freqInc = nextFreq-baseFreq;

                Quantity f0(baseFreq, "Hz");
                Quantity cdelt(freqInc, "Hz");

                ASKAPLOG_INFO_STR(logger,"Base Freq " << f0);
                ASKAPLOG_INFO_STR(logger,"Freq inc (CDELT) " << cdelt);
                if (writePsf) {
                    itsPsfCube.reset(new askap::cp::CubeBuilder<casacore::Float>(comms, 0, subset, nChanCube, f0, cdelt, "psf.image"));
                }
                if (writeResidual) {
                    itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(comms, 0, subset, nChanCube, f0, cdelt, "residual"));
                }
                if (writeModel) {
                    itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(comms, 0, subset, nChanCube, f0, cdelt, "image"));
                }
                if (writeRestored) {
                    itsRestoredCube.reset(new askap::cp::CubeBuilder<casacore::Float>(comms, 0, subset, nChanCube, f0, cdelt, "restored"));
                }
                initialiseBeamList(nChanCube);
                initialiseWeightsList(nChanCube);
            } else {
                // this should work fine as the cubes will exist by the time
                // they are needed.
                if (writePsf) {
                    itsPsfCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,"psf.image"));
                }
                if (writeResidual) {
                    itsResidualCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,"residual"));
                }
                if (writeModel) {
                    itsModelCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,"image"));
                }
                if (writeRestored) {
                    itsRestoredCube.reset(new askap::cp::CubeBuilder<casacore::Float>(subset,"restored"));
                }
            }

            // built the output cubes for this image, all wait until done (for parallel write)
            if (!serialWrite) {
                comms.barrier();
            }


            // create a stats object for this image.
            askap::utils::StatsAndMask statsAndMask(comms);
            if ( calcstats && itsRestoredCube) {
              statsAndMask.set(itsRestoredCube->filename(),itsRestoredCube->imageHandler());
            }


            // What fraction of the full problem does a rank have
            int firstChannel, numChannelsLocal;
            std::tie(firstChannel, numChannelsLocal) = get_channel_allocation(comms, nChanCube);
            ASKAPLOG_INFO_STR(logger,"Rank " << comms.rank() << " - RankAllocation starts at " <<  firstChannel <<
                                     " and is " << numChannelsLocal << " in size");
            bool firstPassForMaster = true;

            casacore::IPosition outShape(4,shape[0],shape[1],1,1);
            casacore::IPosition outBlc(4,0);
            casacore::IPosition outTrc(4,0);
            ASKAPLOG_INFO_STR(logger,"Output image shape " << shape);

            for (int channel = firstChannel; channel < firstChannel + numChannelsLocal; channel++) {

                //FIXME: this is just looping over each channel of the allocation

                ASKAPLOG_INFO_STR(logger,"Processing Channel " << channel);

                if (imagePlaneInput) {
                    psfImage = casacore::Array<casacore::Float>(outShape,0);
                    pcfImage = casacore::Array<casacore::Float>(outShape,0);
                    dirtyImage = casacore::Array<casacore::Float>(outShape,0);
                } else {
                    psfGrid = casacore::Array<casacore::Complex>(outShape,0);
                    pcfGrid = casacore::Array<casacore::Complex>(outShape,0);
                    visGrid = casacore::Array<casacore::Complex>(outShape,0);
                }

                // accumulate multiple inputs (if nCubes>1)
                for (uint i = 0; i < nCubes; i++) {
                    casacore::IPosition inblc(shapes(i).nelements(),0); // input bottom left corner of this allocation
                    casacore::IPosition intrc(shapes(i)-1); // get the top right

                    // assumes pol, chan are axis 2 and 3
                    inblc[3] = channel;
                    intrc[3] = channel;
                    ASKAPCHECK(intrc[2]==0,"Cannot handle >1 polarisation plane in the cubes");
                    outBlc[0] = shape[0]/2 - shapes(i)[0]/2; 
                    outTrc[0] = shape[0]/2 + shapes(i)[0]/2 - 1; 
                    outBlc[1] = shape[1]/2 - shapes(i)[1]/2; 
                    outTrc[1] = shape[1]/2 + shapes(i)[1]/2 - 1; 
                    if (imagePlaneInput) {
                        psfImage += weights[i] * iaccF->read(psfGridCubeNames[i], inblc, intrc);
                        pcfImage += weights[i] * iaccF->read(pcfGridCubeNames[i], inblc, intrc);
                        dirtyImage += weights[i] * iaccF->read(visGridCubeNames[i], inblc, intrc);
                    } else {
                        casacore::Array<casacore::Complex> psfGridRef = psfGrid(outBlc,outTrc);
                        casacore::Array<casacore::Complex> pcfGridRef = pcfGrid(outBlc,outTrc);
                        casacore::Array<casacore::Complex> visGridRef = visGrid(outBlc,outTrc);
                        if (combineRealImag) {
                            psfGridRef += casacore::makeComplex(weights[i] * iaccF->read(psfGridCubeNames[i]+".real",inblc, intrc),
                                weights[i] * iaccF->read(psfGridCubeNames[i]+".imag",inblc,intrc));
                            pcfGridRef += casacore::makeComplex(weights[i] * iaccF->read(pcfGridCubeNames[i]+".real",inblc, intrc),
                                weights[i] * iaccF->read(pcfGridCubeNames[i]+".imag",inblc,intrc));
                            visGridRef += casacore::makeComplex(weights[i] * iaccF->read(visGridCubeNames[i]+".real",inblc, intrc),
                                weights[i] * iaccF->read(visGridCubeNames[i]+".imag",inblc, intrc));
                        } else {
                            psfGridRef += casacore::Complex(weights[i]) * iaccC->read(psfGridCubeNames[i], inblc, intrc);
                            pcfGridRef += casacore::Complex(weights[i]) * iaccC->read(pcfGridCubeNames[i], inblc, intrc);
                            visGridRef += casacore::Complex(weights[i]) * iaccC->read(visGridCubeNames[i], inblc, intrc);
                        }
                    }
                }

                // do the work
                casacore::IPosition subShape = (imagePlaneInput ? dirtyImage.shape() : visGrid.shape());
                imagemath::MultiDimArrayPlaneIter planeIter(subShape);

                for ( ; planeIter.hasMore(); planeIter.next()) {
                    /// FIXME: this is supposed to loop over the polarisations as well as channels
                    /// FIXME: but i have not sorted out the output indexes for this to work

                    casacore::IPosition curpos = planeIter.position();
                    ASKAPLOG_INFO_STR(logger, "Processing from position: " << curpos);

                    // the inputs
                    casacore::Matrix<casacore::Float> psfIn;
                    casacore::Matrix<casacore::Float> pcfIn;
                    casacore::Matrix<casacore::Float> dirtyIn;

                    if (!imagePlaneInput) {
                        casacore::Matrix<casacore::Complex> psfPlane = planeIter.getPlane(psfGrid, curpos).nonDegenerate();
                        casacore::Matrix<casacore::Complex> pcfPlane = planeIter.getPlane(pcfGrid, curpos).nonDegenerate();
                        casacore::Matrix<casacore::Complex> visPlane = planeIter.getPlane(visGrid, curpos).nonDegenerate();
                        getRealFFT(psfIn,psfPlane);
                        getRealFFT(pcfIn,pcfPlane);
                        getRealFFT(dirtyIn,visPlane);
                    } else {
                        psfIn = planeIter.getPlane(psfImage, curpos);
                        pcfIn = planeIter.getPlane(pcfImage, curpos);
                        dirtyIn = planeIter.getPlane(dirtyImage, curpos);
                    }
                    float maxPsf = max(psfIn);
                    itsWeightsList[channel] = maxPsf;
                    ASKAPLOG_INFO_STR(logger,"Max PSF array:" << maxPsf);

                    // the outputs
                    casacore::Array<casacore::Float> psfOut;
                    casacore::Array<casacore::Float> dirty;
                    casacore::Array<casacore::Float> model;
                    casacore::Array<casacore::Float> restored;
                    casacore::Vector<casacore::Quantum<double>> beam(3);
                    doTheWork(subset, dirtyIn, psfIn, pcfIn, psfOut, dirty, model, restored, writeRestored, oversampling, fov, beam);
                    itsBeamList[channel] = beam;

                    if (serialWrite) {
                        if (comms.isMaster() && firstPassForMaster) {
                          ASKAPLOG_INFO_STR(logger, "Ensuring serial access to cubes");
                          firstPassForMaster = false;
                        }
                        else { // this is essentially a serializer - it is required for CASA image types
                        // but not FITS
                          int buf;
                          int from = (comms.rank() > 0 ? comms.rank() - 1 : comms.nProcs()-1);
                          ASKAPLOG_INFO_STR(logger, "Waiting for trigger from rank " << from);
                          comms.receive((void *) &buf,sizeof(int),from);
                        }
                    }
                    // write out the slice
                    // FIXME: THis is the issue with npol I need to use some position

                    // Image based cleaning means oversampling already done - no need to use flexible here
                    if (writePsf) {
                        itsPsfCube->writeRigidSlice(psfOut, channel);
                    }
                    if (writeResidual) {
                        itsResidualCube->writeRigidSlice(dirty, channel);
                    }
                    if (writeModel) {
                        itsModelCube->writeRigidSlice(model, channel);
                    }
                    if (writeRestored) {
                        itsRestoredCube->writeRigidSlice(restored, channel);
                    }
                    if (serialWrite) {
                      int buf = 0;
                      const int to = (comms.rank() == comms.nProcs()-1 ? 0 : comms.rank()+1);
                      comms.send((void *) &buf,sizeof(int),to);
                    }
                    if ( calcstats && writeRestored) {
                      statsAndMask.calculate(channel,restored);
                    }
                }
            }

            if (comms.isMaster() && serialWrite) {
                // one last message to rank 0 to read before we start receiving stats
                int buf;
                int from = comms.nProcs()-1;
                comms.receive((void *) &buf,sizeof(int),from);
            }

            writeBeamInfo(comms);
            writeWeightsInfo(comms);

            if ( calcstats && writeRestored ) {
              // Since the processing of the image channels is distributed among the MPI ranks,
              // the master has to collect all the stats from the worker ranks prior to writing
              // the stats to the image table
              if (comms.isMaster()) {
                statsAndMask.receiveStats();
                statsAndMask.writeStatsToImageTable(itsRestoredCube->filename());
                if ( outputStats != "" ) {
                    statsAndMask.writeStatsToFile(outputStats);
                }
              } else {
                statsAndMask.sendStats();
              }
            }

            stats.logSummary();
            comms.barrier();
            return 0;
        }

    private:
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

/// get the real part of the FFT of the input
void CdeconvolverApp::getRealFFT(casacore::Matrix<casacore::Float> &fArray,
                                casacore::Matrix<casacore::Complex> &cArray) {

    fArray.resize(cArray.shape());
    // Limit number of fft threads to 8 (more is slower for our fft sizes)
    scimath::FFT2DWrapper<imtypeComplex> fft2d(true,8);
    #ifdef ASKAP_FLOAT_IMAGE_PARAMS
    fft2d(cArray,false);
    casacore::real(fArray,cArray);
    #else
    casacore::Matrix<casacore::DComplex> scratch(cArray.shape());
    casacore::convertArray<casacore::DComplex,casacore::Complex>(scratch, cArray);
    fft2d(scratch, false);
    casacore::convertArray<casacore::Float, casacore::Double>(fArray,real(scratch));
    #endif
    fArray *= static_cast<casacore::Float>(fArray.nelements());
}

void CdeconvolverApp::doTheWork(const LOFAR::ParameterSet subset,
                                casacore::Array<casacore::Float> &dirtyIn,
                                casacore::Array<casacore::Float> &psfIn,
                                casacore::Array<casacore::Float> &pcfIn,
                                casacore::Array<casacore::Float> &psfOut,
                                casacore::Array<casacore::Float> &dirtyOut,
                                casacore::Array<casacore::Float> &model,
                                casacore::Array<casacore::Float> &restored,
                                bool writeRestored,
                                float oversampling,
                                casacore::Vector<casacore::Double>& fov,
                                casacore::Vector<casacore::Quantum<double>> & beam)
{

    ASKAPLOG_DEBUG_STR(logger,"Array Shape: " << dirtyIn.shape());

    // Convolution correction probably should pull the support from the PARSET
    // alpha is usually 1
    // support is ususally 3
    int alpha = 1;
    int support = 3;
    scimath::SpheroidalFunction sf(casacore::C::pi*support, alpha);
    const double cutoff = subset.getDouble("restore.beam.cutoff",0.5);
    const int maxsupport = subset.getInt("restore.beam.maxsupport",101);


    #ifdef ASKAP_FLOAT_IMAGE_PARAMS
    SphFuncVisGridder::correctConvolution(dirtyIn,sf,support,true);
    SphFuncVisGridder::correctConvolution(psfIn,sf,support,true);
    #else
    {
        casacore::Array<casacore::Double> dBuffer(dirtyIn.shape());
        casacore::convertArray<casacore::Double, casacore::Float> (dBuffer,dirtyIn);
        SphFuncVisGridder::correctConvolution(dBuffer,sf,support,true);
        casacore::convertArray<casacore::Float, casacore::Double>(dirtyIn,dBuffer);

        casacore::convertArray<casacore::Double, casacore::Float> (dBuffer,psfIn);
        SphFuncVisGridder::correctConvolution(dBuffer,sf,support,true);
        casacore::convertArray<casacore::Float, casacore::Double>(psfIn,dBuffer);
    }
    #endif
    // *** Preconditioning ***
    if (subset.isDefined("preconditioner.Names")) {
        ASKAPLOG_INFO_STR(logger,"Preparing for preconditioning");

        const vector<string> preconditioners=subset.getStringVector("preconditioner.Names",vector<string>());

        // could follow ImageSolverFactory and use addPreconditioner to add each to an ImageSolver.
        // but just keep separate for now

        if ( preconditioners.size() == 0 ) {
            ASKAPLOG_WARN_STR(logger," - no preconditioners given. Deconvolving unfiltered images.");
        }

        for (vector<string>::const_iterator pc = preconditioners.begin(); pc != preconditioners.end(); ++pc) {
            if ( (*pc)=="Wiener" ) {
                ASKAPLOG_INFO_STR(logger," - using a Wiener filter");
            }
            else if ( (*pc) == "GaussianTaper" ) {
                ASKAPLOG_INFO_STR(logger," - using a GaussianTaper");
            }
            else if ( (*pc) == "None" ) {
                ASKAPLOG_INFO_STR(logger," - no preconditioning specified. Deconvolving unfiltered images.");
            }
            else {
                ASKAPTHROW(AskapError, "Unknown preconditioner "<<*pc);
            }
        }

        for (vector<string>::const_iterator pc = preconditioners.begin(); pc != preconditioners.end(); ++pc) {

            // The preconditioner assumes that the PCF is accumulated in the image domain so FFT it.
            // I'm not normalising it as I dont think I need to.

            if ( (*pc)=="Wiener" ) {
                ASKAPLOG_INFO_STR(logger,"Applying Wiener filter");

                // Preconditioning Assuming they only want Wiener.
                boost::shared_ptr<WienerPreconditioner>
                    wp = WienerPreconditioner::createPreconditioner(subset.makeSubset("preconditioner.Wiener."));
                ASKAPASSERT(wp);

                wp->doPreconditioning(psfIn,dirtyIn,pcfIn);

            } else if ( (*pc) == "GaussianTaper") {
                ASKAPLOG_INFO_STR(logger,"Applying GaussianTaper");

                // Copied from ImageSolverFactory::configurePreconditioners. See that for comments and updates.
                // Should really have a function like WienerPreconditioner::createPreconditioner() to set all this up.
                ASKAPCHECK(subset.isDefined("preconditioner.GaussianTaper"),
                           "preconditioner.GaussianTaper parameter is required to use GaussianTaper");
                const vector<double> taper = SynthesisParamsHelper::convertQuantity(
                                                 subset.getStringVector("preconditioner.GaussianTaper"),"rad");
                ASKAPCHECK((taper.size() == 3) || (taper.size() == 1),
                           "preconditioner.GaussianTaper can have either single element or "
                           " a vector of 3 elements. You supplied a vector of "<<taper.size()<<" elements");

                if (fov(0)==0 || fov(1)==0) {
                    ASKAPCHECK(subset.isDefined("Images.shape") && subset.isDefined("Images.cellsize"),
                              "Imager.shape and Imager.cellsize should be defined to convert the taper fwhm "
                              "specified in angular units in the image plane into uv cells");
                    const std::vector<double> cellsize = SynthesisParamsHelper::convertQuantity(
                                                            subset.getStringVector("Images.cellsize"),"rad");
                    const std::vector<int> shape = subset.getInt32Vector("Images.shape");
                    ASKAPCHECK((cellsize.size() == 2) && (shape.size() == 2),
                              "Images.cellsize and Images.shape parameters should have exactly two values");
                    fov(0) = abs(cellsize[0]) * shape[0];
                    fov(1) = abs(cellsize[1]) * shape[1];
                    if (fov(0)>0) {
                        ASKAPLOG_INFO_STR(logger,"Image Field of View "<<fov(0)<<" radians"<<" or "<<fov(0)*57.3<<" degrees");
                    }

                }

                const bool isPsfSize = subset.getBool("preconditioner.GaussianTaper.isPsfSize",False);
                const double tol = subset.getDouble("preconditioner.GaussianTaper.tolerance",0.005);

                /*
                 * leave this unless needed
                 *
                // additional scaling factor due to padding. by default - no padding
                const boost::shared_ptr<ImageCleaningSolver>
                    ics = boost::dynamic_pointer_cast<ImageCleaningSolver>(solver);
                const double paddingFactor = ics ? ics->paddingFactor() : 1.;
                */
                const double paddingFactor = 1.;

                // factors which appear in nominator are effectively half sizes in radians
                const double xFactor = 4. * log(2.) * fov(0) * paddingFactor / casacore::C::pi;
                const double yFactor = 4. * log(2.) * fov(1) * paddingFactor / casacore::C::pi;

                boost::shared_ptr<GaussianTaperPreconditioner> gp;

                if (taper.size() == 3) {
                  ASKAPDEBUGASSERT((taper[0]!=0) && (taper[1]!=0));
                  gp.reset(new GaussianTaperPreconditioner(xFactor/taper[0],yFactor/taper[1],taper[2],
                                                           isPsfSize,cutoff,maxsupport,tol));
                  //solver->addPreconditioner(IImagePreconditioner::ShPtr(
                  //    new GaussianTaperPreconditioner(xFactor/taper[0],yFactor/taper[1],taper[2],
                  //                                    isPsfSize,cutoff,tol)));
                } else {
                  ASKAPDEBUGASSERT(taper[0]!=0);
                  if (std::abs(xFactor-yFactor)<4e-15) {
                    // the image is square, can use the short cut
                    gp.reset(new GaussianTaperPreconditioner(xFactor/taper[0],isPsfSize,cutoff,maxsupport,tol));
                    //solver->addPreconditioner(IImagePreconditioner::ShPtr(
                    //    new GaussianTaperPreconditioner(xFactor/taper[0],isPsfSize,cutoff,tol)));
                  } else {
                    // the image is rectangular. Although the gaussian taper is symmetric in
                    // angular coordinates, it will be elongated along the vertical axis in
                    // the uv-coordinates.
                    gp.reset(new GaussianTaperPreconditioner(xFactor/taper[0], yFactor/taper[0],
                                                             0., isPsfSize, cutoff, maxsupport, tol));
                    //solver->addPreconditioner(IImagePreconditioner::ShPtr(
                    //    new GaussianTaperPreconditioner(xFactor/taper[0],yFactor/taper[0],0.,isPsfSize,cutoff,tol)));
                  } // xFactor!=yFactor
                } // else: taper.size() == 3

                gp->doPreconditioning(psfIn,dirtyIn,pcfIn);

            }

        } // loop over all preconditioners

    } else {
        ASKAPLOG_INFO_STR(logger,"No preconditioning");
    }

    // *** Normalisation ***
    const float unnormalisedMaxPSF = max(psfIn);
    if (unnormalisedMaxPSF<=0.) {
        ASKAPTHROW(AskapError, "PSF Error. Peak is " << unnormalisedMaxPSF);
    }

    ASKAPLOG_INFO_STR(logger, "Normalising PSF");
    psfIn /= unnormalisedMaxPSF;
    ASKAPLOG_INFO_STR(logger, "Peak of PSF before normalisation = " << unnormalisedMaxPSF<< ", after normalisation = " <<max(psfIn));
    ASKAPLOG_INFO_STR(logger, "Normalising Dirty image");
    dirtyIn /= unnormalisedMaxPSF;
    ASKAPLOG_INFO_STR(logger, "Peak of Dirty Image after normalisation = " << max(dirtyIn));

    // Now oversample the dirty image and psf if requested
    if (oversampling > 1) {
        ASKAPLOG_INFO_STR(logger,"Oversampling image and PSF by factor "<<oversampling<<" before Clean");
        SynthesisParamsHelper::oversample(dirtyIn, oversampling);
        SynthesisParamsHelper::oversample(psfIn, oversampling);
    }

    // *** Deconvolution ***
    // We have a normalised corrected preconditioned image
    DeconvolverBase<Float, Complex>::ShPtr deconvolver;
    string algorithm = subset.getString("solver.Clean.algorithm", "BasisfunctionMFS");

    if (algorithm == "Basisfunction") {
        ASKAPLOG_INFO_STR(logger, "Constructing Basisfunction Clean solver");
        deconvolver.reset(new DeconvolverBasisFunction<Float, Complex>(dirtyIn, psfIn));
        ASKAPASSERT(deconvolver);
    } else if (algorithm == "BasisfunctionMFS") {
        ASKAPLOG_INFO_STR(logger, "Constructing MultiTermBasisfunction Clean solver");
        deconvolver.reset(new DeconvolverMultiTermBasisFunction<Float, Complex>(dirtyIn, psfIn));
        ASKAPASSERT(deconvolver);
    } else if (algorithm == "Hogbom") {
        ASKAPLOG_INFO_STR(logger, "Constructing Hogbom Clean deconvolver");
        deconvolver.reset(new DeconvolverHogbom<Float, Complex>(dirtyIn, psfIn));
        ASKAPASSERT(deconvolver);
    } else {
        ASKAPTHROW(AskapError, "Unknown Clean algorithm " << algorithm);
    }

    // copy cleaning parameters and add any extra stopping criteria
    LOFAR::ParameterSet cleanset = subset.makeSubset("solver.Clean.");

    Float sigmaValue(0);
    Matrix<imtype> madMap;
    // could make the following a function that returns the updated parset and add to configure line
    const std::string parName = "threshold.minorcycle";
    if (subset.isDefined(parName)) {
        const std::vector<std::string> thresholds = subset.getStringVector(parName);
        ASKAPCHECK(thresholds.size() && (thresholds.size()<4), "Parameter "<<parName<<
                   " must contain either 1 element or a vector of 2-3 elements, you have "<< thresholds.size());
        bool absoluteThresholdDefined = false;
        bool relativeThresholdDefined = false;
        bool absoluteThreshold2Defined = false;
        bool noiseThresholdDefined = false;
        bool noiseThreshold2Defined = false;
        // do we want spatially variant sigma thresholds?
        const uInt boxSize = subset.getUint("solver.Clean.noiseboxsize",0);
        for (const string& t : thresholds) {
            casacore::Quantity cThreshold;
            // check for noise thresholds
            const string sigma("sigma");
            const size_t pos = t.rfind(sigma);
            if (pos != std::string::npos && pos == t.size() - sigma.size()) {
              ASKAPCHECK(!noiseThreshold2Defined, "Parameter "<<parName<<
                         " defines noise threshold thrice ("<<t<<")");
              ASKAPCHECK(!absoluteThresholdDefined,"Cannot mix "<<
              "absolute and noise thresholds");
              casacore::Quantity::read(cThreshold, t.substr(0,pos));
              if (noiseThresholdDefined) {
                  noiseThreshold2Defined = true;
                  std::ostringstream pstr;
                  pstr<<cThreshold.getValue("") * sigmaValue;
                  cleanset.add("targetobjective2", pstr.str().c_str());
                  ASKAPLOG_INFO_STR(logger, "Will stop deep minor cycle at the noise threshold of "<<
                                      cThreshold.getValue("")<<" sigma");
              } else {
                  noiseThresholdDefined = true;
                // get noise for thresholds if needed
                // get mad estimate for sigma
                // may need to take mask into account?
                Float mad = casacore::madfm(dirtyIn);
                sigmaValue = 1.48f * mad;
                boost::shared_ptr<DeconvolverMultiTermBasisFunction<Float,Complex>> dcmtbf = 
                    boost::dynamic_pointer_cast<DeconvolverMultiTermBasisFunction<Float,Complex>>(deconvolver);
                if (dcmtbf && boxSize > 0) {
                    // get mad map for position dependent threshold
                    madMap = casacore::boxedArrayMath(dirtyIn.nonDegenerate(),
                        IPosition(2,boxSize),MadfmFunc<imtype>());
                    //normalise madMap to overall mad and send it to cleaner
                    if (mad > 0) {
                        madMap /= mad;
                        // do we want to enforce madMap >= 1 ?
                    }
                   dcmtbf->setNoiseMap(madMap, boxSize);
                }
                std::ostringstream pstr;
                pstr<<cThreshold.getValue("") * sigmaValue;
                cleanset.add("targetobjective", pstr.str().c_str());
                ASKAPLOG_INFO_STR(logger, "Will stop minor cycle at the noise threshold of "<<
                                    cThreshold.getValue("")<<" sigma");
                if (boxSize > 0) {
                ASKAPLOG_INFO_STR(logger, "Will use a spatially variant noise threshold with box size of "<<
                                    boxSize<<" pixels");
                }
              }
            } else {
                casacore::Quantity::read(cThreshold, t);
                cThreshold.convert();
                if (cThreshold.isConform("Jy")) {
                    ASKAPCHECK(!absoluteThreshold2Defined, "Parameter "<<parName<<
                            " defines absolute threshold thrice ("<<t<<")");
                    ASKAPCHECK(!noiseThresholdDefined,"Cannot mix "<<
                    "absolute and noise thresholds");

                    if (absoluteThresholdDefined) {
                    absoluteThreshold2Defined = true;
                        std::ostringstream pstr;
                        pstr<<cThreshold.getValue("Jy");
                        cleanset.add("targetobjective2", pstr.str().c_str());
                        ASKAPLOG_INFO_STR(logger, "Will stop deep minor cycle at the absolute threshold of "<<
                                        cThreshold.getValue("mJy")<<" mJy");
                    } else {
                        absoluteThresholdDefined = true;
                        std::ostringstream pstr;
                        pstr<<cThreshold.getValue("Jy");
                        cleanset.add("targetobjective", pstr.str().c_str());
                        ASKAPLOG_INFO_STR(logger, "Will stop the minor cycle at the absolute threshold of "<<
                                        pstr.str().c_str()<<" Jy");
                    }
                } else if (cThreshold.isConform("")) {
                    ASKAPCHECK(!relativeThresholdDefined, "Parameter "<<parName<<
                            " defines relative threshold twice ("<<t<<")");
                    relativeThresholdDefined = true;
                    std::ostringstream pstr;
                    pstr<<cThreshold.getValue();
                    cleanset.add("fractionalthreshold", pstr.str().c_str());
                    ASKAPLOG_INFO_STR(logger, "Will stop minor cycle at the relative threshold of "<<
                                    cThreshold.getValue()*100.<<"\%");
                } else {
                    ASKAPTHROW(AskapError, "Unable to convert units in the quantity "<<
                            cThreshold<<" to either Jy or a dimensionless quantity");
                }
            }
        }
    }
    // tortured way to get the (oversampled) cellsize of the output cubes (only one has to exist)
    casacore::Vector<casacore::Double> increments;
    if (itsRestoredCube) {
        increments = itsRestoredCube->imageHandler()->coordSys(itsRestoredCube->filename()).directionCoordinate().increment();
    } else if (itsResidualCube) {
        increments = itsResidualCube->imageHandler()->coordSys(itsResidualCube->filename()).directionCoordinate().increment();
    } else if (itsModelCube) {
        increments = itsModelCube->imageHandler()->coordSys(itsModelCube->filename()).directionCoordinate().increment();
    } else if (itsPsfCube) {
        increments = itsPsfCube->imageHandler()->coordSys(itsPsfCube->filename()).directionCoordinate().increment();
    } else {
        ASKAPTHROW(AskapError, "Cannot determine output image increments");
    }

    // Fit the PSF and set restoring beam in Deconvolver
    if (writeRestored) {
        const vector<string> beampar = subset.getStringVector("restore.beam");
        if (beampar.size() == 1) {
            ASKAPCHECK(beampar[0] == "fit",
                "beam parameter should be either equal to 'fit' or contain 3 elements defining the beam size."<<
                " You have "<<beampar[0]);
            ASKAPLOG_INFO_STR(logger, "Fitting restoring beam");
            casacore::Vector<double> result = SynthesisParamsHelper::fitBeam(psfIn, cutoff, maxsupport);
            ASKAPDEBUGASSERT(result.size() == 3);
            // Deconvolver wants pixels and degrees for beam
            cleanset.replace("beam","["+std::to_string(result[0])+","+std::to_string(result[1])+","+std::to_string(result[2]/C::pi * 180.0)+"]");
            // now convert to arcsec & sky PA in deg (code from SynthesisParamsHelper)
            beam[0] = casacore::Quantum<double>(fabs(increments[0])*result[0],"rad").get("arcsec");
            beam[1] = casacore::Quantum<double>(fabs(increments[1])*result[1],"rad").get("arcsec");
            double pa = increments[0]<0 ? result[2] : -result[2];
            if (pa < -casacore::C::pi/2) {
                pa += casacore::C::pi;
            }
            beam[2] = casacore::Quantum<double>(pa,"rad").get("deg");

        } else {
            ASKAPCHECK(beampar.size() == 3, "Need three elements for beam or a single word 'fit'. You have "<<subset.getString("restore.beam"));
            for (int i=0; i<3; ++i) {
                casacore::Quantity::read(beam(i), beampar[i]);
            }
            // convert to pixels for deconvolver
            cleanset.replace("beam","["+std::to_string(beam[0].get("rad").getValue()/fabs(increments[0]))+","+
                std::to_string(beam[1].get("rad").getValue()/increments[1])+","+std::to_string(beam[2].get("rad").getValue())+"]");

        }
        ASKAPLOG_INFO_STR(logger, "Restore will convolve with the 2D gaussian: " << beam);
    } else {
        // deconvolver wants a beam even if we're not restoring
        cleanset.replace("beam","[0,0,0]");
    }

    ASKAPLOG_INFO_STR(logger,"Configure deconvolver");
    cleanset.replace("configurebase","true");
    deconvolver->configure(cleanset);


    ASKAPLOG_INFO_STR(logger,"Do the deconvolution");
    deconvolver->deconvolve();

    dirtyOut = deconvolver->dirty();
    model = deconvolver->model();
    ASKAPLOG_INFO_STR(logger,"After deconvolution: Peak of model    : "<<max(model)<<", Peak of residuals: "<<max(dirtyOut));

    psfOut = psfIn;

    if (writeRestored) {
        // this was set by "name" for facets, but surely it should be equal for all facets
        casacore::Vector< casacore::Array<casacore::Float> > restored_vec(1); // sometimes there is more than one restored image
        ASKAPLOG_INFO_STR(logger,"Restore the model");
        if(deconvolver->restore(restored_vec)) {
            restored = restored_vec(0);
        }
    }
}
void CdeconvolverApp::writeBeamInfo(askap::askapparallel::AskapParallel &comms)
{
    if (itsRestoredCube) {
        askap::accessors::BeamLogger beamlog;
        beamlog.beamlist() = itsBeamList;
        if (comms.nProcs() > 1) {
          beamlog.gather(comms, 0, true);
        }

        if (comms.isMaster()) {
            ASKAPLOG_INFO_STR(logger, "Channel-dependent restoring beams will be written to image " << itsRestoredCube->filename());
            itsRestoredCube->addBeamList(beamlog.beamlist());

            if (itsParset.getString("imagetype") == "fits") {
              // can't write ref beam to casa image if per channel beams are stored
              ASKAPLOG_DEBUG_STR(logger, "Writing reference restoring beam to header of restored cube");
              casa::Vector<casa::Quantum<double> > refbeam = beamlog.beam(itsBeamReferenceChannel);
              itsRestoredCube->addBeam(refbeam);
            } else {
              itsRestoredCube->setUnits("Jy/beam");
            }
        }
    }
}
void CdeconvolverApp::writeWeightsInfo(askap::askapparallel::AskapParallel &comms)
{
    askap::accessors::WeightsLog weightsLog;
    weightsLog.weightslist() = itsWeightsList;
    if (comms.nProcs() > 1) {
      weightsLog.gather(comms, 0, true);
    }

    if (comms.isMaster()) {
        casacore::Record wtInfo = weightsLog.toRecord();
        if (itsRestoredCube) {
            ASKAPLOG_INFO_STR(logger, "Channel-dependent weights will be written to image " << itsRestoredCube->filename());
            itsRestoredCube->setInfo(wtInfo);
        }
        if (itsResidualCube) {
            ASKAPLOG_INFO_STR(logger, "Channel-dependent weights will be written to image " << itsResidualCube->filename());
            itsResidualCube->setInfo(wtInfo);
        }
    }
}

void CdeconvolverApp::initialiseBeamList(const unsigned int numChannels)
{

    casa::Vector<casa::Quantum<double> > beamVec(3);
    beamVec[0] = casa::Quantum<double>(0., "rad");
    beamVec[1] = casa::Quantum<double>(0., "rad");
    beamVec[2] = casa::Quantum<double>(0., "deg");

    for(unsigned int i=0;i<numChannels;i++) {
        itsBeamList[i] = beamVec;
    }

}
void CdeconvolverApp::initialiseWeightsList(const unsigned int numChannels)
{
  for(unsigned int i=0;i<numChannels;i++) {
      itsWeightsList[i] = 0.0;
  }
}


// Main function
int main(int argc, char* argv[])
{
    CdeconvolverApp app;
    return app.main(argc, argv);
}
