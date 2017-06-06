
///local includes
#include <askap_synthesis.h>
#include <utils/LinmosUtils.h>
#include <measurementequation/SynthesisParamsHelper.h>
///ASKAP includes
#include <askap/Application.h>
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/StatReporter.h>
#include <askapparallel/AskapParallel.h>

///CASA includes
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/images/Images/ImageProxy.h>

///3rd party
#include <Common/ParameterSet.h>


ASKAP_LOGGER(logger, ".linmos");


using namespace askap;
using namespace askap::synthesis;


// @brief do the merge
/// @param[in] parset subset with parameters
static void mergeMPI(const LOFAR::ParameterSet &parset, askap::askapparallel::AskapParallel &comms) {

    ASKAPLOG_INFO_STR(logger, "ASKAP linear (parallel) mosaic task " << ASKAP_PACKAGE_VERSION);
    ASKAPLOG_INFO_STR(logger, "Parset parameters:\n" << parset);

    imagemath::LinmosAccumulator<float> accumulator;

    // Original shape
    int originalNchan = -1;
    // load the parset
    if ( !accumulator.loadParset(parset) ) return;

    // initialise an image accessor
    accessors::IImageAccess& iacc = SynthesisParamsHelper::imageHandler();
    // if we have Taylor terms and we need to correct them for the beam spectral
    // index - do it now ...


    // loop over the mosaics, reading each in an adding to the output pixel arrays
    vector<string> inImgNames, inWgtNames, inSenNames;
    string outImgName, outWgtName, outSenName;
    map<string,string> outWgtNames = accumulator.outWgtNames();
    for(map<string,string>::iterator ii=outWgtNames.begin(); ii!=outWgtNames.end(); ++ii) {

        // get output files for this mosaic
        outImgName = (*ii).first;
        outWgtName = accumulator.outWgtNames()[outImgName];
        ASKAPLOG_INFO_STR(logger, "++++++++++++++++++++++++++++++++++++++++++");
        ASKAPLOG_INFO_STR(logger, "Preparing mosaic " << outImgName);
        if (!accumulator.outWgtDuplicates()[outImgName]) {
            ASKAPLOG_INFO_STR(logger, " - also weights image " << outWgtName);
        }
        accumulator.doSensitivity(false);
        if (accumulator.genSensitivityImage()[outImgName]) {
            outSenName = accumulator.outSenNames()[outImgName];
            accumulator.doSensitivity(true);
            ASKAPLOG_INFO_STR(logger, " - also sensitivity image " << outSenName);
        }

        // get input files for this mosaic
        inImgNames = accumulator.inImgNameVecs()[outImgName];
        ASKAPLOG_INFO_STR(logger, " - input images: "<<inImgNames);

        if (accumulator.weightType() == FROM_WEIGHT_IMAGES) {
            inWgtNames = accumulator.inWgtNameVecs()[outImgName];
            ASKAPLOG_INFO_STR(logger, " - input weights images: " << inWgtNames);
        }
        else if (accumulator.weightType() == FROM_BP_MODEL) {
            accumulator.beamCentres(loadBeamCentres(parset,iacc,inImgNames));
        }

        if (accumulator.doSensitivity()) {
            inSenNames = accumulator.inSenNameVecs()[outImgName];
            ASKAPLOG_INFO_STR(logger, " - input sensitivity images: " << inSenNames);
        }




    // set the output coordinate system and shape, based on the overlap of input images

        vector<IPosition> inShapeVec;
        vector<CoordinateSystem> inCoordSysVec;
        int myAllocationSize = 0;
        int myAllocationStart = 0;

        for (vector<string>::iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {


            const casa::IPosition shape = iacc.shape(*it);


            ASKAPCHECK(shape.nelements()>=3,"Work with at least 3D cubes!");

            ASKAPLOG_INFO_STR(logger," - ImageAccess Shape " << shape);


            casa::IPosition blc(shape.nelements(),0);
            casa::IPosition trc(shape);
            if (comms.rank() >= trc[3]) {
                ASKAPLOG_WARN_STR(logger,"Rank " << comms.rank() << " has no work to merge");
                return;
            }
            if (trc[3] % comms.nProcs() != 0) {
                ASKAPLOG_WARN_STR(logger,"Unbalanced allocation: num of ranks:" << comms.nProcs() << " not a factor of number of channels: "<< trc[3]);
            }
            if (comms.nProcs() >= trc[3]) {
                myAllocationSize = 1;
            }
            else {
                myAllocationSize = trc[3]/comms.nProcs();
            }

            myAllocationStart = comms.rank()*myAllocationSize;

                // unless last rank
            if (comms.rank() == comms.nProcs()-1) {
                myAllocationSize = trc[3] - myAllocationStart; // we are using End is Last
            }


            blc[3] = myAllocationStart;
            trc[0] = trc[0]-1;
            trc[1] = trc[1]-1;
            trc[2] = trc[2]-1;
            trc[3] = myAllocationStart + myAllocationSize-1;


            ASKAPLOG_INFO_STR(logger,"Allocation starts at " << myAllocationStart << " and is " << myAllocationSize << " in size");

            ASKAPCHECK(blc[3]>=0 && blc[3]<shape[3], "Start channel is outside the number of channels or negative, shape: "<<shape);
            ASKAPCHECK(trc[3]<=shape[3], "Subcube extends beyond the original cube, shape:"<<shape);

            ASKAPLOG_INFO_STR(logger, " - Corners " << "blc  = " << blc << ", trc = " << trc << "\n");
            inCoordSysVec.push_back(iacc.coordSysSlice(*it,blc,trc));
            // reset the shape to be the size ...
            trc[0] = trc[0]+1;
            trc[1] = trc[1]+1;
            trc[2] = trc[2]+1;
            trc[3] = myAllocationSize;
            const casa::IPosition shape3(trc);
            ASKAPLOG_INFO_STR(logger, " - Calculated Shape " << shape3);
            inShapeVec.push_back(shape3);


        }

        accumulator.setOutputParameters(inShapeVec, inCoordSysVec);
        ASKAPLOG_INFO_STR(logger, " - Output Shape " << accumulator.outShape());
        // set up the output pixel arrays
        Array<float> outPix(accumulator.outShape(),0.);
        Array<float> outWgtPix(accumulator.outShape(),0.);
        Array<bool> outMask(accumulator.outShape(),0.);

        Array<float> outSenPix;
        if (accumulator.doSensitivity()) {
            outSenPix = Array<float>(accumulator.outShape(),0.);
        }

        // set up an indexing vector for the arrays
        casa::IPosition curpos(outPix.shape());
        ASKAPASSERT(curpos.nelements()>=2);
        for (uInt dim=0; dim<curpos.nelements(); ++dim) {
            curpos[dim] = 0;
        }


        // loop over the input images, reading each in an adding to the output pixel arrays
        for (uInt img = 0; img < inImgNames.size(); ++img ) {

            // short cuts
            string inImgName = inImgNames[img];
            string inWgtName, inSenName;

            ASKAPLOG_INFO_STR(logger, "Processing input image " << inImgName);
            if (accumulator.weightType() == FROM_WEIGHT_IMAGES) {
                inWgtName = inWgtNames[img];
                ASKAPLOG_INFO_STR(logger, " - and input weight image " << inWgtName);
            }
            if (accumulator.doSensitivity()) {
                inSenName = inSenNames[img];
                ASKAPLOG_INFO_STR(logger, " - and input sensitivity image " << inSenName);
            }

            //casa::PagedImage<casa::Float> inImg(inImgName);
            const casa::IPosition shape = iacc.shape(inImgName);
            casa::IPosition blc(shape.nelements(),0);
            casa::IPosition trc(shape);


            if (originalNchan < 0) {
                originalNchan = trc[3];
            }
            else {
                ASKAPCHECK(originalNchan == trc[3],"Nchan missmatch in merge" );
            }
            // this assumes all allocations
            blc[3] = myAllocationStart;
            trc[0] = trc[0]-1;
            trc[1] = trc[1]-1;
            trc[2] = trc[2]-1;
            trc[3] = myAllocationStart + myAllocationSize-1;



            accumulator.setInputParameters(inShapeVec[img], inCoordSysVec[img], img);

            Array<float> inPix = iacc.read(inImgName,blc,trc);

            Array<float> inWgtPix;
            Array<float> inSenPix;


            if (accumulator.weightType() == FROM_WEIGHT_IMAGES) {

                //casa::PagedImage<casa::Float> inImg(inWgtName);
                const casa::IPosition shape = iacc.shape(inWgtName);
                casa::IPosition blc(shape.nelements(),0);
                casa::IPosition trc(shape);

                blc[3] = myAllocationStart;
                trc[0] = trc[0]-1;
                trc[1] = trc[1]-1;
                trc[2] = trc[2]-1;
                trc[3] = myAllocationStart + myAllocationSize-1;

                inWgtPix = iacc.read(inWgtName,blc,trc);

                ASKAPASSERT(inPix.shape() == inWgtPix.shape());
            }
            if (accumulator.doSensitivity()) {

                // casa::PagedImage<casa::Float> inImg(inSenName);
                const casa::IPosition shape = iacc.shape(inSenName);
                casa::IPosition blc(shape.nelements(),0);
                casa::IPosition trc(shape);

                blc[3] = myAllocationStart;
                trc[0] = trc[0]-1;
                trc[1] = trc[1]-1;
                trc[2] = trc[2]-1;
                trc[3] = myAllocationStart + myAllocationSize-1;


                inSenPix = iacc.read(inSenName,blc,trc);

                ASKAPASSERT(inPix.shape() == inSenPix.shape());
            }

            // set up an iterator for all directionCoordinate planes in the input image
            scimath::MultiDimArrayPlaneIter planeIter(accumulator.inShape());

            // test whether to simply add weighted pixels, or whether a regrid is required
            bool regridRequired = !accumulator.coordinatesAreEqual();

            // if regridding is required, set up buffer some images
            if ( regridRequired ) {

                ASKAPLOG_INFO_STR(logger, " - regridding -- input pixel grid is different from the output");
                // currently all output planes have full-size, so only initialise once
                // would be faster if this was reduced to the size of the current input image
                if ( accumulator.outputBufferSetupRequired() ) {
                    ASKAPLOG_INFO_STR(logger, " - initialising output buffers and the regridder");
                    // set up temp images required for regridding
                    //accumulator.initialiseOutputBuffers();
                    // set up regridder
                    accumulator.initialiseRegridder();
                }

                // set up temp images required for regridding
                // need to do this here if some do and some do not have sensitivity images
                accumulator.initialiseOutputBuffers();

                // set up temp images required for regridding
                // are those of the previous iteration correctly freed?
                accumulator.initialiseInputBuffers();

            }

            else {
                ASKAPLOG_INFO_STR(logger, " - not regridding -- input pixel grid is the same as the output");
            }
            // iterator over planes (e.g. freq & polarisation), regridding and accumulating weights and weighted images
            for (; planeIter.hasMore(); planeIter.next()) {

                // set the indices of any higher-order dimensions for this slice
                curpos = planeIter.position();

                ASKAPLOG_INFO_STR(logger, " - slice " << curpos);

                if ( regridRequired ) {

                    // load input buffer for the current plane
                    accumulator.loadInputBuffers(planeIter, inPix, inWgtPix, inSenPix);
                    // call regrid for any buffered images
                    accumulator.regrid();
                    // update the accululation arrays for this plane
                    accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, curpos);

                } else {

                    // Update the accululation arrays for this plane.
                    accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, inPix, inWgtPix, inSenPix, curpos);

                }

            }
        }
        //build the mask
        //use the outWgtPix to define the mask
        //i dont care about planes etc ... just going to run through

        float itsCutoff = 0.01;

        if (parset.isDefined("cutoff")) itsCutoff = parset.getFloat("cutoff");
        Array<bool>::iterator iterMask = outMask.begin();
        Array<float>::iterator iterWgt = outWgtPix.begin();



        /// This logic is in addition to the mask in the accumulator
        /// which works on an individual beam weight
        /// this is masking the output mosaick - it therefore has slightly
        /// different criteria. In this case the final weight has to be equal to
        /// or bigger than the cutoff.
        /// There is a possible failure mode where due to rounding a pixel maybe
        /// have been masked by the accumulator but missed here.
        /// The first time I implemented this I just used a looser conditional:
        /// > instead of >= - but in the second attempt I decided to replace all
        /// masked pixels by NaN - which has the nice secondary effect of implementing
        /// the FITS mask.

        float wgtCutoff = itsCutoff * itsCutoff;
        for( ; iterWgt != outWgtPix.end() ; iterWgt++ ) {
            if (*iterWgt >= wgtCutoff) {
                *iterMask = casa::True;
            } else {
                *iterMask = casa::False;
                setNaN(*iterWgt);
            }
            iterMask++;
        }


        // deweight the image pixels
        // use another iterator to loop over planes
        ASKAPLOG_INFO_STR(logger, "Deweighting accumulated images");
        scimath::MultiDimArrayPlaneIter deweightIter(accumulator.outShape());
        for (; deweightIter.hasMore(); deweightIter.next()) {
            curpos = deweightIter.position();
            accumulator.deweightPlane(outPix, outWgtPix, outSenPix, curpos);
        }
        // set one of the input images as a reference for metadata (the first by default)
        uint psfref = 0;
        if (parset.isDefined("psfref")) psfref = parset.getUint("psfref");

        ASKAPLOG_INFO_STR(logger, "Getting brightness info for the output image from input number " << psfref);
        // get pixel units from the selected reference image

        string units = iacc.getUnits(inImgNames[psfref]);
        ASKAPLOG_INFO_STR(logger, "Got units as " << units);

        ASKAPLOG_INFO_STR(logger, "Getting PSF beam info for the output image from input number " << psfref);
        // get psf beam information from the selected reference image
        Vector<Quantum<double> > psf;
        Vector<Quantum<double> > psftmp = iacc.beamInfo(inImgNames[psfref]);
        if (psftmp.nelements()<3) {
            ASKAPLOG_WARN_STR(logger, inImgNames[psfref] <<
                ": beamInfo needs at least 3 elements. Not writing PSF");
        }
        else if ((psftmp[0].getValue("rad")==0) || (psftmp[1].getValue("rad")==0)) {
            ASKAPLOG_WARN_STR(logger, inImgNames[psfref] <<
                ": beamInfo invalid. Not writing PSF");
        }
        else {
            psf = psftmp;
        }
        // write accumulated images and weight images
        ASKAPLOG_INFO_STR(logger, "Writing accumulated image to " << outImgName);
        casa::IPosition outShape = accumulator.outShape();



        if (comms.isMaster()) {
            outShape[3] = originalNchan;

            ASKAPLOG_INFO_STR(logger, " Creating output file - Shape " << outShape << " OriginalNchan " << originalNchan);
            iacc.create(outImgName, outShape, accumulator.outCoordSys());
            iacc.makeDefaultMask(outImgName);

        }
        else {
            int buf;
            int from = comms.rank() - 1;
            comms.receive((void *) &buf,sizeof(int),from);
        }
        casa::IPosition loc(outShape.nelements(),0);
        loc[3] = myAllocationStart;
        ASKAPLOG_INFO_STR(logger, " - location " << loc);
        iacc.write(outImgName,outPix,loc);
        iacc.writeMask(outImgName,outMask,loc);
        iacc.setUnits(outImgName,units);

        if (psf.nelements()>=3)
            iacc.setBeamInfo(outImgName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));

        if (accumulator.outWgtDuplicates()[outImgName]) {
            ASKAPLOG_INFO_STR(logger, "Accumulated weight image " << outWgtName << " already written");
        } else {
            if (comms.isMaster()) {
                ASKAPLOG_INFO_STR(logger, "Writing accumulated weight image to " << outWgtName);
                iacc.create(outWgtName, outShape, accumulator.outCoordSys());
                iacc.makeDefaultMask(outWgtName);
            }
            iacc.write(outWgtName,outWgtPix,loc);
            iacc.writeMask(outWgtName,outMask,loc);
            iacc.setUnits(outWgtName,units);
            if (psf.nelements()>=3)
                iacc.setBeamInfo(outWgtName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
        }

        if (accumulator.doSensitivity()) {
            if (comms.isMaster()) {
                ASKAPLOG_INFO_STR(logger, "Writing accumulated sensitivity image to " << outSenName);
                iacc.create(outSenName, outShape, accumulator.outCoordSys());
                iacc.makeDefaultMask(outSenName);
            }
            iacc.write(outSenName,outSenPix,loc);
            iacc.writeMask(outSenName,outMask,loc);
            iacc.setUnits(outSenName,units);
            if (psf.nelements()>=3)
                iacc.setBeamInfo(outSenName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
        }

        if (comms.rank() < comms.nProcs()-1) {
            int buf;
            int to = comms.rank()+1;
            comms.send((void *) &buf,sizeof(int),to);
        }

    }
}
class linmosMPIApp : public askap::Application
{
    public:

        virtual int run(int argc, char* argv[]) {

            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;
                LOFAR::ParameterSet subset(config().makeSubset("linmos."));
                SynthesisParamsHelper::setUpImageHandler(subset);
                mergeMPI(subset, comms);
                stats.logSummary();
                return 0;
            }
            catch (const askap::AskapError& e) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << e.what());
                std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
                return 1;
            } catch (const std::exception& e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << e.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what()
                    << std::endl;
                return 1;
            }


    };

};

int main(int argc, char *argv[])
{
    linmosMPIApp app;
    return app.main(argc, argv);
}
