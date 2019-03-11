
/// local includes
#include <askap_synthesis.h>
#include <utils/LinmosUtils.h>
#include <measurementequation/SynthesisParamsHelper.h>
/// ASKAP includes
#include <askap/Application.h>
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/StatReporter.h>
#include <askapparallel/AskapParallel.h>

/// CASA includes
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/images/Images/ImageProxy.h>

/// 3rd party
#include <Common/ParameterSet.h>


ASKAP_LOGGER(logger, ".linmos");


using namespace askap::synthesis;

namespace askap {

// @brief do the merge
/// @param[in] parset subset with parameters
static void mergeMPI(const LOFAR::ParameterSet &parset, askap::askapparallel::AskapParallel &comms) {

  ASKAPLOG_INFO_STR(logger, "ASKAP linear (parallel) mosaic task (MPI+LOWMEM)" << ASKAP_PACKAGE_VERSION);
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

  // Ahh this loops over the output mosaicks first

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
    ASKAPLOG_INFO_STR(logger, "output mosaic " <<outImgName << " has input images: "<<inImgNames);

    if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED ) {
      inWgtNames = accumulator.inWgtNameVecs()[outImgName];
      ASKAPLOG_INFO_STR(logger, " - input weights images: " << inWgtNames);
    }

    if (accumulator.weightType() == FROM_BP_MODEL) {
      accumulator.beamCentres(loadBeamCentres(parset,iacc,inImgNames));
    }

    if (accumulator.doSensitivity()) {
      inSenNames = accumulator.inSenNameVecs()[outImgName];
      ASKAPLOG_INFO_STR(logger, " - input sensitivity images: " << inSenNames);
    }

    // set the output coordinate system and shape, based on the overlap of input images
    // for this output image

    vector<IPosition> inShapeVec;
    vector<CoordinateSystem> inCoordSysVec;
    // What fraction of the full problem does this rank
    // THe units here are channels - not polarisations.
    int myFullAllocationSize = 0;
    int myFullAllocationStart = 0;
    int myFullAllocationStop = 0;
    // Where a rank is in its allocation
    int myAllocationSize = 1;
    int myAllocationStart = 0;
    int myAllocationStop = myAllocationStart + myAllocationSize;
    //
    // this calculates the allocations for the images but stupidly does it for every image ...
    // anyway worry about that later
    //

    // I need to change this so that this is done for the current channel ....
    //
    // Lets get an example infile

    string testImage = inImgNames[0];

    //
    // lets get its shape

    const casa::IPosition shape = iacc.shape(testImage);

    ASKAPCHECK(shape.nelements()>=3,"Work with at least 3D cubes!");
    ASKAPLOG_INFO_STR(logger," - ImageAccess Shape " << shape);

    // lets calculate the allocations ...

    casa::IPosition blc(shape.nelements(),0);
    casa::IPosition trc(shape);
    originalNchan = trc[3];

    if (comms.rank() >= trc[3]) {
      ASKAPLOG_WARN_STR(logger,"Rank " << comms.rank() << " has no work to merge");
      return;
    }
    if (trc[3] % comms.nProcs() != 0) {
      ASKAPLOG_WARN_STR(logger,"Unbalanced allocation: num of ranks:" << comms.nProcs() << " not a factor of number of channels: "<< trc[3]);
    }
    if (comms.nProcs() >= trc[3]) {
      myFullAllocationSize = 1;
    }
    else {
      myFullAllocationSize = trc[3]/comms.nProcs();
    }

    myFullAllocationStart = comms.rank()*myFullAllocationSize;
    myFullAllocationStop = myFullAllocationStart + myFullAllocationSize;

    // unless last rank
    if (comms.rank() == comms.nProcs()-1) {
      myFullAllocationSize = trc[3] - myFullAllocationStart; // we are using End is Last
    }

    ASKAPLOG_INFO_STR(logger,"FullAllocation starts at " << myFullAllocationStart << " and is " << myFullAllocationSize << " in size");


    // So the plan is to iterate over each channel ...
    // Calculate the inShapes for each channel and file ....

    for (myAllocationStart = myFullAllocationStart; myAllocationStart < myFullAllocationStop; myAllocationStart = myAllocationStart +  myAllocationSize) {




      for (vector<string>::iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {

        ASKAPLOG_INFO_STR(logger,"Processing Channel " << myAllocationStart << " of input image " << *it << " which is part of output mosaick " << outImgName);


        const casa::IPosition shape = iacc.shape(*it);

        ASKAPCHECK(shape.nelements()>=3,"Work with at least 3D cubes!");

        ASKAPLOG_INFO_STR(logger," - ImageAccess Shape " << shape);

        casa::IPosition inblc(shape.nelements(),0); // input bottom left corner of this allocation
        casa::IPosition intrc(shape);
        originalNchan = intrc[3];
        myAllocationStop = myAllocationStart + myAllocationSize;
        inblc[3] = myAllocationStart;
        // change the indexing?? .... not sure I understand this ....
        // ahh maybe the shape in the number of elements but it uses 0 indexing
        intrc[0] = intrc[0]-1;
        intrc[1] = intrc[1]-1;
        intrc[2] = intrc[2]-1;
        intrc[3] = myAllocationStart + myAllocationSize-1;



        ASKAPCHECK(inblc[3]>=0 && inblc[3]<shape[3], "Start channel is outside the number of channels or negative, shape: "<<shape);
        ASKAPCHECK(trc[3]<=shape[3], "Subcube extends beyond the original cube, shape:"<<shape);

        ASKAPLOG_INFO_STR(logger, " - Corners " << "input bottom lc  = " << inblc << ", input top rc = " << intrc << "\n");
        inCoordSysVec.push_back(iacc.coordSysSlice(*it,inblc,intrc));
        // reset the shape to be the size ...
        intrc[0] = intrc[0]+1;
        intrc[1] = intrc[1]+1;
        intrc[2] = intrc[2]+1;
        intrc[3] = myAllocationSize;
        const casa::IPosition shape3(intrc);
        ASKAPLOG_INFO_STR(logger, " - Calculated Shape for this accumulator and this image is" << shape3);
        inShapeVec.push_back(shape3);


      } // got the input shapes for this output image


      // I wonder if we can re-use the accumulator ... lets find out

      casa::IPosition example = inShapeVec[0];
      ASKAPLOG_INFO_STR(logger, "Number of channels in allocations is " << example[3]);

      accumulator.setOutputParameters(inShapeVec, inCoordSysVec);
      ASKAPLOG_INFO_STR(logger, " - Output Shape " << accumulator.outShape());

      casa::IPosition sliceShape = accumulator.outShape();



      // Build the full output cube here:
      // test for master .... (and first channel)
      // this loop uses the accumulator.outShape() method - but only the spatial dimensions
      //


      if (comms.isMaster() && myAllocationStart == myFullAllocationStart) { // build this cube - does not need to loop over the outWgtNames

        ASKAPLOG_INFO_STR(logger, "++++++++++++++++++++++++++++++++++++++++++");
        ASKAPLOG_INFO_STR(logger, "Building output mosaic " << outImgName);
        ASKAPLOG_INFO_STR(logger, "++++++++++++++++++++++++++++++++++++++++++");
        casa::IPosition outShape  = accumulator.outShape();
        // has the channel dimension of the allocation - so lets fix that.
        outShape[3] = originalNchan;

        ASKAPLOG_INFO_STR(logger, " Creating output file - Shape " << outShape << " OriginalNchan " << originalNchan);
        iacc.create(outImgName, outShape, accumulator.outCoordSys());
        iacc.makeDefaultMask(outImgName);

        if (accumulator.outWgtDuplicates()[outImgName]) {
          ASKAPLOG_INFO_STR(logger, "Accumulated weight image " << outWgtName << " already written");
        } else {
          outWgtName = accumulator.outWgtNames()[outImgName];
          ASKAPLOG_INFO_STR(logger, "Writing accumulated weight image to " << outWgtName);
          iacc.create(outWgtName, outShape, accumulator.outCoordSys());
          iacc.makeDefaultMask(outWgtName);

        }
        if (accumulator.doSensitivity()) {
          outSenName = accumulator.outSenNames()[outImgName];
          ASKAPLOG_INFO_STR(logger, "Writing accumulated sensitivity image to " << outSenName);
          iacc.create(outSenName, outShape, accumulator.outCoordSys());
          iacc.makeDefaultMask(outSenName);

        }
      } // built the output cube for this image

      // here we have to loop over the channels and let everything else take care of the
      // other planes ...
      //
      // set up the output pixel arrays
      // lets do this one channel at a time


      Array<float> outPix(sliceShape,0.);
      Array<float> outWgtPix(sliceShape,0.);
      Array<bool> outMask(sliceShape,0.);
      Array<float> outSenPix;
      //
      if (accumulator.doSensitivity()) {
        outSenPix = Array<float>(sliceShape,0.);
      }

      // set up an indexing vector for the arrays
      casa::IPosition curpos(outPix.shape());
      ASKAPASSERT(curpos.nelements()>=2);
      for (uInt dim=0; dim<curpos.nelements(); ++dim) {
        curpos[dim] = 0;
      }

      // iterator over planes (e.g. freq & polarisation), regridding and accumulating weights and weighted images

      scimath::MultiDimArrayPlaneIter planeIter(accumulator.inShape());
      // loop over the input images, reading each in an adding to the output pixel arrays
      // remember this is for the current output mosaick

      for (; planeIter.hasMore(); planeIter.next()) { // this is a loop over the polarisations as well as channels
        for (uInt img = 0; img < inImgNames.size(); ++img ) {
        // set up an iterator for all directionCoordinate planes in the input images



        // short cuts
          string inImgName = inImgNames[img];
          string inWgtName, inSenName;

          ASKAPLOG_INFO_STR(logger, "Processing input image " << inImgName);
          if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED) {
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


          ASKAPLOG_INFO_STR(logger, "Shapes " << shape << " blc " << blc << " trc " << trc << " inpix " << inPix.shape());

          Array<float> inWgtPix;
          Array<float> inSenPix;

          if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED) {

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

          // test whether to simply add weighted pixels, or whether a regrid is required
          bool regridRequired = (!accumulator.coordinatesAreEqual()) ;

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

          // set the indices of any higher-order dimensions for this slice
          curpos = planeIter.position();
          // i think we have to edit the plane Iter ...
          ASKAPLOG_INFO_STR(logger, " - input slice " << curpos);

          if ( regridRequired ) {

            // load input buffer for the current plane
            accumulator.loadInputBuffers(planeIter, inPix, inWgtPix, inSenPix);
            // call regrid for any buffered images
            accumulator.regrid();
            // update the accululation arrays for this plane
            accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, curpos);

          } else {

            // Update the accumulation arrays for this plane.
            accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, inPix, inWgtPix, inSenPix, curpos);

          }


        } // over the input images for this
      } // iterated over the polarisation - the accumulator is FULL for this CHANNEL

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

      /// Since I added the scheme to incorporate a variance weight in the mosaicking
      /// I can no longer assume max weight is 1.

      float minVal, maxVal;
      IPosition minPos, maxPos;

      minMax(minVal,maxVal,minPos,maxPos,outWgtPix);
      ASKAPLOG_INFO_STR(logger, "Maximum pixel weight is " << maxVal);
      ASKAPLOG_INFO_STR(logger, "Power fraction cutoff is " << itsCutoff*itsCutoff);


      float wgtCutoff = itsCutoff * itsCutoff * maxVal;
      for( ; iterWgt != outWgtPix.end() ; iterWgt++ ) {
        if (*iterWgt >= wgtCutoff) {
          *iterMask = casa::True;
        }
        else {
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


        ASKAPLOG_INFO_STR(logger, "Ensuring serial access to cubes");


      }
      else { // this is essentially a serializer - it is required for CASA image types
      // but not FITS
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
      }
      else {

        iacc.write(outWgtName,outWgtPix,loc);
        iacc.writeMask(outWgtName,outMask,loc);
        iacc.setUnits(outWgtName,units);
        if (psf.nelements()>=3)
          iacc.setBeamInfo(outWgtName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
      }

      if (accumulator.doSensitivity()) {

        iacc.write(outSenName,outSenPix,loc);
        iacc.writeMask(outSenName,outMask,loc);
        iacc.setUnits(outSenName,units);

      if (psf.nelements()>=3)
        iacc.setBeamInfo(outSenName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
      }

      if (comms.rank() < comms.nProcs()-1) { // last rank doesnot use this method
        int buf;
        int to = comms.rank()+1;
        comms.send((void *) &buf,sizeof(int),to);
      }
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

} // end namespace askap

int main(int argc, char *argv[])
{
    askap::linmosMPIApp app;
    return app.main(argc, argv);
}
