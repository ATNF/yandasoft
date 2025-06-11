/// @file linmos.cc
///
/// @brief combine a number of images as a linear mosaic
/// @details This is a utility to merge images into a mosaic. Images can be set
/// explicitly or found automatically based on input tags.
///
/// @copyright (c) 2012,2014,2015,2021 CSIRO
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>
/// @author Daniel Mitchell <daniel.mitchell@csiro.au>

// Package level header file
#include "askap/askap_synthesis.h"

// ASKAP includes
#include <askap/utils/LinmosUtils.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/imageaccess/WeightsLog.h>
#include <askap/imagemath/linmos/LinmosAccumulator.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>

// 3rd party
#include <Common/ParameterSet.h>


ASKAP_LOGGER(logger, ".linmos");

//using namespace casa;
using namespace askap::synthesis;

namespace askap {

/// @brief do the merge
/// @param[in] parset subset with parameters
static void merge(const LOFAR::ParameterSet &parset) {

    ASKAPLOG_INFO_STR(linmoslogger, "ASKAP linear mosaic task " << ASKAP_PACKAGE_VERSION);
    ASKAPLOG_INFO_STR(linmoslogger, "Parset parameters:\n" << parset);

    // initialise an image accumulator
    imagemath::LinmosAccumulator<float> accumulator;

    // load the parset
    if ( !accumulator.loadParset(parset) ) return;

    // initialise an image accessor
    accessors::IImageAccess<casacore::Float>& iacc = SynthesisParamsHelper::imageHandler();

    // get the imageHistory
    const std::vector<std::string> historyLines = parset.getStringVector("imageHistory",{},false);

    // get the list of keywords to copy from the input
    // Set some defaults for simple beam mosaic, won't be correct for more complex cases
    const std::vector<std::string> keywordsToCopy = parset.getStringVector("keywordsToCopy",
        {"TELESCOP","PROJECT","SBID","DATE-OBS","DURATION"});

    // loop over the mosaics, reading each in and adding to the output pixel arrays
    vector<string> inImgNames, inWgtNames, inSenNames, inStokesINames;
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

        if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED) {
            inWgtNames = accumulator.inWgtNameVecs()[outImgName];
            if (inWgtNames.size()==0) {
              // want weights, but have not specified any files - try to read from image headers
              bool useWtFromHdr = true;
              for (auto i : inImgNames) {
                useWtFromHdr = readWeightsTable(i).size()>0;
                if (!useWtFromHdr) {
                  break;
                }
              }
              ASKAPCHECK(useWtFromHdr, "Image weighting requested, but no weights found in headers or files");
              ASKAPLOG_INFO_STR(linmoslogger,"Reading weights from input images");
              accumulator.setUseWtFromHdr(useWtFromHdr);
            } else if (accumulator.useWeightsLog()) {
                ASKAPLOG_INFO_STR(logger, " - input weightslog files: " << inWgtNames);
            } else {
                ASKAPLOG_INFO_STR(logger, " - input weights images: " << inWgtNames);
            }
        }

        if (accumulator.weightType() == FROM_BP_MODEL || accumulator.weightType() == COMBINED) {
            accumulator.setBeamCentres(loadBeamCentres(parset,iacc,inImgNames));
        }
        if (accumulator.doSensitivity()) {
            inSenNames = accumulator.inSenNameVecs()[outImgName];
            ASKAPLOG_INFO_STR(logger, " - input sensitivity images: " << inSenNames);
        }

        if (accumulator.doLeakage()) {
            inStokesINames = accumulator.inStokesINameVecs()[outImgName];
            if (inStokesINames.size()>0) {
                ASKAPLOG_INFO_STR(logger, " - and using input Stokes I images: " << inStokesINames);
            }
        }

        // set the output coordinate system and shape, based on the overlap of input images
        vector<IPosition> inShapeVec;
        vector<CoordinateSystem> inCoordSysVec;
        for (vector<string>::iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {
            inShapeVec.push_back(iacc.shape(*it));
            inCoordSysVec.push_back(iacc.coordSys(*it));
        }
        accumulator.setOutputParameters(inShapeVec, inCoordSysVec);

        // set up the output pixel arrays
        Array<float> outPix(accumulator.outShape(),0.f);
        Array<float> outWgtPix(accumulator.outShape(),0.f);
        Array<bool> outMask(accumulator.outShape(),false);

        Array<float> outSenPix;
        if (accumulator.doSensitivity()) {
            outSenPix = Array<float>(accumulator.outShape(),0.f);
        }

        // set up an indexing vector for the arrays
        IPosition curpos(outPix.ndim(),0);
        ASKAPASSERT(curpos.nelements()>=2);

        // loop over the input images, reading each in an adding to the output pixel arrays
        for (uInt img = 0; img < inImgNames.size(); ++img ) {

            // short cuts
            string inImgName = inImgNames[img];
            string inWgtName, inSenName, inStokesIName;

            ASKAPLOG_INFO_STR(logger, "Processing input image " << inImgName);
            if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED) {
              if (!accumulator.useWtFromHdr()) {
                inWgtName = inWgtNames[img];
                if (accumulator.useWeightsLog()) {
                    ASKAPLOG_INFO_STR(logger, " - and input weightslog " << inWgtName);
                } else {
                    ASKAPLOG_INFO_STR(logger, " - and input weight image " << inWgtName);
                }
              }
            }
            if (accumulator.doSensitivity()) {
                inSenName = inSenNames[img];
                ASKAPLOG_INFO_STR(logger, " - and input sensitivity image " << inSenName);
            }

            if (accumulator.doLeakage() && inStokesINames.size()>img) {
                inStokesIName = inStokesINames[img];
                ASKAPLOG_INFO_STR(logger, " - and input Stokes I image " << inStokesIName);
            }

            // set the input coordinate system and shape
            accumulator.setInputParameters(iacc.shape(inImgName), iacc.coordSys(inImgName), img);

            Array<float> inPix = iacc.read(inImgName);

            if (parset.getBool("removebeam",false)) {

                Array<float> taylor0;
                Array<float> taylor1;
                Array<float> taylor2;

                ASKAPLOG_INFO_STR(linmoslogger, "Scaling Taylor terms -- inImage = " << inImgNames[img]);
                // need to get all the taylor terms for this image
                string ImgName = inImgName;
                int inPixIsTaylor = 0;
                for (int n = 0; n < accumulator.numTaylorTerms(); ++n) {
                    const string taylorN = "taylor." + boost::lexical_cast<string>(n);
                    // find the taylor.0 image for this image
                    size_t pos0 = ImgName.find(taylorN);
                    if (pos0!=string::npos) {
                        ImgName.replace(pos0, taylorN.length(), accumulator.taylorTag());
                        ASKAPLOG_INFO_STR(linmoslogger, "This is a Taylor " << n << " image");
                        inPixIsTaylor = n;
                        break;
                    }

                    // now go through each taylor term

                }


                ASKAPLOG_INFO_STR(linmoslogger, "To avoid altering images on disk re-reading the Taylor terms");
                for (int n = 0; n < accumulator.numTaylorTerms(); ++n) {

                    size_t pos0 = ImgName.find(accumulator.taylorTag());
                    if (pos0!=string::npos) {
                        const string taylorN = "taylor." + boost::lexical_cast<string>(n);
                        ImgName.replace(pos0, taylorN.length(), taylorN);


                        switch (n)
                        {
                        case 0:
                            ASKAPLOG_INFO_STR(linmoslogger, "Reading -- Taylor0");
                            ASKAPLOG_INFO_STR(linmoslogger, "Reading -- inImage = " << ImgName);
                            taylor0 = iacc.read(ImgName);
                            ASKAPLOG_INFO_STR(linmoslogger, "Shape -- " << taylor0.shape());
                            break;
                        case 1:
                            ASKAPLOG_INFO_STR(linmoslogger, "Reading -- Taylor1");
                            ASKAPLOG_INFO_STR(linmoslogger, "Reading -- inImage = " << ImgName);
                            taylor1 = iacc.read(ImgName);
                            ASKAPLOG_INFO_STR(linmoslogger, "Shape -- " << taylor1.shape());
                            break;
                        case 2:
                            ASKAPLOG_INFO_STR(linmoslogger, "Reading -- Taylor2");
                            ASKAPLOG_INFO_STR(linmoslogger, "Reading -- inImage = " << ImgName);
                            taylor2 = iacc.read(ImgName);
                            ASKAPLOG_INFO_STR(linmoslogger, "Shape -- " << taylor2.shape());
                            break;

                        }
                        ImgName.replace(pos0, accumulator.taylorTag().length(), accumulator.taylorTag());
                    }

                    // now go through each taylor term

                }


                casa::IPosition thispos(taylor0.shape().nelements(),0);
                ASKAPLOG_INFO_STR(logger, " removing Beam for Taylor terms - slice " << thispos);
                accumulator.removeBeamFromTaylorTerms(taylor0,taylor1,taylor2,thispos,iacc.coordSys(inImgName));


                // now we need to set the inPix to be the scaled version
                // Note this means we are reading the Taylor terms 3 times for every
                // read. But I'm not sure this matters.

                switch (inPixIsTaylor)
                {
                case 0:
                    inPix = taylor0;
                    break;
                case 1:
                    inPix = taylor1;
                    break;
                case 2:
                    inPix = taylor2;
                    break;
                }

            }

            if (parset.getBool("removeleakage",false)) {
                ASKAPCHECK(inPix.shape()[2]==1,"Pol axis should have size 1 for removeleakage");
                // only do this if we're processing a Q, U or V image
                int pol = 0;
                size_t pos = inImgName.find(".q.");
                bool found = false;
                if (pos != string::npos) {
                    found = true;
                    pol = 1;
                }
                if (!found) {
                    pos = inImgName.find(".u.");
                    if (pos != string::npos) {
                        found = true;
                        pol = 2;
                    }
                }
                if (!found) {
                    pos = inImgName.find(".v.");
                    if (pos !=string::npos) {
                        found = true;
                        pol = 3;
                    }
                }
                if (found) {
                    // find corresponding Stokes I image
                    if (inStokesIName=="") {
                        //try to find it
                        inStokesIName = inImgName;
                        inStokesIName.replace(pos,3,".i.");
                    }
                    Array<float> stokesI = iacc.read(inStokesIName);
                    ASKAPCHECK(stokesI.shape()==inPix.shape(),"Stokes I and Pol image shapes don't match");

                    // do leakage correction
                    ASKAPLOG_INFO_STR(logger," removing Stokes I leakage using "<<inStokesIName);

                    // iterator over planes (e.g. freq & polarisation)
                    for (imagemath::MultiDimArrayPlaneIter planeIter(accumulator.inShape()); planeIter.hasMore(); planeIter.next()) {
                        // set the indices of any higher-order dimensions for this slice
                        curpos = planeIter.position();
                        // removeLeakage works on single frequency planes
                        Array<float> inPlane = planeIter.getPlane(inPix);
                        Array<float> stokesIplane = planeIter.getPlane(stokesI);
                        accumulator.removeLeakage(inPlane,stokesIplane,pol,curpos,iacc.coordSys(inImgName));
                    }
                } else {
                    ASKAPLOG_WARN_STR(logger,"Skipping removeLeakage - cannot determine polarisation of input");
                }
            }

            Array<float> inWgtPix;
            Array<float> inSenPix;
            if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED) {
                if (accumulator.useWeightsLog()) {
                    ASKAPLOG_INFO_STR(logger,"Reading weights log file :"<< inWgtName);
                    accessors::WeightsLog wtlog(inWgtName);
                    wtlog.read();

                    inWgtPix.resize(inPix.shape());
                    ASKAPCHECK(inWgtPix.ndim()==4,"Cannot use weightslog for non-standard cubes - need axes ra,dec,pol,freq");
                    // iterator over planes (e.g. freq & polarisation)
                    for (imagemath::MultiDimArrayPlaneIter planeIter(accumulator.inShape()); planeIter.hasMore(); planeIter.next()) {
                        // set the indices of any higher-order dimensions for this slice
                        curpos = planeIter.position();
                        // set weight value for this frequency plane
                        planeIter.getPlane(inWgtPix) = wtlog.weight(curpos[3]);
                    }
                } else if (accumulator.useWtFromHdr()) {
                    ASKAPLOG_INFO_STR(logger,"Reading weights from input image file :"<< inImgName);
                    const casacore::Vector<casacore::Float> wts = readWeightsTable(inImgName);
                    const size_t size = wts.size();
                    if (size > 1) {
                      inWgtPix.resize(inPix.shape());
                      ASKAPCHECK(inWgtPix.ndim()==4,"Cannot use weightslog for non-standard cubes - need axes ra,dec,pol,freq");
                      for (imagemath::MultiDimArrayPlaneIter planeIter(accumulator.inShape()); planeIter.hasMore(); planeIter.next()) {
                        // set the indices of any higher-order dimensions for this slice
                        curpos = planeIter.position();
                        // set weight value for this frequency plane
                        ASKAPCHECK(curpos[3] < size,"Not enough channels in weights table");
                        planeIter.getPlane(inWgtPix) = wts(curpos[3]);
                      }
                    } else if (size == 1) {
                        inWgtPix.resize(inPix.shape());
                        inWgtPix = wts(0);
                    }
                    ASKAPCHECK(inWgtPix.nelements()>0,"No weights found in image header or extension for image: "<<inImgName);
                } else {
                    // use weights images
                    inWgtPix = iacc.read(inWgtName);
                    ASKAPASSERT(inPix.shape() == inWgtPix.shape());
                }
            }
            if (accumulator.doSensitivity()) {
                inSenPix = iacc.read(inSenName);
                ASKAPASSERT(inPix.shape() == inSenPix.shape());
            }

            // set up an iterator for all directionCoordinate planes in the input image
            imagemath::MultiDimArrayPlaneIter planeIter(accumulator.inShape());

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
                // are those of the previous iteration correctly freed?
                accumulator.initialiseOutputBuffers();
                accumulator.initialiseInputBuffers();
            }
            else {
                ASKAPLOG_INFO_STR(logger, " - not regridding -- input pixel grid is the same as the output");
                // not regridding so point output image buffers at the input buffers
                accumulator.initialiseInputBuffers();
                accumulator.redirectOutputBuffers();
            }

            // iterator over planes (e.g. freq & polarisation), regridding and accumulating weights and weighted images
            for (; planeIter.hasMore(); planeIter.next()) {

                // set the indices of any higher-order dimensions for this slice
                curpos = planeIter.position();

                ASKAPLOG_INFO_STR(logger, " - slice " << curpos);

                // load input buffer for the current plane
                accumulator.loadAndWeightInputBuffers(curpos, inPix, inWgtPix, inSenPix);

                if ( regridRequired ) {
                    // call regrid for any buffered images
                    accumulator.regrid();
                }

                // update the accululation arrays for this plane
                accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, curpos);

            }

        } // img loop (over input images)

        //build the mask
        //use the outWgtPix to define the mask

        float itsCutoff = 0.01;

        if (parset.isDefined("cutoff")) itsCutoff = parset.getFloat("cutoff");
        if (parset.isDefined("finalcutoff")) itsCutoff = parset.getFloat("finalcutoff");


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
        /// I can no longer assume max weight is 1

        // Need to do this plane by plane, since weight varies by channel

        for (imagemath::MultiDimArrayPlaneIter outIter(accumulator.outShape()); outIter.hasMore(); outIter.next()) {

            float minVal, maxVal;
            IPosition minPos, maxPos;
            Array<float> outWgtPlane = outIter.getPlane(outWgtPix);
            Array<bool> outMaskPlane = outIter.getPlane(outMask);
            minMax(minVal,maxVal,minPos,maxPos,outWgtPlane);
            ASKAPLOG_INFO_STR(logger, "Maximum pixel weight is " << maxVal << ", slice "<<outIter.position());

            float wgtCutoff = itsCutoff * itsCutoff * maxVal;

            for(size_t i=0;i<outMaskPlane.size();i++){
                if (outWgtPlane.data()[i] >= wgtCutoff) {
                    outMaskPlane.data()[i] = casa::True;
                } else {
                    outMaskPlane.data()[i] = casa::False;
                    setNaN(outWgtPlane.data()[i]);
                }
            }
        }

        ASKAPLOG_INFO_STR(logger, "Power fraction cutoff is " << itsCutoff*itsCutoff);

        // deweight the image pixels
        // use another iterator to loop over planes
        ASKAPLOG_INFO_STR(logger, "Deweighting accumulated images");
        imagemath::MultiDimArrayPlaneIter deweightIter(accumulator.outShape());
        for (; deweightIter.hasMore(); deweightIter.next()) {
            curpos = deweightIter.position();
            accumulator.deweightPlane(outPix, outWgtPix, outSenPix, curpos);
        }

        // set one of the input images as a reference for metadata (the first by default)
        uint psfref = 0;
        if (parset.isDefined("psfref")) psfref = parset.getUint("psfref");
        ASKAPLOG_INFO_STR(logger, "Getting PSF beam info for the output image from input number " << psfref);
        // get pixel units from the selected reference image

        string units = iacc.getUnits(inImgNames[psfref]);
        ASKAPLOG_INFO_STR(logger, "Got units as " << units);

        // get psf beam information from the selected reference image
        Vector<Quantum<double>> psf = iacc.beamInfo(inImgNames[psfref]);
        bool psfValid = (psf.nelements()==3) && (psf[0].getValue("rad")>0) && (psf[1].getValue("rad")>0);
        accessors::BeamList refBeamList = iacc.beamList(inImgNames[psfref]);

        // write accumulated images and weight images
        ASKAPLOG_INFO_STR(logger, "Writing accumulated image to " << outImgName);
        iacc.create(outImgName, accumulator.outShape(), accumulator.outCoordSys());
        copyKeywords(outImgName, accumulator.getReference(outImgName), keywordsToCopy);
        iacc.addHistory(outImgName,historyLines);
        iacc.write(outImgName,outPix);
        iacc.writeMask(outImgName,outMask);
        iacc.setUnits(outImgName,units);
        if (psfValid) {
            iacc.setBeamInfo(outImgName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
        }
        if (!refBeamList.empty()) {
            iacc.setBeamInfo(outImgName,refBeamList);
        }

        if (accumulator.outWgtDuplicates()[outImgName]) {
            ASKAPLOG_INFO_STR(logger, "Accumulated weight image " << outWgtName << " already written");
        } else {
            ASKAPLOG_INFO_STR(logger, "Writing accumulated weight image to " << outWgtName);
            iacc.create(outWgtName, accumulator.outShape(), accumulator.outCoordSys());
            copyKeywords(outWgtName, accumulator.getReference(outImgName), keywordsToCopy);
            iacc.addHistory(outWgtName,historyLines);
            iacc.write(outWgtName,outWgtPix);
            iacc.writeMask(outWgtName,outMask);
            iacc.setUnits(outWgtName,units);
            if (psfValid) {
                iacc.setBeamInfo(outWgtName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
            }
            if (!refBeamList.empty()) {
                iacc.setBeamInfo(outWgtName,refBeamList);
            }
        }

        if (accumulator.doSensitivity()) {
            ASKAPLOG_INFO_STR(logger, "Writing accumulated sensitivity image to " << outSenName);
            iacc.create(outSenName, accumulator.outShape(), accumulator.outCoordSys());
            copyKeywords(outSenName, accumulator.getReference(outImgName), keywordsToCopy);
            iacc.addHistory(outSenName,historyLines);
            iacc.write(outSenName,outSenPix);
            iacc.writeMask(outSenName,outMask);
            iacc.setUnits(outSenName,units);
            if (psfValid) {
                iacc.setBeamInfo(outSenName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
            }
            if (!refBeamList.empty()) {
                iacc.setBeamInfo(outSenName,refBeamList);
            }
        }
        // Save table of mosaic pointing centres & beamsizes?
        saveMosaicTable(outImgName,inImgNames,accumulator.getBeamCentres());

    } // ii loop (separate mosaics for different image types)

};

class LinmosApp : public askap::Application
{
    public:
        int run(int argc, char* argv[]) final
        {
            StatReporter stats;
            LOFAR::ParameterSet subset(config().makeSubset("linmos."));
            SynthesisParamsHelper::setUpImageHandler(subset);
            merge(subset);
            stats.logSummary();
            return 0;
        }

    private:
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

} // end namespace askap

int main(int argc, char *argv[])
{
    askap::LinmosApp app;
    return app.main(argc, argv);
}
