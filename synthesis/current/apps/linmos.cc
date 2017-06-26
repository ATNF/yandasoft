/// @file linmos.cc
///
/// @brief combine a number of images as a linear mosaic
/// @details This is a utility to merge images into a mosaic. Images can be set
/// explicitly or found automatically based on input tags.
///
/// @copyright (c) 2012,2014,2015 CSIRO
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
#include "askap_synthesis.h"

// System includes
#include <sstream>
#include <typeinfo>
#include <iostream>

// other 3rd party
#include <Common/ParameterSet.h>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/images/Images/ImageRegrid.h>

// Local packages includes
#include <measurementequation/SynthesisParamsHelper.h>
#include <linmos/LinmosAccumulator.h>
#include <utils/LinmosUtils.h>

ASKAP_LOGGER(logger, ".linmos");

using namespace casa;
using namespace askap;
using namespace askap::synthesis;


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
    accessors::IImageAccess& iacc = SynthesisParamsHelper::imageHandler();

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
        for (vector<string>::iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {
            inShapeVec.push_back(iacc.shape(*it));
            inCoordSysVec.push_back(iacc.coordSys(*it));
        }
        accumulator.setOutputParameters(inShapeVec, inCoordSysVec);

        // set up the output pixel arrays
        Array<float> outPix(accumulator.outShape(),0.);
        Array<float> outWgtPix(accumulator.outShape(),0.);
        Array<bool> outMask(accumulator.outShape(),0.);

        Array<float> outSenPix;
        if (accumulator.doSensitivity()) {
            outSenPix = Array<float>(accumulator.outShape(),0.);
        }

        // set up an indexing vector for the arrays
        IPosition curpos(outPix.shape());
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

            // set the input coordinate system and shape
            accumulator.setInputParameters(iacc.shape(inImgName), iacc.coordSys(inImgName), img);

            Array<float> inPix = iacc.read(inImgName);
            Array<float> inWgtPix;
            Array<float> inSenPix;
            if (accumulator.weightType() == FROM_WEIGHT_IMAGES) {
                inWgtPix = iacc.read(inWgtName);
                ASKAPASSERT(inPix.shape() == inWgtPix.shape());
            }
            if (accumulator.doSensitivity()) {
                inSenPix = iacc.read(inSenName);
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

            } else {
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

        } // img loop (over input images)
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
        ASKAPLOG_INFO_STR(logger, "Getting PSF beam info for the output image from input number " << psfref);
        // get pixel units from the selected reference image

        string units = iacc.getUnits(inImgNames[psfref]);
        ASKAPLOG_INFO_STR(logger, "Got units as " << units);

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
        iacc.create(outImgName, accumulator.outShape(), accumulator.outCoordSys());
        iacc.write(outImgName,outPix);
        iacc.writeMask(outImgName,outMask);
        iacc.setUnits(outImgName,units);
        if (psf.nelements()>=3)
            iacc.setBeamInfo(outImgName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));

        if (accumulator.outWgtDuplicates()[outImgName]) {
            ASKAPLOG_INFO_STR(logger, "Accumulated weight image " << outWgtName << " already written");
        } else {
            ASKAPLOG_INFO_STR(logger, "Writing accumulated weight image to " << outWgtName);
            iacc.create(outWgtName, accumulator.outShape(), accumulator.outCoordSys());
            iacc.write(outWgtName,outWgtPix);
            iacc.writeMask(outWgtName,outMask);

            iacc.setUnits(outWgtName,units);
            if (psf.nelements()>=3)
                iacc.setBeamInfo(outWgtName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
        }

        if (accumulator.doSensitivity()) {
            ASKAPLOG_INFO_STR(logger, "Writing accumulated sensitivity image to " << outSenName);
            iacc.create(outSenName, accumulator.outShape(), accumulator.outCoordSys());
            iacc.write(outSenName,outSenPix);
            iacc.writeMask(outSenName,outMask);
            iacc.setUnits(outSenName,units);
            if (psf.nelements()>=3)
                iacc.setBeamInfo(outSenName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
        }

    } // ii loop (separate mosaics for different image types)

};

class LinmosApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            StatReporter stats;
            LOFAR::ParameterSet subset(config().makeSubset("linmos."));
            SynthesisParamsHelper::setUpImageHandler(subset);
            merge(subset);
            stats.logSummary();
            return 0;
        }
};

int main(int argc, char *argv[])
{
    LinmosApp app;
    return app.main(argc, argv);
}
