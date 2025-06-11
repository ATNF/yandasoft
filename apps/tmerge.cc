/// @file tmerge.cc
///
/// @brief test/development program for updates to normal equation merging.
///
/// @copyright (c) 2015 CSIRO
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
/// @author Daniel Mitchell <daniel.mitchell@csiro.au>

// Package level header file
#include "askap/askap_synthesis.h"

//# Includes
#include <askap/askap/Application.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/StatReporter.h>

#include <askap/scimath/fitting/ImagingNormalEquations.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/scimath/utils/ImageUtils.h>
#include <linmos/LinmosAccumulator.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>

ASKAP_LOGGER(logger, ".tmerge");

using namespace askap;
using namespace askap::scimath;
using namespace askap::synthesis;

class TmergeApp : public askap::Application
{
    public:

        int run(int argc, char* argv[]) final
        {
            StatReporter stats;
            SynthesisParamsHelper::setUpImageHandler(config());

            ASKAPLOG_INFO_STR(linmoslogger,
                "ASKAP program to test normal equation merging " << ASKAP_PACKAGE_VERSION);

            // initialise an image accumulator
            imagemath::LinmosAccumulator<imtype> accumulator;

            // initialise an image accessor
            accessors::IImageAccess<casacore::Float>& iacc = SynthesisParamsHelper::imageHandler();

            // initialise a shared pointer for the normal equation
            ImagingNormalEquations::ShPtr ne(new ImagingNormalEquations());
            ASKAPCHECK(ne, "normal equation is not correctly initialized");

            // give all parameters the same name
            const std::string paramName = "merged_image";

            ASKAPLOG_INFO_STR(logger, "testing normal equation merge function, tmerge");

            // get input image names
            const std::vector<std::string> images = config().getStringVector("Names", true);

            // set the output coordinate system and shape, based on the overlap of input images
            vector<IPosition> inShapeVec;
            vector<CoordinateSystem> inCoordSysVec;
            for (vector<string>::const_iterator it = images.begin(); it != images.end(); ++it) {
                inShapeVec.push_back(iacc.shape(*it));
                inCoordSysVec.push_back(iacc.coordSys(*it));
            }
            accumulator.setOutputParameters(inShapeVec, inCoordSysVec);
            IPosition outShape = accumulator.outShape();
            CoordinateSystem outCoordSys = accumulator.outCoordSys();

            // initialise the first parameter named paramName. All images will be merged onto this.
            casa::IPosition reference(outShape.nelements());
            for (uInt dim=0; dim<outShape.nelements(); ++dim) {
                reference[dim] = outCoordSys.referencePixel()[dim];
            }
            casa::Vector<imtype> outPSFVec;
            casa::Vector<imtype> outPreconVec;
            casa::Vector<imtype> outWeightVec(IPosition(1, outShape.product()),0.);
            casa::Vector<imtype> outDataVec(IPosition(1, outShape.product()),0.);
            ne->addSlice(paramName, outPSFVec, outWeightVec, outPreconVec,
                         outDataVec, outShape, reference, outCoordSys);

            // load images into the model
            for (std::vector<std::string>::const_iterator ci = images.begin(); ci != images.end(); ++ci) {

                const casa::CoordinateSystem coordSys = iacc.coordSys(*ci);

                // initialise a model and load the current image
                Params::ShPtr model;
                model.reset(new Params(true));
                SynthesisParamsHelper::loadImageParameter(*model,paramName,*ci);
                ASKAPCHECK(model->size()==1, "Expect single images");
                ASKAPLOG_INFO_STR(logger, "loaded image with " <<
                    model->valueT(paramName).nelements() << " nelements");

                // generate a normal equation and add the model
                ASKAPLOG_INFO_STR(logger, "merging as a normal eq. data vector");
                ImagingNormalEquations imageNe = ImagingNormalEquations(*model);
                const casa::IPosition imageShape(model->shape(paramName));
                for (uInt dim=0; dim<imageShape.nelements(); ++dim) {
                    reference[dim] = coordSys.referencePixel()[dim];
                }
                casa::IPosition vecShape(1, model->valueT(paramName).nelements());
                casa::Vector<imtype> inPSFVec;
                casa::Vector<imtype> inPreconVec;
                casa::Vector<imtype> inWeightVec(vecShape,1.); // overwritten if linmos is used
                casa::Vector<imtype> inDataVec(model->valueT(paramName).reform(vecShape));
                imageNe.addSlice(paramName, inPSFVec, inWeightVec, inPreconVec,
                                 inDataVec, imageShape, reference, coordSys);

                // merge with the normal equation mosaic
                ne->merge(imageNe);

            }

            // recover mosaics for output
            casa::Array<imtype> outPix(ne->dataVectorT(paramName).reform(outShape));
            casa::Array<imtype> outWgtPix(ne->normalMatrixDiagonal(paramName).reform(outShape));
            casa::Array<imtype> outSenPix(outShape,static_cast<imtype>(1.));

            std::string outName;
            casa::Array<float> scratch(outShape);

            // get pixel units and psf beam information from a reference image
            Table tmpTable(images[0]);
            string units = tmpTable.keywordSet().asString("units");
            Vector<Quantum<double> > psf;
            Vector<Quantum<double> > psftmp = iacc.beamInfo(images[0]);
            if (psftmp.nelements()<3) {
                ASKAPLOG_WARN_STR(logger, images[0] <<
                    ": beamInfo needs at least 3 elements. Not writing PSF");
            }
            else if ((psftmp[0].getValue("rad")==0) || (psftmp[1].getValue("rad")==0)) {
                ASKAPLOG_WARN_STR(logger, images[0] <<
                    ": beamInfo invalid. Not writing PSF");
            }
            else {
                psf = psftmp;
            }

            // write accumulated weight image
            outName = paramName+".wgt";
            #ifdef ASKAP_FLOAT_IMAGE_PARAMS
            scratch = outWgtPix;
            #else
            casa::convertArray<float, double>(scratch, outWgtPix);
            #endif
            iacc.create(outName, outShape, outCoordSys);
            iacc.write(outName,scratch);
            iacc.setUnits(outName,units);
            if (psf.nelements()>=3) {
                iacc.setBeamInfo(outName, psf[0].getValue("rad"),
                                 psf[1].getValue("rad"), psf[2].getValue("rad"));
            }

            // write accumulated image
            outName = paramName+".img";
            #ifdef ASKAP_FLOAT_IMAGE_PARAMS
            scratch = outPix;
            #else
            casa::convertArray<float, double>(scratch, outPix);
            #endif
            iacc.create(outName, outShape, outCoordSys);
            iacc.write(outName,scratch);
            iacc.setUnits(outName,units);
            if (psf.nelements()>=3) {
                iacc.setBeamInfo(outName, psf[0].getValue("rad"),
                                 psf[1].getValue("rad"), psf[2].getValue("rad"));
            }

            // de-weight image
            IPosition curpos(outShape);
            ASKAPASSERT(curpos.nelements()>=2);
            for (uInt dim=0; dim<curpos.nelements(); ++dim) {
                curpos[dim] = 0;
            }
            ASKAPLOG_INFO_STR(logger, "Deweighting accumulated images");
            imagemath::MultiDimArrayPlaneIter deweightIter(outShape);
            for (; deweightIter.hasMore(); deweightIter.next()) {
                curpos = deweightIter.position();
                accumulator.doSensitivity(false);
                accumulator.deweightPlane(outPix, outWgtPix, outSenPix, curpos);
            }

            // write de-weighted image
            outName = paramName+"_deweighted.img";
            #ifdef ASKAP_FLOAT_IMAGE_PARAMS
            scratch = outPix;
            #else
            casa::convertArray<float, double>(scratch, outPix);
            #endif
            iacc.create(outName, outShape, outCoordSys);
            iacc.write(outName,scratch);
            iacc.setUnits(outName,units);
            if (psf.nelements()>=3) {
                iacc.setBeamInfo(outName, psf[0].getValue("rad"),
                                 psf[1].getValue("rad"), psf[2].getValue("rad"));
            }

            stats.logSummary();
            return 0;
        }
    private:
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

int main(int argc, char *argv[])
{
    TmergeApp app;
    return app.main(argc, argv);
}
