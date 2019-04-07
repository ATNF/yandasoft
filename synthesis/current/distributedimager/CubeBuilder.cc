/// @file CubeBuilder.cc
///
/// @copyright (c) 2013 CSIRO
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>

// Include own header file first
#include <distributedimager/CubeBuilder.h>

// Include package level header file
#include <askap_synthesis.h>

// System includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

// ASKAPsoft includes
#include <askap/AskapError.h>
#include <askap/AskapLogging.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <imageaccess/ImageAccessFactory.h>
#include <Common/ParameterSet.h>
#include <utils/PolConverter.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/measures/Measures/Stokes.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/casa/Quanta/Unit.h>
#include <casacore/casa/Quanta/QC.h>

ASKAP_LOGGER(logger, ".CubeBuilder");

using namespace askap::cp;
using namespace casa;
using namespace std;
using namespace askap::synthesis;
CubeBuilder::CubeBuilder(const LOFAR::ParameterSet& parset,const std::string& name) {
    // as long as the cube exists all should be fine
    ASKAPLOG_INFO_STR(logger, "Instantiating Cube Builder with existing cube ");
    itsCube = accessors::imageAccessFactory(parset);

    vector<string> filenames;
    if (parset.isDefined("Images.Names")) {
        filenames = parset.getStringVector("Images.Names", true);
        itsFilename = filenames[0];
    }
    else if(parset.isDefined("Images.name")) {
        itsFilename = parset.getString("Images.name");
    }
    else {
        ASKAPLOG_ERROR_STR(logger, "Could not find the image name(s) ");
    }
    ASKAPCHECK(itsFilename.substr(0,5)=="image",
               "Images.name (Names) must start with 'image' starts with " << itsFilename.substr(0,5));

    // If necessary, replace "image" with _name_ (e.g. "psf", "weights")
    // unless name='restored', in which case we append ".restored"
    if (!name.empty()) {
        if (name == "restored") {
            itsFilename = itsFilename + ".restored";
        } else {
            const string orig = "image";
            const size_t f = itsFilename.find(orig);
            itsFilename.replace(f, orig.length(), name);
        }
    }
    
    ASKAPLOG_INFO_STR(logger, "Instantiated Cube Builder with existing cube " << itsFilename);
}
CubeBuilder::CubeBuilder(const LOFAR::ParameterSet& parset,
                         const casacore::uInt nchan,
                         const casacore::Quantity& f0,
                         const casacore::Quantity& inc,
                         const std::string& name)
{
    ASKAPLOG_INFO_STR(logger, "Instantiating Cube Builder by creating cube ");
    itsCube = accessors::imageAccessFactory(parset);

    vector<string> filenames;
    if (parset.isDefined("Images.Names")) {
        filenames = parset.getStringVector("Images.Names", true);
        itsFilename = filenames[0];
    }
    else if(parset.isDefined("Images.name")) {
        itsFilename = parset.getString("Images.name");
    }
    else {
        ASKAPLOG_ERROR_STR(logger, "Could not find the image name(s) ");
    }
    ASKAPCHECK(itsFilename.substr(0,5)=="image",
               "Images.name (Names) must start with 'image' starts with " << itsFilename.substr(0,5));

    // If necessary, replace "image" with _name_ (e.g. "psf", "weights")
    // unless name='restored', in which case we append ".restored"
    if (!name.empty()) {
        if (name == "restored") {
            itsFilename = itsFilename + ".restored";
        } else {
            const string orig = "image";
            const size_t f = itsFilename.find(orig);
            itsFilename.replace(f, orig.length(), name);
        }
    }

    const std::string restFreqString = parset.getString("Images.restFrequency", "-1.");
    if (restFreqString == "HI") {
#ifdef HAVE_CASACORE3
        itsRestFrequency = casacore::QC::HI();
#else
        itsRestFrequency = casacore::QC::HI;
#endif // HAVE_CASACORE3
    } else {
        itsRestFrequency = SynthesisParamsHelper::convertQuantity(restFreqString, "Hz");
    }

    // Polarisation
    const std::vector<std::string>
        stokesVec = parset.getStringVector("Images.polarisation", std::vector<std::string>(1,"I"));
    // there could be many ways to define stokes, e.g. ["XX YY"] or ["XX","YY"] or "XX,YY"
    // to allow some flexibility we have to concatenate all elements first and then
    // allow the parser from PolConverter to take care of extracting the products.
    std::string stokesStr;
    for (size_t i=0; i<stokesVec.size(); ++i) {
        stokesStr += stokesVec[i];
    }
    itsStokes = scimath::PolConverter::fromString(stokesStr);
    const casacore::uInt npol=itsStokes.size();

    // Get the image shape
    const vector<casacore::uInt> imageShapeVector = parset.getUintVector("Images.shape");
    const casacore::uInt nx = imageShapeVector[0];
    const casacore::uInt ny = imageShapeVector[1];
    const casacore::IPosition cubeShape(4, nx, ny, npol, nchan);

    // Use a tile shape appropriate for plane-by-plane access
    casacore::IPosition tileShape(cubeShape.nelements(), 1);
    tileShape(0) = 256;
    tileShape(1) = 256;

    const casacore::CoordinateSystem csys = createCoordinateSystem(parset, nx, ny, f0, inc);

    ASKAPLOG_INFO_STR(logger, "Creating Cube " << itsFilename <<
                       " with shape [xsize:" << nx << " ysize:" << ny <<
                       " npol:" << npol << " nchan:" << nchan <<
                       "], f0: " << f0.getValue("MHz") << " MHz, finc: " <<
                       inc.getValue("kHz") << " kHz");

    itsCube->create(itsFilename, cubeShape, csys);

    // default flux units are Jy/pixel. If we set the restoring beam
    // later on, can set to Jy/beam
    itsCube->setUnits(itsFilename,"Jy/pixel");

    ASKAPLOG_INFO_STR(logger, "Instantiated Cube Builder by creating cube " << itsFilename);
}

CubeBuilder::~CubeBuilder()
{
}

void CubeBuilder::writeSlice(const casacore::Array<float>& arr, const casacore::uInt chan)
{
    casacore::IPosition where(4, 0, 0, 0, chan);
    itsCube->write(itsFilename,arr, where);
}

casacore::CoordinateSystem
CubeBuilder::createCoordinateSystem(const LOFAR::ParameterSet& parset,
                                    const casacore::uInt nx,
                                    const casacore::uInt ny,
                                    const casacore::Quantity& f0,
                                    const casacore::Quantity& inc)
{
    CoordinateSystem coordsys;
    const vector<string> dirVector = parset.getStringVector("Images.direction");
    const vector<string> cellSizeVector = parset.getStringVector("Images.cellsize");


    // Direction Coordinate
    {
        Matrix<Double> xform(2, 2);
        xform = 0.0;
        xform.diagonal() = 1.0;
        const Quantum<Double> ra = asQuantity(dirVector.at(0), "deg");
        const Quantum<Double> dec = asQuantity(dirVector.at(1), "deg");
        ASKAPLOG_DEBUG_STR(logger, "Direction: " << ra.getValue() << " degrees, "
                           << dec.getValue() << " degrees");

        const Quantum<Double> xcellsize = asQuantity(cellSizeVector.at(0), "arcsec") * -1.0;
        const Quantum<Double> ycellsize = asQuantity(cellSizeVector.at(1), "arcsec");
        ASKAPLOG_DEBUG_STR(logger, "Cellsize: " << xcellsize.getValue()
                           << " arcsec, " << ycellsize.getValue() << " arcsec");

        casacore::MDirection::Types type;
        casacore::MDirection::getType(type, dirVector.at(2));
        const DirectionCoordinate radec(type, Projection(Projection::SIN),
                                        ra, dec, xcellsize, ycellsize,
                                        xform, nx / 2, ny / 2);

        coordsys.addCoordinate(radec);
    }

    // Stokes Coordinate
    {

        // To make a StokesCoordinate, need to convert the StokesTypes
        // into integers explicitly
        casacore::Vector<casacore::Int> stokes(itsStokes.size());
        for(unsigned int i=0;i<stokes.size();i++){
            stokes[i] = itsStokes[i];
        }
        const StokesCoordinate stokescoord(stokes);
        coordsys.addCoordinate(stokescoord);

    }
    // Spectral Coordinate
    {
        const Double refPix = 0.0;  // is the reference pixel

        MFrequency::Types freqRef=MFrequency::TOPO;
        // setup frequency frame
        const std::string freqFrame = parset.getString("freqframe","topo");
        if (freqFrame == "topo") {
            ASKAPLOG_INFO_STR(logger, "Image cube frequencies will be treated as topocentric");
            freqRef = casacore::MFrequency::TOPO;
        } else if (freqFrame == "lsrk") {
            ASKAPLOG_INFO_STR(logger, "Image cube frequencies will be treated as lsrk");
            freqRef = casacore::MFrequency::LSRK;
        } else if (freqFrame == "bary") {
        ASKAPLOG_INFO_STR(logger, "Image cube frequencies will be treated as barycentric");
            freqRef = casacore::MFrequency::BARY;
        } else {
            ASKAPTHROW(AskapError, "Unsupported frequency frame "<<freqFrame);
        }
    
        
        SpectralCoordinate sc(freqRef, f0, inc, refPix);

        // add rest frequency, but only if requested, and only for
        // image.blah, residual.blah, image.blah.restored
        if (itsRestFrequency.getValue("Hz") > 0.) {
            if ((itsFilename.find("image.") != string::npos) ||
                    (itsFilename.find("residual.") != string::npos)) {

                if (!sc.setRestFrequency(itsRestFrequency.getValue("Hz"))) {
                    ASKAPLOG_ERROR_STR(logger, "Could not set the rest frequency to " <<
                                       itsRestFrequency.getValue("Hz") << "Hz");
                }
            }
        }

        coordsys.addCoordinate(sc);
    }

    return coordsys;
}

void CubeBuilder::addBeam(casacore::Vector<casacore::Quantum<double> > &beam)
{
        itsCube->setBeamInfo(itsFilename,beam[0].getValue("rad"),beam[1].getValue("rad"),beam[2].getValue("rad"));
        setUnits("Jy/beam");
}

void CubeBuilder::setUnits(const std::string &units)
{
    itsCube->setUnits(itsFilename,units);
}
