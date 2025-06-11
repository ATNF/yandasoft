/// @file CubeBuilder.tcc
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
#include <askap/distributedimager/CubeBuilder.h>

// Include package level header file
// #include <askap/askap_synthesis.h>

// System includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

// ASKAPsoft includes
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapLogging.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/imageaccess/ImageAccessFactory.h>
#include <askap/imageaccess/CasaImageAccess.h>
#include <Common/ParameterSet.h>
#include <askap/scimath/utils/PolConverter.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/LinearCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/fits/FITS/FITSDateUtil.h>
#include <casacore/measures/Measures/Stokes.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/Quanta/Unit.h>
#include <casacore/casa/Quanta/QC.h>

ASKAP_LOGGER(CubeBuilderLogger, ".CubeBuilder");

using namespace casa;
using namespace std;
using namespace askap::synthesis;

namespace askap {
namespace cp {

template <> inline
casacore::CoordinateSystem
CubeBuilder<casacore::Complex>::createCoordinateSystem(const LOFAR::ParameterSet& parset,
                                    const casacore::uInt nx,
                                    const casacore::uInt ny,
                                    const casacore::Quantity& f0,
                                    const casacore::Quantity& inc)
{
  // This specialisation will be just for the grid cubes. After all who else wants a complex cube
  // The coordinates are therefore those of the UV grid.
  return createUVCoordinateSystem(parset, nx, ny, f0, inc);
}

template <class T> inline
casacore::CoordinateSystem
CubeBuilder<T>::createUVCoordinateSystem(const LOFAR::ParameterSet& parset,
                                    const casacore::uInt nx,
                                    const casacore::uInt ny,
                                    const casacore::Quantity& f0,
                                    const casacore::Quantity& inc)
{
    // This specialisation will be just for the grid cubes. After all who else wants a complex cube
    // The coordinates are therefore those of the UV grid.
    CoordinateSystem coordsys;
    const vector<string> dirVector = parset.getStringVector("Images.direction");
    const vector<string> cellSizeVector = parset.getStringVector("Images.cellsize");
    // Get the image shape (I get nx and ny passed to me so I probably dont need this)
    const vector<casacore::uInt> imageShapeVector = parset.getUintVector("Images.shape");
    ASKAPASSERT(imageShapeVector.size() >= 2u);


    const Quantum<Double> xcellsize = asQuantity(cellSizeVector.at(0), "arcsec");
    const Quantum<Double> ycellsize = asQuantity(cellSizeVector.at(1), "arcsec");
    // UV cellsize is probably
    casacore::Vector<casacore::Double> UVCellSize(2);
    ASKAPDEBUGASSERT(imageShapeVector.size()>=2);
    UVCellSize[0] = 1./(xcellsize.getValue("rad")*(imageShapeVector[0]));
    UVCellSize[1] = 1./(ycellsize.getValue("rad")*(imageShapeVector[1]));

    // Direction Coordinate - which is now UV in metres - and now Linear
    {
        casacore::Vector<casacore::String> units(2, casacore::String("m"));

        casacore::Vector<casacore::String> names(2);
        names[0] = "u";
        names[1] = "v";

        casacore::Matrix<casacore::Double> xform(2, 2, 0.);
        xform.diagonal() = 1.0;

        casacore::Vector<casacore::Double> crpix(2);
        casacore::Vector<casacore::Double> crval(2, 0.);

        for (size_t dim = 0; dim < 2; ++dim) {
            // note from MV: this way it will be an integer division. Do we want that?
            crpix[dim] = imageShapeVector[dim]/2;
        }

        LinearCoordinate lc(names, units, crval, UVCellSize,xform,crpix);

        coordsys.addCoordinate(lc);
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
            ASKAPLOG_INFO_STR(CubeBuilderLogger, "Image cube frequencies will be treated as topocentric");
            freqRef = casacore::MFrequency::TOPO;
        } else if (freqFrame == "lsrk") {
            ASKAPLOG_INFO_STR(CubeBuilderLogger, "Image cube frequencies will be treated as lsrk");
            freqRef = casacore::MFrequency::LSRK;
        } else if (freqFrame == "bary") {
        ASKAPLOG_INFO_STR(CubeBuilderLogger, "Image cube frequencies will be treated as barycentric");
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
                    ASKAPLOG_ERROR_STR(CubeBuilderLogger, "Could not set the rest frequency to " <<
                                       itsRestFrequency.getValue("Hz") << "Hz");
                }
            }
        }

        coordsys.addCoordinate(sc);
    }

    return coordsys;
}

template <> inline
CubeBuilder<casacore::Complex>::CubeBuilder(const LOFAR::ParameterSet& parset,const std::string& name) {
    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiating Cube Builder by co-opting existing Complex cube");
    itsFilename = makeImageName(parset, name);
    boost::shared_ptr<CasaImageAccess<casacore::Complex> > iaCASA(new CasaImageAccess<casacore::Complex>());
    itsCube = iaCASA;
}

template <class T>
CubeBuilder<T>::CubeBuilder(const LOFAR::ParameterSet& parset,const std::string& name) {
    // as long as the cube exists all should be fine
    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiating Cube Builder co-opting existing cube");
    itsFilename = makeImageName(parset, name);
    itsCube = accessors::imageAccessFactory(parset);
    // Check whether image param is stored at a lower resolution
    if (parset.isDefined("Images.extraoversampling")) {
        itsExtraOversamplingFactor = parset.getFloat("Images.extraoversampling");
        // The parameter should only be defined if has a legitimate value (is set by the code). Check anyway.
        ASKAPDEBUGASSERT(*itsExtraOversamplingFactor > 1.);
        ASKAPLOG_INFO_STR(CubeBuilderLogger, "Using extraoversampling " << *itsExtraOversamplingFactor);
    }
}

template <> inline
CubeBuilder<casacore::Complex>::CubeBuilder(const LOFAR::ParameterSet& parset,
                         const casacore::uInt nchan,
                         const casacore::Quantity& f0,
                         const casacore::Quantity& inc,
                         const std::string& name,
                         const bool uvcoord)
{
    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiating Cube Builder by creating Complex cube");
    itsFilename = makeImageName(parset, name);
    boost::shared_ptr<CasaImageAccess<casacore::Complex> > iaCASA(new CasaImageAccess<casacore::Complex>());
    itsCube = iaCASA;

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

    // Check whether image param is stored at a lower resolution
    if (parset.isDefined("Images.extraoversampling")) {
        itsExtraOversamplingFactor = parset.getFloat("Images.extraoversampling");
        // The parameter should only be defined if has a legitimate value (is set by the code). Check anyway.
        ASKAPDEBUGASSERT(*itsExtraOversamplingFactor > 1.);
        ASKAPLOG_INFO_STR(CubeBuilderLogger, "Using extraoversampling " << *itsExtraOversamplingFactor);
    }

    // Get the image shape
    const vector<casacore::uInt> imageShapeVector = parset.getUintVector("Images.shape");
    casacore::uInt nx = imageShapeVector[0];
    casacore::uInt ny = imageShapeVector[1];
    if (itsExtraOversamplingFactor) {
        const casacore::IPosition fullShape =
            scimath::PaddingUtils::paddedShape(casacore::IPosition(2,nx,ny),*itsExtraOversamplingFactor);
        nx = fullShape[0];
        ny = fullShape[1];
    }
    const casacore::IPosition cubeShape(4, nx, ny, npol, nchan);

    const casacore::CoordinateSystem csys = createCoordinateSystem(parset, nx, ny, f0, inc);

    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Creating Cube " << itsFilename <<
                       " with shape [xsize:" << nx << " ysize:" << ny <<
                       " npol:" << npol << " nchan:" << nchan <<
                       "], f0: " << f0.getValue("MHz") << " MHz, finc: " <<
                       inc.getValue("kHz") << " kHz");

    itsCube->create(itsFilename, cubeShape, csys);

    // default flux units are Jy/pixel. If we set the restoring beam
    // later on, can set to Jy/beam
    itsCube->setUnits(itsFilename,"Jy/pixel");

    // set the header keywords
    itsCube->setMetadataKeywords(itsFilename,parset.makeSubset("header."));

    // set the image HISTORY keywords
    const std::vector<std::string> historyLines = parset.getStringVector("imageHistory",std::vector<std::string> {});
    itsCube->addHistory(itsFilename,historyLines);

    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiated Cube Builder by creating cube " << itsFilename);
}

template <class T>
CubeBuilder<T>::CubeBuilder(const LOFAR::ParameterSet& parset,
                         const casacore::uInt nchan,
                         const casacore::Quantity& f0,
                         const casacore::Quantity& inc,
                         const std::string& name,
                         const bool uvcoord)
{
    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiating Cube Builder by creating cube");
    itsFilename = makeImageName(parset, name);
    itsCube = accessors::imageAccessFactory(parset);

    setupCube(parset, nchan, f0, inc, uvcoord);
}

template <class T>
CubeBuilder<T>::CubeBuilder(askapparallel::AskapParallel& comm,
                         size_t comm_index,
                         const LOFAR::ParameterSet& parset,
                         const casacore::uInt nchan,
                         const casacore::Quantity& f0,
                         const casacore::Quantity& inc,
                         const std::string& name,
                         const bool uvcoord)
{
    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiating Cube Builder by creating cube");
    itsFilename = makeImageName(parset, name);

    itsCube = accessors::imageAccessFactory(parset, comm, comm_index);

    setupCube(parset, nchan, f0, inc, uvcoord);
}

template <class T>
CubeBuilder<T>::CubeBuilder(const LOFAR::ParameterSet& parset,
                            const casacore::IPosition& shape,
                            const casacore::CoordinateSystem & coordSys,
                            const std::string& name)
{
    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiating Cube Builder by creating cube");
    itsFilename = makeImageName(parset, name);
    itsCube = accessors::imageAccessFactory(parset);

    // Check whether input image is stored at a lower resolution
    if (parset.isDefined("Images.extraoversampling")) {
        itsExtraOversamplingFactor = parset.getFloat("Images.extraoversampling");
        // The parameter should only be defined if has a legitimate value (is set by the code). Check anyway.
        ASKAPDEBUGASSERT(*itsExtraOversamplingFactor > 1.);
        ASKAPLOG_INFO_STR(CubeBuilderLogger, "Using extraoversampling " << *itsExtraOversamplingFactor);
    }

    // Get the image shape
    ASKAPDEBUGASSERT(shape.size()==4);
    casacore::uInt nx = shape[0];
    casacore::uInt ny = shape[1];
    const casacore::uInt npol = shape[2];
    const casacore::uInt nchan = shape[3];
    if (itsExtraOversamplingFactor) {
        const casacore::IPosition fullShape =
            scimath::PaddingUtils::paddedShape(casacore::IPosition(2,nx,ny),*itsExtraOversamplingFactor);
        nx = fullShape[0];
        ny = fullShape[1];
    }
    const casacore::IPosition cubeShape(4, nx, ny, npol, nchan);

    casacore::CoordinateSystem newCoordSys = coordSys;
    casacore::DirectionCoordinate dirCoord = coordSys.directionCoordinate();
    if (itsExtraOversamplingFactor) {
        // fix increments and ref pixel
        casacore::Vector<casacore::Double> inc = dirCoord.increment();
        inc(0) /= *itsExtraOversamplingFactor;
        inc(1) /= *itsExtraOversamplingFactor;
        dirCoord.setIncrement(inc);
        casacore::Vector<casacore::Double> refPix = dirCoord.referencePixel();
        refPix(0) = nx / 2;
        refPix(1) = ny / 2;
        dirCoord.setReferencePixel(refPix);
        newCoordSys.replaceCoordinate(dirCoord,coordSys.findCoordinate(casacore::Coordinate::DIRECTION));
    }

    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Creating Cube " << itsFilename <<
                       " with shape [xsize:" << nx << " ysize:" << ny <<
                       " npol:" << npol << " nchan:" << nchan <<"]");

    itsCube->create(itsFilename, cubeShape, newCoordSys);

    // default flux units are Jy/pixel. If we set the restoring beam
    // later on, can set to Jy/beam
    itsCube->setUnits(itsFilename,"Jy/pixel");

    // set the header keywords
    itsCube->setMetadataKeywords(itsFilename,parset.makeSubset("header."));

    // set the image HISTORY keywords
    const std::vector<std::string> historyLines = parset.getStringVector("imageHistory",std::vector<std::string> {});
    itsCube->addHistory(itsFilename,historyLines);

    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiated Cube Builder by creating cube " << itsFilename);
}

template < class T >
CubeBuilder<T>::~CubeBuilder()
{
}

template <class T>
void CubeBuilder<T>::setupCube(const LOFAR::ParameterSet& parset, 
                               const casacore::uInt nchan,
                               const casacore::Quantity& f0,
                               const casacore::Quantity& inc,
                               const bool uvcoord) {

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

    // Check whether image param is stored at a lower resolution
    if (parset.isDefined("Images.extraoversampling")) {
        itsExtraOversamplingFactor = parset.getFloat("Images.extraoversampling");
        // The parameter should only be defined if has a legitimate value (is set by the code). Check anyway.
        ASKAPDEBUGASSERT(*itsExtraOversamplingFactor > 1.);
        ASKAPLOG_INFO_STR(CubeBuilderLogger, "Using extraoversampling " << *itsExtraOversamplingFactor);
    }

    // Get the image shape
    const vector<casacore::uInt> imageShapeVector = parset.getUintVector("Images.shape");
    casacore::uInt nx = imageShapeVector[0];
    casacore::uInt ny = imageShapeVector[1];
    if (itsExtraOversamplingFactor) {
        const casacore::IPosition fullShape =
            scimath::PaddingUtils::paddedShape(casacore::IPosition(2,nx,ny),*itsExtraOversamplingFactor);
        nx = fullShape[0];
        ny = fullShape[1];
    }
    const casacore::IPosition cubeShape(4, nx, ny, npol, nchan);

    const casacore::CoordinateSystem csys =
    (uvcoord ? createUVCoordinateSystem(parset, nx, ny, f0, inc) :
               createCoordinateSystem(parset, nx, ny, f0, inc));

    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Creating Cube " << itsFilename <<
                       " with shape [xsize:" << nx << " ysize:" << ny <<
                       " npol:" << npol << " nchan:" << nchan <<
                       "], f0: " << f0.getValue("MHz") << " MHz, finc: " <<
                       inc.getValue("kHz") << " kHz");

    itsCube->create(itsFilename, cubeShape, csys);

    // default flux units are Jy/pixel. If we set the restoring beam
    // later on, can set to Jy/beam
    itsCube->setUnits(itsFilename,"Jy/pixel");

    // set the header keywords
    itsCube->setMetadataKeywords(itsFilename,parset.makeSubset("header."));

    // set the image HISTORY keywords
    const std::vector<std::string> historyLines = parset.getStringVector("imageHistory",std::vector<std::string> {});
    itsCube->addHistory(itsFilename,historyLines);


    ASKAPLOG_INFO_STR(CubeBuilderLogger, "Instantiated Cube Builder by creating cube " << itsFilename);
}

template <class T>
std::string CubeBuilder<T>::makeImageName(const LOFAR::ParameterSet& parset, const std:: string& name) {
    vector<string> filenames;
    std::string filename;
    if (parset.isDefined("Images.Names")) {
        filenames = parset.getStringVector("Images.Names", true);
        filename = filenames[0];
    }
    else if(parset.isDefined("Images.name")) {
        filename = parset.getString("Images.name");
    }
    else {
        ASKAPLOG_ERROR_STR(CubeBuilderLogger, "Could not find the image name(s) ");
    }
    ASKAPCHECK(filename.substr(0,5)=="image",
               "Images.name (Names) must start with 'image' starts with " << filename.substr(0,5));

    // If necessary, replace "image" with _name_ (e.g. "psf", "weights")
    // unless name='restored', in which case we append ".restored"
    if (!name.empty()) {
        if (name == "restored") {
            filename = filename + ".restored";
        } else {
            const string orig = "image";
            const size_t f = filename.find(orig);
            // if name has .real or .imag we want to move that to the end
            const string real = ".real";
            const string imag = ".imag";
            const size_t realpos = name.find(real);
            size_t imagpos = name.find(imag);
            // we want to find .imag but not .image
            if (name.length() > imagpos+5) imagpos = string::npos;
            if (realpos == string::npos && imagpos == string::npos) {
                filename.replace(f, orig.length(), name);
            } else if (realpos != string::npos) {
                string basename = name.substr(0,realpos);
                filename.replace(f, orig.length(), basename);
                filename.append(".real");
            } else if (imagpos != string::npos) {
                string basename = name.substr(0,imagpos);
                filename.replace(f, orig.length(), basename);
                filename.append(".imag");
            }
        }
    }
    return filename;
}

template < class T >
void CubeBuilder<T>::writeRigidSlice(const casacore::Array<T>& arr, const casacore::uInt chan)
{
    casacore::IPosition where(4, 0, 0, 0, chan);
    itsCube->write(itsFilename, arr, where);
}

template < class T >
const casacore::Array<float> CubeBuilder<T>::createFlexibleSlice(const casacore::Array<float>& arr)
{
    if (itsExtraOversamplingFactor) {
        // Image param is stored at a lower resolution, so increase to desired resolution before writing
        casacore::Array<float> fullresarr(scimath::PaddingUtils::paddedShape(arr.shape(),*itsExtraOversamplingFactor));
        scimath::PaddingUtils::fftPad(arr,fullresarr);
        return fullresarr;
    }
    else {
        return arr;
    }
}


template < class T >
const casacore::Array<float> CubeBuilder<T>::writeFlexibleSlice(const casacore::Array<float>& arr, const casacore::uInt chan)
{
    const Array<float> farr = createFlexibleSlice(arr);
    casacore::IPosition where(4, 0, 0, 0, chan);
    itsCube->write(itsFilename, farr, where);
    return farr;
}

template < class T >
const casacore::Array<float> CubeBuilder<T>::readRigidSlice(const casacore::uInt chan)
{
    const casacore::IPosition blc(4, 0, 0, 0, chan);
    casacore::IPosition trc = itsCube->shape(itsFilename) - 1;
    ASKAPDEBUGASSERT(chan >=0 && chan <=trc(3));
    trc(3) = chan;
    casacore::Array<float> arr = itsCube->read(itsFilename, blc, trc);
    return arr;
}

template < class T >
casacore::CoordinateSystem
CubeBuilder<T>::createCoordinateSystem(const LOFAR::ParameterSet& parset,
                                    const casacore::uInt nx,
                                    const casacore::uInt ny,
                                    const casacore::Quantity& f0,
                                    const casacore::Quantity& inc)
{
    CoordinateSystem coordsys;
    vector<string> dirVector = parset.getStringVector("Images.direction");
    if (parset.getBool("updatedirection",False)) {
        // override with image specific direction if present - for mosaic case
        vector<string> names = parset.getStringVector("Images.Names",{},false);
        if (names.size()>0) {
            if (parset.isDefined("Images."+names[0]+".direction")) {
                dirVector = parset.getStringVector("Images."+names[0]+".direction");
            }
        }
    }
    const vector<string> cellSizeVector = parset.getStringVector("Images.cellsize");


    // Direction Coordinate
    {
        Matrix<Double> xform(2, 2);
        xform = 0.0;
        xform.diagonal() = 1.0;
        const MDirection dir = asMDirection(dirVector);
        const Quantum<Double> ra = dir.getValue().getLong("deg");
        const Quantum<Double> dec = dir.getValue().getLat("deg");
        ASKAPLOG_DEBUG_STR(CubeBuilderLogger, "Direction: " << ra.getValue() << " degrees, "
                           << dec.getValue() << " degrees");

        Quantum<Double> xcellsize = asQuantity(cellSizeVector.at(0), "arcsec") * -1.;
        Quantum<Double> ycellsize = asQuantity(cellSizeVector.at(1), "arcsec");

        if (itsExtraOversamplingFactor) {
            // have already checked this, but may as well check again
            ASKAPDEBUGASSERT(*itsExtraOversamplingFactor > 1.);
            xcellsize /= *itsExtraOversamplingFactor;
            ycellsize /= *itsExtraOversamplingFactor;
        }
        ASKAPLOG_DEBUG_STR(CubeBuilderLogger, "Cellsize: " << xcellsize.getValue()
                           << " arcsec, " << ycellsize.getValue() << " arcsec");

        casacore::MDirection::Types type;
        casacore::MDirection::getType(type, dirVector.at(2));
        const DirectionCoordinate radec(type, Projection(Projection::SIN), ra, dec, xcellsize, ycellsize, xform,
                                        nx / 2, ny / 2);

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
            ASKAPLOG_INFO_STR(CubeBuilderLogger, "Image cube frequencies will be treated as topocentric");
            freqRef = casacore::MFrequency::TOPO;
        } else if (freqFrame == "lsrk") {
            ASKAPLOG_INFO_STR(CubeBuilderLogger, "Image cube frequencies will be treated as lsrk");
            freqRef = casacore::MFrequency::LSRK;
        } else if (freqFrame == "bary") {
        ASKAPLOG_INFO_STR(CubeBuilderLogger, "Image cube frequencies will be treated as barycentric");
            freqRef = casacore::MFrequency::BARY;
        } else {
            ASKAPTHROW(AskapError, "Unsupported frequency frame "<<freqFrame);
        }


        SpectralCoordinate sc(freqRef, f0, inc, refPix);

        // add rest frequency, but only if requested, and only for
        // image.blah, residual.blah, image.blah.restored
        // Prefer to have them all consistent and wcs compliant. (MHW)
        if (itsRestFrequency.getValue("Hz") > 0.) {
            //if ((itsFilename.find("image.") != string::npos) ||
            //        (itsFilename.find("residual.") != string::npos)) {

                if (!sc.setRestFrequency(itsRestFrequency.getValue("Hz"))) {
                    ASKAPLOG_ERROR_STR(CubeBuilderLogger, "Could not set the rest frequency to " <<
                                       itsRestFrequency.getValue("Hz") << "Hz");
                }
            //}
        }

        coordsys.addCoordinate(sc);
    }

    return coordsys;
}

template <class T>
void CubeBuilder<T>::addBeam(casacore::Vector<casacore::Quantum<double> > &beam)
{
        itsCube->setBeamInfo(itsFilename,beam[0].getValue("rad"),beam[1].getValue("rad"),beam[2].getValue("rad"));
        itsCube->setUnits(itsFilename,"Jy/beam");
}

template <class T>
void CubeBuilder<T>::addBeamList(const BeamList & beamList)
{
    itsCube->setBeamInfo(itsFilename, beamList);
}

template <class T>
void CubeBuilder<T>::setInfo(const casacore::Record & info)
{
    itsCube->setInfo(itsFilename, info);
}

template <class T>
void CubeBuilder<T>::writeImageHistory(const std::vector<std::string>& historyLines)
{
    if ( ! historyLines.empty() ) {
        if ( itsCube ) {
            itsCube->addHistory(this->itsFilename,historyLines);
        }
    }
}

template <class T>
boost::shared_ptr<accessors::IImageAccess<T>>  CubeBuilder<T>::imageHandler()
{
    return itsCube;
}

template <class T>
boost::optional<float> CubeBuilder<T>::oversamplingFactor()
{
    return itsExtraOversamplingFactor;
}
} // namespace cp
} // namespace askap
