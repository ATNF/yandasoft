/// @file CubeBuilder.h
///
/// Class to run the creation of a new cube
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
/// @author Ben Humphreys <Ben.Humphreys@csiro.au>
///
#ifndef ASKAP_CP_SIMAGER_CUBEBUILDER_H
#define ASKAP_CP_SIMAGER_CUBEBUILDER_H

// System includes
#include <string>

// ASKAPsoft includes
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <Common/ParameterSet.h>
#include <askap/imageaccess/ImageAccessFactory.h>

#include <casacore/images/Images/PagedImage.h>
#include <casacore/lattices/Lattices/PagedArray.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/casa/Quanta.h>
#include <askap/distributedimager/CubeComms.h>

namespace askap {
namespace cp {

template <class T>
class CubeBuilder {
    public:
        /// @brief Construct from parset and some extra parameters
        /// @details Construct CubeBuilder from parset, frequency parameters and name.
        /// @param[in] parset - ParameterSet, used to get the image parameters
        /// @param[in] nchan - uInt Number of channels in the cube
        /// @param[in] f0 - Quantity Start frequency of the cube
        /// @param[in] inc - Quantity Frequency increments of the cube
        /// @param[in] name - string The type of image (e.g., "restored","psf","visgrid")
        /// @param[in] uvcoord - bool label the cube with UV coordinates, this is for
        /// complex casa type cubes with gridded visibilities
        CubeBuilder(const LOFAR::ParameterSet& parset,
                    const casacore::uInt nchan,
                    const casacore::Quantity& f0,
                    const casacore::Quantity& inc,
                    const std::string& name = "",
                    const bool uvcoord = false);

        /// @brief Construct from parset and some extra parameters
        /// @details Construct CubeBuilder from parset, shape, coordinate system and name.
        /// @param[in] parset - ParameterSet, used to get the image parameters
        /// @param[in] shape - IPosition giving the image shape
        /// @param[in] coordsys - A casacore coordinate system for the cube
        /// @param[in] name - string The type of image (e.g., "restored","psf","visgrid")
        CubeBuilder(const LOFAR::ParameterSet& parset,
                    const casacore::IPosition & shape,
                    const casacore::CoordinateSystem & coordSys,
                    const std::string& name = "");

        CubeBuilder(askapparallel::AskapParallel& comm,
                    size_t comm_index,
                    const LOFAR::ParameterSet& parset,
                    const casacore::uInt nchan,
                    const casacore::Quantity& f0,
                    const casacore::Quantity& inc,
                    const std::string& name = "",
                    const bool uvcoord = false);

        /// @brief Construct from parset and name
        /// @details Construct CubeBuilder from parset and name. This version is for
        /// the ranks that don't actually create the cube, but just use it to write their output to
        /// @param[in] parset - ParameterSet, used to get the image parameters
        /// @param[in] name - string The type of image (e.g., "restored","psf","visgrid")
        CubeBuilder(const LOFAR::ParameterSet& parset,const std:: string& name);

        /// Destructor
        ~CubeBuilder();

        /// @brief Create the image filename
        /// @details This create the filename for the image using the base name from the Parset and
        /// the type of image specified by the name parameter
        /// @param[in] parset - ParameterSet, used to get the base image name
        /// @param[in] name - The type of image (e.g., "restored","psf","visgrid")
        static std::string makeImageName(const LOFAR::ParameterSet& parset, const std:: string& name);

        /// @brief Write a channel to the cube
        /// @details This writes a single channel (slice) to the cube without changing
        /// the shape of the array given.
        /// @param[in] arr - Array channel image data to write to the cube
        /// @param[in] chan - uInt channel in the cube to write to
        void writeRigidSlice(const casacore::Array<T>& arr, const casacore::uInt chan);

        /// @brief Oversample the array if needed
        /// @details This will oversample the input array if needed and return it, 
        /// @params[in] arr - Array channel image data
        const casacore::Array<float> createFlexibleSlice(const casacore::Array<float>& arr);

        /// @brief Write a channel to the cube, oversampled if needed
        /// @details This writes a single channel (slice) to the cube, oversampling the
        /// array given if needed (i.e., Nyquist gridding was used to reduce the array size).
        /// The array written out is also returned, so it is available for statistics calculation.
        /// @param[in] arr - Array channel image data to write to the cube
        /// @param[in] chan - uInt channel in the cube to write to
        /// @return a reference copy of the possibly oversampled array written to the cube
        const casacore::Array<float> writeFlexibleSlice(const casacore::Array<float>& arr, const casacore::uInt chan);

        /// @brief Read a channel from the cube
        /// @details This reads a single channel (slice) to the cube
        /// @param[in] chan - uInt channel in the cube to read
        /// @return The array read from the cube
        const casacore::Array<float> readRigidSlice(const casacore::uInt chan);

        /// @brief create a coordinate system
        /// @details This creates a coordinate system for a cube
        /// @param[in] parset - ParameterSet, used to get the image parameters
        /// @param[in] nx - uInt grid size in x direction
        /// @param[in] ny - uInt grid size in y direction
        /// @param[in] f0 - Quantity Start frequency of the cube
        /// @param[in] inc - Quantity Frequency increments of the cube
        casacore::CoordinateSystem
        createCoordinateSystem(const LOFAR::ParameterSet& parset,
                               const casacore::uInt nx,
                               const casacore::uInt ny,
                               const casacore::Quantity& f0,
                               const casacore::Quantity& inc);

        /// @brief create a coordinate system for a UV grid
        /// @details This creates a coordinate system for a cube containing gridded uv data
        /// @param[in] parset - ParameterSet, used to get the image parameters
        /// @param[in] nx - uInt grid size in x direction
        /// @param[in] ny - uInt grid size in y direction
        /// @param[in] f0 - Quantity Start frequency of the cube
        /// @param[in] inc - Quantity Frequency increments of the cube
        casacore::CoordinateSystem
        createUVCoordinateSystem(const LOFAR::ParameterSet& parset,
                                 const casacore::uInt nx,
                                 const casacore::uInt ny,
                                 const casacore::Quantity& f0,
                                 const casacore::Quantity& inc);

        /// @brief set the restoring beam size for this image
        /// @details This sets the 'overall' beam size for the cube
        /// @param[in] beam - Vector<Quantum<Double>>, the beam parameters (major, minor, PA)
        void addBeam(casacore::Vector<casacore::Quantum<double> > &beam);

        /// @brief set units for this image
        /// @details This sets the brightness units for the cube (e.g., "Jy/beam")
        /// @param[in] units - string, the units
        void setUnits(const std::string &units)
        { itsCube->setUnits(itsFilename, units); }

        /// @brief set observing date for this image
        /// @details This sets the observing date for the cube, usually the starting date
        /// @param[in] dateObs - MVEpoch, the observation date
        void setDateObs(const casacore::MVEpoch &dateObs);

        /// @brief write the image history
        /// @details This sets the image history for the cube, a number of lines
        /// describing processing history and other details
        /// @param[in] historyLines - vector<string> the history to be recorded with the image
        void writeImageHistory(const std::vector<std::string>& historyLines);

        /// @brief set the restoring beam size for each channel in this image
        /// @details This sets the beam sizes for the cube, one set for each channel
        /// @param[in] beamList - BeamList, the list of all beam parameters
        void addBeamList(const accessors::BeamList & beamList);

        /// @brief set the image info from a record
        /// @details This sets the casacore image Info
        /// @param[in] info - Record, extra image information to be recorded
        void setInfo(const casacore::Record & info);

        /// @brief the filename of this cube
        /// @details Access to the name of the cube on disk
        /// @return string, filename of the cube
        std::string filename() const{return itsFilename;};

        /// @brief the image accessor for this cube
        /// @details Generic access to the cube using standard interface
        /// @return shared pointer to image accessor
        boost::shared_ptr<accessors::IImageAccess<T>> imageHandler();

        /// @brief the oversampling factor for this cube
        /// @details If Nyquist gridding is used this will contain the oversampling factor,
        /// it may be needed to oversample some of the arrays before writing them to the cube
        /// @return optional<float> the oversampling factor
        boost::optional<float> oversamplingFactor();


    private:

        void setupCube(const LOFAR::ParameterSet& parset,
                       const casacore::uInt nchan,
                       const casacore::Quantity& f0,
                       const casacore::Quantity& inc,
                       const bool uvcoord);

        /// Shared pointer to the Image accessor
        boost::shared_ptr<accessors::IImageAccess<T> > itsCube;

        /// Image name from parset - must start with "image."
        std::string itsFilename;

        /// Rest frequency to be written to the cubes
        casacore::Quantum<double> itsRestFrequency;

        /// Description of the polarisation properties of the output cubes
        casacore::Vector<casacore::Stokes::StokesTypes> itsStokes;

        /// @brief extra oversampling factor to use when building cubes
        boost::optional<float> itsExtraOversamplingFactor;
};

}
}
#include "CubeBuilder.tcc"

#endif
