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

namespace askap {
namespace cp {

template <class T>
class CubeBuilder {
    public:
        /// Constructor
        CubeBuilder(const LOFAR::ParameterSet& parset,
                    const casacore::uInt nchan,
                    const casacore::Quantity& f0,
                    const casacore::Quantity& inc,
                    const std::string& name = "");

        CubeBuilder(const LOFAR::ParameterSet& parset,const std:: string& name);

        /// Destructor
        ~CubeBuilder();

        void writeRigidSlice(const casacore::Array<T>& arr, const casacore::uInt chan);
        void writeFlexibleSlice(const casacore::Array<float>& arr, const casacore::uInt chan);

        casacore::CoordinateSystem
        createCoordinateSystem(const LOFAR::ParameterSet& parset,
                               const casacore::uInt nx,
                               const casacore::uInt ny,
                               const casacore::Quantity& f0,
                               const casacore::Quantity& inc);

        void addBeam(casacore::Vector<casacore::Quantum<double> > &beam);
        void setUnits(const std::string &units);
        void setDateObs(const casacore::MVEpoch &dateObs);

        std::string filename() const{return itsFilename;};

    private:


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
