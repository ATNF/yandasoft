/// @file MSGroupInfo.h
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
#ifndef ASKAP_CP_SIMAGER_MSGROUPINFO_H
#define ASKAP_CP_SIMAGER_MSGROUPINFO_H

// System includes
#include <string>
#include <vector>

// ASKAPsoft includes
#include <casacore/casa/Quanta.h>

namespace askap {
namespace cp {

class MSGroupInfo {
    public:
        /// Constructor
        MSGroupInfo();

        /// Constructor
        MSGroupInfo(const std::vector<std::string>& ms);

        /// Destructor
        ~MSGroupInfo();

        casa::uInt getNumChannels(const casa::uInt n) const;

        casa::uInt getTotalNumChannels() const;

        casa::Quantity getFirstFreq() const;

        casa::Quantity getFreqInc() const;

    private:
        std::vector<casa::uInt> itsNumChannels;
        casa::uInt itsTotalNumChannels;
        casa::Quantity itsFirstFreq;
        casa::Quantity itsFreqInc;
};

}
}

#endif
