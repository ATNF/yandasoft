/// @file SpectralLineWorkRequest.cc
///
/// @copyright (c) 2009 CSIRO
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
#include <messages/SpectralLineWorkRequest.h>

// System includes
#include <limits>

// ASKAPsoft includes
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>
#include <casacore/casa/Arrays/Array.h>
#include <Blob/BlobArray.h>
#include <Blob/BlobSTL.h>
#include <fitting/Params.h>

// Using
using namespace askap::cp;

const unsigned int SpectralLineWorkRequest::CHANNEL_UNINITIALISED
    = std::numeric_limits<unsigned int>::max();

SpectralLineWorkRequest::SpectralLineWorkRequest()
    : itsGlobalChannel(CHANNEL_UNINITIALISED)
{
}

SpectralLineWorkRequest::~SpectralLineWorkRequest()
{
}

IMessage::MessageType SpectralLineWorkRequest::getMessageType(void) const
{
    return IMessage::SPECTRALLINE_WORKREQUEST;
}

/////////////////////////////////////////////////////////////////////
// Setters
/////////////////////////////////////////////////////////////////////
void SpectralLineWorkRequest::set_globalChannel(unsigned int chan)
{
    itsGlobalChannel = chan;
}

void SpectralLineWorkRequest::set_params(askap::scimath::Params::ShPtr params)
{
    itsParams = params;
}

/////////////////////////////////////////////////////////////////////
// Getters
/////////////////////////////////////////////////////////////////////
unsigned int SpectralLineWorkRequest::get_globalChannel(void) const
{
    return itsGlobalChannel;
}

askap::scimath::Params::ShPtr SpectralLineWorkRequest::get_params(void)
{
    return itsParams;
}

/////////////////////////////////////////////////////////////////////
// Serializers
/////////////////////////////////////////////////////////////////////
void SpectralLineWorkRequest::writeToBlob(LOFAR::BlobOStream& os) const
{
    os << itsGlobalChannel;
    const bool hasParams = itsParams.get() == 0 ? false : true;
    os << hasParams;
    if (hasParams) {
        os << *itsParams;
    }
}

void SpectralLineWorkRequest::readFromBlob(LOFAR::BlobIStream& is)
{
    is >> itsGlobalChannel;
    bool hasParams;
    is >> hasParams;
    if (hasParams) {
        itsParams.reset(new askap::scimath::Params());
        is >> *itsParams;
    } else {
        itsParams.reset();
    }
}
