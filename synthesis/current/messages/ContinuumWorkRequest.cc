/// @file ContinuumWorkRequest.cc
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
#include <messages/ContinuumWorkRequest.h>

// System includes
#include <limits>

// ASKAPsoft includes
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>
#include <casacore/casa/Arrays/Array.h>
#include <Blob/BlobArray.h>
#include <Blob/BlobSTL.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobIBufVector.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobOBufVector.h>

#include <fitting/Params.h>

// Using
using namespace askap::cp;

const unsigned int ContinuumWorkRequest::CHANNEL_UNINITIALISED
    = std::numeric_limits<unsigned int>::max();

ContinuumWorkRequest::ContinuumWorkRequest()
    : itsGlobalChannel(CHANNEL_UNINITIALISED)
{
}

ContinuumWorkRequest::~ContinuumWorkRequest()
{
}

IMessage::MessageType ContinuumWorkRequest::getMessageType(void) const
{
    return IMessage::SPECTRALLINE_WORKREQUEST;
}

/////////////////////////////////////////////////////////////////////
// Setters
/////////////////////////////////////////////////////////////////////
void ContinuumWorkRequest::set_globalChannel(unsigned int chan)
{
    itsGlobalChannel = chan;
}

void ContinuumWorkRequest::set_params(askap::scimath::Params::ShPtr params)
{
    itsParams = params;
}

/////////////////////////////////////////////////////////////////////
// Getters
/////////////////////////////////////////////////////////////////////
unsigned int ContinuumWorkRequest::get_globalChannel(void) const
{
    return itsGlobalChannel;
}

askap::scimath::Params::ShPtr ContinuumWorkRequest::get_params(void)
{
    return itsParams;
}

/////////////////////////////////////////////////////////////////////
// Serializers
/////////////////////////////////////////////////////////////////////
void ContinuumWorkRequest::writeToBlob(LOFAR::BlobOStream& os) const
{
    os << itsGlobalChannel;
    const bool hasParams = itsParams.get() == 0 ? false : true;
    os << hasParams;
    if (hasParams) {
        os << *itsParams;
    }
}

void ContinuumWorkRequest::readFromBlob(LOFAR::BlobIStream& is)
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
/////////////////////////////////////////////////////////////////////
// Communicators
/////////////////////////////////////////////////////////////////////

void ContinuumWorkRequest::sendRequest(int master, askapparallel::AskapParallel& comm)
{

    size_t communicator = 0; //default
    std::vector<int8_t> buf;
    LOFAR::BlobOBufVector<int8_t> bv(buf);
    LOFAR::BlobOStream out(bv);
    out.putStart("Message", 1);
    this->writeToBlob(out);
    out.putEnd();

    int messageType = this->getMessageType();

    // First send the size of the buffer
    const unsigned long size = buf.size();
    comm.send(&size,sizeof(long),master,messageType,communicator);

    // Now send the actual byte stream
    comm.send(&buf[0], size * sizeof(int8_t), master, messageType,communicator);


}
void ContinuumWorkRequest::receiveRequest(int& id, askapparallel::AskapParallel& comm) {

    unsigned long size = 0;
    size_t communicator = 0; //default
    id = comm.receiveAnySrc(&size, sizeof(long), this->getMessageType(), communicator);


 // Receive the byte stream
    std::vector<int8_t> buf;
    buf.resize(size);

    ASKAPCHECK(buf.size() == size,
               "ContinuumWorkRequest::receiveRequest() buf is too small");

    comm.receive(&buf[0], size * sizeof(char), id, this->getMessageType(), communicator);

    // Decode
    LOFAR::BlobIBufVector<int8_t> bv(buf);
    LOFAR::BlobIStream in(bv);
    int version = in.getStart("Message");
    ASKAPASSERT(version == 1);

    this->readFromBlob(in);
}
