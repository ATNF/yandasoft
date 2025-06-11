/// @file ContinuumWorkUnit.cc
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
/// @author Stephen Ord <stephen.ord@csiro.au>

// Include own header file first
#include <askap/messages/ContinuumWorkUnit.h>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>
#include <casacore/casa/Arrays/Array.h>
#include <Blob/BlobArray.h>
#include <Blob/BlobSTL.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobIBufVector.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobOBufVector.h>
#include <askap/askapparallel/AskapParallel.h>

// Using
using namespace askap::cp;

ContinuumWorkUnit::ContinuumWorkUnit()
    : itsGlobalChannel(-1), itsLocalChannel(-1), itsNumContiguousChannels(1u), itsChannelFrequency(0.)
{
}

ContinuumWorkUnit::~ContinuumWorkUnit()
{
}

IMessage::MessageType ContinuumWorkUnit::getMessageType(void) const
{
    return IMessage::SPECTRALLINE_WORKUNIT;
}

/////////////////////////////////////////////////////////////////////
// Setters
/////////////////////////////////////////////////////////////////////
void ContinuumWorkUnit::set_payloadType(PayloadType type)
{
    itsPayloadType = type;
}

void ContinuumWorkUnit::set_dataset(std::string dataset)
{
        itsDataset = dataset;
}

void ContinuumWorkUnit::set_globalChannel(unsigned int chan)
{
    itsGlobalChannel = chan;
}

void ContinuumWorkUnit::set_localChannel(unsigned int chan)
{
    itsLocalChannel = chan;
}

void ContinuumWorkUnit::set_nchan(unsigned int nchan) {
    itsNumContiguousChannels = nchan;
}

void ContinuumWorkUnit::set_channelFrequency(double freq)
{
    itsChannelFrequency = freq;
}

void ContinuumWorkUnit::set_beam(unsigned int beam)
{
    itsBeam = beam;
}
void ContinuumWorkUnit::set_channelWidth(double width) {
    itsChannelWidth = width;
}
void ContinuumWorkUnit::set_writer(unsigned int writer) {
    itsWriter = writer;
}
/////////////////////////////////////////////////////////////////////
// Getters
/////////////////////////////////////////////////////////////////////
ContinuumWorkUnit::PayloadType ContinuumWorkUnit::get_payloadType(void) const
{
    return itsPayloadType;
}

std::string ContinuumWorkUnit::get_dataset(void) const
{
        return itsDataset;
}

unsigned int ContinuumWorkUnit::get_globalChannel(void) const
{
    return itsGlobalChannel;
}

unsigned int ContinuumWorkUnit::get_localChannel(void) const
{
    return itsLocalChannel;
}

double ContinuumWorkUnit::get_channelFrequency(void) const
{
    return itsChannelFrequency;
}
unsigned int ContinuumWorkUnit::get_beam(void) const
{
    return itsBeam;
}
double ContinuumWorkUnit::get_channelWidth(void) const
{
    return itsChannelWidth;
}
unsigned int ContinuumWorkUnit::get_writer(void) const
{
    return itsWriter;
}

unsigned int ContinuumWorkUnit::get_nchan(void) const
{
    return itsNumContiguousChannels;
}

/////////////////////////////////////////////////////////////////////
// Serializers
/////////////////////////////////////////////////////////////////////
void ContinuumWorkUnit::writeToBlob(LOFAR::BlobOStream& os) const
{
    os << static_cast<int>(itsPayloadType);
    os << itsDataset;
    os << itsGlobalChannel;
    os << itsChannelFrequency;
    os << itsChannelWidth;
    os << itsLocalChannel;
    os << itsBeam;
    os << itsWriter;
}

void ContinuumWorkUnit::readFromBlob(LOFAR::BlobIStream& is)
{
    int payloadType;

    is >> payloadType;
    is >> itsDataset;
    is >> itsGlobalChannel;
    is >> itsChannelFrequency;
    is >> itsChannelWidth;
    is >> itsLocalChannel;
    is >> itsBeam;
    is >> itsWriter;

    itsPayloadType = static_cast<PayloadType>(payloadType);
}
/////////////////////////////////////////////////////////////////////
// Communicators
/////////////////////////////////////////////////////////////////////

void ContinuumWorkUnit::sendUnit(int target, askapparallel::AskapParallel& comm)
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
    comm.send(&size,sizeof(long),target,messageType,communicator);

    // Now send the actual byte stream
    comm.send(&buf[0], size * sizeof(int8_t), target, messageType,communicator);


}
void ContinuumWorkUnit::receiveUnit(int& id, askapparallel::AskapParallel& comm) {

    unsigned long size = 0;
    size_t communicator = 0; //default
    id = comm.receiveAnySrc(&size, sizeof(long), this->getMessageType(), communicator);


    // Receive the byte stream
    std::vector<int8_t> buf;
    buf.resize(size);

    ASKAPCHECK(buf.size() == size,
               "ContinuumWorkUnit::receiveUnit() buf is too small");

    comm.receive(&buf[0], size * sizeof(char), id, this->getMessageType(), communicator);

    // Decode
    LOFAR::BlobIBufVector<int8_t> bv(buf);
    LOFAR::BlobIStream in(bv);
    int version = in.getStart("Message");
    ASKAPASSERT(version == 1);

    this->readFromBlob(in);
}
void ContinuumWorkUnit::receiveUnitFrom(const int sender, askapparallel::AskapParallel& comm) {
    unsigned long size = 0;
    size_t communicator = 0; //default
    comm.receive(&size, sizeof(long), sender, this->getMessageType(),communicator);
    // Receive the byte stream
    std::vector<int8_t> buf;
    buf.resize(size);

    ASKAPCHECK(buf.size() == size,
               "ContinuumWorkUnit::receiveUnit() buf is too small");

    comm.receive(&buf[0], size * sizeof(char), sender, this->getMessageType(), communicator);

    // Decode
    LOFAR::BlobIBufVector<int8_t> bv(buf);
    LOFAR::BlobIStream in(bv);
    int version = in.getStart("Message");
    ASKAPASSERT(version == 1);

    this->readFromBlob(in);
}
