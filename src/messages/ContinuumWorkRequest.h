/// @file ContinuumWorkRequest.h
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

#ifndef ASKAP_CP_SIMAGER_CONTINUUMWORKREQUEST_H
#define ASKAP_CP_SIMAGER_CONTINUUMWORKREQUEST_H

// ASKAPsoft includes
#include <messages/IMessage.h>
#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>
#include <fitting/Params.h>
#include <askapparallel/AskapParallel.h>

namespace askap {
    namespace cp {

        class ContinuumWorkRequest : public IMessage
        {
            public:
                // Default channel number (max value of "unsigned int") used to
                // indicate the channel number has not been set.
                static const unsigned int CHANNEL_UNINITIALISED;

                /// @brief Constructor.
                ContinuumWorkRequest();

                /// @brief Destructor.
                virtual ~ContinuumWorkRequest();

                /// @brief Messages must be self-identifying and must
                /// their type via this interface.
                ///
                /// @note Messages must be self-identifying and must return
                /// their type via this interface. While they can also be
                /// identified by their class type, this method easily translates
                /// to an int which can be used to tag messags (eg. MPI tags).
                virtual MessageType getMessageType(void) const;

                // Setters
                void set_globalChannel(unsigned int chan);
                void set_params(askap::scimath::Params::ShPtr params);

                // Getters
                unsigned int get_globalChannel(void) const;
                askap::scimath::Params::ShPtr get_params(void);

                // Serializer functions

                /// @brief write the object to a blob stream
                /// @param[in] os the output stream
                virtual void writeToBlob(LOFAR::BlobOStream& os) const;

                /// @brief read the object from a blob stream
                /// @param[in] is the input stream
                virtual void readFromBlob(LOFAR::BlobIStream& is);

                /// @brief Send this request to the master
                /// @param[in] master is the id of the node to which the request is sent
                /// @param[in] comm is the communicator to be used
                void sendRequest(int master, askap::askapparallel::AskapParallel& comm);

                /// @brief Receive this request from anyone
                /// @param[out] id becomes the id of the node sending the request
                /// @param[in] comm is the communicator to be used
                void receiveRequest(int& id, askap::askapparallel::AskapParallel& comm);


            private:

                unsigned int itsGlobalChannel;
                askap::scimath::Params::ShPtr itsParams;


            };

    };
};

#endif
