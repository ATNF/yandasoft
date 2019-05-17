/// @file IMessage.h
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

#ifndef ASKAP_CP_SIMAGER_IMESSAGE_H
#define ASKAP_CP_SIMAGER_IMESSAGE_H

// Boost includes
#include <boost/shared_ptr.hpp>

// ASKAPsoft includes
#include <fitting/ISerializable.h>

namespace askap {
    namespace cp {

        class IMessage : public ISerializable
        {
            public:
                enum MessageType {
                    SPECTRALLINE_WORKUNIT,
                    SPECTRALLINE_WORKREQUEST,
                    SPECTRALLINE_WORKRESULT
                };

                /// @brief Destructor.
                virtual ~IMessage();

                /// @brief Messages must be self-identifying and must
                /// their type via this interface.
                ///
                /// @note Messages must be self-identifying and must return
                /// their type via this interface. While they can also be
                /// identified by their class type, this method easily translates
                /// to an int which can be used to tag messags (eg. MPI tags).
                virtual MessageType getMessageType(void) const = 0;

        };

        /// Short cut for shared pointer to an IMessage
        typedef boost::shared_ptr<IMessage> IMessageSharedPtr;
    };
};

#endif
