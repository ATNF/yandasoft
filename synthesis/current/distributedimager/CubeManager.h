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
/// @author Stephen Ord <Stephen.Ord@csiro.au>
///

#ifndef ASKAP_CP_IMAGER_CUBEMANAGER_H
#define ASKAP_CP_IMAGER_CUBEMANAGER_H


#include <map.h>

namespace askap {
namespace cp {

    /// a utility class that contains all the channels a worker is responsible for and its rank
    
    class cubeClient {
    
    public:
    /// what is the status of this channel
        int getChannelStatus(int channel) const {return channelMap[channel];};
        int myRank() { return Rank;};
    private:
    /// a map of the processing status of each channel 0 for not processed and 1 for processed
        std::map<int,int> channelMap;
    /// rank of the client
        int Rank;
    };
    
    /// a utility class that holds all the clients for a given writer
    
    class cubeWriter {
        
    public:
        
        std::vector<cubeClient>& myClients() const { return this->clientList;};
        int myRank() {return this->Rank;} ;
        
    private:
        
        std::vector<cubeClient> clientList;
        int Rank;
    };
    
    /// a utility class that holds all the writers for a given instance
    
    class cubeManager {

    
    public:
        
        std::vector<vector<cubeWriters> >& getWriters() const {return theWriters;};
        
    private:
        
        std::vector<vector<cubeWriter> > theWriters;
    
};
}
}
#endif
