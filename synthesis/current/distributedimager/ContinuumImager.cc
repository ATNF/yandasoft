/// @file ContinuumImager.cc
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
/// @author Stephen Ord <stephen.ord@csiro.au>

// Include own header file first
#include "ContinuumImager.h"

// Include package level header file
#include <askap_synthesis.h>

// System includes
#include <string>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <Common/ParameterSet.h>

// Local includes

#include "distributedimager/ContinuumMaster.h"
#include "distributedimager/ContinuumWorker.h"

ASKAP_LOGGER(logger, ".ContinuumImager");

using namespace askap::cp;
using namespace askap;

ContinuumImager::ContinuumImager(LOFAR::ParameterSet& parset,
                                       CubeComms& comms) :
    itsParset(parset), itsComms(comms)
{
    if (isMaster()) {
        ASKAPLOG_INFO_STR(logger,
                          "ASKAP Distributed Continuum Imager - " << ASKAP_PACKAGE_VERSION);
    }
    itsComms.buildCommIndex();
}

ContinuumImager::~ContinuumImager()
{
}

void ContinuumImager::run(void) {

    if (isMaster()) {
        ContinuumMaster master(itsParset,itsComms);
        master.run();
    } else {
        ContinuumWorker worker(itsParset,itsComms);
        worker.run();
    }
    itsComms.barrier(0);
}

bool ContinuumImager::isMaster(void)
{
    return (itsComms.isMaster());
}
