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
#include <askap/askap_synthesis.h>

// System includes
#include <string>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>
#include <Common/ParameterSet.h>

// Local includes

#include "askap/distributedimager/ContinuumMaster.h"
#include "askap/distributedimager/ContinuumWorker.h"

ASKAP_LOGGER(logger, ".ContinuumImager");

using namespace askap::cp;
using namespace askap;

ContinuumImager::ContinuumImager(LOFAR::ParameterSet& parset,
                                       CubeComms& comms, StatReporter& stats) :
    itsParset(parset), itsComms(comms), itsStats(stats)
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
    ASKAPTRACE("ContinuumImager::run");

    if (isMaster()) {
        ContinuumMaster master(itsParset,itsComms,itsStats);
        master.run();
    } else {
        ContinuumWorker worker(itsParset,itsComms,itsStats);
        worker.run();
    }
    ASKAPLOG_INFO_STR(logger,"Reached world barrier");
    itsComms.barrier(0);
    ASKAPLOG_INFO_STR(logger,"Passed world barrier");
}

bool ContinuumImager::isMaster(void)
{
    return (itsComms.isMaster());
}
