/// @file MEParallel.cc
///
/// @brief Base class for parallel synthesis applications
/// @details
/// Supports parallel algorithms by providing methods for initialization
/// of MPI connections, sending normal equations and models around.
/// There is assumed to be one solver and many prediffers.
///
/// @copyright (c) 2007 CSIRO
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///

// Include own header file first
#include <parallel/MEParallel.h>

// Package level header file
#include <askap_synthesis.h>

// System includes
#include <cmath>

// Askapsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askapparallel/AskapParallel.h>
#include <askapparallel/BlobIBufMW.h>
#include <askapparallel/BlobOBufMW.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobOStream.h>
#include <Common/ParameterSet.h>
#include <fitting/Equation.h>
#include <fitting/Solver.h>
#include <fitting/INormalEquations.h>
#include <fitting/ImagingNormalEquations.h>
#include <fitting/GenericNormalEquations.h>
#include <profile/AskapProfiler.h>
#include <casa/OS/Timer.h>

ASKAP_LOGGER(logger, ".parallel");

using namespace askap;
using namespace askap::askapparallel;
using namespace askap::scimath;

namespace askap {
namespace synthesis {

MEParallel::MEParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset) :
        SynParallel(comms, parset)
{
    itsSolver = Solver::ShPtr(new Solver);
    itsNe = ImagingNormalEquations::ShPtr(new ImagingNormalEquations(*itsModel));
}

MEParallel::~MEParallel()
{
}

// Send the normal equations from this worker to the master as a blob
void MEParallel::sendNE()
{
    ASKAPTRACE("MEParallel::sendNE");

    if (itsComms.isParallel() && itsComms.isWorker()) {
        casa::Timer timer;
        timer.mark();
        ASKAPLOG_DEBUG_STR(logger, "Reducing normal equations to the solver");
        reduceNE(itsNe);
        ASKAPLOG_DEBUG_STR(logger, "Reduced normal equations to the solver in "
                              << timer.real() << " seconds ");
    }
}

// Receive the normal equations as a blob
void MEParallel::receiveNE()
{
    ASKAPTRACE("MEParallel::receiveNE");

    ASKAPCHECK(itsSolver, "Solver not yet defined");

    if (itsComms.isParallel() && itsComms.isMaster()) {
        ASKAPLOG_INFO_STR(logger, "Initialising solver");
        itsSolver->init();

        ASKAPLOG_INFO_STR(logger, "Waiting for normal equation reduction to complete");
        casa::Timer timer;
        timer.mark();

        reduceNE(itsNe);
        itsSolver->addNormalEquations(*itsNe);

        ASKAPLOG_INFO_STR(logger, "Received normal equations from all prediffers in "
                              << timer.real() << " seconds");
    }
}

/*
 * This method performs a graph reduction (using a binary tree topology)
 * from all processes to rank zero. The sequence of workers is mapped to
 * a binary tree like so:
 * - Rank 0 is the root
 * - To find the parent of a node: floor((rank - 1) / 2)
 * - To find the left child of a node: (2 * rank) + 1
 * - To find the right child of a node: (2 * rank) + 2
 *
 * The height of the tree (and hence the number of parallel reductions is
 * given by floor(log2(nProcs)) where nProcs is the number of processes
 * participating in the reduction.
 *
 * Below is an example of the tree for a 7 process reduction. This tree has a
 * height of 2.
 *              0
 *           /     \
 *          1       2
 *        /  \     /  \
 *       3    4   5    6
 *
 * In the above case the result is a perfect binary tree (all leaves are at the
 * same depth) however this method will also handle the imperfect case.
 */ 
void MEParallel::reduceNE(askap::scimath::INormalEquations::ShPtr ne)
{
    // Number of processes in the reduction
    const int nProcs = itsComms.nProcs();

    // Rank (zero-based)
    const int rank = itsComms.rank();

    // This is the height of the binary tree. For example:
    // floor(log2(1)) == 0
    // floor(log2(3)) == 1
    // floor(log2(4)) == 2
    const int height = int(floor(log2(nProcs)));

    // The depth of a node n is the length of the path from the root to the node.
    // The depth of the root node is zero. For example:
    // floor(log2(0 + 1, 2)) == 0
    // floor(log2(1 + 1, 2)) == 1
    // floor(log2(2 + 1, 2)) == 1
    // floor(log2(3 + 1, 2)) == 2
    const int depth = int(floor(log2(rank + 1)));
    ASKAPCHECK(depth <= height, "Depth exceeds height");

    // One reduction step is executed for each level of the tree
    // (except the top level)
    for (int level = height; level > 0; --level) {

        // Only two levels participate in each step, upper level are receivers
        // and lower level are senders
        if (depth == level) {
            // This round I am a sender
            const int parent = int(floor((rank - 1) / 2));
            sendNormalEquations(ne, parent);

        } else if (depth == level - 1) {
            // This round I am a receiver

            // Receive from the left child if it exists
            const int left = (2 * rank) + 1;
            if (left < nProcs) {
                ne->merge(*receiveNormalEquations(left));
            }

            // Receive from the right child if it exists
            const int right = (2 * rank) + 2;
            if (right < nProcs) {
                ne->merge(*receiveNormalEquations(right));
            }
        } else {
            // This round I am a non-participant
        }
    }
}

void MEParallel::sendNormalEquations(const askap::scimath::INormalEquations::ShPtr ne, int dest)
{
    ASKAPDEBUGTRACE("MEParallel::sendNormalEquations");

    casa::Timer timer;
    timer.mark();
    ASKAPLOG_DEBUG_STR(logger, "Sending normal equations to rank " << dest);

    BlobOBufMW bobmw(itsComms, dest);
    LOFAR::BlobOStream out(bobmw);
    out.putStart("ne", 1);
    out << itsComms.rank() << *ne;
    out.putEnd();
    bobmw.flush();
    ASKAPLOG_DEBUG_STR(logger, "Sent normal equations to rank " << dest << " in "
            << timer.real() << " seconds ");
}

askap::scimath::INormalEquations::ShPtr MEParallel::receiveNormalEquations(int source)
{
    ASKAPDEBUGTRACE("MEParallel::receiveNormalEquations");

    askap::scimath::INormalEquations::ShPtr ne;

    // This code needs to create an empty/pristine normal equations
    // instance. To do this the type is needed, information which is not directly
    // available. The clone/reset method below functionally works but is prone
    // to significant memory bloat, especially in the case of huge ImagingNormalEquations.
    // So (a bit of a hack) this code deals with those types is knows, and just clones
    // in the case it doesn't.
    if (dynamic_cast<ImagingNormalEquations*>(itsNe.get())) {
        ne = ImagingNormalEquations::ShPtr(new ImagingNormalEquations());
    } else if (dynamic_cast<GenericNormalEquations*>(itsNe.get())) {
        ne = GenericNormalEquations::ShPtr(new GenericNormalEquations());
    } else {
        ne = itsNe->clone(); 
        ne->reset(); // Reset the normal equation, not the pointer!
    }

    ASKAPLOG_DEBUG_STR(logger, "Waiting to receive normal equations from rank " << source);
    casa::Timer timer;
    timer.mark();

    BlobIBufMW bibmw(itsComms, source);
    LOFAR::BlobIStream in(bibmw);
    const int version = in.getStart("ne");
    ASKAPASSERT(version == 1);
    int rank;
    in >> rank >> *ne;
    in.getEnd();
    ASKAPCHECK(rank == source, "Received normal equations are from an unexpected source");
    ASKAPLOG_DEBUG_STR(logger, "Received normal equations from rank " << source
            << " after " << timer.real() << " seconds");

    return ne;
}

void MEParallel::writeModel(const std::string &)
{
}

}
}
