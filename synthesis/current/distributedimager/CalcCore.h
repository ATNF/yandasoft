/// @file CalcCore.h
///
/// @copyright (c) 2016 CSIRO
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
/// @author Stephen Ord <stephen.ord@csiro.au>
/// based upon SolverCore by
/// Ben Humphreys <ben.humphreys@csiro.au>

#ifndef ASKAP_CP_ASKAP_IMAGER_CALCCORE_H
#define ASKAP_CP_ASKAP_IMAGER_CALCCORE_H

// System includes
#include <string>

// ASKAPsoft includes
#include <Common/ParameterSet.h>
#include <fitting/INormalEquations.h>
#include <fitting/Params.h>
#include <fitting/Solver.h>
#include <fitting/Quality.h>
#include <parallel/ImagerParallel.h>
#include <gridding/IVisGridder.h>
#include <gridding/VisGridderFactory.h>

#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Arrays/Vector.h>
#include <dataaccess/TableDataSource.h>

// Local includes


namespace askap {
namespace cp {

/// @brief Core Normal Equation Calculation functionality required by the imager.
    class CalcCore : public synthesis::ImagerParallel

    {
    public:
        /// @brief Constructor
        CalcCore(LOFAR::ParameterSet& parset,
                   askap::askapparallel::AskapParallel& comms,
                   accessors::TableDataSource ds, int localChannel=1);

        /// @brief Constructor that maintains the gridder
        CalcCore(LOFAR::ParameterSet& parset,
                askap::askapparallel::AskapParallel& comms,
                accessors::TableDataSource ds, askap::synthesis::IVisGridder::ShPtr gdr,
                 int localChannel=1);

        virtual ~CalcCore();

        /// @brief Calc the normal equations
        /// @detail Overrides the virtual function in the ImagerParallel base
        void calcNE();

        void solveNE();

        void doCalc();

        void init();

        void reset();

        void zero();


        void check();

        void restoreImage();

        void writeLocalModel(const std::string& postfix);

        askap::synthesis::IVisGridder::ShPtr gridder() { return itsGridder_p;};

    private:

        // Parameter set
        LOFAR::ParameterSet& itsParset;

        // Communications class
        askap::askapparallel::AskapParallel& itsComms;

        // Solver
        askap::scimath::Solver::ShPtr itsSolver;

        // Restore Solver
        bool itsRestore;

        // Data: vector of the stored datasources
        accessors::TableDataSource itsData;


        // Pointer to the gridder
        askap::synthesis::IVisGridder::ShPtr itsGridder_p;

        // Its channel in the dataset
        int itsChannel;



        // No support for assignment
        CalcCore& operator=(const CalcCore& rhs);

        // No support for copy constructor
        CalcCore(const CalcCore& src);
};
};
};

#endif
