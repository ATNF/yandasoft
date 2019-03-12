/// @file lsqrtest.cc
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
/// @author Vitaliy Ogarko <vogarko@gmail.com>
///
/// This is an application for testing parallel functionality of the lsqr solver.
/// It should be merged with unit tests in LSQRSolverTest.h once we have the parallel unit tests support.
///

#include <iostream>
#include <cassert>

#include <lsqr_solver/LSQRSolver.h>
#include <lsqr_solver/SparseMatrix.h>

#include <askapparallel/AskapParallel.h>

using namespace askap;

//---------------------------------------------------------------------
// Testing overdetermined system.
//---------------------------------------------------------------------
// Consider a regression with constant, linear and quadratic terms:
// f(x) = b1 + b2 * x + b3 * x^2
//
// We build a matrix 3 x N:
//
//   1  x1  x1^2
//   1  x2  x2^2
//   ...
//   1  xn  xn^2,
//
// with x_i = i/N, i = 1, ..., N.
//
// Apply LSQR method and compare to the solution x = (b1, b2, b3).
// An example from [1].
// [1] Least Squares Estimation, Sara A. van de Geer, Volume 2, pp. 1041-1045,
//     in Encyclopedia of Statistics in Behavioral Science, 2005.
//---------------------------------------------------------------------
void test_overdetermined(int myrank, int nbproc)
{
    if (nbproc != 1 && nbproc != 3)
    {
        std::cout << "WARNING: Test test_overdetermined can only be run on 1 and 3 CPUs!" << std::endl;
        return;
    }

    size_t nelements_total = 3;
    size_t nrows = 1000;
    double rmin = 1.e-14;
    size_t niter = 100;

    size_t nelements = nelements_total / nbproc;

    lsqr::SparseMatrix matrix(nrows, nelements * nrows);

    lsqr::Vector b_RHS(nrows, 0.0);

    lsqr::Vector b(3);
    b[0] = 1.0;
    b[1] = - 3.0;
    b[2] = 0.0;

    // Building the matrix with right hand side.
    for (size_t i = 0; i < nrows; ++i)
    {
        matrix.NewRow();

        double xi = double(i) / double(nrows);

        if (nbproc == 1)
        {
            matrix.Add(1.0, 0);
            matrix.Add(xi, 1);
            matrix.Add(xi * xi, 2);
        }
        else if (nbproc == 3)
        {
            matrix.Add(pow(xi, double(myrank)), 0);
        }

        b_RHS[i] = b[0] + b[1] * xi + b[2] * xi * xi;
    }
    matrix.Finalize(nelements);

    lsqr::LSQRSolver solver(nrows, nelements);

    lsqr::Vector x(nelements, 0.0);
    solver.Solve(niter, rmin, matrix, b_RHS, x, myrank, nbproc);

    double epsilon = 1.e-14;

    // Check the solution.
    if (nbproc == 1)
    {
        for (int i = 0; i < 3; i++)
        {
            assert(std::abs(b[i] - x[i]) < epsilon);
        }
    }
    else if (nbproc == 3)
    {
        assert(std::abs(b[myrank] - x[0]) < epsilon);
    }
}

// This application should be run using 1 and 3 CPUs only.
int main(int argc, char *argv[])
{
    askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

    int nbproc = comms.nProcs();
    int myrank = comms.rank();

    test_overdetermined(myrank, nbproc);
}

