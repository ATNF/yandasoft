/// @file
/// 
/// @brief interface to eigen problem solver
/// @details
///
/// An interface module for low-level routines to solve eigenproblem or
/// a singular value decomposition problem
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

#include <calweightsolver/EigenSolve.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <fitting/GSLSVDReplacement.h>

#include <askap/AskapError.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

using namespace casa;
using namespace std;

using namespace askap;

// return results
const casa::Matrix<casa::Complex>& EigenSolver::getEigenVectors() const throw()
{
 return vec;
}

const casa::Vector<casa::Double>& EigenSolver::getEigenValues() const throw()
{
 return val;
}

const casa::Matrix<casa::Complex>& EigenSolver::getV() const throw()
{
 return vecV;
}

// solve eigenproblem and fill vec and val data members
void EigenSolver::solveEigen(const casa::Matrix<casa::Complex> &in)
                               throw(casa::AipsError)
{
  // temporary code which works for real matrices only
  // extra copying is happening here
  ASKAPASSERT(in.nrow() == in.ncolumn());
  casa::uInt size = in.nrow();
  gsl_matrix *V = gsl_matrix_alloc(size,size);
  ASKAPDEBUGASSERT(V!=NULL);
  
  gsl_vector *S = gsl_vector_alloc(size);
  gsl_matrix *A = gsl_matrix_alloc(size,size);
  
  for (size_t row=0; row<in.nrow();++row) {
       for (size_t column=0; column<in.ncolumn();++column) {
            gsl_matrix_set(A,row,column,real(in(row,column)));
       }
  }
  try {
    scimath::SVDecomp(A,V,S);
  
    vec.resize(in.nrow(),in.ncolumn());
    vecV.resize(in.nrow(),in.ncolumn());
    val.resize(size);
     for (size_t row=0; row<in.nrow();++row) {
          for (size_t column=0; column<in.ncolumn();++column) {
               vec(row,column) = gsl_matrix_get(A,row,column);
               vecV(row,column) = gsl_matrix_get(V,row,column);
          }
          val[row] = gsl_vector_get(S,row);
     }
  }
  catch (...) {
    gsl_vector_free(S);
    gsl_matrix_free(V);
    gsl_matrix_free(A);
     throw;
  }
  
  gsl_vector_free(S);
  gsl_matrix_free(V);
  gsl_matrix_free(A);
}
