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

// Include own header file first
#include "GridKernel.h"

/// Use pointers instead of casacore::Matrix operators to grid
//#define ASKAP_GRID_WITH_POINTERS 1

/// Use BLAS
//#define ASKAP_GRID_WITH_BLAS 1

#ifdef ASKAP_GRID_WITH_BLAS
#ifdef __APPLE_CC__
#include <vecLib/cblas.h>
#else
#include <cblas.h>
#endif
#endif

namespace askap {
namespace synthesis {

std::string GridKernel::info() {
#ifdef ASKAP_GRID_WITH_BLAS
	return std::string("Gridding with BLAS");
#else
#ifdef ASKAP_GRID_WITH_POINTERS
	return std::string("Gridding with casacore::Matrix pointers");
#else
	return std::string("Standard gridding/degridding with casacore::Matrix");
#endif
#endif
}

/// Totally selfcontained gridding
void GridKernel::grid(casacore::Matrix<casacore::Complex>& grid,
		casacore::Matrix<casacore::Complex>& convFunc, const casacore::Complex& cVis,
		const int iu, const int iv, const int support) {

#if defined ( ASKAP_GRID_WITH_POINTERS ) || defined ( ASKAP_GRID_WITH_BLAS )
#if defined ( ASKAP_GRID_WITH_POINTERS )
    casacore::Float rVis = cVis.real();
    casacore::Float iVis = cVis.imag();
#endif
	for (int suppv = -support; suppv < +support; suppv++) {
		const int voff = suppv + support;
		const int uoff = -support + support;
#ifdef ASKAP_GRID_WITH_BLAS
        casacore::Complex *wtPtr = &convFunc(uoff, voff);
        casacore::Complex *gridPtr = &(grid(iu - support, iv + suppv));
		cblas_caxpy(2*support+1, &cVis, wtPtr, 1, gridPtr, 1);
#else
        // Writing the multiply in real/imag is twice as fast with gcc
        casacore::Float *wtPtrF = reinterpret_cast<casacore::Float *> (&convFunc(uoff, voff));
        casacore::Float *gridPtrF = reinterpret_cast<casacore::Float *> (&grid(iu - support, iv + suppv));
        for (int suppu = -support; suppu < +support; suppu++, wtPtrF+=2, gridPtrF+=2) {
            gridPtrF[0] += rVis * wtPtrF[0] - iVis * wtPtrF[1];
            gridPtrF[1] += rVis * wtPtrF[1] + iVis * wtPtrF[0];
    	}
#endif
	}
#else
	for (int suppv=-support; suppv<+support; suppv++)
	{
		const int voff=suppv+support;
		for (int suppu=-support; suppu<+support; suppu++)
		{
			const int uoff=suppu+support;
			casacore::Complex wt=convFunc(uoff, voff);
			grid(iu+suppu, iv+suppv)+=cVis*wt;
		}
	}
#endif
}

/// Totally selfcontained degridding
void GridKernel::degrid(casacore::Complex& cVis,
		const casacore::Matrix<casacore::Complex>& convFunc,
		const casacore::Matrix<casacore::Complex>& grid,
        const int iu, const int iv, const int support) {
	/// Degridding from grid to visibility. Here we just take a weighted sum of the visibility
	/// data using the convolution function as the weighting function.
	cVis = 0.0;
#if defined ( ASKAP_GRID_WITH_POINTERS ) || defined ( ASKAP_GRID_WITH_BLAS )
	for (int suppv = -support; suppv < +support; suppv++) {
		const int voff = suppv + support;
		const int uoff = -support + support;
        const casacore::Complex *wtPtr = &convFunc(uoff, voff);
        const casacore::Complex *gridPtr = &(grid(iu - support, iv + suppv));
#ifdef ASKAP_GRID_WITH_BLAS
		casacore::Complex dot;
		cblas_cdotc_sub(2*support+1, gridPtr, 1, wtPtr, 1, &dot);
		cVis+=dot;
#else
        // Writing the multiply in real/imag is twice as fast with gcc
        // Doing the 'reinterpret_cast' thing to avoid complex like in grid, doesn't help here.
        for (int suppu = -support; suppu < +support; suppu++, wtPtr++, gridPtr++) {
            cVis += casacore::Complex( (*wtPtr).real()*(*gridPtr).real()+(*wtPtr).imag()*(*gridPtr).imag(),
                                  -(*wtPtr).real()*(*gridPtr).imag()+(*wtPtr).imag()*(*gridPtr).real());
		}
#endif
	}
#else
	for (int suppv=-support; suppv<+support; suppv++)
	{
		const int voff=suppv+support;
		for (int suppu=-support; suppu<+support; suppu++)
		{
			const int uoff=suppu+support;
			casacore::Complex wt=convFunc(uoff, voff);
			cVis+=wt*conj(grid(iu+suppu, iv+suppv));
		}
	}
#endif
}

}
}
