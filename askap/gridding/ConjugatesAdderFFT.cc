/// @file
/// @brief Class adding conjugates to the plane of accumulated uv-weights
/// @details This is an implementation of UVWeight calculator interface which adds conjugates to
/// ensure the FT of it is real. It is not necessary on its own (we could've just ignored the imaginary
/// part of the FT instead), but allows us to rely on an assumption in other filters that accumulated weights are 
/// symmetric and simplify their code. This code follows the original approach by Daniel Mitchell via FFT. 
/// Note, it is possible to write a weight builder class which ensures conjugate symmetry automatically, but 
/// it will likely take more operations (unless observation is very short). It seems possible, however, to do it
/// without FFT as as only relative weight matters. For now, do it as a separate step. This will allow to 
/// change algorithm without introducing much technical debt.
/// 
///
/// @copyright (c) 2023 CSIRO
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

// own includes
#include <askap/gridding/ConjugatesAdderFFT.h>
#include <askap/scimath/fft/FFT2DWrapper.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <askap/profile/AskapProfiler.h>

#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.conjugatesadderfft");

namespace askap {

namespace synthesis {

/// @brief perform processing for the given weight (single grid slice along the 3rd axis)
/// @details For performance reasons, slices along the 3rd axis are taken inside finalise method
/// of the builder (this can be changed if we ever had any effect where frequency dependence matter).
/// At this stage, we can guarantee that supplied matrix has contiguous storage.
/// @param[in] wt weight to work with (it is modified in situ).
/// @note The shape is supposed to stay intact.
void ConjugatesAdderFFT::process(casacore::Matrix<float> &wt) const
{
   ASKAPTRACE("ConjugatesAdderFFT::process");
   ASKAPLOG_DEBUG_STR(logger, "Before ConjugatesAdderFFT sum of grid = " << sum(wt));

   // Note, we could've moved wrapper to the class body. It would allow us to reuse the plan on subsequent calls.
   // However, we expect to build and finalise weight only once and then reuse the result. So if the high-level code
   // does the right thing, this caching won't be necessary.
   // Also, for this particular use it doesn't matter where the centre is, so bypass shifting harmonics
   scimath::FFT2DWrapper<casacore::Complex> wrapper(false);

   casacore::Matrix<casacore::Complex> buffer(wt.shape(), casacore::Complex(0.f));
   casacore::setReal(buffer, wt);
   wrapper(buffer, false);
   buffer += casacore::conj(buffer);
   wrapper(buffer, true);
   casacore::real(wt, buffer);
   ASKAPLOG_DEBUG_STR(logger, "After ConjugatesAdderFFT sum of grid = " << sum(wt));
}

} // namespace synthesis

} // namespace askap

