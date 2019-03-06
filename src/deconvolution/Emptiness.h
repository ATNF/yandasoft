/// @file Emptiness.h
/// @brief Base class for a deconvolver
/// @details This interface class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
/// @ingroup Deconvolver
///
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

#ifndef ASKAP_SYNTHESIS_EMPTINESS_H
#define ASKAP_SYNTHESIS_EMPTINESS_H

#include <deconvolution/EntropyBase.h>

namespace askap {

    namespace synthesis {

        // <summary> Maximum Emptiness measure used by MEM
        // </summary>

        // Emptiness measure
        template<class T>
        class Emptiness : public EntropyBase<T> {
            public:

                typedef boost::shared_ptr<Emptiness<T> > ShPtr;

                enum GRADTYPE {H = 0, C, F, J };

                // calculate the entropy for the whole image
                virtual T entropy(const casa::Array<T>& model);

                // calculate the entropy for the whole image
                virtual T entropy(const casa::Array<T>& model, const casa::Array<T>& mask);

                // calculate the gradient entropy for the whole image
                virtual void gradEntropy(casa::Array<T>& gradH, casa::Array<T>& rHess,
                                         const casa::Array<T>& model);

                // calculate the gradient entropy for the whole image
                virtual void gradEntropy(casa::Array<T>& gradH, casa::Array<T>& rHess,
                                         const casa::Array<T>& model, const Array<T>& mask);
        };

    } // namespace synthesis

} // namespace askap

#include <deconvolution/Emptiness.tcc>

#endif
