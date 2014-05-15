/// @file DeconvolverFactory.cc
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

#include <askap_synthesis.h>

#include <askap/AskapError.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".deconvolver.factory");
#include <casa/BasicSL/String.h>   // for downcase
#include <images/Images/PagedImage.h>

#include <deconvolution/DeconvolverHelpers.h>
#include <deconvolution/DeconvolverFactory.h>
#include <deconvolution/DeconvolverBase.h>
#include <deconvolution/DeconvolverBasisFunction.h>
#include <deconvolution/DeconvolverMultiTermBasisFunction.h>
#include <deconvolution/DeconvolverEntropy.h>
#include <deconvolution/DeconvolverFista.h>
#include <deconvolution/DeconvolverHogbom.h>
#include <deconvolution/MultiScaleBasisFunction.h>

namespace askap {
    namespace synthesis {

        DeconvolverFactory::DeconvolverFactory()
        {
        }

        DeconvolverBase<Float, Complex>::ShPtr DeconvolverFactory::make(const LOFAR::ParameterSet &parset)
        {

            DeconvolverBase<Float, Complex>::ShPtr deconvolver;

            if (parset.getString("solver") == "Fista") {
                ASKAPLOG_INFO_STR(logger, "Constructing Fista deconvolver");
                Array<Float> dirty(DeconvolverHelpers::getArrayFromImage("dirty", parset));
                Array<Float> psf(DeconvolverHelpers::getArrayFromImage("psf", parset));
                deconvolver.reset(new DeconvolverFista<Float, Complex>(dirty, psf));
                ASKAPASSERT(deconvolver);

                // Now get the control parameters
                LOFAR::ParameterSet subset(parset.makeSubset("solver.Fista."));
                deconvolver->configure(subset);
            } else if (parset.getString("solver") == "Entropy") {
                ASKAPLOG_INFO_STR(logger, "Constructing Entropy deconvolver");
                Array<Float> dirty(DeconvolverHelpers::getArrayFromImage("dirty", parset));
                Array<Float> psf(DeconvolverHelpers::getArrayFromImage("psf", parset));
                deconvolver.reset(new DeconvolverEntropy<Float, Complex>(dirty, psf));
                ASKAPASSERT(deconvolver);

                deconvolver->configure(parset.makeSubset("solver.Entropy."));
            } else {
                Array<Float> dirty(DeconvolverHelpers::getArrayFromImage("dirty", parset));
                Array<Float> psf(DeconvolverHelpers::getArrayFromImage("psf", parset));

                string algorithm = parset.getString("solver.Clean.algorithm", "Basisfunction");

                if (algorithm == "Basisfunction") {
                    ASKAPLOG_INFO_STR(logger, "Constructing Basisfunction Clean solver");
                    deconvolver.reset(new DeconvolverBasisFunction<Float, Complex>(dirty, psf));
                    ASKAPASSERT(deconvolver);
                } else if (algorithm == "MultiTermBasisfunction") {
                    ASKAPLOG_INFO_STR(logger, "Constructing MultiTermBasisfunction Clean solver");
                    deconvolver.reset(new DeconvolverMultiTermBasisFunction<Float, Complex>(dirty, psf));
                    ASKAPASSERT(deconvolver);
                } else if (algorithm == "Hogbom") {
                    ASKAPLOG_INFO_STR(logger, "Constructing Hogbom Clean deconvolver");
                    deconvolver.reset(new DeconvolverHogbom<Float, Complex>(dirty, psf));
                    ASKAPASSERT(deconvolver);
                } else {
                    ASKAPTHROW(AskapError, "Unknown Clean algorithm " << algorithm);
                }
                deconvolver->configure(parset.makeSubset("solver.Clean."));
            }
            if (parset.getString("weight", "") != "") {
                deconvolver->setWeight(DeconvolverHelpers::getArrayFromImage("weight", parset));
            }
            return deconvolver;

        }
    }
}
