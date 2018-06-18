/// @file AltWProjectVisGridder.cc
///
/// @copyright (c) 2007,2016 CSIRO
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

// Package level header file
#include <askap_synthesis.h>

// System includes
#include <cmath>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <askap/AskapUtil.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <fft/FFTWrapper.h>
#include <profile/AskapProfiler.h>
//ASKAPSoft package includes
#include <gridding/WProjectVisGridder.h>
#include <gridding/SupportSearcher.h>

// Local package includes
#include "AltWProjectVisGridder.h"



ASKAP_LOGGER(logger, ".gridding.AltWProjectVisGridder");

namespace askap {
namespace synthesis {

AltWProjectVisGridder::AltWProjectVisGridder(const double wmax,
                                       const int nwplanes,
                                       const double cutoff,
                                       const int overSample,
                                       const int maxSupport,
                                       const int limitSupport,
                                       const std::string& name,
                                       const float alpha) :
        WProjectVisGridder(wmax, nwplanes, cutoff, overSample,maxSupport,limitSupport,name,alpha)

{

}

AltWProjectVisGridder::~AltWProjectVisGridder()
{
}


/// Clone a copy of this Gridder
IVisGridder::ShPtr AltWProjectVisGridder::clone()
{
    return IVisGridder::ShPtr(new AltWProjectVisGridder(*this));
}


/// @brief static method to create gridder
/// @details Each gridder should have a static factory method, which is
/// able to create a particular type of the gridder and initialise it with
/// the parameters taken form the given parset. It is assumed that the
/// method receives a subset of parameters where the gridder name is already
/// taken out.
/// @param[in] parset input parset file
/// @return a shared pointer to the gridder instance
IVisGridder::ShPtr AltWProjectVisGridder::createGridder(const LOFAR::ParameterSet& parset)
{
    const double wmax = parset.getDouble("wmax", 35000.0);
    const int nwplanes = parset.getInt32("nwplanes", 65);
    const double cutoff = parset.getDouble("cutoff", 1e-3);
    const int oversample = parset.getInt32("oversample", 8);
    const int maxSupport = parset.getInt32("maxsupport", 256);
    const int limitSupport = parset.getInt32("limitsupport", 0);
    const string tablename = parset.getString("tablename", "");
    const float alpha=parset.getFloat("alpha", 1.);

    ASKAPLOG_INFO_STR(logger, "Gridding using W projection with " << nwplanes << " w-planes");
    boost::shared_ptr<AltWProjectVisGridder> gridder(new AltWProjectVisGridder(wmax, nwplanes,
            cutoff, oversample, maxSupport, limitSupport, tablename, alpha));
    gridder->configureGridder(parset);
    gridder->configureWSampling(parset);
    return gridder;
}



/// @brief assignment operator
/// @details Defined as private, so it can't be called (to enforce usage of the
/// copy constructor
/// @param[in] other input object
/// @return reference to itself
AltWProjectVisGridder& AltWProjectVisGridder::operator=(const AltWProjectVisGridder &)
{
    ASKAPTHROW(AskapError, "This method is not supposed to be called!");
    return *this;
}


} // namespace askap
} // namespace synthesis
