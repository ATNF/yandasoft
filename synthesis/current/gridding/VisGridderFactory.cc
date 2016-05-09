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

// Package level header file
#include <askap_synthesis.h>

// ASKAPsoft includes
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.visgridderfactory");
#include <askap/AskapError.h>
#include <casacore/casa/OS/DynLib.h>        // for dynamic library loading
#include <casacore/casa/BasicSL/String.h>   // for downcase
#include <casacore/scimath/Mathematics/Interpolate2D.h>

// Local package includes
#include <gridding/VisGridderFactory.h>
#include <gridding/VisGridderWithPadding.h>
#include <gridding/BoxVisGridder.h>
#include <gridding/SphFuncVisGridder.h>
#include <gridding/WProjectVisGridder.h>
#include <gridding/AWProjectVisGridder.h>
#include <gridding/WStackVisGridder.h>
#include <gridding/AProjectWStackVisGridder.h>
#include <gridding/SnapShotImagingGridderAdapter.h>
#include <gridding/SmearingGridderAdapter.h>
#include <gridding/VisWeightsMultiFrequency.h>
#include <measurementequation/SynthesisParamsHelper.h>

namespace askap {
namespace synthesis {

  // Define the static registry.
  std::map<std::string, VisGridderFactory::GridderCreator*>
  VisGridderFactory::theirRegistry;


  VisGridderFactory::VisGridderFactory() {
  }

  void VisGridderFactory::registerGridder (const std::string& name,
                                           VisGridderFactory::GridderCreator* creatorFunc)
  {
    ASKAPLOG_DEBUG_STR(logger, "     - Adding "<<name<<" gridder to the registry");
    theirRegistry[name] = creatorFunc;
  }

  IVisGridder::ShPtr VisGridderFactory::createGridder (const std::string& name,
                                    const LOFAR::ParameterSet& parset)
  {
    std::map<std::string,GridderCreator*>::const_iterator it = theirRegistry.find (name);
    if (it == theirRegistry.end()) {
      // Unknown gridder. Try to load the data manager from a dynamic library
      // with that lowercase name (without possible template extension).
      std::string libname(toLower(name));
      const std::string::size_type pos = libname.find_first_of (".<");
      if (pos != std::string::npos) {
        libname = libname.substr (0, pos);      // only take before . or <
      }
      // Try to load the dynamic library and execute its register function.
      // Do not dlclose the library.
      ASKAPLOG_INFO_STR(logger, "Gridder "<<name<<
                 " is not in the registry, attempting to load it dynamically");
      casa::DynLib dl(libname, string("libaskap_"), "register_"+libname, false);
      if (dl.getHandle()) {
        // Successfully loaded. Get the creator function.
        ASKAPLOG_INFO_STR(logger, "Dynamically loaded gridder " << name);
        // the first thing the gridder in the shared library is supposed to do is to
        // register itself. Therefore, its name will appear in the registry.
        it = theirRegistry.find (name);
      }
    }
    if (it == theirRegistry.end()) {
      ASKAPTHROW(AskapError, "Unknown gridder " << name);
    }
    // Execute the registered function.
    return it->second(parset);
  }

  // Make the gridder object for the gridder given in the parset file.
  // Currently the standard gridders are still handled by this function.
  // In the (near) future it should be done by putting creator functions
  // for these gridders in the registry and use that.
IVisGridder::ShPtr VisGridderFactory::make(const LOFAR::ParameterSet &parset) {
    if (theirRegistry.size() == 0) {
        // this is the first call of the method, we need to fill the registry with
        // all pre-defined gridders
        ASKAPLOG_DEBUG_STR(logger, "Filling the gridder registry with pre-defined gridders");
        addPreDefinedGridder<WProjectVisGridder>();
        addPreDefinedGridder<WStackVisGridder>();
        addPreDefinedGridder<AProjectWStackVisGridder>();
        addPreDefinedGridder<AWProjectVisGridder>();
        addPreDefinedGridder<BoxVisGridder>();
        addPreDefinedGridder<SphFuncVisGridder>();        
    }

    // buffer for the result
    IVisGridder::ShPtr gridder;
    /// @todo Better handling of string case
    std::string prefix("gridder");	
    const string gridderName = parset.getString(prefix);
    prefix += "." + gridderName + ".";
    ASKAPLOG_DEBUG_STR(logger, "Attempting to greate gridder "<<gridderName);
    gridder = createGridder (gridderName, parset.makeSubset(prefix));

    ASKAPASSERT(gridder);
    if (parset.isDefined("gridder.padding")) {
        const float padding =parset.getFloat("gridder.padding");
        ASKAPLOG_INFO_STR(logger, "Use padding at the gridder level, padding factor = "<<padding);
        boost::shared_ptr<VisGridderWithPadding> vg = 
            boost::dynamic_pointer_cast<VisGridderWithPadding>(gridder);
        ASKAPCHECK(vg, "Gridder type ("<<parset.getString("gridder")<<
                ") is incompatible with the padding option");
        vg->setPaddingFactor(padding);           
    } else {
        ASKAPLOG_INFO_STR(logger,"No padding at the gridder level");
    }

    if (parset.isDefined("gridder.MaxPointingSeparation")) {
        const double threshold = SynthesisParamsHelper::convertQuantity(
                    parset.getString("gridder.MaxPointingSeparation","-1rad"),"rad");
        ASKAPLOG_INFO_STR(logger,"MaxPointingSeparation is used, data from pointing centres further than "<<
                threshold*180./casa::C::pi<<" deg from the image centre will be rejected");
        boost::shared_ptr<TableVisGridder> tvg = 
	        boost::dynamic_pointer_cast<TableVisGridder>(gridder);
	    ASKAPCHECK(tvg, "Gridder type ("<<gridderName<<") is incompatible with the MaxPointingSeparation option");
	    tvg->maxPointingSeparation(threshold);        
    } else {
            ASKAPLOG_INFO_STR(logger,"MaxPointingSeparation is not used for gridder: "<<gridderName<<
                                     ", all data will be used");
    }
    
    if (parset.isDefined("gridder.alldatapsf")) {
        const bool useAll = parset.getBool("gridder.alldatapsf");
        if (useAll) {
            ASKAPLOG_INFO_STR(logger, "Use all data for PSF calculations instead of the representative feed and field");
        } else {
            ASKAPLOG_INFO_STR(logger, "Use representative feed and field for PSF calculation");
        }
        boost::shared_ptr<TableVisGridder> tvg = 
            boost::dynamic_pointer_cast<TableVisGridder>(gridder);
        ASKAPCHECK(tvg, "Gridder type ("<<parset.getString("gridder")<<
                ") is incompatible with the alldatapsf option");
        tvg->useAllDataForPSF(useAll);
    } else {
        ASKAPLOG_INFO_STR(logger, "gridder.alldatapsf option is not used, default to representative feed and field for PSF calculation");
    }	

    {
        const bool osWeight = parset.getBool("gridder.oversampleweight",false);
        if (osWeight) {
            ASKAPLOG_INFO_STR(logger, "Weight will be tracked per oversampling plane");
        } else {
            ASKAPLOG_INFO_STR(logger, "First oversampling plane will always be used for weight computation");
        }
        boost::shared_ptr<TableVisGridder> tvg = 
            boost::dynamic_pointer_cast<TableVisGridder>(gridder);
        if (tvg) {
            tvg->trackWeightPerPlane(osWeight);
        } else {
            ASKAPLOG_WARN_STR(logger,"Gridder type ("<<parset.getString("gridder")<<
                    ") is incompatible with the oversampleweight option (trying to set it to "<<osWeight<<")");
        }
    }	

    // Initialize the Visibility Weights
    if (parset.getString("visweights","")=="MFS")
    {
        double reffreq=parset.getDouble("visweights.MFS.reffreq", 1.405e+09);
        ASKAPLOG_INFO_STR(logger, "Initialising for MFS with reference frequency " << reffreq << " Hz");
        ASKAPLOG_INFO_STR(logger, "Assuming that the number of Taylor terms is greater than 1");

        gridder->initVisWeights(IVisWeights::ShPtr(new VisWeightsMultiFrequency(reffreq)));
    }
    else // Null....
    {
        //		gridder->initVisWeights(IVisWeights::ShPtr(new VisWeightsMultiFrequency()));
    }

    if (parset.getBool("gridder.snapshotimaging",false)) {
        ASKAPLOG_INFO_STR(logger, "A gridder adapter will be set up to do snap-shot imaging");
        const double wtolerance = parset.getDouble("gridder.snapshotimaging.wtolerance"); 
        ASKAPLOG_INFO_STR(logger, "w-coordinate tolerance is "<<wtolerance<<" wavelengths");
        const casa::uInt decimate = parset.getUint("gridder.snapshotimaging.coorddecimation", 3);
        const casa::Interpolate2D::Method method = interpMethod(
                parset.getString("gridder.snapshotimaging.interpmethod", "cubic"));

        const bool PredictWPlane = parset.getBool("gridder.snapshotimaging.longtrack",false);
        if (PredictWPlane) {    
            ASKAPLOG_INFO_STR(logger, "W-plane fitting via prediction will be used in snapshot mode");
        }
        boost::shared_ptr<SnapShotImagingGridderAdapter> adapter(
                new SnapShotImagingGridderAdapter(gridder, wtolerance, decimate, method, PredictWPlane));
        const double clippingFactor = parset.getDouble("gridder.snapshotimaging.clipping", 0.);
        ASKAPLOG_INFO_STR(logger, "Clipping factor is "<<clippingFactor);
        adapter->setClippingFactor(float(clippingFactor));
        const bool doPSFReprojection = parset.getBool("gridder.snapshotimaging.reprojectpsf", false);
        adapter->setPSFReprojection(doPSFReprojection);
        // possible additional configuration comes here
        gridder = adapter;
    }

    if (parset.getBool("gridder.bwsmearing",false)) {
        ASKAPLOG_INFO_STR(logger, "A gridder adapter will be set up to simulate bandwidth smearing");
        const double chanBW = parset.getDouble("gridder.bwsmearing.chanbw",1e6);
        const casa::Int nSteps = parset.getInt32("gridder.bwsmearing.nsteps", 10);
        ASKAPCHECK(nSteps > 0, "Number of steps is supposed to be positive");
        ASKAPLOG_INFO_STR(logger, "  assumed channel bandwidth = "<<chanBW/1e6<<" MHz, number of integration steps = "<<nSteps); 
        boost::shared_ptr<SmearingGridderAdapter> adapter(new SmearingGridderAdapter(gridder, chanBW, casa::uInt(nSteps)));
        ASKAPDEBUGASSERT(adapter);
        gridder = adapter;        
    } 

    return gridder;
}

casa::Interpolate2D::Method VisGridderFactory::interpMethod(casa::String str) {
    str.downcase();
    if (str.compare("nearest") == 0) {
        ASKAPLOG_INFO_STR(logger, "Regridding interpolation method: NEAREST");
        return casa::Interpolate2D::NEAREST;
    } else if (str.compare("linear") == 0) {
        ASKAPLOG_INFO_STR(logger, "Regridding interpolation method: LINEAR");
        return casa::Interpolate2D::LINEAR;
    } else if (str.compare("cubic") == 0) {
        ASKAPLOG_INFO_STR(logger, "Regridding interpolation method: CUBIC");
        return casa::Interpolate2D::CUBIC;
    } else if (str.compare("lanczos") == 0) {
        ASKAPLOG_INFO_STR(logger, "Regridding interpolation method: LANCZOS");
        return casa::Interpolate2D::LANCZOS;
    } else {
        ASKAPTHROW(AskapError, "Unknown interpolation method: " << str);
    }
}

} // namespace synthesis

} // namespace askap

