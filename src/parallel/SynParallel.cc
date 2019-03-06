/// @file
///
/// @brief Base class for parallel applications
/// @details
/// Supports algorithms by providing methods for initialization
/// of MPI connections, sending and models around.
/// There is assumed to be one master and many workers.
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
#include <parallel/SynParallel.h>

#include <measurementequation/SynthesisParamsHelper.h>
#include <measurementequation/ImageParamsHelper.h>
#include <gridding/VisGridderFactory.h>
#include <gridding/TableVisGridder.h>
#include <gridding/IVisGridder.h>
#include <profile/AskapProfiler.h>


#include <sstream>
#include <set>
#include <vector>
#include <string>
#include <algorithm>

#include <Blob/BlobString.h>
#include <Blob/BlobIBufString.h>
#include <Blob/BlobOBufString.h>
#include <Blob/BlobIStream.h>
#include <Blob/BlobOStream.h>

#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Utilities/Regex.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/OS/Path.h>

#include <askap_synthesis.h>
#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".parallel");

#include <askap/AskapError.h>
#include <askap_synthesis.h>

using namespace std;
using namespace askap;
using namespace askap::askapparallel;
using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    SynParallel::SynParallel(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset) : 
                         itsComms(comms), itsParset(parset)
    {
      itsModel.reset(new Params());
      ASKAPCHECK(itsModel, "Model not defined correctly");

      // setup frequency frame
      const std::string freqFrame = parset.getString("freqframe","topo");
      if (freqFrame == "topo") {
          ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as topocentric");
          itsFreqRefFrame = casa::MFrequency::Ref(casa::MFrequency::TOPO);
      } else if (freqFrame == "lsrk") {
          ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as lsrk");
          itsFreqRefFrame = casa::MFrequency::Ref(casa::MFrequency::LSRK);
      } else if (freqFrame == "bary") {
          ASKAPLOG_INFO_STR(logger, "Parset frequencies will be treated as barycentric");
          itsFreqRefFrame = casa::MFrequency::Ref(casa::MFrequency::BARY);
      } else {
          ASKAPTHROW(AskapError, "Unsupported frequency frame "<<freqFrame);
      }    
      const std::string parString = itsComms.isParallel() ? "parallel" : "serial";
      std::string mwString = itsComms.isWorker() ? "worker" : "master";
      if (itsComms.isWorker() && itsComms.isMaster()) {
          mwString += "&master";
      } 
      ASKAPLOG_INFO_STR(logger, "SynParallel in "<<parString<<" mode("<<mwString<<"), rank = "<<itsComms.rank()<<
                        " nProcs="<<itsComms.nProcs());
    }

    SynParallel::~SynParallel()
    {
    }

    askap::scimath::Params::ShPtr& SynParallel::params()
    {
      return itsModel;
    }

    // Send the model to all workers
    void SynParallel::broadcastModel()
    {
      ASKAPTRACE("SynParallel::broadcastModel");

      if (itsComms.isParallel() && itsComms.isMaster())
      {
        ASKAPCHECK(itsModel, "Model not defined prior to broadcast")
        ASKAPLOG_DEBUG_STR(logger, "Current model held by the master: "<<*itsModel);
        casa::Timer timer;
        timer.mark();

        const std::vector<std::string> names = parametersToBroadcast();
        if (itsComms.nGroups() == 1) {
            ASKAPLOG_INFO_STR(logger, "Sending the whole model to all workers");
            if (names.size() == itsModel->names().size()) {
                ASKAPLOG_INFO_STR(logger, "About to broadcast all model parameters: "<<names);
                broadcastModelImpl(*itsModel);
            } else {
                ASKAPLOG_INFO_STR(logger, "About to broadcast the following model parameters: "<<names);
                scimath::Params buffer;
                buffer.makeSlice(*itsModel, names);
                broadcastModelImpl(buffer);
            }
        } else {
            ASKAPLOG_INFO_STR(logger, "Distribute model between "<<itsComms.nGroups()<<
                  " groups of workers");
            // build two lists of parameters: parameters to distribute and parameters to send to all groups
            std::vector<std::string> names2distribute;
            std::vector<std::string> names2keep;
            names2distribute.reserve(names.size());
            names2keep.reserve(names.size());            
            for (std::vector<std::string>::const_iterator ci = names.begin(); ci!=names.end(); ++ci) {
                 // distribute only parameters starting with "image" for now
                 if (ci->find("image") == 0) {
                     names2distribute.push_back(*ci);
                 } else {
                     names2keep.push_back(*ci);
                 }
            }
            //
            ASKAPDEBUGASSERT(itsComms.nGroups() > 1);
            // number of parameters per group (note the last group can have more)
            const size_t nPerGroup = names2distribute.size() / itsComms.nGroups();
            // this check is not relevant if all parameters are in names2keep
            if (names2distribute.size() > 0) {
                ASKAPCHECK(nPerGroup > 0, "The model has too few parameters ("<<
                      names2distribute.size()<<") to distribute between "<< itsComms.nGroups()<<" groups");
            } else {
                ASKAPCHECK(names2keep.size() > 0, "The model has too few parameters ("<<
                      names2keep.size()<<")");
            }
            
            std::vector<std::string> currentNames;
            currentNames.reserve(itsComms.nGroups() + nPerGroup - 1 + names2keep.size());
            scimath::Params buffer;
            for (size_t group = 0, index = 0; group<itsComms.nGroups(); ++group, index+=nPerGroup) {
                 const size_t nPerCurrentGroup = (group + 1 < itsComms.nGroups()) ? 
                          nPerGroup : names2distribute.size() - index; 
                 ASKAPDEBUGASSERT((names2distribute.size() > index) || (names2distribute.size() == 0));
                 if (nPerCurrentGroup != nPerGroup) {
                     ASKAPLOG_WARN_STR(logger, "An unbalanced distribution of the model has been detected. "
                                       " the last group ("<<group<<") will have "<<nPerCurrentGroup<<
                                       " parameters vs. "<<nPerGroup<<" for other groups");
                 }
                 currentNames.resize(nPerCurrentGroup + names2keep.size());
                 for (size_t i = 0; i<nPerCurrentGroup; ++i) {
                      currentNames[i] = names2distribute[index + i];
                 }
                 for (size_t i = 0; i<names2keep.size(); ++i) {
                      currentNames[nPerCurrentGroup + i] = names2keep[i];
                 }
                 ASKAPLOG_INFO_STR(logger, "Group "<<group<<
                        " will get the following parameters: "<<currentNames);
                 buffer.makeSlice(*itsModel, currentNames);
                 ASKAPLOG_INFO_STR(logger, "Sending the model to appropriate workers (group "<<
                                 group<<") ");
                 itsComms.useGroupOfWorkers(group);
                 try {
                    broadcastModelImpl(buffer);
                 }
                 catch (const std::exception &) {
                    itsComms.useAllWorkers();
                    throw;
                 }
                 itsComms.useAllWorkers();
            }
        }

        ASKAPLOG_INFO_STR(logger, "Broadcast model to the workers in "<< timer.real()
                           << " seconds ");
      }
    }
    // drop a new model in place
    void SynParallel::replaceModel(scimath::Params::ShPtr Model)
    {
        ASKAPTRACE("SynParallel::replaceModel");
        if (itsComms.isParallel() && itsComms.isWorker())
        {
            ASKAPCHECK(itsModel, "Model not defined prior to receiving")
            // copy over the model - not just the reference
            *itsModel = *Model;
        }
    }
    // Receive the model from the master
    void SynParallel::receiveModel()
    {
      ASKAPTRACE("SynParallel::receiveModel");

      if (itsComms.isParallel() && itsComms.isWorker())
      {
        ASKAPCHECK(itsModel, "Model not defined prior to receiving")
        casa::Timer timer;
        timer.mark();

        if (itsComms.nGroups() == 1) {
            ASKAPLOG_INFO_STR(logger, "Wait to receive the whole model from the master");
            receiveModelImpl(*itsModel);
        } else {
            size_t currentGroup = itsComms.nGroups(); // just a flag that group index is not found
            for (size_t group = 0; group< itsComms.nGroups(); ++group) {
                 if (itsComms.inGroup(group)) {
                     ASKAPCHECK(currentGroup == itsComms.nGroups(), 
                           "Each worker can belong to one and only one group! "
                           "For some reason it belongs to groups "<<currentGroup<<" and "<<group);
                     currentGroup = group;
                 }
            }
            ASKAPCHECK(currentGroup < itsComms.nGroups(), "The worker at rank="<<itsComms.rank()<<
                       "does not seem to belong to any group!");
            ASKAPLOG_INFO_STR(logger, 
                 "Wait to receive from the master a part of the model appropriate for the group "<<
                 currentGroup);
            itsComms.useGroupOfWorkers(currentGroup);
            try {
               receiveModelImpl(*itsModel);
            }
            catch (const std::exception &) {
               itsComms.useAllWorkers();
               throw;
            }
            itsComms.useAllWorkers();
        }

        ASKAPLOG_INFO_STR(logger, "Received model from the master in "<< timer.real()
                           << " seconds ");
        ASKAPLOG_DEBUG_STR(logger, "Current model held by the worker: "<<*itsModel);
      }
    }
      
    /// @brief actual implementation of the model broadcast
    /// @details This method is only supposed to be called from the master.
    /// @param[in] model the model to send
    void SynParallel::broadcastModelImpl(const scimath::Params &model)
    {
        ASKAPDEBUGTRACE("SynParallel::broadcastModelImpl");

        ASKAPDEBUGASSERT(itsComms.isParallel() && itsComms.isMaster());
        LOFAR::BlobString bs;
        bs.resize(0);
        LOFAR::BlobOBufString bob(bs);
        LOFAR::BlobOStream out(bob);
        out.putStart("model", 1);
        out << model;
        out.putEnd();
        itsComms.broadcastBlob(bs ,0);
    }

    /// @brief actual implementation of the model receive
    /// @details This method is only supposed to be called from workers. 
    /// There should be one to one match between the number of calls to 
    /// broadcastModelImpl and receiveModelImpl.
    /// @param[in] model the model to fill
    void SynParallel::receiveModelImpl(scimath::Params &model)
    {
        ASKAPDEBUGTRACE("SynParallel::receiveModelImpl");

        ASKAPDEBUGASSERT(itsComms.isParallel() && itsComms.isWorker());
        LOFAR::BlobString bs;
        bs.resize(0);
        itsComms.broadcastBlob(bs, 0);
        LOFAR::BlobIBufString bib(bs);
        LOFAR::BlobIStream in(bib);
        int version=in.getStart("model");
        ASKAPASSERT(version==1);
        in >> model;
        in.getEnd();
    }
    
    /// @brief helper method to identify model parameters to broadcast
    /// @details We use itsModel to buffer some derived images like psf, weights, etc
    /// which are not required for prediffers. It just wastes memory and CPU time if
    /// we broadcast them. At the same time, some auxilliary parameters like peak
    /// residual value need to be broadcast (so the major cycle can terminate in workers).
    /// This method returns the vector with all parameters to be broadcast. By default
    /// it returns all parameter names. This method is supposed to be overridden in
    /// derived classes (e.g. ImagerParallel) where a different behavior is needed.
    /// @return a vector with parameters to broadcast
    std::vector<std::string> SynParallel::parametersToBroadcast() const
    {
       ASKAPDEBUGASSERT(itsModel);
       return itsModel->names();
    }
    

    std::string SynParallel::substitute(const std::string& s) const
    {
       return itsComms.substitute(s);
    }

    /// @brief helper method to create and configure gridder
    /// @details It is expected to be called from the constructor of derived classes
    /// @param[in] comms communications object
    /// @param[in] parset parameter set      
    IVisGridder::ShPtr SynParallel::createGridder(const askap::askapparallel::AskapParallel& comms, 
                           const LOFAR::ParameterSet& parset)
    {
       // Create the gridder using a factory acting on a parameterset
       IVisGridder::ShPtr gridder = VisGridderFactory::make(parset);
       ASKAPCHECK(gridder, "Gridder is not defined correctly");              
       if (comms.isParallel()) {
           const int rankStoringCF = parset.getInt32("rankstoringcf", 1);
           if (comms.rank() == rankStoringCF) {
               ASKAPLOG_INFO_STR(logger, "Rank "<<rankStoringCF<<
                       " will attempt to export convolution functions (if export is requested)");
           } else {
               boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(gridder);
               if (tvg) {
                   ASKAPLOG_INFO_STR(logger, "Export of convolution functions will be inhibited for rank "<<
                         comms.rank());
                   tvg->setTableName("");
               } else {
                   ASKAPLOG_INFO_STR(logger, "Unable to inhibit export of CFs for rank "<<
                         comms.rank()<<" - operation not supported by the chosen type of the gridder");
               }
           }
       }
       return gridder;
    }
    
    /// @brief read the models from parset file to the given params object
    /// @details The model can be composed from both images and components. This
    /// method populates Params object by adding model data read from the parset file.
    /// The model is given by shared pointer because the same method can be used for both
    /// simulations and calibration (the former populates itsModel, the latter populates
    /// itsPerfectModel) 
    /// @param[in] pModel shared pointer to the params object (must exist)
    void SynParallel::readModels(const scimath::Params::ShPtr &pModel) const
    {
      ASKAPTRACE("SynParallel::readModels");

      ASKAPCHECK(pModel, "model is not initialised prior to call to SynParallel::readModels");
      
      LOFAR::ParameterSet parset(itsParset);
  
      if (itsParset.isDefined("sources.definition")) {
          parset = LOFAR::ParameterSet(substitute(itsParset.getString("sources.definition")));
      }
      
      const std::vector<std::string> sources = parset.getStringVector("sources.names");
      std::set<std::string> loadedImageModels;
      for (size_t i=0; i<sources.size(); ++i) {
	       const std::string modelPar = std::string("sources.")+sources[i]+".model";
	       const std::string compPar = std::string("sources.")+sources[i]+".components";
	       // check that only one is defined
	       ASKAPCHECK(parset.isDefined(compPar) != parset.isDefined(modelPar),
	            "The model should be defined with either image (via "<<modelPar<<") or components (via "<<
	             compPar<<"), not both");
	       // 
           if (parset.isDefined(modelPar)) {
               const std::vector<std::string> vecModels = parset.getStringVector(modelPar);
               const int nTaylorTerms = parset.getInt32(std::string("sources.")+sources[i]+".nterms",1);                                                      
               ASKAPCHECK(nTaylorTerms>0, "Number of Taylor terms is supposed to be a positive number, you gave "<<
                         nTaylorTerms);
               if (nTaylorTerms>1) {
                   ASKAPLOG_INFO_STR(logger,"Simulation from model presented by Taylor series (a.k.a. MFS-model) with "<<
                               nTaylorTerms<<" terms");
               }              
               ASKAPCHECK((vecModels.size() == 1) || (int(vecModels.size()) == nTaylorTerms), 
                    "Number of model images given by "<<modelPar<<" should be either 1 or one per taylor term, you gave "<<
                    vecModels.size()<<" nTaylorTerms="<<nTaylorTerms);
               ImageParamsHelper iph("image."+sources[i]);
               // for simulations we don't need cross-terms 
               for (int order = 0; order<nTaylorTerms; ++order) {
                    if (nTaylorTerms > 1) {
                        // this is an MFS case, setup Taylor terms
                        iph.makeTaylorTerm(order);
                        ASKAPLOG_INFO_STR(logger,"Processing Taylor term "<<order);
                    }
                    std::string model = substitute(vecModels[vecModels.size() == 1 ? 0 : order]);
                    if (vecModels.size() == 1) {
                        // only base name is given, need to add taylor suffix
                        model += iph.suffix();
                    }
                                        
                    if (std::find(loadedImageModels.begin(),loadedImageModels.end(),model) != loadedImageModels.end()) {
                        ASKAPLOG_INFO_STR(logger, "Model " << model << " has already been loaded, reusing it for "<< sources[i]);
                        if (vecModels.size()!=1) {
                            ASKAPLOG_WARN_STR(logger, "MFS simulation will not work correctly if you specified the same model "<<
                                 model<<" for multiple Taylor terms");
                        }
                    } else {
                        ASKAPLOG_INFO_STR(logger, "Adding image " << model << " as model for "<< sources[i]
                                           << ", parameter name: "<<iph.paramName() );
                        // need to patch model to append taylor suffix
                        SynthesisParamsHelper::loadImageParameter(*pModel, iph.paramName(), model);
                        loadedImageModels.insert(model);
                    }
               }
           } else {
               // loop through components
               ASKAPLOG_INFO_STR(logger, "Adding components as model for "<< sources[i] );
               const vector<string> compList = parset.getStringVector(compPar);
               for (vector<string>::const_iterator cmp = compList.begin(); cmp != compList.end(); ++cmp) {
                    ASKAPLOG_INFO_STR(logger, "Loading component " << *cmp << " as part of the model for " << sources[i]);
                    SynthesisParamsHelper::copyComponent(pModel, parset,*cmp,"sources.");
                }
           }
      }
      ASKAPLOG_INFO_STR(logger, "Successfully read models");      
    }
  }
}
