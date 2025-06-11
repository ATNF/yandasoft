/// @file tMajorCycleComms.cc
///
/// @breif A tool to simulate major cycle communication pattern for debugging / comms tests
/// @details
/// It fakes the model and normal equations (with the desired size) and performs model broadcast
/// followed by the tree reduction of normal equations. It is possible to set delays specific to
/// master and worker as well as override them for a chosen rank or multiple ranks.
///
/// Control parameters are passed in from a LOFAR ParameterSet file in a standard way for all yandasoft apps.
/// This tool doesn't write or read anything but does understand some generic parameters because it reuses
/// low level imager code.
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

// Package level header file
#include <askap/askap_synthesis.h>

// System includes
#include <stdexcept>
#include <iostream>
#include <csignal>

// boost includes
#include <boost/scoped_ptr.hpp>

// ASKAPsoft includes
#include <askap/askap/Application.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/StatReporter.h>
#include <askap/measurementequation/MEParsetInterface.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/profile/AskapProfiler.h>
#include <askap/parallel/MEParallel.h>
#include <askap/scimath/fitting/ImagingNormalEquations.h>


ASKAP_LOGGER(logger, ".tMajorCycleComms");

using namespace askap;
using namespace askap::synthesis;
using namespace askap::scimath;

/// @brief class pretending to be imager but doing nothing except initialising the structures
/// @details It is derived from MEParallel rather than MEParallelApp to avoid the need to define
/// gridders and other stuff.
class FakeImager : public askap::synthesis::MEParallel {
public:
    FakeImager(askap::askapparallel::AskapParallel& comms, const LOFAR::ParameterSet& parset) :
      MEParallel(comms,parset,true), itsThisRankDelay(0)
    {
      if (itsComms.isMaster()) {
          ASKAPLOG_INFO_STR(logger, "Initializing the model images");

          /// Create the specified images from the definition in the
          /// parameter set. We can solve for any number of images
          /// at once (but you may/will run out of memory!)
          SynthesisParamsHelper::setUpImages(itsModel, parset.makeSubset("Images."));
      }

      const unsigned int defaultDelay = parset.getUint32("delay.default",10);
      const std::string delayParamStr = "delay.rank"+utility::toString(comms.rank());
      itsThisRankDelay = parset.isDefined(delayParamStr) ? parset.getUint32(delayParamStr) : defaultDelay;
    }

    /// @brief stub for model writer
    virtual void writeModel(const std::string &) {}

    /// @brief stub for calculation of normal equations
    virtual void calcNE() {
       ASKAPTRACE("FakeImager::calcNE");
       boost::shared_ptr<ImagingNormalEquations> ne(new ImagingNormalEquations(*itsModel));
       itsNe = ne;

       if (itsComms.isWorker()) {
          // fill fake NE here
          const std::vector<std::string> completions(itsModel->completions("image"));
          for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();++it) {
               const string imageName("image"+(*it));
               const casacore::IPosition imageShape(itsModel->shape(imageName));
               ASKAPLOG_INFO_STR(logger,"Name: " << imageName << " Shape: " << imageShape);
               casacore::Array<imtype> imageDeriv(imageShape);
               casacore::Array<imtype> imagePrecon(imageShape);
               casacore::Array<imtype> imagePSF(imageShape);
               casacore::Array<imtype> imageWeight(imageShape);
               imageWeight.set(static_cast<imtype>(1.));
               imagePSF.set(static_cast<imtype>(0.));
               imageDeriv.set(static_cast<imtype>(-1.));
               imagePrecon.set(static_cast<imtype>(-1.));

               {
                 casacore::IPosition vecShape(1, imagePSF.nelements());
                 casacore::IPosition reference(4, imageShape(0)/2, imageShape(1)/2, 0, 0);
                 casacore::Vector<imtype> imagePSFVec(imagePSF.reform(vecShape));
                 casacore::Vector<imtype> imageWeightVec(imageWeight.reform(vecShape));
                 casacore::Vector<imtype> imageDerivVec(imageDeriv.reform(vecShape));
                 casacore::Vector<imtype> imagePreconVec(imagePrecon.reform(vecShape));
                 ne->addSlice(imageName, imagePSFVec, imageWeightVec, imagePreconVec,
                        imageDerivVec, imageShape, reference,SynthesisParamsHelper::coordinateSystem(*itsModel,imageName));
               }

          }
          if (itsThisRankDelay > 0u) {
              ASKAPLOG_INFO_STR(logger, "Pausing for "<<itsThisRankDelay<<" seconds");
              sleep(itsThisRankDelay);
          }
          if (itsComms.isParallel()) {
              sendNE();
          }
       }
    }

    /// @brief stub for solving normal equations
    virtual void solveNE() {
       ASKAPTRACE("FakeImager::solveNE");
       if (itsComms.isMaster()) {
           // Receive the normal equations
           if (itsComms.isParallel()) {
               receiveNE();
           }
           ASKAPLOG_INFO_STR(logger, "Received normal equations");
       }
       const std::vector<std::string> completions(itsModel->completions("image"));
       for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();++it) {
            casacore::Array<imtype>& model = itsModel->valueT("image"+*it);
            casacore::IPosition shape = model.shape();
            if (shape.product() > 0u) {
                for (int dim = 0; dim < shape.nelements(); ++dim) {
                     shape[dim] /= 2;
                }
                model(shape) += static_cast<imtype>(1.);
            }
       }
       // in this test we don't use peak_residual, but update it anyway as it is sent around in the comms pattern
       if (itsModel->has("peak_residual")) {
           const double peak = itsModel->scalarValue("peak_residual") + 1.;
           itsModel->update("peak_residual",peak);
       } else {
           itsModel->add("peak_residual",1.);
       }
       itsModel->fix("peak_residual");
       if (itsThisRankDelay > 0u) {
           ASKAPLOG_INFO_STR(logger, "Pausing for "<<itsThisRankDelay<<" seconds");
           sleep(itsThisRankDelay);
       }
    }

protected:
    /// @brief helper method to indentify model parameters to broadcast
    /// @details We use itsModel to buffer some derived images like psf, weights, etc
    /// which are not required for prediffers. It just wastes memory and CPU time if
    /// we broadcast them. At the same time, some auxilliary parameters like peak
    /// residual value need to be broadcast (so the major cycle can terminate in workers).
    /// This method returns the vector with all parameters to be broadcast. By default
    /// it returns all parameter names, so it is overridden here to broadcast only
    /// model images and the peak_residual metadata.
    /// @return a vector with parameters to broadcast
    virtual std::vector<std::string> parametersToBroadcast() const
    {
       ASKAPDEBUGASSERT(itsModel);
       const std::vector<std::string> names = itsModel->names();
       std::vector<std::string> result;
       result.reserve(names.size());
       for (std::vector<std::string>::const_iterator ci=names.begin(); ci!=names.end(); ++ci) {
            if ((ci->find("image") == 0) || (ci->find("peak_residual") == 0)) {
                result.push_back(*ci);
            }
       }
       return result;
    }

private:
    /// @brief optional delay
    unsigned int itsThisRankDelay;
};

/// @brief application class, does the cycles
class MajorCycleCommsTestApp : public askap::Application
{
    public:
        int run(int argc, char* argv[]) final
        {
            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;
                // in principle, we can remove Cimager prefix, but having it simplifies copying
                // from real imager parsets
                LOFAR::ParameterSet subset(config().makeSubset("Cimager."));

                boost::scoped_ptr<askap::ProfileSingleton::Initialiser> profiler;
                if (parameterExists("profile")) {
                    std::string profileFileName("profile.cimager");
                    if (subset.isDefined("Images.Names")){
                        profileFileName += "."+subset.getStringVector("Images.Names")[0];
                    }
                    if (comms.isParallel()) {
                        profileFileName += ".rank"+utility::toString(comms.rank());
                    }
                    profiler.reset(new askap::ProfileSingleton::Initialiser(profileFileName));
                }

                // Put everything in scope to ensure that all destructors are called
                // before the final message
                {
                    // imager-specific configuration of the master/worker to allow groups of workers
                    const int nWorkerGroups = subset.getInt32("nworkergroups", 1);
                    ASKAPCHECK(nWorkerGroups > 0, "nworkergroups is supposed to be greater than 0");
                    if (nWorkerGroups > 1) {
                        ASKAPLOG_INFO_STR(logger, "Model parameters will be distributed between "<<nWorkerGroups<<
                                " groups of workers");
                        ASKAPCHECK(comms.isParallel(), "This option is only allowed in the parallel mode");
                        comms.defineGroups(nWorkerGroups);
                    } else {
                        ASKAPLOG_INFO_STR(logger, "All workers are treated as identical");
                    }

                    // Perform %w substitutions for all keys.
                    // NOTE: This MUST happen after AskapParallel::defineGroups() is called
                    for (LOFAR::ParameterSet::iterator it = subset.begin();
                            it != subset.end(); ++it) {
                        it->second = LOFAR::ParameterValue(comms.substitute(it->second));
                    }

                    FakeImager tester(comms, subset);
                    ASKAPLOG_INFO_STR(logger, "ASKAP tMajorCycleComms test utility " << ASKAP_PACKAGE_VERSION);

                    const int nCycles = subset.getInt32("ncycles", 0);
                    if (nCycles == 0) {
                        /// No cycling - just make a dirty image
                        tester.broadcastModel();
                        tester.receiveModel();
                        tester.calcNE();
                        tester.receiveNE();
                    } else {

                        // Distribute initial model
                        tester.broadcastModel();
                        tester.receiveModel();

                        /// Perform multiple major cycles
                        for (int cycle = 0; cycle < nCycles; ++cycle) {

                            ASKAPLOG_INFO_STR(logger, "*** Starting major cycle " << cycle << " ***");
                            tester.calcNE();
                            tester.solveNE();

                            stats.logSummary();

                            if (cycle + 1 >= nCycles) {
                                ASKAPLOG_INFO_STR(logger, "Reached " << nCycles <<
                                    " cycle(s), the maximum number of major cycles. Stopping.");
                            }

                            // Distribute current model
                            tester.broadcastModel();
                            tester.receiveModel();
                        }

                        ASKAPLOG_INFO_STR(logger, "*** Finished major cycles ***");
                        tester.calcNE();
                        tester.receiveNE();
                    }
                }
                stats.logSummary();
            } catch (const askap::AskapError& x) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << x.what());
                std::cerr << "Askap error in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            } catch (const std::exception& x) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << x.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << x.what() << std::endl;
                exit(1);
            }

            return 0;
        }

    private:
        std::string getVersion() const final {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

int main(int argc, char *argv[])
{
    MajorCycleCommsTestApp app;
    app.addParameter("profile", "p", "Write profiling output files", false);
    return app.main(argc, argv);
}
