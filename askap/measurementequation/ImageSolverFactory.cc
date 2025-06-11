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

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.imagesolverfactory");

#include <askap/askap/AskapError.h>

#include <askap/measurementequation/ImageSolverFactory.h>

#include <askap/measurementequation/ImageSolver.h>
#include <askap/measurementequation/ImageMultiScaleSolver.h>
#include <askap/measurementequation/ImageBasisFunctionSolver.h>
#include <askap/measurementequation/ImageFistaSolver.h>
#include <askap/measurementequation/ImageMSMFSolver.h>
#include <askap/measurementequation/ImageAMSMFSolver.h>
#include <askap/measurementequation/IImagePreconditioner.h>
#include <askap/measurementequation/WienerPreconditioner.h>
#include <askap/measurementequation/NormWienerPreconditioner.h>
#include <askap/measurementequation/RobustPreconditioner.h>
#include <askap/measurementequation/GaussianTaperPreconditioner.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>

using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    ImageSolverFactory::ImageSolverFactory()
    {
    }

    ImageSolverFactory::~ImageSolverFactory()
    {
    }

    /// @brief Helper method to configure minor cycle threshold(s)
    /// @details This method parses threshold.minorcycle parameter
    /// of the parset file. The parameter can be either a single string
    /// or a vector of two strings. A number without units is interpreted
    /// as a fractional stopping threshold (w.r.t the peak residual),
    /// so does the number with the percentage sign. An absolute flux given
    /// in Jy or related units is interpreted as an absolute threshold.
    /// Either one or both of these thresholds can be given in the same
    /// time.
    /// Optionally a second absolute threshold can be given for some Cleaners
    /// that support a 'deep' cleaning mode using the existing model as a mask
    /// @param[in] parset parameter set to extract the input from
    /// @param[in] solver shared pointer to the solver to be configured
    void ImageSolverFactory::configureThresholds(const LOFAR::ParameterSet &parset,
                                                 const boost::shared_ptr<ImageSolver> &solver)
    {
      ASKAPDEBUGASSERT(solver);
      const std::string parName = "threshold.minorcycle";
      if (parset.isDefined(parName)) {
        const std::vector<std::string> thresholds = parset.getStringVector(parName);
        ASKAPCHECK(thresholds.size() && (thresholds.size()<4), "Parameter "<<parName<<
                   " must contain either 1 element or a vector of 2-3 elements, you have "<<
                   thresholds.size()<<" elements");
        bool absoluteThresholdDefined = false;
        bool relativeThresholdDefined = false;
        bool absoluteThreshold2Defined = false;
        bool noiseThresholdDefined = false;
        bool noiseThreshold2Defined = false;
        for (const string& t : thresholds) {
          casacore::Quantity cThreshold;
          // check for noise thresholds
          const string sigma("sigma");
          const size_t pos = t.rfind(sigma);
          // do we want spatially variant sigma thresholds?
           const uInt boxSize = parset.getUint("solver.Clean.noiseboxsize",0);

          if (pos != std::string::npos && pos == t.size() - sigma.size()) {
              ASKAPCHECK(!noiseThreshold2Defined, "Parameter "<<parName<<
                         " defines noise threshold thrice ("<<t<<")");
              ASKAPCHECK(!absoluteThresholdDefined,"Cannot mix "<<
              "absolute and noise thresholds");
              casacore::Quantity::read(cThreshold, t.substr(0,pos));
              if (noiseThresholdDefined) {
                  noiseThreshold2Defined = true;
                  boost::shared_ptr<ImageCleaningSolver> ics =
                      boost::dynamic_pointer_cast<ImageAMSMFSolver>(solver);
                  if (ics) {
                    ics->setDeepNoiseThreshold(cThreshold.getValue());
                    ASKAPLOG_INFO_STR(logger, "Will stop deep minor cycle at the noise threshold of "<<
                                      cThreshold.getValue("")<<" sigma");
                  } else {
                    ASKAPLOG_INFO_STR(logger, "The type of the image solver used does not allow to specify "
                                      "a deep noise threshold, ignoring "<<t<<" in "<<parName);
                  }
              } else {
                  noiseThresholdDefined = true;
                  boost::shared_ptr<ImageCleaningSolver> ics =
                      boost::dynamic_pointer_cast<ImageAMSMFSolver>(solver);
                  if (ics) {
                    ics->setNoiseThreshold(cThreshold.getValue());
                    ASKAPLOG_INFO_STR(logger, "Will stop minor cycle at the noise threshold of "<<
                                      cThreshold.getValue("")<<" sigma");
                    if (boxSize > 0) {
                      ics->setNoiseBoxSize(boxSize);
                      ASKAPLOG_INFO_STR(logger, "Will use a spatially variant noise threshold with box size of "<<
                                      boxSize<<" pixels");
                    }
                  } else {
                    ASKAPLOG_INFO_STR(logger, "The type of the image solver used does not allow to specify "
                                      "a noise threshold, ignoring "<<t<<" in "<<parName);
                  }
              }
          } else {
              casacore::Quantity::read(cThreshold, t);
              cThreshold.convert();
              if (cThreshold.isConform("Jy")) {
                ASKAPCHECK(!absoluteThreshold2Defined, "Parameter "<<parName<<
                           " defines absolute threshold thrice ("<<t<<")");
                ASKAPCHECK(!noiseThresholdDefined,"Cannot mix "<<
                "absolute and noise thresholds");

                if (absoluteThresholdDefined) {
                  absoluteThreshold2Defined = true;
                  boost::shared_ptr<ImageCleaningSolver> ics =
                      boost::dynamic_pointer_cast<ImageAMSMFSolver>(solver);
                  if (ics) {
                    ics->setDeepThreshold(cThreshold.getValue("Jy"));
                    ASKAPLOG_INFO_STR(logger, "Will stop deep minor cycle at the absolute threshold of "<<
                                      cThreshold.getValue("mJy")<<" mJy");
                  } else {
                    ASKAPLOG_INFO_STR(logger, "The type of the image solver used does not allow to specify "
                                      "a deep threshold, ignoring "<<t<<" in "<<parName);
                  }
                } else {
                  absoluteThresholdDefined = true;
                  solver->setThreshold(cThreshold);
                  ASKAPLOG_INFO_STR(logger, "Will stop the minor cycle at the absolute threshold of "<<
                                    cThreshold.getValue("mJy")<<" mJy");
                }
              } else if (cThreshold.isConform("")) {
                ASKAPCHECK(!relativeThresholdDefined, "Parameter "<<parName<<
                           " defines relative threshold twice ("<<t<<")");
                relativeThresholdDefined = true;

                boost::shared_ptr<ImageCleaningSolver> ics =
                  boost::dynamic_pointer_cast<ImageCleaningSolver>(solver);
                if (ics) {
                  ics->setFractionalThreshold(cThreshold.getValue());
                  ASKAPLOG_INFO_STR(logger, "Will stop minor cycle at the relative threshold of "<<
                                    cThreshold.getValue()*100.<<"\%");
                } else {
                  ASKAPLOG_INFO_STR(logger, "The type of the image solver used does not allow to specify "
                                    "a fractional threshold, ignoring "<<t<<" in "<<parName);
                }
              } else {
                ASKAPTHROW(AskapError, "Unable to convert units in the quantity "<<
                           cThreshold<<" to either Jy or a dimensionless quantity");
              } // if - isConform
          } // if - noise threshold
        } // for - parameter loop
      } // if - parameter defined

      const std::string parName2 = "threshold.masking";
      if (parset.isDefined(parName2)) {
        boost::shared_ptr<ImageCleaningSolver> ics =
          boost::dynamic_pointer_cast<ImageCleaningSolver>(solver);
        if (ics) {
          ics->setMaskingThreshold(parset.getFloat(parName2, -1));
        } else {
          ASKAPLOG_INFO_STR(logger, "The type of the image solver used does not allow to specify "
                            "masking threshold, ignoring "<<parName2);
        }
      } // if - parameter defined

      const std::string parName3 = "threshold.firstimage";
      boost::shared_ptr<ImageCleaningSolver> ics =
            boost::dynamic_pointer_cast<ImageAMSMFSolver>(solver);
      if (ics) {
        ASKAPLOG_INFO_STR(logger, "Using the first image to set thresholds for the other images");
        ics->setUseFirstImageForThresholds(parset.getBool(parName3, true));
      } else {
        ASKAPLOG_INFO_STR(logger, "The type of the image solver used does not allow to specify "
            "to use thresholds from first image, ignoring "<<parName3);
      }
    } // method


    /// @brief Helper method to configure preconditioners
    /// @details This method parses reconditioner.Names and sets up
    /// given preconditioners. Parameters necessary for this setup are also
    /// extracted from parset
    /// @param[in] parset parameter set to extract the input from
    /// @param[in] solver shared pointer to the solver to be configured
    /// @note This method is made public (temporary) because the same
    /// functionality is required to setup the restore solver. In future, we
    /// expect to set up the restore solver through the same factory
    void ImageSolverFactory::configurePreconditioners(const LOFAR::ParameterSet &parset,
                                                      const boost::shared_ptr<ImageSolver> &solver)
    {
      ASKAPDEBUGASSERT(solver);
      // Set Up the Preconditioners - a whole list of 'em
      const vector<string> preconditioners=parset.getStringVector("preconditioner.Names",std::vector<std::string>());
      for (vector<string>::const_iterator pc = preconditioners.begin(); pc != preconditioners.end(); ++pc) {
        if ( (*pc)=="Wiener" ) {
          // Use the cached preconditioner function (if it exists & size matches).
          solver->addPreconditioner(WienerPreconditioner::createPreconditioner(
              parset.makeSubset("preconditioner.Wiener."), true));
        } else if ( (*pc)=="NormWiener" ) {
          const float robustness = parset.getFloat("preconditioner.NormWiener.robustness",0.0);
          solver->addPreconditioner(IImagePreconditioner::ShPtr(new NormWienerPreconditioner(robustness)));
        } else if ( (*pc)=="Robust" ) {
          const float robustness = parset.getFloat("preconditioner.Robust.robustness",0.0);
          solver->addPreconditioner(IImagePreconditioner::ShPtr(new RobustPreconditioner(robustness)));
        } else if ( (*pc) == "GaussianTaper") {
          // at this stage we have to define tapers in uv-cells, rather than in klambda
          // because the physical cell size is unknown to solver factory.
          // Theoretically we could parse the parameters here and extract the cell size and
          // shape, but it can be defined separately for each image. We need to find
          // the way of dealing with this complication.
          ASKAPCHECK(parset.isDefined("preconditioner.GaussianTaper"),
                     "preconditioner.GaussianTaper showing the taper size should be defined to use GaussianTaper");
          const vector<double> taper =
              SynthesisParamsHelper::convertQuantity(parset.getStringVector("preconditioner.GaussianTaper"),"rad");
          ASKAPCHECK((taper.size() == 3) || (taper.size() == 1),
                     "preconditioner.GaussianTaper can have either single element or "
                     " a vector of 3 elements. You supplied a vector of "<<taper.size()<<" elements");
          ASKAPCHECK(parset.isDefined("Images.shape") && parset.isDefined("Images.cellsize"),
                     "Images.shape and Images.cellsize should be defined to convert the taper fwhm specified in "
                     "angular units in the image plane into uv cells");
          const std::vector<double> cellsize =
              SynthesisParamsHelper::convertQuantity(parset.getStringVector("Images.cellsize"),"rad");
          const std::vector<int> shape = parset.getInt32Vector("Images.shape");
          ASKAPCHECK((cellsize.size() == 2) && (shape.size() == 2),
                     "Images.cellsize and Images.shape parameters should have exactly two values");
          bool isPsfSize = parset.getBool("preconditioner.GaussianTaper.isPsfSize",False);
          double tol = parset.getDouble("preconditioner.GaussianTaper.tolerance",0.005);
          // Try to use the same cutoff & support as used by the restore solver so final beams will match
          double cutoff = parset.getDouble("restore.beam.cutoff",0.5);
          int support = parset.getInt("restore.beam.maxsupport",101);

          // additional scaling factor due to padding. by default - no padding
          const boost::shared_ptr<ImageCleaningSolver> ics = boost::dynamic_pointer_cast<ImageCleaningSolver>(solver);
          const double paddingFactor = ics ? ics->paddingFactor() : 1.;

          // factors which appear in nominator are effectively half sizes in radians
          const double xFactor = 4. * log(2.) * cellsize[0]*double(shape[0])*paddingFactor / casacore::C::pi;
          const double yFactor = 4. * log(2.) * cellsize[1]*double(shape[1])*paddingFactor / casacore::C::pi;

          if (taper.size() == 3) {

            ASKAPDEBUGASSERT((taper[0]!=0) && (taper[1]!=0));
            solver->addPreconditioner(IImagePreconditioner::ShPtr(
            new GaussianTaperPreconditioner(xFactor/taper[0],yFactor/taper[1],taper[2],
                                            isPsfSize,cutoff,support,tol)));
          } else {
            ASKAPDEBUGASSERT(taper[0]!=0);
            if (std::abs(xFactor-yFactor)<4e-15) {
              // the image is square, can use the short cut
              solver->addPreconditioner(IImagePreconditioner::ShPtr(
              new GaussianTaperPreconditioner(xFactor/taper[0],isPsfSize,cutoff,support,tol)));
            } else {
              // the image is rectangular. Although the gaussian taper is symmetric in
              // angular coordinates, it will be elongated along the vertical axis in
              // the uv-coordinates.
              solver->addPreconditioner(IImagePreconditioner::ShPtr(
              new GaussianTaperPreconditioner(xFactor/taper[0],yFactor/taper[0],0.,isPsfSize,cutoff,support,tol)));
            } // xFactor!=yFactor
          } // else: taper.size() == 3

        //} else if ( (*pc)=="ApproxPsf" ) { //later, add the option of specifying a beam, or fitting for it.
        //  solver->addPreconditioner(IImagePreconditioner::ShPtr(new ApproxPsfPreconditioner()));
        } else if ( (*pc)!="None") {
          ASKAPTHROW(AskapError, "Unknown preconditioner "<<*pc);
        }
      } // loop over all preconditioners
    }  // method

    Solver::ShPtr ImageSolverFactory::make(const LOFAR::ParameterSet &parset) {
      ImageSolver::ShPtr solver;
      if(parset.getString("solver")=="Clean") {
        std::vector<float> defaultScales(3);
        defaultScales[0]=0.0;
        defaultScales[1]=10.0;
        defaultScales[2]=30.0;

        // temporary
        ASKAPCHECK(!parset.isDefined("solver.Clean.threshold"),
                   "The use of the parameter solver.Clean.threshold is deprecated, use threshold.minorcycle instead");

        string algorithm=parset.getString("solver.Clean.algorithm","MultiScale");
        // MV: a bit of technical debt highlighted by casacore's interface change. In principle, we could've
        // had scales as std::vector in the interface to avoid the explicit construction (in this particular case,
        // there is no benefit of using casacore::Vector) + we have a similar code in a number of places, could've achieved a better reuse
        const casacore::Vector<float> scales = casacore::Vector<float>(parset.getFloatVector("solver.Clean.scales", defaultScales));

        if (algorithm=="MultiScale") {
          ASKAPLOG_INFO_STR(logger, "Constructing MultiScale Clean solver (ASKAP version)" );
          solver.reset(new ImageMultiScaleSolver(scales));
          solver->setAlgorithm(algorithm);
        }
        else if(algorithm=="Basisfunction") {
         ASKAPLOG_INFO_STR(logger, "Constructing Basisfunction Clean solver ");
         solver.reset(new ImageBasisFunctionSolver(scales));
        }
        else if(algorithm=="Hogbom") {
          ASKAPLOG_INFO_STR(logger, "Constructing Hogbom Clean solver");
          solver.reset(new ImageMultiScaleSolver());
        }
        else if ((algorithm=="MSMFS")||(algorithm=="MultiScaleMFS")) {
          ASKAPCHECK(!parset.isDefined("solver.nterms"), "Specify nterms for each image instead of using solver.nterms");
          ASKAPCHECK(!parset.isDefined("solver.Clean.nterms"),
              "Specify nterms for each image instead of using solver.Clean.nterms");
          solver.reset(new ImageMSMFSolver(casacore::Vector<float>(scales)));
          ASKAPLOG_INFO_STR(logger, "Constructed image multiscale multi-frequency solver (CASA version)" );
        }
        else if ((algorithm=="BFMFS")||(algorithm=="BasisfunctionMFS")) {
          ASKAPCHECK(!parset.isDefined("solver.nterms"), "Specify nterms for each image instead of using solver.nterms");
          ASKAPCHECK(!parset.isDefined("solver.Clean.nterms"),
              "Specify nterms for each image instead of using solver.Clean.nterms");
          solver.reset(new ImageAMSMFSolver());
          ASKAPLOG_INFO_STR(logger, "Constructed basis function multi-frequency solver" );
        }
        else {
          ASKAPTHROW(AskapError, "Unknown Clean algorithm " << algorithm);
        }
        solver->setAlgorithm(algorithm);
        solver->configure(parset.makeSubset("solver.Clean."));
        solver->setGain(parset.getFloat("solver.Clean.gain", 0.7));
        solver->setNiter(parset.getInt32("solver.Clean.niter", 100));

        if (parset.isDefined("solver.Clean.speedup")) {
          boost::shared_ptr<ImageMultiScaleSolver> mss = boost::dynamic_pointer_cast<ImageMultiScaleSolver>(solver);
          boost::shared_ptr<ImageMSMFSolver> mfmss = boost::dynamic_pointer_cast<ImageMSMFSolver>(solver);
          ASKAPCHECK(mss || mfmss, "speed up can currently be used with multi scale solvers only");
          const float factor = parset.getFloat("solver.Clean.speedup");
          ASKAPLOG_INFO_STR(logger,"Using speed up factor of "<<factor<<" in lattice clean");
          if (mss) {
            mss->setSpeedUp(factor);
          }
          if (mfmss) {
            mfmss->setSpeedUp(factor);
          }
        }

        if (parset.isDefined("Images.extraoversampling")) {
          const boost::shared_ptr<ImageMultiScaleSolver>
              mss = boost::dynamic_pointer_cast<ImageMultiScaleSolver>(solver);
          const boost::shared_ptr<ImageAMSMFSolver>
              amsmfs = boost::dynamic_pointer_cast<ImageAMSMFSolver>(solver);
          ASKAPCHECK(mss || amsmfs, "Nyquist gridding only supported with Hogbom or BasisfunctionMFS solvers");
          const float factor = parset.getFloat("Images.extraoversampling");
          ASKAPCHECK(factor>1., "extra oversampling factor should be greater than 1");
          ASKAPLOG_INFO_STR(logger,"Using an extra oversampling factor of "<<factor<<" in lattice clean");
          if (mss) {
            mss->setExtraOversampling(factor);
          }
          if (amsmfs) {
            amsmfs->setExtraOversampling(factor);
          }
        }

      } else if(parset.getString("solver")=="Fista") {
        ASKAPLOG_INFO_STR(logger, "Constructing FISTA solver" );
        solver.reset(new ImageFistaSolver());
        solver->configure(parset.makeSubset("solver.Fista."));
      } else {
        ASKAPLOG_INFO_STR(logger, "Constructing dirty image solver" );
        solver = ImageSolver::ShPtr(new ImageSolver);
        solver->configure(parset.makeSubset("solver.Dirty."));
      }
      configureThresholds(parset, solver);
      configurePreconditioners(parset, solver);
      return solver;
    }
  }
}
