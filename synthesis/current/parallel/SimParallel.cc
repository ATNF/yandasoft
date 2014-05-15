/// @file SimParallel.cc
///
/// @brief Class for parallel simulation using CASA NewMSSimulator
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
#include <parallel/SimParallel.h>

// System includes
#include <iostream>
#include <sstream>

// ASKAPsoft includes
#include <casa/OS/Timer.h>
#include <askap_synthesis.h>
#include <askap/AskapError.h>
#include <askap/AskapLogging.h>
#include <askap/AskapUtil.h>
#include <askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
ASKAP_LOGGER(logger, ".parallel");

// Local package includes
#include <dataaccess/DataAccessError.h>
#include <dataaccess/TableDataSource.h>
#include <dataaccess/ParsetInterface.h>
#include <measurementequation/ImageFFTEquation.h>
#include <measurementequation/SynthesisParamsHelper.h>
#include <measurementequation/ImageRestoreSolver.h>
#include <measurementequation/CalibrationME.h>
#include <measurementequation/NoXPolGain.h>
#include <measurementequation/Product.h>
#include <measurementequation/LeakageTerm.h>
#include <measurementequation/ImagingEquationAdapter.h>
#include <measurementequation/SumOfTwoMEs.h>
#include <measurementequation/GaussianNoiseME.h>
#include <measurementequation/ComponentEquation.h>
#include <measurementequation/ImageSolverFactory.h>
#include <calibaccess/CalibAccessFactory.h>
#include <measurementequation/CalibParamsMEAdapter.h>
#include <gridding/VisGridderFactory.h>
#include <parallel/ParallelWriteIterator.h>

using namespace std;
using namespace askap;
using namespace askap::askapparallel;
using namespace LOFAR;
using namespace askap::accessors;

namespace askap
{
namespace synthesis
{

SimParallel::SimParallel(askap::askapparallel::AskapParallel& comms,
                         const LOFAR::ParameterSet& parset) :
        SynParallel(comms,parset), itsModelReadByMaster(true), itsMSWrittenByMaster(false), itsNoiseVariance(-1.), 
        itsDoChecksForNoise(false)
{
  itsModelReadByMaster = parset.getBool("modelReadByMaster", true);
  itsMSWrittenByMaster = parset.getBool("msWrittenByMaster", false);
  ASKAPCHECK(getFreqRefFrame().getType() == casa::MFrequency::Ref(casa::MFrequency::TOPO).getType(), 
             "Only topocentric reference frame is currently understood by the simulator");
  if (itsMSWrittenByMaster) {
      ASKAPCHECK(comms.isParallel(), "msWrittenByMaster can only be used in the parallel case");
      ASKAPLOG_INFO_STR(logger, "Master will receive data from workers and write a single measurement set");      
  } else if (comms.isParallel()) {
      ASKAPLOG_INFO_STR(logger, "Each worker will write its own measurement set");        
  }
}

void SimParallel::init()
{
    // set up image handler
    SynthesisParamsHelper::setUpImageHandler(parset());
    if (itsComms.isMaster()) {
        if (itsModelReadByMaster) {
            readModels();
            broadcastModel();
        }
    }
    if (itsComms.isWorker()) {
        if (itsModelReadByMaster) {
            receiveModel();
        } else {
            readModels();
        }
    }

    // the following code can be executed both in master and in workers depending who is
    // writing the measurement set
    if (itsComms.isWorker() != itsMSWrittenByMaster) {
        string msname(substitute(parset().getString("dataset",
                                 "test%w.ms")));
        int bucketSize  = parset().getInt32("stman.bucketsize", 32768);
        int tileNcorr = parset().getInt32("stman.tilencorr", 4);
        int tileNchan = parset().getInt32("stman.tilenchan", 32);
        itsSim.reset(new Simulator(msname, bucketSize, tileNcorr, tileNchan));

        itsMs.reset(new casa::MeasurementSet(msname, casa::Table::Update));

        // The antenna info is kept in a separate parset file
        readAntennas();

        // Get the source definitions and get the model from the master
        readSources();

        // Get the feed definitions
        readFeeds();

        // Get the spectral window definitions. Not all of these need to be
        // used.
        readSpws();

        // Get miscellaneous information about the simulation
        readSimulation();

        // extract noise figure if needed
        itsDoChecksForNoise = false;
        if (parset().getBool("noise", false)) {
            itsNoiseVariance = getNoise(parset().makeSubset("noise."));
            ASKAPCHECK(itsNoiseVariance>0., 
               "Noise variance is supposed to be positive, you have "<<itsNoiseVariance);
            const double rms = sqrt(itsNoiseVariance);
            ASKAPLOG_INFO_STR(logger, 
               "SIGMA column will be scaled to account for simulated Gaussian noise (variance=" << 
                itsNoiseVariance << " Jy^2 or sigma="<<rms<<" Jy)");
            itsSim->setNoiseRMS(rms);
        }
    }
    // calibration is only dealt with in workers
    if (itsComms.isWorker()) {    
        // initialise calibration solution source, if the visibilities are going to be corrupted
        if (parset().getBool("corrupt", false)) {
            if (parset().isDefined("corrupt.gainsfile")) {
                ASKAPCHECK(!parset().isDefined("calibaccess") && !parset().isDefined("calibaccess.parset"),
                      "corrupt.gainsfile option conflicts with calibaccess, the former is deprecated - please correct your parset file");
                ASKAPLOG_WARN_STR(logger, "The parset has deprecated corrupt.gainsfile keyword, use calibaccess.parset instead");
                LOFAR::ParameterSet tmpParset(parset());
                tmpParset.add("calibaccess.parset", parset().getString("corrupt.gainsfile"));
                tmpParset.add("calibaccess", "parset");
                // setup solution source from the temporary parset
                itsSolutionSource = CalibAccessFactory::roCalSolutionSource(tmpParset);
            } else {
                itsSolutionSource = CalibAccessFactory::roCalSolutionSource(parset());
            }
            ASKAPASSERT(itsSolutionSource);                      
        } else {
            ASKAPLOG_INFO_STR(logger, "Calibration effects will not be simulated");
            itsSolutionSource.reset();
        }        
    }
}

SimParallel::~SimParallel()
{
    if ((itsComms.isWorker() != itsMSWrittenByMaster) && itsMs) {
        itsMs->flush();
    }
}

void SimParallel::readAntennas()
{
    ParameterSet parset(SynParallel::parset());

    if (SynParallel::parset().isDefined("antennas.definition")) {
        parset = ParameterSet(substitute(SynParallel::parset().getString("antennas.definition")));
    }

    /// Csimulator.name = ASKAP
    const std::string telName = parset.getString("antennas.telescope");
    ASKAPLOG_INFO_STR(logger, "Simulating " << telName);
    ostringstream oos;
    oos << "antennas." << telName << ".";
    ParameterSet antParset(parset.makeSubset(oos.str()));

    /// Csimulator.ASKAP.number=45
    ASKAPCHECK(antParset.isDefined("names"), "Subset (antennas."<<telName<<") of the antenna definition parset does not have 'names' keyword.");
    vector<string> antNames(antParset.getStringVector("names"));
    int nAnt = antNames.size();
    ASKAPCHECK(nAnt > 0, "No antennas defined in parset file");

    /// Csimulator.ASKAP.mount=equatorial
    string mount = antParset.getString("mount", "equatorial");
    ASKAPCHECK((mount == "equatorial") || (mount == "alt-az"), "Antenna mount unknown");

    /// Csimulator.ASKAP.mount=equatorial
    double diameter = asQuantity(antParset.getString("diameter", "12m")).getValue("m");
    ASKAPCHECK(diameter > 0.0, "Antenna diameter not positive");

    /// Csimulator.ASKAP.coordinates=local
    string coordinates = antParset.getString("coordinates", "local");
    ASKAPCHECK((coordinates == "local") || (coordinates == "global"), "Coordinates type unknown");

    /// Csimulator.ASKAP.scale=0.333
    const double scale = antParset.getDouble("scale", 1.0);

    /// Now we get the coordinates for each antenna in turn
    casa::Vector<double> x(nAnt);
    casa::Vector<double> y(nAnt);
    casa::Vector<double> z(nAnt);

    casa::Vector<double> dishDiameter(nAnt);

    casa::Vector<double> offset(nAnt);
    offset.set(0.0);
    casa::Vector<casa::String> mounts(nAnt);
    casa::Vector<casa::String> name(nAnt);

    /// Antenna information in the form:
    /// antennas.ASKAP.antenna0=[x,y,z]
    /// ...
    for (int iant = 0; iant < nAnt; iant++) {
        const vector<double> xyz = antParset.getDoubleVector(antNames[iant]);
        x[iant] = xyz[0] * scale;
        y[iant] = xyz[1] * scale;
        z[iant] = xyz[2] * scale;
        mounts[iant] = mount;
        dishDiameter[iant] = diameter;
        name[iant] = antNames[iant];
    }

    /// Csimulator.ASKAP.location=[+115deg, -26deg, 192km, WGS84]
    casa::MPosition location;
    if(coordinates == "local") {
        location = asMPosition(antParset.getStringVector("location"));
    }
    itsSim->initAnt(telName, x, y, z, dishDiameter, offset, mounts, name,
                    casa::String(coordinates), location);
    ASKAPLOG_INFO_STR(logger, "Successfully defined " << nAnt
            << " antennas of " << telName);
}

void SimParallel::readFeeds()
{
    ParameterSet parset(SynParallel::parset());

    if (SynParallel::parset().isDefined("feeds.definition")) {
        parset = ParameterSet(substitute(SynParallel::parset().getString("feeds.definition")));
    }

    vector<string> feedNames(parset.getStringVector("feeds.names"));
    int nFeeds = feedNames.size();
    ASKAPCHECK(nFeeds > 0, "No feeds specified");

    casa::Vector<double> x(nFeeds);
    casa::Vector<double> y(nFeeds);
    casa::Vector<casa::String> pol(nFeeds);

    casa::String mode = parset.getString("feeds.mode", "perfect X Y");

    for (int feed = 0; feed < nFeeds; feed++) {
        ostringstream os;
        os << "feeds." << feedNames[feed];
        string subFeed(substitute(os.str()));
        vector<double> xy(parset.getDoubleVector(subFeed));
        x[feed] = xy[0];
        y[feed] = xy[1];
        pol[feed] = "X Y";
    }

    if (parset.isDefined("feeds.spacing")) {
        casa::Quantity qspacing = asQuantity(parset.getString("feeds.spacing"));
        double spacing = qspacing.getValue("rad");
        ASKAPLOG_INFO_STR(logger, "Scaling feed specifications by " << qspacing);
        x *= spacing;
        y *= spacing;
    }

    itsSim->initFeeds(mode, x, y, pol);
    ASKAPLOG_INFO_STR(logger, "Successfully defined " << nFeeds << " feeds");
}

/// Csimulator.sources.names = [3C273, 1934-638]
/// Csimulator.sources.3C273.direction =
/// Csimulator.sources.1934-638.direction =
void SimParallel::readSources()
{
    ParameterSet parset(SynParallel::parset());

    if (SynParallel::parset().isDefined("sources.definition")) {
        parset = ParameterSet(substitute(SynParallel::parset().getString("sources.definition")));
    }

    const vector<string> sources = parset.getStringVector("sources.names");

    for (size_t i = 0; i < sources.size(); ++i) {
        {
            ostringstream oos;
            oos << "sources." << sources[i] << ".direction";
            ASKAPLOG_INFO_STR(logger, "Simulating source " << sources[i]);
            casa::MDirection direction(asMDirection(parset.getStringVector(oos.str())));
            itsSim->initFields(casa::String(sources[i]), direction, casa::String(""));
        }
    }

    ASKAPLOG_INFO_STR(logger, "Successfully defined sources");
}

void SimParallel::readSpws()
{
    ParameterSet parset(SynParallel::parset());

    if (SynParallel::parset().isDefined("spws.definition")) {
        parset = ParameterSet(substitute(SynParallel::parset().getString("spws.definition")));
    }

    vector<string> names(parset.getStringVector("spws.names"));
    const int nSpw = names.size();
    ASKAPCHECK(nSpw > 0, "No spectral windows defined");

    for (int spw = 0; spw < nSpw; spw++) {
        ostringstream os;
        os << "spws." << names[spw];
        vector<string> line = parset.getStringVector(os.str());
        ASKAPASSERT(line.size() >= 4);
        const casa::Quantity startFreq = asQuantity(line[1]);
        const casa::Quantity freqInc = asQuantity(line[2]);
        ASKAPCHECK(startFreq.isConform("Hz"), "start frequency for spectral window " << names[spw] << " is supposed to be in units convertible to Hz, you gave " <<
                   line[1]);
        ASKAPCHECK(freqInc.isConform("Hz"), "frequency increment for spectral window " << names[spw] << " is supposed to be in units convertible to Hz, you gave " <<
                   line[1]);
        itsSim->initSpWindows(names[spw], askap::utility::fromString<int>(line[0]),
                              startFreq, freqInc, freqInc, line[3]);
    }

    ASKAPLOG_INFO_STR(logger, "Successfully defined " << nSpw << " spectral windows");
}

void SimParallel::readSimulation()
{
    ParameterSet parset(SynParallel::parset());

    if (SynParallel::parset().isDefined("simulation.definition")) {
        parset = ParameterSet(substitute(SynParallel::parset().getString("simulation.definition")));
    }

    /// Csimulator.simulate.blockage=0.1
    itsSim->setFractionBlockageLimit(parset.getDouble("simulation.blockage", 0.0));
    /// Csimulator.simulate.elevationlimit=8deg
    itsSim->setElevationLimit(asQuantity(parset.getString(
                                  "simulation.elevationlimit", "8deg")));
    /// Csimulator.simulate.autocorrwt=0.0
    itsSim->setAutoCorrelationWt(parset.getFloat("simulation.autocorrwt", 0.0));

    /// Csimulator.simulate.integrationtime=10s
    casa::Quantity
    integrationTime(asQuantity(parset.getString(
                        "simulation.integrationtime", "10s")));
    /// Csimulator.simulate.usehourangles=true
    bool useHourAngles(parset.getBool("simulation.usehourangles", true));
    /// Csimulator.simulate.referencetime=2007Mar07
    vector<string> refTimeString(parset.getStringVector("simulation.referencetime"));
    casa::MEpoch refTime(asMEpoch(refTimeString));
    itsSim->settimes(integrationTime, useHourAngles, refTime);
    ASKAPLOG_INFO_STR(logger, "Successfully set simulation parameters");
}

/// Csimulator.observe.number=2
/// Csimulator.scan1=[1934-638, LBand1, 0s, 120s]
/// Csimulator.scan2=[3C273, LBand1, 120s, 1920s]
/// etc.
void SimParallel::simulate()
{

    if (itsComms.isWorker() != itsMSWrittenByMaster) {
        /// Now that the simulator is defined, we can observe each scan
        ParameterSet parset(SynParallel::parset());

        if (SynParallel::parset().isDefined("observe.definition")) {
            parset = ParameterSet(substitute(SynParallel::parset().getString("observe.definition")));
        }

        int nScans = parset.getInt32("observe.number", 0);
        ASKAPCHECK(nScans > 0, "No scans defined");

        for (int scan = 0; scan < nScans; scan++) {
            ostringstream oos;
            oos << "observe.scan" << scan;
            vector<string> line = parset.getStringVector(oos.str());
            ASKAPCHECK(line.size()>=4, "Expect at least 4-elements in observe.scanN, you have "<<line);
            string source = substitute(line[0]);
            string spw = substitute(line[1]);
            ASKAPLOG_INFO_STR(logger, "Observing scan " << scan << " on source " << source
                                  << " at band " << spw << " from "
                                  << asQuantity(line[2]) << " to "
                                  << asQuantity(line[3]));
            itsSim->observe(source, spw,
                            asQuantity(line[2]),
                            asQuantity(line[3]));
            if (itsDoChecksForNoise) {
                // integration time, spectral resolution and array sizes can change for every scan (in principle)
                // the following checks consistency and throws an exception in the case of mismatch
                // (this avoids nasty problems when wrong value is used; it seems impractical to spend time implementing
                // unconstrained case)
                itsSim->areaTimesSqrtBT();
            }
        }

        ASKAPLOG_INFO_STR(logger, "Successfully simulated " << nScans << " scans");
        ASKAPDEBUGASSERT(itsMs);
        itsMs->flush();
    }
    if (itsComms.isWorker() != itsMSWrittenByMaster) {
        predict(itsMs->tableName());
    } else if (itsMSWrittenByMaster && itsComms.isWorker()) {
        // to run clients
        predict("");
    }
}

void SimParallel::predict(const string& ms)
{
    if (itsComms.isWorker() != itsMSWrittenByMaster) {
        ASKAPDEBUGASSERT(ms != "");
        casa::Timer timer;
        timer.mark();
        ASKAPLOG_INFO_STR(logger, "Simulating data for " << ms);
        ASKAPDEBUGASSERT(itsModel);
        ASKAPLOG_INFO_STR(logger, "Model is " << *itsModel);
        TableDataSource ds(ms, TableDataSource::WRITE_PERMITTED);
        IDataSelectorPtr sel = ds.createSelector();
        sel << parset();
        IDataConverterPtr conv = ds.createConverter();
        conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO), "Hz");
        conv->setDirectionFrame(casa::MDirection::Ref(casa::MDirection::J2000));
        // ensure that time is counted in seconds since 0 MJD
        conv->setEpochFrame(); 
        IDataSharedIter it = ds.createIterator(sel, conv);
        if (itsComms.isWorker()) {
            // default case, workers write their own measurement sets
            predict(it);
        } 
        if (itsComms.isMaster() && itsMSWrittenByMaster) {
            // server code
            ParallelWriteIterator::masterIteration(itsComms, it);
        }        
        ASKAPLOG_INFO_STR(logger,  "Predicted data for " << ms << " in " << timer.real() << " seconds ");
    }
    if (itsComms.isWorker() && itsMSWrittenByMaster) {
        // client code
        ASKAPDEBUGASSERT(ms == "");
        casa::Timer timer;
        timer.mark();
        IDataSharedIter it(new ParallelWriteIterator(itsComms));
        predict(it);
        ASKAPLOG_INFO_STR(logger,  "Finished prediction in worker at rank " << itsComms.rank() << " in " << timer.real() << " seconds ");
    }
}       

/// Predict data for current model
/// @param it data iterator to store the result to
void SimParallel::predict(IDataSharedIter &it)
{
    if (itsComms.isWorker()) { 
        ASKAPDEBUGASSERT(it);      
        /// Create the gridder using a factory acting on a
        /// parameterset
        IVisGridder::ShPtr gridder = createGridder(itsComms, parset());
        ASKAPCHECK(gridder, "Gridder not defined correctly");

        // a part of the equation defined via image
        askap::scimath::Equation::ShPtr imgEquation;

        if (SynthesisParamsHelper::hasImage(itsModel)) {
            ASKAPLOG_INFO_STR(logger, "Sky model contains at least one image, building an image-specific equation");
            // it should ignore parameters which are not applicable (e.g. components)
            imgEquation.reset(new ImageFFTEquation(*itsModel, it, gridder));
        }

        // a part of the equation defined via components
        boost::shared_ptr<ComponentEquation> compEquation;

        if (SynthesisParamsHelper::hasComponent(itsModel)) {
            // model is a number of components
            ASKAPLOG_INFO_STR(logger, "Sky model contains at least one component, building a component-specific equation");
            // it doesn't matter which iterator is passed below. It is not used
            // it should ignore parameters which are not applicable (e.g. images)
            compEquation.reset(new ComponentEquation(*itsModel, it));
        }

        // the measurement equation used for prediction
        // actual type depends on what we are simulating
        // therefore it is uninitialized at the moment
        askap::scimath::Equation::ShPtr equation;

        if (imgEquation && !compEquation) {
            ASKAPLOG_INFO_STR(logger, "Pure image-based model (no components defined)");
            equation = imgEquation;
        } else if (compEquation && !imgEquation) {
            ASKAPLOG_INFO_STR(logger, "Pure component-based model (no images defined)");
            equation = compEquation;
        } else if (imgEquation && compEquation) {
            ASKAPLOG_INFO_STR(logger, "Making a sum of image-based and component-based equations");
            equation = imgEquation;
            addEquation(equation, compEquation, it);
        } else {
            ASKAPTHROW(AskapError, "No sky models are defined");
        }

        if (itsSolutionSource) {
            corruptEquation(equation, it);
        } else {
            ASKAPLOG_INFO_STR(logger, "Calibration effects are not simulated");
        }

        ASKAPCHECK(equation, "Equation is not defined correctly");

        if (itsNoiseVariance > 0.) {
            ASKAPDEBUGASSERT(itsSim);

            const casa::Int seed1 = getSeed("noise.seed1","time");                  
            const casa::Int seed2 = getSeed("noise.seed2","%w");                  

            ASKAPLOG_INFO_STR(logger, "Set seed1 to " << seed1);
            ASKAPLOG_INFO_STR(logger, "Set seed2 to " << seed2);

            if (itsSim->relAntennaWeight().nelements()) {
                ASKAPLOG_INFO_STR(logger, "Gaussian noise (matching SIGMA column of the measurement set) will be added to visibilities");
                boost::shared_ptr<GaussianNoiseME const> noiseME(new GaussianNoiseME(seed1, seed2));
                addEquation(equation, noiseME, it);
            } else {
                ASKAPLOG_INFO_STR(logger, "Gaussian noise (variance=" << itsNoiseVariance <<
                         " Jy^2 or sigma="<<sqrt(itsNoiseVariance)<<" Jy) will be added to visibilities");

                boost::shared_ptr<GaussianNoiseME const> noiseME(new
                       GaussianNoiseME(itsNoiseVariance, seed1, seed2));
                addEquation(equation, noiseME, it);
           }
        }

        equation->predict();
    }
}

/// @brief helper method to obtain the noise per visibility
/// @details Depending on the parameters, the noise is either read
/// directly from the parset (parameters should not have any prefix, e.g.
/// just variance, rms, etc) or calculated using the simulator settings
/// and specified parameters such as Tsys and efficiency (should also be
/// defined without any prefix).
/// @param[in] parset ParameterSet for inputs
/// @return noise variance per visibility in Jy
double SimParallel::getNoise(const LOFAR::ParameterSet& parset) const
{
   ASKAPCHECK(parset.isDefined("Tsys") == parset.isDefined("efficiency"), 
      "Tsys and efficiency parset parameters should either be both defined (for automatic noise calculation) or not "
      "if noise magnitude is overridden with an explicit value given by rms or variance");
   if (parset.isDefined("Tsys")) {
       // automatic noise estimate
       itsDoChecksForNoise = true;
       ASKAPCHECK(!parset.isDefined("rms") && !parset.isDefined("variance"), 
          "If an automatic noise estimate is used, neither 'rms', nor 'variance' parset parameters should be given");
   
       const casa::Vector<double> tSysVector = parset.getDoubleVector("Tsys");
       const casa::Vector<double> effVector = parset.getDoubleVector("efficiency");
       ASKAPCHECK(tSysVector.nelements()>=1, "At least one Tsys has to be defined");
       ASKAPCHECK(effVector.nelements()>=1, "At least one efficiency has to be defined");
       ASKAPLOG_INFO_STR(logger, "Noise level is estimated automatically using Tsys="<<tSysVector<<
                                 " K and efficiency="<<effVector);
       for (casa::uInt ant = 0; ant<tSysVector.nelements(); ++ant) {
            ASKAPCHECK(tSysVector[ant]>0, "Tsys is supposed to be positive, you have "<<tSysVector<<" ant="<<ant);
       }
       for (casa::uInt ant = 0; ant<effVector.nelements(); ++ant) {
            ASKAPCHECK((effVector[ant] > 0) && (effVector[ant] <= 1.), "Efficiency is supposed to be from (0,1] interval, you have "<<
                        effVector<<" ant="<<ant);
       }
       double TsysOverEff = tSysVector[0] / effVector[0];
       ASKAPCHECK((tSysVector.nelements() == effVector.nelements()) || (tSysVector.nelements() == 1) || (effVector.nelements() == 1),
                  "If multiple Tsys and efficiencies are given, their numbers should be equal");
       ASKAPASSERT(itsSim);                          
       if (tSysVector.nelements() * effVector.nelements() > 1) {
           casa::Vector<double> relWeights(tSysVector.nelements() != 1 ? tSysVector.nelements() : effVector.nelements(), 1.);
           for (casa::uInt ant=0; ant<relWeights.nelements(); ++ant) {
                relWeights[ant] = tSysVector[ant < tSysVector.nelements() ? ant : 0] / 
                                  effVector[ant < effVector.nelements() ? ant : 0];
                if (relWeights[ant] > TsysOverEff) {
                    TsysOverEff = relWeights[ant];
                }
           }
           for (casa::uInt ant=0; ant<relWeights.nelements(); ++ant) {
                relWeights[ant] /= TsysOverEff;
           }
           itsSim->setRelAntennaWeight(relWeights);
       }
       const double rms = 1e26*sqrt(2.)*1.38e-23*TsysOverEff/itsSim->areaTimesSqrtBT();
       ASKAPLOG_INFO_STR(logger, " resulting in peak rms of "<<rms<<" Jy");
       return rms*rms;
   } 
    
   ASKAPCHECK(parset.isDefined("rms") || parset.isDefined("variance"), 
      "If the noise level is explicitly given, either 'rms', or 'variance' parset parameters should be defined");
                    
   ASKAPCHECK(parset.isDefined("rms") != parset.isDefined("variance"), 
      "Please give either 'rms' or 'variance' parset parameter, but not both!");
   if (parset.isDefined("rms")) {
       const double rms = parset.getDouble("rms");
       ASKAPLOG_INFO_STR(logger, "Noise level is given explicitly as rms of "<<rms<<" Jy");
       return rms*rms;
   }
   const double variance = parset.getDouble("variance");
   ASKAPLOG_INFO_STR(logger, "Noise level is given explicitly as a variance of "<<variance<<" Jy^2");
   return variance;
}

/// @brief read seed for the random generator
/// @details This is a helper method to read in seed used to set up random number generator.
/// It applies nesessary substitution rules.
/// @param[in] parname name of the parameter
/// @param[in] defval default value (as string)
/// @return seed
casa::Int SimParallel::getSeed(const std::string &parname,const std::string &defval) const
{
   const std::string seedStr(parset().getString(parname,defval));
   if (seedStr == "time") {
       return casa::Int(time(0));
   }    
   return utility::fromString<casa::Int>(substitute(seedStr));
}

/// @brief a helper method to corrupt the data (opposite to calibration)
/// @details Applying gains require different operations depending on
/// the type of the measurement equation (accessor-based or iterator-based).
/// It is encapsulated in this method. The method accesses itsParset to
/// extract the information about calibration model.
/// @param[in] equation a non-const reference to the shared pointer holding
/// an equation to update
/// @param[in] it iterator over the dataset (this is a legacy of the current
/// design of the imaging code, when equation requires an iterator. It should
/// get away at some stage)
void SimParallel::corruptEquation(boost::shared_ptr<scimath::Equation> &equation,
                                  const IDataSharedIter &it)
{
    ASKAPLOG_INFO_STR(logger, "Making equation to simulate calibration effects");
    boost::shared_ptr<IMeasurementEquation> accessorBasedEquation =
        boost::dynamic_pointer_cast<IMeasurementEquation>(equation);

    if (!accessorBasedEquation) {
        // initialize an adapter
        // to use imaging equation with the calibration framework
        // form a replacement equation first
        const boost::shared_ptr<ImagingEquationAdapter>
        new_equation(new ImagingEquationAdapter);

        // the actual equation is locked inside ImagingEquationAdapter
        // in a shared pointer. We can change equation variable
        // after the next line
        new_equation->assign(equation);
        accessorBasedEquation = new_equation;
    }

    ASKAPDEBUGASSERT(accessorBasedEquation);
    
    boost::shared_ptr<CalibrationMEBase> calME;
    const bool polLeakage = parset().getBool("corrupt.leakage",false);
    if (polLeakage) {
        ASKAPLOG_INFO_STR(logger, "Polarisation leakage will be simulated");
        calME.reset(new CalibrationME<Product<NoXPolGain, LeakageTerm> >(scimath::Params(), it, accessorBasedEquation));
    } else {
        ASKAPLOG_INFO_STR(logger, "Only parallel-hand gains will be simulated. Polarisation leakage will not be simulated.");
        calME.reset(new CalibrationME<NoXPolGain>(scimath::Params(), it, accessorBasedEquation));
    }
    ASKAPDEBUGASSERT(calME);
    // set up the adapter
    equation.reset(new CalibParamsMEAdapter(calME, itsSolutionSource, it));
}

/// @brief a helper method to add up an equation
/// @details Some times it is necessary to replace a measurement equation
/// with a sum of two equations. Typical use cases are adding noise to
/// the visibility data and simulating using a composite model containing
/// both components and images. This method replaces the input equation
/// with the sum of the input equation and some other equation also passed
/// as a parameter. It takes care of equation types and instantiates
/// adapters if necessary.
/// @param[in] equation a non-const reference to the shared pointer holding
/// an equation to update
/// @param[in] other a const reference to the shared pointer holding
/// an equation to be added
/// @param[in] it iterator over the dataset (this is a legacy of the current
/// design of the imaging code, when equation requires an iterator. It should
/// get away at some stage)
/// @note This method can be moved somewhere else, as it may be needed in
/// some other places as well
void SimParallel::addEquation(boost::shared_ptr<scimath::Equation> &equation,
                              const boost::shared_ptr<IMeasurementEquation const> &other,
                              const IDataSharedIter &it)
{
    boost::shared_ptr<IMeasurementEquation> accessorBasedEquation =
        boost::dynamic_pointer_cast<IMeasurementEquation>(equation);

    if (!accessorBasedEquation) {
        // form a replacement equation first
        const boost::shared_ptr<ImagingEquationAdapter>
        new_equation(new ImagingEquationAdapter);
        // the actual equation is locked inside ImagingEquationAdapter
        // in a shared pointer. We can change equation variable
        // after the next line
        new_equation->assign(equation);
        accessorBasedEquation = new_equation;
    }

    // we need to instantiate a new variable and then assign it to avoid
    // assigning a pointer to itself behind the scene, which would happen if
    /// everything is accessor-based up front (although
    // such situations are probably dealt with correctly inside the shared pointer)
    const boost::shared_ptr<scimath::Equation> result(new
            SumOfTwoMEs(accessorBasedEquation, other, it));
    equation = result;
}

}
}
