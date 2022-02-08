/// @file imcontsub.cc
///
/// @brief Image based Continuum subtraction
///
/// @copyright (c) 2019 CSIRO
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
/// @author Mark Wieringa

// Include package level header file
#include <askap/askap_synthesis.h>


// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/Application.h>
#include <askap/askap/StatReporter.h>
#include <askap/askapparallel/AskapParallel.h>
#include <askap/imageaccess/FitsImageAccessParallel.h>
#include <Common/ParameterSet.h>

// casacore includes
#include <casacore/scimath/Fitting/LinearFitSVD.h>
#include <casacore/scimath/Functionals/Polynomial.h>
#include <casacore/casa/OS/CanonicalConversion.h>
#include <casacore/fits/FITS/FITSDateUtil.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/casa/Quanta/MVTime.h>
// robust contsub C++ version

ASKAP_LOGGER(logger, ".imcontsub");

//using namespace std;
using namespace askap;
using namespace askap::accessors;
using namespace casacore;
//using namespace askap::synthesis;

class ImContSubApp : public askap::Application
{
public:
    virtual int run(int argc, char* argv[])
    {
        // This class must have scope outside the main try/catch block
        askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

        try {
            StatReporter stats;
            LOFAR::ParameterSet subset(config().makeSubset("imcontsub."));
            // we only deal with fits files
            String imagetype = subset.getString("imcontsub.imagetype","fits");
            ASKAPCHECK(imagetype=="fits","imcontsub can only process fits files");

            ASKAPLOG_INFO_STR(logger, "ASKAP image based continuum subtraction application " << ASKAP_PACKAGE_VERSION);

            String infile = subset.getString("inputfitscube","");
            String outfile = subset.getString("outputfitscube","");

            // create output file if empty
            if (outfile=="") outfile = infile + ".contsub";

            float threshold = subset.getFloat("threshold",2.0);
            int order = subset.getInt("order",2);
            int blocksize = subset.getInt("blocksize",0);
            int shift = subset.getInt("shift",0);
            bool interleave = subset.getBool("interleave",false);


            FitsImageAccessParallel accessor(comms);

            // work out channel shift due to bary/lsrk to topo correction
            Vector<int> channelShift = dopplerCorrection(infile, accessor, subset);

            if (comms.isMaster()) {
                ASKAPLOG_INFO_STR(logger,"In = "<<infile <<", Out = "<<
                                      outfile <<", threshold = "<<threshold << ", order = "<< order <<
                                      ", blocksize = " << blocksize << ", shift = "<< shift <<
                                      ", interleave = "<< interleave);
                uint nchan = channelShift.nelements();
                if (channelShift(0)!=0 or channelShift(nchan-1)!=0) {
                    ASKAPLOG_INFO_STR(logger,"Channel shift at start and end of spectrum: "<<channelShift(0)<<", "<<channelShift(nchan-1));
                }
                ASKAPLOG_INFO_STR(logger,"master creates the new output file and copies header");
                ASKAPLOG_INFO_STR(logger,"master creates the new output file and copies header "<<infile<<", "<<outfile);
                accessor.copy_header(infile, outfile);
            }

            // All wait for header to be written
            comms.barrier();

            // Now process the rest of the file in parallel
            // Specify axis of cube to distribute over: 1=y -> array dimension returned: (nx,n,nchan)
            const int iax = 1;
            Array<Float> arr = accessor.read_all(infile, iax);
            // remove degenerate 3rd or 4th axis - cube constructor will fail if there isn't one
            arr.removeDegenerate();
            ASKAPCHECK(arr.shape().size()==3,"imcontsub can only deal with 3D data cubes");
            Cube<Float> cube(arr);
            int nz = cube.shape()(2);
            // Are we processing in blocks of channels (to match beamforming intervals)?
            if (blocksize==0) {
                blocksize = nz;
                shift = 0;
            }

            // Process spectrum by spectrum
            ASKAPLOG_INFO_STR(logger,"Process the spectra");
            // subtract in blocks
            // bary/lsrk complication: channels shifted wrt topo observing frame, shift is freq dependent
            // need to count blocks and channels in topo frame and subtract corresponding bary/lsrk channels
            int step = blocksize;
            if (interleave) step = step / 2;
            Vector<Float> workvec(nz);
            Vector<Float> spec(step);
            for (uint y = 0; y< cube.shape()(1); y++ ) {
                for (uint x = 0; x < cube.shape()(0); x++ ) {
                    // get a reference to the current spectrum from the cube
                    Vector<Float> refvec(cube(Slice(x,1),Slice(y,1),Slice()));
                    // make a working copy
                    workvec = refvec;
                    int ic = 0;
                    int stop = 0;
                    int lastStopsub = 0;
                    while (stop < nz) {
                        int start = -shift + ic * step;
                        stop = min(start + blocksize, nz);
                        // now find start and stop channel before bary/lsrk correction
                        int topostart = start + channelShift(min(max(0,start),nz-1));
                        int topostop = stop + channelShift(stop-1);
                        // the interval we're going to subtract (smaller when interleaving)
                        int startsub = topostart;
                        int stopsub = topostop;
                        if (interleave) {
                            // all but first & last interval - use central 50%
                            if (ic > 0) {
                                startsub = topostart + step/2;
                            }
                            if (stop < nz) {
                                stopsub = topostop - step/2 - (step/2)%2;
                            }
                        }
                        if (lastStopsub > 0) {
                            // make sure we don't skip a channel in bary mode
                            startsub = lastStopsub;
                        }
                        lastStopsub = stopsub;
                        startsub = max(0, startsub);
                        stopsub =  min(stopsub, nz);
                        topostart = max(0, topostart);
                        topostop = min(topostop, nz);
                        // size can change, spec will resize if needed
                        spec.assign(workvec(Slice(topostart, topostop - topostart)));
                        process_spectrum(spec, threshold, order);
                        const size_t length = stopsub - startsub;
                        refvec(Slice(startsub,length)) = spec(Slice(startsub-topostart,length));
                        ic += 1;
                    }
                }
            }


            // Write results to output file - make sure we use the same axis as for reading
            accessor.write_all(outfile,arr,iax);
            ASKAPLOG_INFO_STR(logger,"Done");
            // Done
            stats.logSummary();
            } catch (const AskapError& x) {
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


        Vector<int> dopplerCorrection(String& infile, FitsImageAccessParallel& accessor, LOFAR::ParameterSet& parset)
        {
            IPosition shape = accessor.shape(infile);
            ASKAPCHECK(shape.size()>2,"imcontsub needs at least 3 dimensions in the image");

            // get some header info we need
            String specsys = accessor.getMetadataKeyword(infile, "SPECSYS");

            // default is TOPO
            if (specsys=="") specsys = "TOPOCENT";
            CoordinateSystem cs = accessor.coordSys(infile);
            int spectralAxis = cs.spectralAxisNumber();
            // if we can't find the spectral axis, assume largest of 3rd or 4th axis for topo only
            if (spectralAxis < 0) {
                ASKAPCHECK(specsys.startsWith("TOPO"),"Cannot find spectral axis in cube");
                spectralAxis = 2;
                if (shape.size()>3 && shape(3) > shape(2)) spectralAxis = 3;
            }
            int nFreq = shape(spectralAxis);
            // Fill with zero shift for TOPO
            Vector<int> channelShift(nFreq,0);
            if (!specsys.startsWith("TOPO")) {
                MFrequency::Ref topoFrame = MFrequency::Ref(MFrequency::TOPO);
                MFrequency::Ref freqRefFrame = topoFrame;
                if (specsys.startsWith("BARY")) freqRefFrame = MFrequency::Ref(MFrequency::BARY);
                if (specsys.startsWith("LSR")) freqRefFrame = MFrequency::Ref(MFrequency::LSRK);
                // we need to calculate doppler shift: need direction, location, time, frequency info
                CoordinateSystem cs = accessor.coordSys(infile);
                DirectionCoordinate dc = cs.directionCoordinate();
                MDirection dir(MVDirection(dc.referenceValue()),dc.directionType());
                SpectralCoordinate sc = cs.spectralCoordinate();
                // ASKAP antenna 0
                MPosition pos(MVPosition(-2556088.476234,  5097405.971301, -2848428.398018),MPosition::ITRF);
                string dateObs = accessor.getMetadataKeyword(infile,"DATE-OBS");
                if (dateObs == "") {
                    dateObs = parset.getString("date-obs","");
                    ASKAPCHECK(dateObs!="","Please specify the observing date (DATE-OBS is missing in FITS)");
                }
                string timeSys = accessor.getMetadataKeyword(infile,"TIMESYS");
                if (timeSys=="") timeSys = "UTC";
                MVTime mvtime;
                MEpoch::Types system;
                FITSDateUtil::fromFITS(mvtime, system, dateObs, timeSys);
                MEpoch epoch(MVEpoch(mvtime.get()),system);
                // convert channel frequencies back to TOPO frame
                Vector<double> freqs(nFreq);
                Vector<double> topoFreqs(nFreq);
                for (int chan=0; chan<nFreq; chan++) {
                    sc.toWorld(freqs(chan),double(chan));
                }
                sc.setReferenceConversion(MFrequency::TOPO,epoch,pos,dir);
                for (int chan=0; chan<nFreq; chan++) {
                    sc.toWorld(topoFreqs(chan),double(chan));
                }

                for (int chan=0; chan<nFreq; chan++) {
                    channelShift(chan) = round((freqs(chan) - topoFreqs(chan))/sc.increment()(0));
                }
            }
            return channelShift;
        }

        void process_spectrum(Vector<Float>& vec, float threshold, int order, bool log=false) {
            size_t n = vec.nelements();

            // first remove spectral index

            // use every 10th point, more for shorter spectra
            const int inc = max(1,min(n/100,10));
            const size_t n1 = n/inc;
            Vector<Float> y1(n1);
            for (size_t i=0; i<n1; i++) y1(i) = vec(i*inc);

            // mask out NaNs and extreme outliers and count how many points we have left
            // Note: fractile partially sorts the array!
            Float q5 = fractile(y1,0.05f);
            Float q95 = fractile(y1,0.95f);
            int ngood = 0;
            for (size_t i=0; i<n1; i++) {
                if (!isNaN(y1(i)) && (y1(i) >= q5) && (y1(i) <= q95)) ngood++;
            }
            if (log) {
                ASKAPLOG_INFO_STR(logger,"q5="<<q5<<", q95="<<q95<< ", y1="<< y1);
                ASKAPLOG_INFO_STR(logger,"n1="<<n1<<", ngood="<<ngood);
            }
            Float xmean = 0;
            Float offset = 0;
            Float slope = 0;
            if (ngood > 0) {
                Vector<Float> x(ngood), y(ngood);
                for (size_t i=0, j=0; i<n1; i++) {
                    Float s = vec(i*inc);
                    if (!isNaN(s) && (s >= q5) && (s <= q95)) {
                        y(j) = s;
                        x(j++) = i*inc;
                    }
                }
                if (log) ASKAPLOG_INFO_STR(logger,"y="<<y);


                if (ngood > 1) {
                    xmean = mean(x);
                    x -= xmean;
                    // fit line to spectrum
                    slope = sum(x * y) / sum(x * x);
                    //offset = sum(y) / ngood;
                    offset = fractile(y,0.50f);
                }
            }

            // subtract linear fit and filter out NaNs
            Vector<Float> x_full(n), y_full(n), tmp(n);
            Vector<Bool> mask(n);
            indgen(x_full);
            x_full -= xmean;
            ngood = 0;
            for (size_t i = 0; i<n; i++) {
                y_full(i) = vec(i) - (offset + slope * x_full(i));
                mask(i) = !isNaN(y_full(i));
                if (mask(i)) {
                    tmp(ngood++) = y_full(i);
                }
            }
            if (log) ASKAPLOG_INFO_STR(logger,"offset = "<<offset<<", slope = "<<slope);

            // now select data within threshold of median
            if (ngood==0) return;
            Float q15 = fractile(tmp(Slice(0,ngood)), 0.15f);
            Float q50 = fractile(tmp(Slice(0,ngood)), 0.50f);
            Float iqr = (1/1.35)*2*(q50-q15); // just copying python version - iqr should really be 2*(q50-q25)
            for (size_t i = 0; i<n; i++) mask(i) &= (y_full(i) >= q50 - threshold * iqr) &
                                                    (y_full(i) <= q50 + threshold * iqr);
            if (log) {
                ASKAPLOG_INFO_STR(logger,"q15="<<q15<<", q50="<<q50<< ", iqr="<< iqr);
                for (size_t i = 0; i<n; i++) {
                    if (!mask(i)) ASKAPLOG_INFO_STR(logger,"discarded y("<<i<<")="<<y_full(i)<<", vec("<<i<<")="<<vec(i));
                }
            }

            // fit polynomial of given order
            //LinearFitSVD<Float> fitter;

            // using AutoDiff - slow
            // Polynomial<AutoDiff<Float> > func(order+1);
            // for (int i=0; i<order+1; i++) {
            //     func.setCoefficient( i, 1.0);
            // }
            // fitter.setFunction(func);
            // ASKAPLOG_INFO_STR(logger,"Fitting spectrum");
            // Vector<Float> solution = fitter.fit(x_full, vec, &mask);
            // subtract continuum model from data
            // Vector<Float> model(n);
            // ASKAPLOG_INFO_STR(logger,"Subtracting model");
            // bool ok = fitter.residual(model, x_full, True);
            // vec -= model;

            // using LSQaips - need invert before solve
            // Note: fitter uses double internally, giving float inputs just results in internal copying
            //       and is slower
            LSQaips fitter(order+1);
            Vector<Double> xx(order+1);
            VectorSTLIterator<Double> it(xx);
            for (size_t i=0; i<n; i++) {
                if (mask(i)) {
                    xx(0) = 1;
                    for (int j=1; j<order+1; j++) {
                        xx(j) = xx(j-1) * i;
                    }
                    //ASKAPLOG_INFO_STR(logger,"makeNorm "<<i<<" "<<vec(i));
                    fitter.makeNorm(it,1.0,Double(vec(i)));
                }
            }
            uInt nr1;
            //cout << "Invert = " << fitter.invert(nr1);
            //cout << ", rank=" << nr1 << endl;
            Vector<Double> solution(order+1);
            if (fitter.invert(nr1)) {
                fitter.solve(solution);

                Vector<Float> model(n);
                for (size_t i=0; i<n; i++) {
                    model(i)=solution(order);
                    for (int j=order-1; j>=0; j--) {
                        model(i) = i * model(i) + solution(j);
                    }
                }
                if (log) ASKAPLOG_INFO_STR(logger,"solution="<<solution<<", model(10)="<<model(10));
                vec -= model;
            }
        }

    private:
        void getVersion() const override
        {
            std::cerr << ASKAP_PACKAGE_VERSION << std::endl;
        }

};

// Main function
int main(int argc, char* argv[])
{
    ImContSubApp app;
    return app.main(argc, argv);
}
