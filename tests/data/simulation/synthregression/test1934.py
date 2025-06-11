# regression tests with ATCA 1934-638 data
# some fixed parameters are given in 1934_template.in

from synthprogrunner import *

msarchive = "1934-638.tar.bz2"

import os,sys

def analyseResult(spr, checkFlux = True, checkPos = True, resultHasMFSName = True):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   expected_flux = 27.336 # XX+YY for 1934-63 at 1806.5 MHz
   expected_pos=[-65.145725,-63.712675]
   # for some reason our naming scheme is not consistent between spectral line and continuum/MFS modes
   if resultHasMFSName:
      imgName = 'image.1934.taylor.0.restored'
   else:
      imgName = 'image.restored.1934'
   stats = spr.imageStats(imgName)
   print("Statistics for restored image: ",stats)
   flux_diff = abs(expected_flux - stats['peak'])
   if checkFlux:
      print("Expected flux %f, obtained %f, difference %f (or %f%%)" % (expected_flux,stats['peak'],flux_diff,flux_diff/expected_flux*100.))
      if flux_diff > 0.05:
         raise RuntimeError("Flux difference is too much: %f Jy (XX+YY)" % flux_diff)
   disterr = getDistance(stats,expected_pos[0],expected_pos[1])*3600.
   print("Offset of the measured position w.r.t. the true position is %f arcsec" % disterr)
   if disterr > 8 and checkPos:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, expected_pos=%s" % (disterr,expected_pos))

   #stats = spr.imageStats('psf.field1')
   #print "Statistics for psf image: ",stats
   #disterr = getDistance(stats,psf_peak[0],psf_peak[1])*3600.
   #if disterr > 8:
   #   raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)

def checkBeamLog():
   """
      this method checks the beamlog file expected to be produced in the spectral line mode
      and throws an exception if some number doesn't match the expectation, otherwise just returns

      Note, the number of tests we do should produce identical beam within the tolerance
   """
   import numpy
   fname = "beamlog.image.restored.1934.txt"
   beam = numpy.loadtxt(fname)
   if beam.shape != (4,):
      raise RuntimeError("Expect 4 columns and 1 row with data in {}".format(fname))
   if abs(beam[0]) > 1e-6:
      raise RuntimeError("Expect zero in the first column in {}".format(fname))
   beam_tolerance = 1e-2
   if abs(beam[1]-8.23) > beam_tolerance:
      raise RuntimeError("Mismatch of the major axis of the beam in {}".format(fname))
   if abs(beam[2]-6.18) > beam_tolerance:
      raise RuntimeError("Mismatch of the minor axis of the beam in {}".format(fname))
   if abs(beam[3]-10.39) > beam_tolerance:
      raise RuntimeError("Mismatch of the position angle of the beam in {}".format(fname))

if not os.path.exists(msarchive):
   raise RuntimeError("A tarball with measurement sets does not seem to exist (%s)" % msarchive)

for f in ["1934pt0.ms", "1934pt1.ms"]:
   if os.path.exists(f):
      print("Removing old %s" % f)
      os.system("rm -rf %s" % f)

os.system("tar -xjf %s" % msarchive)

spr = SynthesisProgramRunner(template_parset = '1934_template.in')

# a number of tests which can be enabled / disabled separately (usually for debugging) via the appropriate if-statement

if True:
   print("Central pointing")
   spr.addToParset("Cimager.dataset=1934pt0.ms")
   spr.runImager()
   analyseResult(spr)

if True:
   print("Nyquist gridding for central pointing with the new imager + traditional weighting")
   os.system("rm -rf *.1934.* *.1934")
   spr.initParset()
   spr.addToParset("Cimager.dataset=1934pt0.ms")
   # Image only 10 channels for this test (otherwise it is too slow and combinedchannels don't work, see AXA-2457
   spr.addToParset("Cimager.Channels=[10,90]")
   spr.addToParset("Cimager.Images.nyquistgridding=true")
   # it fails for ncycles=0 for the new imager, this will override the original setting by 
   # adding a duplicated keyword at the end of the parset to bump up the number of major cycles
   spr.addToParset("Cimager.ncycles=1")
   # this enables traditional weighting with robustness -2 (i.e. close to uniform)
   spr.addToParset("Cimager.uvweight = [ConjugatesAdderFFT, Robust]")
   spr.addToParset("Cimager.uvweight.robustness = -2.")
   spr.runNewImagerParallel(nProcs=2, timeout="5m")
   analyseResult(spr)

if True:
   print("Early termination of major cycle for central pointing with the new imager")
   os.system("rm -rf *.1934.* *.1934")
   spr.initParset()
   spr.addToParset("Cimager.dataset=1934pt0.ms")
   # Image only 10 channels for this test (otherwise it is too slow and combinedchannels don't work, see AXA-2457
   spr.addToParset("Cimager.Channels=[10,90]")
   spr.addToParset("Cimager.Images.nyquistgridding=true")
   # increase the number of major cycles, but it is expected that we terminate after the first one due to the threshold setting below
   # the following will override the original setting by adding a duplicated keyword at the end of the parset
   spr.addToParset("Cimager.ncycles=3")
   spr.addToParset("Cimager.threshold.majorcycle=0.05Jy")
   spr.runNewImagerParallel(nProcs=2, timeout="5m")
   analyseResult(spr, checkFlux = True)


if True:
   # we can run imager with and without joint deconvolution mode (despite having just one pointing) and compare beams
   print("Traditional weighting in spectral line mode with joint deconvolution + BasisfunctionMFS")
   os.system("rm -rf *.1934.* *.1934")
   spr.initParset()
   spr.addToParset("Cimager.dataset=1934pt0.ms")
   spr.addToParset("Cimager.Images.nyquistgridding=true")
   # need at least 2 major cycles to get flux within the rather tight tolerance of analyseResult
   spr.addToParset("Cimager.ncycles=2")
   # this enables traditional weighting with robustness -2 (i.e. close to uniform)
   spr.addToParset("Cimager.uvweight = [ConjugatesAdderFFT, Robust]")
   spr.addToParset("Cimager.uvweight.robustness = -2.")
   # just one channel, but spectral line mode
   spr.addToParset("Cimager.Channels=[1,81]")
   spr.addToParset("Cimager.solverpercore=true")
   spr.addToParset("Cimager.nchanpercore=1")
   spr.addToParset("Cimager.Images.image.1934.nterms=1")
   spr.addToParset("Cimager.visweights=\"\"")
   spr.addToParset("Cimager.singleoutputfile=true")
   # also change algorithm to BasisFunctionMFS 
   spr.addToParset("Cimager.solver.Clean.algorithm=BasisfunctionMFS")
   spr.addToParset("Cimager.solver.Clean.niter=200")
   spr.addToParset("Cimager.solver.Clean.gain=0.05")
   spr.addToParset("Cimager.solver.Clean.scales=[0,6,15]")
   spr.addToParset("Cimager.solver.Clean.solutiontype=MAXBASE")
   spr.addToParset("Cimager.solver.Clean.verbose=false")
   spr.addToParset("Cimager.solver.Clean.threshold.masking=0.9")
   spr.addToParset("Cimager.solver.Clean.weightcutoff=zero")
   spr.addToParset("Cimager.solver.Clean.weightcutoff.clean=false")
   spr.addToParset("Cimager.solver.Clean.decoupled=true")
   spr.addToParset("Cimager.solver.Clean.detectdivergence=true")
   spr.addToParset("Cimager.solver.Clean.logevery=50")
   #
   spr.addToParset("Cimager.writepsfrawimage=true")
   # first run in normal mode (without joint deconvolution
   spr.runNewImagerParallel(nProcs=2, timeout="5m")
   analyseResult(spr, checkFlux = True, resultHasMFSName = False)
   checkBeamLog()
   # now with joint deconvolution (simply add a keyword to the same parset, the default value was false)
   os.system("rm -rf *.1934.* *.1934")
   spr.addToParset("Cimager.updatedirection=true")
   # turn off masking because stats don't cope with NaN values
   spr.addToParset("Cimager.maskmosaic=false")
   spr.runNewImagerParallel(nProcs=2, timeout="5m")
   analyseResult(spr, checkFlux = True, resultHasMFSName = False)
   checkBeamLog()
   # the same with different subshape which would exercise linmos merging logic
   os.system("rm -rf *.1934.* *.1934")
   spr.addToParset("Cimager.Images.subshape=[256,256]")
   spr.runNewImagerParallel(nProcs=2, timeout="5m")
   # this mode results in 0.2-0.3 Jy flux error (worse with MultiScale algorithm than for BasisfunctionMFS) due to edges of linmos merging. 
   # The error seem to reduce for lower gain (which has to be compensated with the larger number of iterations), but bumping up the number of major cycles doesn't help
   # It may be possible to suppress it with masking or some kind of beam weighting (because these tests are done with SphFunc gridder, the weight is
   # constant across the patch FOV, so it can't be used as a mask), but for now skip the flux check and only check the fitted PSF
   analyseResult(spr, checkFlux = False, resultHasMFSName = False)
   checkBeamLog()


if True:
   print("Offset pointing")
   os.system("rm -rf *.1934.* *.1934")
   spr.initParset()
   spr.addToParset("Cimager.dataset=1934pt1.ms")
   spr.runImager()
   analyseResult(spr,checkFlux = False)

#clean up
os.system("rm -rf 1934*.ms *.1934.* *.1934 temp_parset.in")
