# regression tests with mosaicing gridders doing single field
# imaging (essentially just primary beam correction, hence the name)
# some fixed parameters are given in pbcorrtest_template.in

from synthprogrunner import *
import math

def analyseResult(spr, expected_rms):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   stats = spr.imageStats('image.field1.restored')
   print "Statistics for restored image: ",stats
   if abs(stats['rms']-expected_rms) > 0.1*expected_rms:
      raise RuntimeError, "rms in the image is notably different from %e Jy, rms=%e" % (expected_rms,stats['rms'])
   if abs(stats['median']) > expected_rms/10.:
      raise RuntimeError, "median exceeds expected_rms/10: expected_rms=%e Jy, median=%e" % (expected_rms,stats['median'])

   stats = spr.imageStats('residual.field1')
   print "Statistics for residual image: ",stats
   if abs(stats['rms']-expected_rms) > 0.1*expected_rms:
      raise RuntimeError, "rms in the image is notably different from %e Jy, rms=%e" % (expected_rms,stats['rms'])
   if abs(stats['median'])>expected_rms/10.:
      raise RuntimeError, "median exceeds expected_rms/10: expected_rms=%e Jy, median=%e" % (expected_rms,stats['median'])

   stats = spr.imageStats("sensitivity.field1")
   print "Statistics for the sensitivity image: ",stats
   if abs(stats['peak'] - expected_rms) > 0.1*expected_rms:
      raise RuntimeError, "peak sensitivity figure %e Jy is notably different from %e Jy" % (stats['peak'],expected_rms)
   if abs(stats['rms']-stats['peak'])>0.1*stats['peak']:
      raise RuntimeError, "Peak sensitivity figure %e Jy should be similar to rms of %e Jy in the sensitivity image (for spheroidal function gridder)" % (stats['peak'],stats['rms'])


def doTests(spr, rmsNoise = 1e-3):
   """
      Actual imaging followed by tests of the results
   """
   # we simulate 270 cycles (9 cuts, 5 min each, 10sec cycle) and 630 baselines
   # and add gaussian noise with sigma=1e-3 Jy
   ncycles = 270
   nbaselines = 630
   noisePerPol = rmsNoise/math.sqrt(ncycles*nbaselines)
   print "Image noise per polarisation plane for %i cycles and %i baselines (sigma per vis is %e Jy) is %e Jy" % (ncycles, nbaselines,rmsNoise,noisePerPol)

   spr.addToParset("Cimager.Images.image.field1.polarisation = [\"XX\"]")
   spr.runImager()
   analyseResult(spr,noisePerPol)

   # test noise figures after polarisation conversion
   spr.initParset()
   spr.addToParset("Cimager.Images.image.field1.polarisation = [\"I\"]")
   spr.runImager()
   analyseResult(spr,noisePerPol*math.sqrt(2.))

   spr.initParset()
   spr.addToParset("Cimager.Images.image.field1.polarisation = [\"Q\"]")
   spr.runImager()
   analyseResult(spr,noisePerPol*math.sqrt(2.))



print "Simulating using explicit definition of the noise level"
spr = SynthesisProgramRunner(template_parset = 'noisetest_template.in')
variance = 1e-6
spr.addToParset("Csimulator.noise.variance = %e" % variance)
spr.runSimulator()
doTests(spr, math.sqrt(variance))

print "Simulating uxing Tsys and efficiency"
Tsys = 30
efficiency = 0.7
rmsNoise = 1e26*math.sqrt(2.)*1.38e-23*Tsys/efficiency/(math.pi*144/4.)/math.sqrt(16e6*10)
spr.initParset()
spr.addToParset("Csimulator.noise.Tsys = %f" % Tsys)
spr.addToParset("Csimulator.noise.efficiency = %f" % efficiency)
spr.runSimulator()
doTests(spr,rmsNoise)
