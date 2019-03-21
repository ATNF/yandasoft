# regression tests with gridders taking w-term into account
# some fixed parameters are given in wtermtest_template.in

from synthprogrunner import *

def analyseResult(spr, checkWeights=True):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   src_offset = 0.006/math.pi*180.
   psf_peak=[-172.5,-45]
   true_peak=sinProjection(psf_peak,src_offset,src_offset)
   stats = spr.imageStats('image.restored.wr.1.field1')
   print "Statistics for restored image: ",stats
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)
   if abs(stats['peak']-1.)>0.1:
      raise RuntimeError, "Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak']

   stats = spr.imageStats('image.wr.1.field1')
   print "Statistics for modelimage: ",stats
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)

   stats = spr.imageStats('psf.wr.1.field1')
   print "Statistics for psf image: ",stats
   disterr = getDistance(stats,psf_peak[0],psf_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)

   stats = spr.imageStats('psf.image.wr.1.field1')
   print "Statistics for preconditioned psf image: ",stats
   disterr = getDistance(stats,psf_peak[0],psf_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)
   if abs(stats['peak']-1.)>0.01:
      raise RuntimeError, "Peak flux in the preconditioned psf image is notably different from 1.0, F=%f" % stats['peak']

   if checkWeights:
      stats = spr.imageStats('weights.wr.1.field1')
      print "Statistics for weight image: ",stats
      if abs(stats['rms']-stats['peak'])>0.1 or abs(stats['rms']-stats['median'])>0.1 or abs(stats['peak']-stats['median'])>0.1:
          raise RuntimeError, "Weight image is expected to be constant for WProject and WStack gridders"

   stats = spr.imageStats('residual.wr.1.field1')
   print "Statistics for residual image: ",stats
   if stats['rms']>0.01 or abs(stats['median'])>0.0001:
      raise RuntimeError, "Residual image has too high rms or median. Please verify"


spr = SynthesisProgramRunner(template_parset = 'wtermtest_template.in')
spr.runSimulator()

spr.addToParset("Cimager.gridder = WProject")
spr.addToParset("Cimager.solverpercore = true")

#spr.runNewImager()
#analyseResult(spr)

spr.addToParset("Cimager.freqframe = bary")
spr.addToParset("Cimager.Channels = [1,0]")
spr.addToParset("Cimager.Frequencies = [1,1.42e9,-10000]")

spr.runNewImager()
analyseResult(spr)
