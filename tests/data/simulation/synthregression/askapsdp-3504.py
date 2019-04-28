# regression tests with gridders taking w-term into account
# some fixed parameters are given in wtermtest_template.in

from synthprogrunner import *

def analyseResult(spr, checkWeights=True):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   src_offset = 0.000970/math.pi*180.
   psf_peak=[294.854316,-63.712]
   true_peak=sinProjection(psf_peak,src_offset,0)
   stats = spr.imageStats('image.restored.field1')
   print "Statistics for restored image: ",stats
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)
   if abs(stats['peak']-1.)>0.1:
      raise RuntimeError, "Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak']

   stats = spr.imageStats('image.field1')
   print "Statistics for modelimage: ",stats
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)

   stats = spr.imageStats('psf.field1')
   print "Statistics for psf image: ",stats
   disterr = getDistance(stats,psf_peak[0],psf_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)

   stats = spr.imageStats('psf.image.field1')
   print "Statistics for preconditioned psf image: ",stats
   disterr = getDistance(stats,psf_peak[0],psf_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)
   if abs(stats['peak']-1.)>0.01:
      raise RuntimeError, "Peak flux in the preconditioned psf image is notably different from 1.0, F=%f" % stats['peak']
   stats = spr.imageStats('residual.field1')

   print "Statistics for residual image: ",stats
   if stats['rms']>0.01 or abs(stats['median'])>0.0001:
      raise RuntimeError, "Residual image has too high rms or median. Please verify"


spr = SynthesisProgramRunner(template_parset = 'askapsdp-sim-3504.in')
spr.runSimulator()

spr = SynthesisProgramRunner(template_parset = 'askapsdp-3504.in')

spr.runNewImager(9)
analyseResult(spr)


