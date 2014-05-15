# regression tests of faceting
# some fixed parameters are given in wtermtest_template.in

from synthprogrunner import *

def analyseResult(spr):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   src_offset = 0.006/math.pi*180.
   psf_peak=[-172.5,-45]
   true_peak=sinProjection(psf_peak,src_offset,src_offset)
   stats = spr.imageStats('image.field1.restored')
   print "Statistics for restored image: ",stats
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)
   if abs(stats['peak']-1.)>0.3:
      raise RuntimeError, "Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak']

   stats = spr.imageStats('image.field1.facet.0.1')
   print "Statistics for modelimage: ",stats
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)

   stats = spr.imageStats('psf.field1.facet.0.1')
   print "Statistics for psf image: ",stats
   # 256 is facetstep of 512 divided by 2, cell is 4 arcsec
   facet_offset = 4.*256./3600.
   facet_psf_peak = sinProjection(psf_peak,facet_offset,facet_offset)
   disterr = getDistance(stats,facet_psf_peak[0],facet_psf_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)

   stats = spr.imageStats('psf.image.field1.facet.0.1')
   print "Statistics for preconditioned psf image: ",stats
   disterr = getDistance(stats,facet_psf_peak[0],facet_psf_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)
   if abs(stats['peak']-1.)>0.01:
      raise RuntimeError, "Peak flux in the preconditioned psf image is notably different from 1.0, F=%f" % stats['peak']

   stats = spr.imageStats('weights.field1.facet.0.1')
   print "Statistics for weight image: ",stats
   if abs(stats['rms']-stats['peak'])>0.1 or abs(stats['rms']-stats['median'])>0.1 or abs(stats['peak']-stats['median'])>0.1:
      raise RuntimeError, "Weight image is expected to be constant for SphFunc gridder"

   stats = spr.imageStats('residual.field1.facet.0.1')
   print "Statistics for residual image: ",stats
   if stats['rms']>0.01 or abs(stats['median'])>0.0001:
      raise RuntimeError, "Residual image has too high rms or median. Please verify"


spr = SynthesisProgramRunner(template_parset = 'facetingtest_template.in')
spr.runSimulator()

print "Testing faceting with overlap"
spr.addToParset("Cimager.Images.shape  = [1024,1024]")
spr.addToParset("Cimager.Images.image.field1.facetstep = 512")
spr.runImager()
analyseResult(spr)

print "Testing faceting with padding"
spr.initParset()
spr.addToParset("Cimager.Images.shape  = [512,512]")
spr.addToParset("Cimager.gridder.padding  = 2")
spr.runImager()
analyseResult(spr)
