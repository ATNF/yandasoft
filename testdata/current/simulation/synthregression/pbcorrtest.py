# regression tests with mosaicing gridders doing single field
# imaging (essentially just primary beam correction, hence the name)
# some fixed parameters are given in pbcorrtest_template.in

from synthprogrunner import *

def analyseResult(spr, expected_flux = 1.):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   src_offset = 0.006/math.pi*180.
   feed_offset = [-0.4,0.4]
   psf_peak=[-172.5,-45]
   true_peak=sinProjection(psf_peak,src_offset-feed_offset[0],src_offset+feed_offset[1])
   stats = spr.imageStats('image.field1.restored')
   print "Statistics for restored image: ",stats
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   print "DAM offsets: ",psf_peak[0],psf_peak[1],true_peak[0],true_peak[1]
   # disterr is dominated by an error in the calculus above: the vector sum
   # of l,m offsets on the great circle does not give the true total offset.
   if disterr > 16:
      raise RuntimeError, "Offset between true and expected position exceeds 2 cell sizes (16 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)
   if abs(stats['peak']-expected_flux)>0.02:
      raise RuntimeError, "Peak flux in the image is notably different from %f Jy (pb-corrected value), F=%f" % (expected_flux,stats['peak'])

   stats = spr.imageStats('weights.field1')
   print "Statistics for weights image: ",stats
   weights_peak=offsetDirection(psf_peak,-feed_offset[0],feed_offset[1])
   disterr = getDistance(stats,weights_peak[0],weights_peak[1])*60.
   if disterr > 1:
      raise RuntimeError, "Offset between true and expected peak weights position exceeds 1 arcmin, d=%f, weights_peak=%s" % (disterr,weights_peak)

   stats = spr.imageStats('psf.field1')
   print "Statistics for psf image: ",stats
   disterr = getDistance(stats,psf_peak[0],psf_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)

   stats = spr.imageStats('residual.field1')
   print "Statistics for residual image: ",stats
   if stats['rms']>0.03 or abs(stats['median'])>0.0001:
      raise RuntimeError, "Residual image has too high rms or median. Please verify"



spr = SynthesisProgramRunner(template_parset = 'pbcorrtest_template.in')
spr.runSimulator()

spr.runImager()
analyseResult(spr,3.26/2.)

spr.initParset()
spr.addToParset("Cimager.restore.equalise = True")
spr.runImager()
analyseResult(spr,3.14/2.)

