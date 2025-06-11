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
#   src_offset = 0.000
   psf_peak=[294.854166667,-63.7125]
   true_peak=sinProjection(psf_peak,src_offset,0)
   stats = spr.imageStats('image.cont.taylor.0.restored')
   print("Statistics for taylor-0 restored image: ",stats)
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))
   if abs(stats['peak']-1.)>0.1:
      raise RuntimeError("Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak'])

   stats = spr.imageStats('image.cont.taylor.0')
   print("Statistics for modelimage: ",stats)
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))

   stats = spr.imageStats('residual.cont.taylor.0')
   print("Statistics for residual image: ",stats)
   if stats['rms']>0.01 or abs(stats['median'])>0.0001:
      raise RuntimeError("Residual image has too high rms or median. Please verify")

import os
os.system("rm -rf image.*")
os.system("rm -rf combined.ms")



spr = SynthesisProgramRunner(template_parset = 'simulator-reversed.in')
spr.addToParset("Csimulator.dataset = combined.ms")

spr.runSimulator()

spr2 = SynthesisProgramRunner(template_parset = 'imager.in')
spr2.addToParset("Cimager.combinechannels = true")
spr2.addToParset("Cimager.dataset = combined.ms")
if "CI" in os.environ:
    spr2.runNewImagerParallel(4)
else:
    spr2.runNewImagerParallel(4)

analyseResult(spr2)

os.system("rm -rf combined.ms *.cont.* temp_parset.in")
