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
   stats = spr.imageStats('image.restored.wr.1.cont')
   print("Statistics for taylor-0 restored image: ",stats)
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 5:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (5 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))
   if abs(stats['peak']-1.)>0.1:
      raise RuntimeError("Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak'])
   
   stats = spr.imageStats('image.wr.1.cont')
   print("Statistics for modelimage: ",stats)
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 5:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (5 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))

   stats = spr.imageStats('residual.wr.1.cont')
   print("Statistics for residual image: ",stats)
   if stats['rms']>0.015 or abs(stats['median'])>0.001:
      raise RuntimeError("Residual image has too high rms or median. Please verify")


spr = SynthesisProgramRunner(template_parset = 'simulator-reversed.in')
spr.addToParset("Csimulator.dataset = reversed.ms")
spr.runSimulator()

spr2 = SynthesisProgramRunner(template_parset = 'testspectral.in')

spr.addToParset("Cimager.dataset = reversed.ms")

import os
os.system("rm -rf image.*")
if "CI" in os.environ:
    spr2.runNewImagerParallel(5)
else:
    spr2.runNewImagerParallel(5)

analyseResult(spr2)

