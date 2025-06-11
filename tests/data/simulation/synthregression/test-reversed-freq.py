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
   stats = spr.imageStats('image.restored.cont')
   print("Statistics for taylor-0 restored image: ",stats)
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 5:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (5 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))
   if abs(stats['peak']-1.)>0.1:
      raise RuntimeError("Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak'])

   stats = spr.imageStats('image.cont')
   print("Statistics for modelimage: ",stats)
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 5:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (5 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))

   stats = spr.imageStats('residual.cont')
   print("Statistics for residual image: ",stats)
   if stats['rms']>0.015 or abs(stats['median'])>0.001:
      raise RuntimeError("Residual image has too high rms or median. Please verify")

import os
os.system("rm -rf image.*")
os.system("rm -rf reversed.ms")

spr = SynthesisProgramRunner(template_parset = 'simulator-reversed.in')
spr.addToParset("Csimulator.dataset = reversed.ms")
spr.runSimulator()

spr2 = SynthesisProgramRunner(template_parset = 'testspectral.in')

spr.addToParset("Cimager.dataset = reversed.ms")

# set OMP_NUM_THREADS to 1, otherwise there could be adverse interaction between OpenMP and MPI
# (5 ranks is close to what we could have in CI or ordinary desktop environment, so no room for OpenMP).
# See AXA-2881 for details.
os.environ['OMP_NUM_THREADS']='1'

if "CI" in os.environ:
    spr2.runNewImagerParallel(5)
else:
    spr2.runNewImagerParallel(5)

analyseResult(spr2)

#clean up
os.system("rm -rf reversed.ms *.cont *.cont.txt temp_parset.in")
