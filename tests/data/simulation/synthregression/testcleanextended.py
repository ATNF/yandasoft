# regression test for restart of imager in spectral mode

from synthprogrunner import *

def analyseResult(spr, checkWeights=True):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   stats1 = spr.imageStats('image.restored.ext1.fits')
   print("Statistics for restored image standard clean: ",stats1)
   stats2 = spr.imageStats('image.restored.ext2.fits')
   print("Statistics for restored image pixellist clean: ",stats2)

   stats3 = spr.imageStats('residual.ext1.fits')
   print("Statistics for residual image standard clean: ",stats3)
   stats4 = spr.imageStats('residual.ext2.fits')
   print("Statistics for residual image pixellist clean: ",stats4)

   diff1 = stats1['rms'] - stats2['rms']
   diff2 = stats3['rms'] - stats4['rms']
   
   if abs(diff1)> 2e-4 or abs(diff2)>2e-4:
      raise RuntimeError("Images differ too much")


import os
os.system("rm -rf *.ext?.fits");
os.system("rm -rf spectral-ext.ms");



spr = SynthesisProgramRunner(template_parset = 'simulator-extended.in')
spr.addToParset("Csimulator.dataset = spectral-ext.ms")

spr.runSimulator()

spr2 = SynthesisProgramRunner(template_parset = 'testcleanextended.in')
spr2.addToParset("Cimager.dataset = spectral-ext.ms")

# set OMP_NUM_THREADS to 1, otherwise there could be adverse interaction between OpenMP and MPI
# (5 ranks is close to what we could have in CI or ordinary desktop environment, so no room for OpenMP).
# See AXA-2881 for details.
os.environ['OMP_NUM_THREADS']='1'

# first run 2 major cycles, pixellist clean
if "CI" in os.environ:
    spr2.runNewImagerParallel(2)
else:
    spr2.runNewImagerParallel(2)

# now change normal clean
spr2.addToParset("Cimager.Images.Names = [image.ext1]")
spr2.addToParset("Cimager.solver.Clean.usepixellists = false")

# run 2 major cycles standard clean
if "CI" in os.environ:
    spr2.runNewImagerParallel(2)
else:
    spr2.runNewImagerParallel(2)


analyseResult(spr2)

#clean up
os.system("rm -rf spectral-ext.ms *.ext?.fits temp_parset.in")
