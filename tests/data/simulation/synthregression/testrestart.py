# regression test for restart of imager in spectral mode

from synthprogrunner import *

def analyseResult(spr, checkWeights=True):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   stats1 = spr.imageStats('image.restored.cont2.fits')
   print("Statistics for restored image 2 major cycles: ",stats1)
   stats2 = spr.imageStats('image.restored.cont1.fits')
   print("Statistics for restored image restarted: ",stats2)

   stats3 = spr.imageStats('residual.cont2.fits')
   print("Statistics for residual image 2 major cycles: ",stats3)
   stats4 = spr.imageStats('residual.cont1.fits')
   print("Statistics for residual image restarted: ",stats4)

   diff1 = stats1['rms'] - stats2['rms']
   diff2 = stats3['rms'] - stats4['rms']
   
   if abs(diff1)> 1e-5 or abs(diff2)>1e-5:
      raise RuntimeError("Images differ too much")


import os
os.system("rm -rf *.cont?.fits");
os.system("rm -rf spectral.ms");



spr = SynthesisProgramRunner(template_parset = 'simulator.in')
spr.addToParset("Csimulator.dataset = spectral.ms")

spr.runSimulator()

spr2 = SynthesisProgramRunner(template_parset = 'testrestart.in')
spr2.addToParset("Cimager.dataset = spectral.ms")

# set OMP_NUM_THREADS to 1, otherwise there could be adverse interaction between OpenMP and MPI
# (5 ranks is close to what we could have in CI or ordinary desktop environment, so no room for OpenMP).
# See AXA-2881 for details.
os.environ['OMP_NUM_THREADS']='1'

# first run 2 major cycles
if "CI" in os.environ:
    spr2.runNewImagerParallel(5)
else:
    spr2.runNewImagerParallel(5)

# now change to 1 major cycle, run twice (reading model 2nd time)
spr2.addToParset("Cimager.Images.Names = [image.cont1]")
spr2.addToParset("Cimager.ncycles = 1")

# run 1 major cycle
if "CI" in os.environ:
    spr2.runNewImagerParallel(5)
else:
    spr2.runNewImagerParallel(5)

spr2.addToParset("Cimager.Images.reuse = true")

# run another major cycle, restarting from previous model
if "CI" in os.environ:
    spr2.runNewImagerParallel(5)
else:
    spr2.runNewImagerParallel(5)


analyseResult(spr2)

#clean up
os.system("rm -rf spectral.ms *.cont?.fits temp_parset.in")
