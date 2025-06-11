# test starting spectral imaging using an MFS model

from synthprogrunner import *

def analyseResult(spr,downsample=False):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   ext=""
   if downsample:
      ext=".downsampled"
   stats = spr.imageStats('image.cont.taylor.0.restored.fits')
   print("Statistics for taylor-0 restored image: ",stats)
   if abs(stats['peak']-1.)>0.1:
      raise RuntimeError("Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak'])

   stats = spr.imageStats('residual.cont.taylor.0.fits')
   print("Statistics for taylor-0 residual image: ",stats)
   if stats['rms']>0.01 or abs(stats['median'])>0.0001:
      raise RuntimeError("Residual image has too high rms or median. Please verify")

   stats = spr.imageStats('image.restored.cube-res'+ext+'.fits')
   print("Statistics for restored spectral cube: ",stats)
   if abs(stats['peak']-1.)>0.1:
      raise RuntimeError("Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak'])

   stats = spr.imageStats('residual.cube-res'+ext+'.fits')
   print("Statistics for residual spectral cube: ",stats)
   if abs(stats['peak'])>0.02:
      raise RuntimeError("Peak flux in the residual image is too large, F=%f" % stats['peak'])
   

import os
os.system("rm -rf image.*")
os.system("rm -rf testmfsmodel.ms")

spr = SynthesisProgramRunner(template_parset = 'simulator-si.in')
spr.addToParset("Csimulator.dataset = testmfsmodel.ms")

spr.runSimulator()
spr2 = SynthesisProgramRunner(template_parset = 'testmfsmodel.in')
spr.addToParset("Cimager.dataset = testmfsmodel.ms")
os.environ['OMP_NUM_THREADS']='1'
print("INFO About to Run new Imager")

if "CI" in os.environ:
    spr2.runNewImagerParallel(4)
else:
    spr2.runNewImagerParallel(4)

spr.addToParset("Cimager.solverpercore = true")
spr.addToParset("Cimager.nworkergroups = 1")
spr.addToParset("Cimager.combinechannels = false")
spr.addToParset("Cimager.singleoutputfile = true")
spr.addToParset("Cimager.write.beamlog = false")
spr.addToParset("Cimager.visweights = ")
spr.addToParset("Cimager.ncycles = 1")
spr.addToParset("Cimager.solver.Clean.niter = 1")
spr.addToParset("Cimager.solver.Clean.gain = 0.001")
spr.addToParset("Cimager.mfsstartingmodel = true")
spr.addToParset("Cimager.Images.Names = [image.cube-res]")
spr.addToParset("Cimager.sources.names = [xyz]")
spr.addToParset("Cimager.sources.xyz.model = image.cont")
spr.addToParset("Cimager.sources.xyz.nterms = 2")
spr.addToParset("Cimager.Images.shape = [192,192]")

if "CI" in os.environ:
    spr2.runNewImagerParallel(2)
else:
    spr2.runNewImagerParallel(2)

analyseResult(spr)

# now try with different cellsize
spr.addToParset("Cimager.Images.shape = [160, 160]")
spr.addToParset("Cimager.Images.cellsize = [3.0arcsec, 3.0arcsec]")
spr.addToParset("Cimager.Images.Names = [image.cube-res.downsampled]")

if "CI" in os.environ:
    spr2.runNewImagerParallel(2)
else:
    spr2.runNewImagerParallel(2)

analyseResult(spr,downsample=True)


#clean up
os.system("rm -rf testmfsmodel.ms *.cont.taylor.* *.cube-res.* temp_parset.in")
