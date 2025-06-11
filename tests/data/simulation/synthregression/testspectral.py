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

   # now check grids
   stats = spr.complexImageStats('psfgrid.cont')
   print("Statistics for psf grid: ", stats)
   if abs(stats['peak'] - 37.44) > 0.1:
      raise RuntimeError("psf grid peak value is wrong, please verify")
   if stats['x']!=65 or stats['y']!=67:
      raise RuntimeError("psf grid peak location is wrong, please verify")
   if abs(stats['rms-real'] - 2.456) > 0.01:
      raise RuntimeError("psf grid rms value is wrong, please verify")
   if abs(stats['median-real'] - 0) > 1e-5:
      raise RuntimeError("psf grid median value is wrong, please verify")
   if abs(stats['rms-imag'] - 0.003827) > 0.001:
      raise RuntimeError("psf grid rms value is wrong, please verify")
   if abs(stats['median-imag'] - 0) > 1e-5:
      raise RuntimeError("psf grid median value is wrong, please verify")

   stats = spr.complexImageStats('pcfgrid.cont')
   print("Statistics for pcf grid: ", stats)
   if abs(stats['peak'] - 284.6) > 0.1:
      raise RuntimeError("pcf grid peak value is wrong, please verify")
   if stats['x']!=59 or stats['y']!=67:
      raise RuntimeError("pcf grid peak location is wrong, please verify")
   if abs(stats['rms-real'] - 3.6928) > 0.01:
      raise RuntimeError("pcf grid rms value is wrong, please verify")
   if abs(stats['median-real'] - 0) > 1e-5:
      raise RuntimeError("pcf grid median value is wrong, please verify")
   if abs(stats['rms-imag'] - 11.0782) > 0.01:
      raise RuntimeError("pcf grid rms value is wrong, please verify")
   if abs(stats['median-imag'] - 0) > 1e-5:
      raise RuntimeError("pcf grid median value is wrong, please verify")

   stats = spr.complexImageStats('visgrid.cont')
   print("Statistics for vis grid: ", stats)
   if abs(stats['peak'] - 5.22346) > 0.1:
      raise RuntimeError("vis grid peak value is wrong, please verify")
   if stats['x']!=65 or stats['y']!=67:
      raise RuntimeError("vis grid peak location is wrong, please verify")
   if abs(stats['rms-real'] - 0.1794) > 0.001:
      raise RuntimeError("vis grid rms real value is wrong, please verify")
   if abs(stats['median-real'] - 0) > 1e-5:
      raise RuntimeError("vis grid median value is wrong, please verify")
   if abs(stats['rms-imag'] - 0.2038) > 0.001:
      raise RuntimeError("vis grid rms imag value is wrong, please verify")
   if abs(stats['median-imag'] - 0) > 1e-5:
      raise RuntimeError("vis grid median value is wrong, please verify")

import os
os.system("rm -rf *.cont");
os.system("rm -rf spectral.ms");



spr = SynthesisProgramRunner(template_parset = 'simulator.in')
spr.addToParset("Csimulator.dataset = spectral.ms")

spr.runSimulator()

spr2 = SynthesisProgramRunner(template_parset = 'testspectral.in')
spr2.addToParset("Cimager.dataset = spectral.ms")

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
os.system("rm -rf spectral.ms *.cont *.cont.txt temp_parset.in")
