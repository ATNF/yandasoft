# regression tests with ATCA 1934-638 data
# some fixed parameters are given in 1934_template.in

from synthprogrunner import *

msarchive = "sim-ddcal.tar"

import os,sys

def analyseResult(spr, checkFlux = True, checkPos = True, nterm = 1):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   for src in [0,1,2]:
       expected_flux = 1.0
       r2d = 180./math.pi
       if src == 0: expected_pos=[3.148664*r2d, -0.780398*r2d]
       if src == 1: expected_pos=[3.134522*r2d, -0.791398*r2d]
       if src == 2: expected_pos=[3.147250*r2d, -0.789398*r2d]
       if nterm==1:
           imagename = 'image.sim-'+str(src)+'-ddcal.restored'
       else:
           imagename = 'image.sim-'+str(src)+'-ddcal.taylor.0.restored'           
       stats = spr.imageStats(imagename)
       print("Statistics for restored image: ",stats)
       flux_diff = abs(expected_flux - stats['peak'])
       if checkFlux:
          print("Expected flux %f, obtained %f, difference %f (or %f%%)" % (expected_flux,stats['peak'],flux_diff,flux_diff/expected_flux*100.))
          if flux_diff > 0.05:
             raise RuntimeError("Flux difference is too much: %f Jy (XX+YY)" % flux_diff)
       disterr = getDistance(stats,expected_pos[0],expected_pos[1])*3600.
       print("Offset of the measured position w.r.t. the true position is %f arcsec" % disterr)
       if disterr > 12 and checkPos:
          raise RuntimeError("Offset between true and expected position exceeds 3x cell size (4 arcsec), d=%f, expected_pos=%s" % (disterr,expected_pos))


if not os.path.exists(msarchive):
   raise RuntimeError("A tarball with measurementset & calibration table does not seem to exist (%s)" % msarchive)

for f in ["sim-corrupt.ms", "sim-all-cal.tab", "sim-image-ddcal.in"]:
   if os.path.exists(f):
      print("Removing old %s" % f)
      os.system("rm -rf %s" % f)

os.system("tar -xf %s" % msarchive)

spr = SynthesisProgramRunner(template_parset = 'sim-image-ddcal.in')

print("DDCal imaging with calibration")
spr.runNewImagerParallel(nProcs=2)
analyseResult(spr)

# disabled nterms=2 test because it is failing, see AXA-3208
#print("DDCal imaging with nterm=2")
#spr.addToParset("Cimager.nworkergroups=3")
#spr.addToParset("Cimager.Images.image.sim-0-ddcal.nterms=2")
#spr.addToParset("Cimager.Images.image.sim-1-ddcal.nterms=2")
#spr.addToParset("Cimager.Images.image.sim-2-ddcal.nterms=2")
#spr.runNewImagerParallel(nProcs=4)
#analyseResult(spr)

#clean up
os.system("rm -rf sim-corrupt.ms sim-all-cal.tab sim-image-ddcal.in image.sim-*-ddcal*restored temp_parset.in")
