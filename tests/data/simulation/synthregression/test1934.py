# regression tests with ATCA 1934-638 data
# some fixed parameters are given in 1934_template.in

from synthprogrunner import *

msarchive = "1934-638.tar.bz2"

import os,sys

def analyseResult(spr, checkFlux = True, checkPos = True):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   expected_flux = 27.336 # XX+YY for 1934-63 at 1806.5 MHz
   expected_pos=[-65.145725,-63.712675]
   stats = spr.imageStats('image.1934.taylor.0.restored')
   print "Statistics for restored image: ",stats
   flux_diff = abs(expected_flux - stats['peak'])
   if checkFlux:
      print "Expected flux %f, obtained %f, difference %f (or %f%%)" % (expected_flux,stats['peak'],flux_diff,flux_diff/expected_flux*100.)
      if flux_diff > 0.05:
         raise RuntimeError, "Flux difference is too much: %f Jy (XX+YY)" % flux_diff
   disterr = getDistance(stats,expected_pos[0],expected_pos[1])*3600.
   print "Offset of the measured position w.r.t. the true position is %f arcsec" % disterr
   if disterr > 8 and checkPos:
      raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, expected_pos=%s" % (disterr,expected_pos)

   #stats = spr.imageStats('psf.field1')
   #print "Statistics for psf image: ",stats
   #disterr = getDistance(stats,psf_peak[0],psf_peak[1])*3600.
   #if disterr > 8:
   #   raise RuntimeError, "Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak)


if not os.path.exists(msarchive):
   raise RuntimeError, "A tarball with measurement sets does not seem to exist (%s)" % msarchive

for f in ["1934pt0.ms", "1934pt1.ms"]:
   if os.path.exists(f):
      print "Removing old %s" % f
      os.system("rm -rf %s" % f)

os.system("tar -xjf %s" % msarchive)

spr = SynthesisProgramRunner(template_parset = '1934_template.in')

print "Central pointing"
spr.addToParset("Cimager.dataset=1934pt0.ms")
spr.runImager()
analyseResult(spr)

print "Offset pointing"
spr.initParset()
spr.addToParset("Cimager.dataset=1934pt1.ms")
spr.runImager()
analyseResult(spr,checkFlux = False)
