# regression tests of delaysolver using ATCA 1934-638 data (small size compared to ASKAP dataset)
# some fixed parameters are given in delaysolver_template.in, also see cpingest.in which delaysolver uses
# We execute delaysolver pretending that the current directory is the MS directory and it grabs cpingest.in from there to
# look at FCM parameters it needs

from synthprogrunner import *

msarchive = "1934-638.tar.bz2"

import os,sys

def analyseResult(expected = [], tolerance=0.1):
   '''
      This method reads corrected_fixeddelay.parset and compares results with the expectation
      (expected - list of delays expected in ns, tolerance - tolerance for comparison in ns)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   with open("corrected_fixeddelay.parset") as f:
       for line in f:
           parts = line.strip().split("=")
           if len(parts)!=2:
              raise RuntimeError("Unexpected line in the result: %s" % (line,))
           if not parts[1].endswith("ns"):
              raise RuntimeError("Result should be in ns: <%s>" % (line,))
           val = float(parts[1][:-2])
           parts = parts[0].split(".")
           if len(parts)!=4:
              raise RuntimeError("Unexpected line in the result: %s, parsing %s" % (line,parts))
           keys = ["common","antenna","ant","delay"]
           for i,key in enumerate(keys):
               if not parts[i].startswith(key):
                  raise RuntimeError("Parse error in %s" % (parts,))
           ant = int(parts[2][3:])
           if ant == 0 or ant > len(expected):
              raise RuntimeError("Unexpected antenna name in %s" % (parts,))
           if abs(val - expected[ant - 1]) >= tolerance:
              raise RuntimeError("Estimated delay is different for antenna %i, expected: %f got: %f" % (ant, expected[ant - 1], val))
              


if not os.path.exists(msarchive):
   raise RuntimeError("A tarball with measurement sets does not seem to exist (%s)" % msarchive)

for f in ["1934pt0.ms", "1934pt1.ms", "corrected_fixeddelay.parset"]:
   if os.path.exists(f):
      print("Removing old %s" % f)
      os.system("rm -rf %s" % f)

os.system("tar -xjf %s" % msarchive)

spr = SynthesisProgramRunner(template_parset = 'delaysolver_template.in')

print("Central pointing")
# runCommand will add -c with the appropriate parset name
spr.runCommand("delaysolver -s . -d . -f 1934pt0.ms")
analyseResult([0, 0.000152135245, -7.34991922e-05, 3.11753594e-05, -0.000234126642, 2.65392955e-05])

print("Offset pointing")
spr.runCommand("delaysolver -s . -d . -f 1934pt1.ms")
analyseResult([0, 3.94, 5.91, 8.76, 21.14, 35.56])

#clean up
os.system("rm -rf 1934*.ms temp_parset.in")
