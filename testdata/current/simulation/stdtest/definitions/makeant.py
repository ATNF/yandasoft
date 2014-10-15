#! /usr/bin/env python

import math
import numpy
import csv
import sys
import os
import string

from locations import *
from antennas import *
from configs import *

def usage(config):
	k = config.keys()
	k.sort()
	print "Usage:\n%s array [format]\n" %(sys.argv[0])
	print "  array = %s" %(" | ".join(k))
	print "         This specifies the array or subset of the array to use.\n"
	print "  format = parset | calc | csv | kml"


config=getConfigs()

outputType = ("parset", "calc", "csv", "kml")
if len(sys.argv) < 2:
	usage(config)
	sys.exit(1)

ot = "parset"
if len(sys.argv) == 3:
	ot = string.lower(sys.argv[2])
	if not ot in outputType:
		usage(config)
		sys.exit(1)

antennas = AntennaList(config[string.upper(sys.argv[1])], 50, "south")

offxyz = [-2556743.707 - -2556745.438, 5097440.315 - 5097448.114, -2847749.657 - -2847753.833]

if ot == "parset":
	antennas.dump()
if ot == "calc":
	antennas.dumpcalc(3, 0.0, offxyz)
if ot == "csv":
	antennas.dumplatlong()
if ot == "kml":
        antennas.dumpKML()
