#!/usr/bin/env python

import math
import csv
import sys
import os
import string

from locations import *
from antennas import *

# The following are the reference points for the meerkat array, and hence the rest of SKA1-mid
meerkatRefLong =  21.44389 
meerkatRefLat  = -30.71317
meerkatRefEasting = 542502.836
meerkatRefNorthing = 6602101.300
meerkatRefZone = 34

askapRefZone = 50

def usage():
	print "Usage:\n%s array [format]\n" %(sys.argv[0])
	print "  array = surveyV2 | lowV2 | midV2 | low | survey | mid | lowCore | lowArms | surveyASKAP | surveyCore | surveyArms | midMeerkat | midCore | midArms"
	print "         This specifies the array or subset of the array to use (those without 'V2' are the 2013 versions).\n"
	print "  format = parset | calc | csv | kml"

config={}
config['LOWCORE'] = AntennaConfig('SKA1LOWCORE',range(1,867), 'ConfigurationData/SKA-low_config_baseline_design_core_stations_2013apr30.csv')
config['LOWARMS'] = AntennaConfig('SKA1LOWARMS',range(1,45), 'ConfigurationData/SKA-low_config_baseline_design_arm_stations_2013apr30.csv')
config['LOW'] = AntennaConfig('SKA1LOW',range(1,911), 'ConfigurationData/SKA-low_config_baseline_design_ALL_stations_2013apr30.csv')
config['SURVEYASKAP'] = AntennaConfig('SKA1SURVEYASKAP', range(1,37), 'ConfigurationData/SKA-survey_config_baseline_design_ASKAP_antennas_2013apr30.csv')
config['SURVEYCORE'] = AntennaConfig('SKA1SURVEYCORE', range(1,19), 'ConfigurationData/SKA1-survey_config_baseline_design_core_antennas_2013apr30.csv')
config['SURVEYARMS'] = AntennaConfig('SKA1SURVEYARMS', range(1,43), 'ConfigurationData/SKA1-survey_config_baseline_design_arm_antennas_2013apr30.csv')
config['SURVEY'] = AntennaConfig('SKA1SURVEY', range(1,97), 'ConfigurationData/SKA-survey_config_baseline_design_ALL_antennas_2013apr30.csv')
config['MIDMEERKAT'] = AntennaConfig('SKA1MIDMEERKAT', range(1,65), 'ConfigurationData/SKA-mid_config_baseline_design_MeerKAT_antennas_2013apr30.csv')
config['MIDCORE'] = AntennaConfig('SKA1MIDCORE', range(1,134), 'ConfigurationData/SKA-mid_config_baseline_design_core_antennas_2013apr30.csv')
config['MIDARMS'] = AntennaConfig('SKA1MIDARMS', range(1,58), 'ConfigurationData/SKA-mid_config_baseline_design_arm_antennas_2013apr30.csv')
config['MID'] = AntennaConfig('SKA1MID', range(1,255), 'ConfigurationData/SKA-mid_config_baseline_design_ALL_antennas_2013apr30.csv')
config['SURVEYV2'] = AntennaConfig('SKA1SURVEYv2', range(1,97), 'ConfigurationData/ska1s96v4b.enu.96x3.csv')
config['LOWV2'] = AntennaConfig('SKA1LOWv2', range(1,753), 'ConfigurationData/ska1low.enu.752x3.csv')
config['MIDV2'] = AntennaConfig('SKA1MIDv2', range(1,255), 'ConfigurationData/ska1midv2s.enu.254x4.csv')


class SKAconfig:
    def __init__(self,config):
        self.config=config
        if(config.name[0:7] == 'SKA1MID'):
            self.refEasting = meerkatRefEasting
            self.refNorthing = meerkatRefNorthing
            self.zone = meerkatRefZone
            self.elevation = 1038
        else:
            self.zone=50
            askap =  AntennaList(AntennaConfig("A27CR3P6B", range(1, 37), "ConfigurationData/ASKAP-SEIC-0005_Antenna_Configuration.csv"), self.zone, "south")
            self.refEasting = askap.antennas[1].easting
            self.refNorthing = askap.antennas[1].northing
            self.elevation = 370.

    def read(self):
        filein=open(self.config.filename)
        self.antlist = AntennaList(self.config,self.zone,"south")
        num=0
        for line in filein:
            num += 1
            name='pad%0d'%num
            cols=line.split(',')
            if len(cols)==2:
                ant = Antenna(num,name,self.refEasting+float(cols[0]),self.refNorthing+float(cols[1]),self.elevation,self.zone,"south")
            else:
                ant = Antenna(num,name,self.refEasting+float(cols[0]),self.refNorthing+float(cols[1]),self.elevation+float(cols[2]),self.zone,"south")
            self.antlist.antennas[num]=ant

    def dump(self,outtype):
        if outtype == "parset":
            self.antlist.dump()
        if outtype == "calc":
            self.antlist.dumpcalc(3, 0.0, offxyz)
        if outtype == "csv":
            self.antlist.dumplatlong()
        if outtype == "kml":
            self.antlist.dumpKML()

arrayType = ("LOW", "SURVEY", "MID", "LOWCORE", "LOWARMS", "SURVEYASKAP", "SURVEYCORE", "SURVEYARMS", "MIDMEERKAT", "MIDCORE", "MIDARMS", 'SURVEYV2', 'LOWV2', 'MIDV2')
outputType = ("parset", "calc", "csv", "kml")
if len(sys.argv) < 2:
	usage()
	sys.exit(1)

at = string.upper(sys.argv[1])
if not at in arrayType:
    usage()
    sys.exit(1)
    
ot = "parset"
if len(sys.argv) == 3:
	ot = string.lower(sys.argv[2])
	if not ot in outputType:
		usage()
		sys.exit(1)

skaconf = SKAconfig(config[at])
skaconf.read()
skaconf.dump(ot)

