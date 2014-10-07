#! /usr/bin/env python

import math
import numpy
import csv
import sys
import os
import string

from locations import *
from antennas import *

def usage(config):
	k = config.keys()
	k.sort()
	print "Usage:\n%s array [format]\n" %(sys.argv[0])
	print "  array = %s" %(" | ".join(k))
	print "         This specifies the array or subset of the array to use.\n"
	print "  format = parset | calc | csv | kml"

datafile = "ConfigurationData/ASKAP-SEIC-0005_Antenna_Configuration.csv"

config = {}
config["A27CR3P6B"] = AntennaConfig("A27CR3P6B", range(1, 37), datafile)
config["A27CR3"] = AntennaConfig("A27CR3", range(1, 31), datafile)
config["BETA"] = AntennaConfig("BETA", (1, 3, 6, 8, 9, 29), datafile)
config["BETA2"] = AntennaConfig("BETA2", (1, 3, 6, 8, 9, 2), datafile)
config["BETA4"] = AntennaConfig("BETA4", (1, 3, 6, 8, 9, 4), datafile)
config["BETA5"] = AntennaConfig("BETA5", (1, 3, 6, 8, 9, 5), datafile)
config["BETA7"] = AntennaConfig("BETA7", (1, 3, 6, 8, 9, 7), datafile)
config["BETA10"] = AntennaConfig("BETA10", (1, 3, 6, 8, 9, 10), datafile)
config["BETA11"] = AntennaConfig("BETA11", (1, 3, 6, 8, 9, 11), datafile)
config["BETA12"] = AntennaConfig("BETA12", (1, 3, 6, 8, 9, 12), datafile)
config["BETA13"] = AntennaConfig("BETA13", (1, 3, 6, 8, 9, 13), datafile)
config["BETA14"] = AntennaConfig("BETA14", (1, 3, 6, 8, 9, 14), datafile)
config["BETA15"] = AntennaConfig("BETA15", (1, 3, 6, 8, 9, 15), datafile)
config["BETA16"] = AntennaConfig("BETA16", (1, 3, 6, 8, 9, 16), datafile)
config["BETA17"] = AntennaConfig("BETA17", (1, 3, 6, 8, 9, 17), datafile)
config["BETA18"] = AntennaConfig("BETA18", (1, 3, 6, 8, 9, 18), datafile)
config["BETA19"] = AntennaConfig("BETA19", (1, 3, 6, 8, 9, 19), datafile)
config["BETA20"] = AntennaConfig("BETA20", (1, 3, 6, 8, 9, 20), datafile)
config["BETA21"] = AntennaConfig("BETA21", (1, 3, 6, 8, 9, 21), datafile)
config["BETA22"] = AntennaConfig("BETA22", (1, 3, 6, 8, 9, 22), datafile)
config["BETA23"] = AntennaConfig("BETA23", (1, 3, 6, 8, 9, 23), datafile)
config["BETA24"] = AntennaConfig("BETA24", (1, 3, 6, 8, 9, 24), datafile)
config["BETA25"] = AntennaConfig("BETA25", (1, 3, 6, 8, 9, 25), datafile)
config["BETA26"] = AntennaConfig("BETA26", (1, 3, 6, 8, 9, 26), datafile)
config["BETA27"] = AntennaConfig("BETA27", (1, 3, 6, 8, 9, 27), datafile)
config["BETA28"] = AntennaConfig("BETA28", (1, 3, 6, 8, 9, 28), datafile)
config["BETA29"] = AntennaConfig("BETA29", (1, 3, 6, 8, 9, 29), datafile)
config["BETA29A"]=AntennaConfig("BETA29A", (1, 3, 6, 8,11, 29), datafile)
config["BETA29B"]=AntennaConfig("BETA29B", (1, 3, 8, 9,11, 29), datafile)
config["BETA30"] = AntennaConfig("BETA30", (1, 3, 6, 8, 9, 30), datafile)
config["BETA31"] = AntennaConfig("BETA31", (1, 3, 6, 8, 9, 31), datafile)
config["BETA32"] = AntennaConfig("BETA32", (1, 3, 6, 8, 9, 32), datafile)
config["BETA33"] = AntennaConfig("BETA33", (1, 3, 6, 8, 9, 33), datafile)
config["BETA34"] = AntennaConfig("BETA34", (1, 3, 6, 8, 9, 34), datafile)
config["BETA35"] = AntennaConfig("BETA35", (1, 3, 6, 8, 9, 35), datafile)
config["BETA36"] = AntennaConfig("BETA36", (1, 3, 6, 8, 9, 36), datafile)
config["ADE6_POSS_SJ1"]          = AntennaConfig("ADE6_POSS_SJ1",          ( 4, 5,11,17,25,26), datafile)
config["ADE6_POSS_SJ2"]          = AntennaConfig("ADE6_POSS_SJ2",          ( 4, 5,11,12,19,20), datafile)
config["ADE6_POSS_MJK"]          = AntennaConfig("ADE6_POSS_MJK",          ( 4, 5,10,11,18,19), datafile)
config["ADE6_POSS_JB"]           = AntennaConfig("ADE6_POSS_JB",           (10,11,18,21,22,29), datafile)
config["ADE6_POSS_JB2"]          = AntennaConfig("ADE6_POSS_JB",           (10,11,18,21,22, 7), datafile)
config["ADE6_POSS_JB3"]          = AntennaConfig("ADE6_POSS_JB",           (10,11,18,21,22,19), datafile)
config["ADE12_POSS_EMU1.5"]      = AntennaConfig("ADE12_POSS_EMU1.5",      ( 4, 5, 7,10,11,12,14,16,18,21,23,26), datafile)
config["ADE12_POSS_EMU3"]        = AntennaConfig("ADE12_POSS_EMU3",        ( 4, 5,10,18,19,20,25,27,28,29,35,36), datafile)
config["ADE12_POSS_EMU6"]        = AntennaConfig("ADE12_POSS_EMU6",        ( 4, 5,10,17,19,20,27,29,32,33,35,36), datafile)
config["ADE12_POSS_WALLABY"]     = AntennaConfig("ADE12_POSS_WALLABY",     ( 2, 4, 5, 7,10,11,12,14,16,17,18,19), datafile)
config["ADE6_AK07_PROPOSED"]     = AntennaConfig("ADE6_AK07_PROPOSED",     ( 2, 4, 5, 7,10,12), datafile)
config["ADE6_AK19_PROPOSED"]     = AntennaConfig("ADE6_AK19_PROPOSED",     ( 2, 4, 5,19,10,12), datafile)
config["ADE6_AK11_PROPOSED"]     = AntennaConfig("ADE6_AK11_PROPOSED",     ( 2, 4, 5,11,10,12), datafile)
config["ADE6_AK14_PROPOSED"]     = AntennaConfig("ADE6_AK14_PROPOSED",     ( 2, 4, 5,14,10,12), datafile)
config["ADE6_AK17_PROPOSED"]     = AntennaConfig("ADE6_AK17_PROPOSED",     ( 2, 4, 5,17,10,12), datafile)
config["ADE12_AK07_PROPOSED"]    = AntennaConfig("ADE12_AK07_PROPOSED",    ( 2, 4, 5, 7,10,12,24,27,30,13,16,28), datafile)
config["ADE12_AK19_PROPOSED"]    = AntennaConfig("ADE12_AK19_PROPOSED",    ( 2, 4, 5,19,10,12,24,27,30,13,16,28), datafile)
config["ADE12_AK11_PROPOSED"]    = AntennaConfig("ADE12_AK11_PROPOSED",    ( 2, 4, 5,11,10,12,24,27,30,13,16,28), datafile)
config["ADE12_AK14_PROPOSED"]    = AntennaConfig("ADE12_AK14_PROPOSED",    ( 2, 4, 5,14,10,12,24,27,30,13,16,28), datafile)
config["ADE12_AK17_PROPOSED"]    = AntennaConfig("ADE12_AK17_PROPOSED",    ( 2, 4, 5,17,10,12,24,27,30,13,16,28), datafile)
config["ADE18_AK07_PROPOSED"]    = AntennaConfig("ADE18_AK07_PROPOSED",    ( 2, 4, 5, 7,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
config["ADE18_AK11_PROPOSED"]    = AntennaConfig("ADE18_AK11_PROPOSED",    ( 2, 4, 5,11,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
config["ADE18_AK21_PROPOSED"]    = AntennaConfig("ADE18_AK21_PROPOSED",    ( 2, 4, 5,21,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
config["ADE18_AK23_PROPOSED"]    = AntennaConfig("ADE18_AK23_PROPOSED",    ( 2, 4, 5,23,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
config["ADE18_AK17_PROPOSED"]    = AntennaConfig("ADE18_AK17_PROPOSED",    ( 2, 4, 5,17,10,12,24,27,30,13,16,28,31,33,35,14,19,26), datafile)
config["ADE12_WALLABY_A"]        = AntennaConfig("ADE12_WALLABY_A",        ( 1, 4, 5, 6, 7,11,13,14,15,19,21,25), datafile)
config["ADE12_WALLABY_B"]        = AntennaConfig("ADE12_WALLABY_B",        ( 1, 2, 6, 9,12,15,16,18,23,27,28,29), datafile)
config["ADE12_WALLABY_C"]        = AntennaConfig("ADE12_WALLABY_C",        ( 2, 3, 4, 6, 7,12,15,18,19,22,23,25), datafile)
#
config["ADE6"]                   = AntennaConfig("ADE6",                   ( 2, 4, 5,14,10,12), datafile)
config["ADE12"]                  = AntennaConfig("ADE12",                  ( 2, 4, 5,14,10,12,24,27,30,13,16,28), datafile)
config["ADE18"]                  = AntennaConfig("ADE18",                  ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26), datafile)
#
# The following are those under consideration for the ASKAP community meeting, October 21, 2014
# 16-antenna configs (ADE12 + 4, either 4BETA or 4 nonBETA)
config["ADE12_4BETA_1368"]       = AntennaConfig("ADE12_4BETA_1368",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3, 6, 8), datafile)
config["ADE12_4BETA_1369"]       = AntennaConfig("ADE12_4BETA_1369",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3, 6, 9), datafile)
config["ADE12_4BETA_136F"]       = AntennaConfig("ADE12_4BETA_136F",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3, 6,15), datafile)
config["ADE12_4BETA_139F"]       = AntennaConfig("ADE12_4BETA_139F",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 3, 9,15), datafile)
config["ADE12_4NONBETA_1"]       = AntennaConfig("ADE12_4NONBETA_1",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28,19,23,26,21), datafile)
config["ADE12_4NONBETA_2"]       = AntennaConfig("ADE12_4NONBETA_2",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28,17,18,20,23), datafile)
config["ADE12_4NONBETA_3"]       = AntennaConfig("ADE12_4NONBETA_3",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28,17,11,22,23), datafile)
config["ADE12_4NONBETA_4"]       = AntennaConfig("ADE12_4NONBETA_4",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28,26,11,22, 7), datafile)
# 20-antenna configs (ADE12 + 8, either 6BETA+2, or 8 non-BETA)
config["ADE12_6BETA_AK31_33"]    = AntennaConfig("ADE12_6BETA_AK31_33",    ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 2, 6, 8, 9,15,31,33), datafile)
config["ADE12_6BETA_AK19_23"]    = AntennaConfig("ADE12_6BETA_AK19_23",    ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 2, 6, 8, 9,15,19,23), datafile)
config["ADE12_6BETA_AK26_33"]    = AntennaConfig("ADE12_6BETA_AK26_33",    ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 2, 6, 8, 9,15,26,33), datafile)
config["ADE12_6BETA_AK31_23"]    = AntennaConfig("ADE12_6BETA_AK31_23",    ( 2, 4, 5,14,10,12,24,27,30,13,16,28, 1, 2, 6, 8, 9,15,31,23), datafile)
config["ADE12_8NONBETA_1"]       = AntennaConfig("ADE12_8NONBETA_1",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26,22,34), datafile)
config["ADE12_8NONBETA_2"]       = AntennaConfig("ADE12_8NONBETA_2",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26,17,21), datafile)
config["ADE12_8NONBETA_3"]       = AntennaConfig("ADE12_8NONBETA_3",       ( 2, 4, 5,14,10,12,24,27,30,13,16,28,31,33,35,19,23,26,32,36), datafile)

outputType = ("parset", "calc", "csv", "kml")
if len(sys.argv) < 2:
	usage(config)
	sys.exit(1)

at = string.upper(sys.argv[1])

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
