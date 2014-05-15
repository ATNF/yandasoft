#! /usr/bin/env python

import math
import numpy
import csv
import sys
import os
import string

from locations import *

class Antenna:
	def __init__(self, pad, name, easting, northing, elevation, zone, hemisphere):
		self.pad = pad
		self.name = name
		self.easting = easting
		self.northing = northing
		self.elevation = elevation
		self.zone = zone
		self.hemisphere = hemisphere

	def toITRF(self):
		lat, lon, el = self.toWGS84()
		return WGS84ToITRF(lat, lon, self.elevation)

	def toWGS84(self):
		lat, lon = ToGeographic(self.easting, self.northing, self.zone, self.hemisphere == "south")
		return lat, lon, self.elevation
	
	def hasEl(self):
		if self.elevation < -0.99e9:
			return False
		return True

	def padName(self):
		return "Pad%02d" %(self.pad)

class AntennaList:
    def __init__(self, config, zone, hemisphere):
        self.config = config
        self.antennas = {}
        if config.filename == 'ConfigurationData/ASKAP-SEIC-0005_Antenna_Configuration.csv':
            nhead = 4
            csvreader = csv.reader(open(config.filename, "rU"), dialect="excel")
            for row in csvreader:
                if nhead > 0:
                    nhead -= 1
                    continue
                if row[2]=='':
                    continue
                pad = int(row[0])
                name = row[2]
                self.antennas[pad] = Antenna(pad, name, float(row[6]), float(row[7]), -1.0e9, zone, hemisphere)

        for ant in self.antennas.values():
            if ant.hasEl() == False:
                ant.elevation = 370.0

    def names(self):
        n = []
        for ant in self.config.antennas:
            n.append(self.antennas[ant].padName())
        return n

    def dump(self):
        print "antennas.telescope = %s" %(self.config.name)
        print "antennas.%s.coordinates = global" %(self.config.name)
        print "antennas.%s.names = [%s]" %(self.config.name, ",".join(self.names()))
        print "antennas.%s.diameter = 12m" %(self.config.name)
        print "antennas.%s.scale = 1.0" %(self.config.name)
        print "antennas.%s.mount = equatorial" %(self.config.name)
        for antnum in self.config.antennas:
            print "# easting: %f; northing: %f; el: %f; Zone: %d %s" %(self.antennas[antnum].easting, self.antennas[antnum].northing, self.antennas[antnum].elevation, self.antennas[antnum].zone, self.antennas[antnum].hemisphere)
            lat, lon, el = self.antennas[antnum].toWGS84()
            print "# lat: %16.12f; long: %16.12f; el: %f" %(lat * Rad2Deg, lon * Rad2Deg, el)
            x, y, z = self.antennas[antnum].toITRF()
            print "antennas.%s.%s = [%f, %f, %f]" %(self.config.name, self.antennas[antnum].padName(), x, y, z)

    def dumpcalc(self, axisType, offset, offxyz):
        print "# Offset of %.3f,%.3f,%.3f applied to all antenna positions." %(offxyz[0], offxyz[1], offxyz[2])
        for antnum in self.config.antennas:
            x, y, z = self.antennas[antnum].toITRF()
            print "%12.3f  %12.3f  %12.3f %3d %5.1f  $%s" %(x + offxyz[0], y + offxyz[1], z + offxyz[2], axisType, offset, string.upper(self.antennas[antnum].padName()))

    def dumplatlong(self):
        print "Pad,easting,northing,long,lat"
        for antnum in self.config.antennas:
            lat, lon, el = self.antennas[antnum].toWGS84()
            print "%s,%f,%f,%16.12f,%16.12f" %(self.antennas[antnum].padName(), self.antennas[antnum].easting, self.antennas[antnum].northing, lon * Rad2Deg, lat * Rad2Deg)


    def dumpKML(self):
        print """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
<Style id="whitecirc">
<IconStyle>
<Icon>
<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>
</Icon>
</IconStyle>
</Style>
"""
        for antnum in self.config.antennas:
            lat, lon, el = self.antennas[antnum].toWGS84()
            print """<Placemark>
<styleUrl>#whitecirc</styleUrl>
<name>%d</name>
<Point>
<coordinates>%16.12f, %16.12f</coordinates>
</Point>
</Placemark>"""%(antnum,lon*Rad2Deg,lat*Rad2Deg)

        print """</Document>
</kml>
"""
                

class AntennaConfig:
    def __init__(self, name, antennas, filename):
        self.name = name
        self.antennas = antennas
        self.filename = filename
