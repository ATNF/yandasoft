#! /usr/bin/env python

import math
import numpy
import csv
import sys
import os
import string

# Ellipsoid model constants (actual values here are for WGS84)
sm_a = 6378137.0
invf = 298.257223563
UTMScaleFactor = 0.9996

f = 1.0 / invf
sm_b = sm_a * (1.0 - f)

Deg2Rad = math.pi / 180.0
Rad2Deg = 180.0 / math.pi

#
# UTMCentralMeridian
#
# Determines the central meridian for the given UTM zone.
#
# Inputs:
#     zone - An integer value designating the UTM zone, range [1,60].
#
# Returns:
#   The central meridian for the given UTM zone, in radians, or zero
#   if the UTM zone parameter is outside the range [1,60].
#   Range of the central meridian is the radian equivalent of [-177,+177].
#
#
def UTMCentralMeridian (zone):
    return Deg2Rad * (-183.0 + (zone * 6.0))


#
# FootpointLatitude
#
# Computes the footpoint latitude for use in converting transverse
# Mercator coordinates to ellipsoidal coordinates.
#
# Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
#   GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
#
# Inputs:
#   y - The UTM northing coordinate, in meters.
#
# Returns:
#   The footpoint latitude, in radians.
#
#
def FootpointLatitude (y):
	# Precalculate n (Eq. 10.18)
	n = (sm_a - sm_b) / (sm_a + sm_b)

	# Precalculate alpha_ (Eq. 10.22)
	# (Same as alpha in Eq. 10.17)
	alpha_ = ((sm_a + sm_b) / 2.0) * (1 + (n ** 2.0 / 4) + (n ** 4.0 / 64))

	# Precalculate y_ (Eq. 10.23)
	y_ = y / alpha_

	# Precalculate beta_ (Eq. 10.22)
	beta_ = (3.0 * n / 2.0) + (-27.0 * n ** 3.0 / 32.0) + (269.0 * n ** 5.0 / 512.0)

	# Precalculate gamma_ (Eq. 10.22)
	gamma_ = (21.0 * n ** 2.0 / 16.0) + (-55.0 * n ** 4.0 / 32.0)

	# Precalculate delta_ (Eq. 10.22)
	delta_ = (151.0 * n ** 3.0 / 96.0) + (-417.0 * n ** 5.0 / 128.0)

	# Precalculate epsilon_ (Eq. 10.22)
	epsilon_ = (1097.0 * n ** 4.0 / 512.0)

	# Now calculate the sum of the series (Eq. 10.21)
	result = y_ + (beta_ * math.sin (2.0 * y_)) + (gamma_ * math.sin (4.0 * y_)) + (delta_ * math.sin (6.0 * y_)) + (epsilon_ * math.sin (8.0 * y_))

	return result

#
# MapXYToLatLon
#
# Converts x and y coordinates in the Transverse Mercator projection to
# a latitude/longitude pair.  Note that Transverse Mercator is not
# the same as UTM; a scale factor is required to convert between them.
#
# Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
#   GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
#
# Inputs:
#   x - The easting of the point, in meters.
#   y - The northing of the point, in meters.
#   lambda0 - Longitude of the central meridian to be used, in radians.
#
# Outputs:
#   philambda - A 2-element containing the latitude and longitude
#               in radians.
#
# Returns:
#   The function does not return a value.
#
# Remarks:
#   The local variables Nf, nuf2, tf, and tf2 serve the same purpose as
#   N, nu2, t, and t2 in MapLatLonToXY, but they are computed with respect
#   to the footpoint latitude phif.
#
#   x1frac, x2frac, x2poly, x3poly, etc. are to enhance readability and
#   to optimize computations.
#
#
def MapXYToLatLon (x, y, lambda0):
	# Get the value of phif, the footpoint latitude.
	phif = FootpointLatitude (y)
	
	# Precalculate ep2
	ep2 = (sm_a ** 2.0 - sm_b ** 2.0) / sm_b ** 2.0
	
	# Precalculate cos (phif)
	cf = math.cos (phif)
	
	# Precalculate nuf2
	nuf2 = ep2 * cf ** 2.0
	
	# Precalculate Nf and initialize Nfpow
	Nf = sm_a ** 2.0 / (sm_b * math.sqrt (1.0 + nuf2))
	Nfpow = Nf
	
	# Precalculate tf
	tf = math.tan (phif)
	tf2 = tf * tf
	tf4 = tf2 * tf2

	# Precalculate fractional coefficients for x**n in the equations
	#  below to simplify the expressions for latitude and longitude.
	x1frac = 1.0 / (Nfpow * cf)

	Nfpow *= Nf   # now equals Nf**2)
	x2frac = tf / (2.0 * Nfpow)

	Nfpow *= Nf   # now equals Nf**3)
	x3frac = 1.0 / (6.0 * Nfpow * cf)

	Nfpow *= Nf   # now equals Nf**4)
	x4frac = tf / (24.0 * Nfpow)

	Nfpow *= Nf   # now equals Nf**5)
	x5frac = 1.0 / (120.0 * Nfpow * cf)

	Nfpow *= Nf   # now equals Nf**6)
	x6frac = tf / (720.0 * Nfpow)

	Nfpow *= Nf   # now equals Nf**7)
	x7frac = 1.0 / (5040.0 * Nfpow * cf)

	Nfpow *= Nf   # now equals Nf**8)
	x8frac = tf / (40320.0 * Nfpow)

	# Precalculate polynomial coefficients for x**n.
	#  -- x**1 does not have a polynomial coefficient.
	x2poly = -1.0 - nuf2
	x3poly = -1.0 - 2 * tf2 - nuf2
	x4poly = 5.0 + 3.0 * tf2 + 6.0 * nuf2 - 6.0 * tf2 * nuf2 - 3.0 * (nuf2 *nuf2) - 9.0 * tf2 * (nuf2 * nuf2)
	x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2
	x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2 + 162.0 * tf2 * nuf2
	x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2)
	x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2)
	
	# Calculate latitude
	phi = phif + x2frac * x2poly * x ** 2.0 + x4frac * x4poly * x ** 4.0 + x6frac * x6poly * x ** 6.0 + x8frac * x8poly * x ** 8.0
	
	# Calculate longitude
	lambdaa = lambda0 + x1frac * x + x3frac * x3poly * x ** 3.0 + x5frac * x5poly * x ** 5.0 + x7frac * x7poly * x ** 7.0
	
	return phi, lambdaa

#
# UTMXYToLatLon
#
# Converts x and y coordinates in the Universal Transverse Mercator
# projection to a latitude/longitude pair.
#
# Inputs:
#	x - The easting of the point, in meters.
#	y - The northing of the point, in meters.
#	zone - The UTM zone in which the point lies.
#	southhemi - True if the point is in the southern hemisphere;
#               false otherwise.
#
# Outputs:
#	latlon - A 2-element array containing the latitude and
#            longitude of the point, in radians.
#
# Returns:
#	The function does not return a value.
#
#
def UTMXYToLatLon (x, y, zone, southhemi):
    x -= 500000.0
    x /= UTMScaleFactor
    	
    # If in southern hemisphere, adjust y accordingly.
    if (southhemi):
    	y -= 10000000.0;
    		
    y /= UTMScaleFactor
    
    cmeridian = UTMCentralMeridian (zone)
    return MapXYToLatLon (x, y, cmeridian)

def ToGeographic (x, y, zone, southhemi):
    if (zone < 1) or (60 < zone):
        print "The UTM zone you entered is out of range.  Please enter a number in the range [1, 60]."
        return None, None
    lat, lon =  UTMXYToLatLon (x, y, zone, southhemi)
    return lat, lon

def WGS84ToITRF (lat, lon, h): # WGS-84 to ITRF
	SINK = math.sin(lat)
	COSK = math.cos(lat)
	e2 = 2.0 * f - f * f
	v = sm_a / math.sqrt(1.0 - e2 * SINK * SINK)
	x = (v + h) * COSK * math.cos(lon)
	y = (v + h) * COSK * math.sin(lon)
	z = ((1 - e2) * v + h) * SINK
	return x, y, z

def ITRFToWGS84 (x, y, z):
	e2 = 2.0 * f - f * f
	E = e2 / (1.0 - e2)
	b = sm_a * (1.0 - f)
	p = math.sqrt(x * x + y * y)
	q = math.atan2(z * sm_a, (p * b))
	lat = math.atan2((z + E * b * math.sin(q) * math.sin(q) * math.sin(q)), (p - e2 * sm_a * math.cos(q) * math.cos(q) * math.cos(q)))
	v = sm_a / math.sqrt(1.0 - e2 * math.sin(lat) * math.sin(lat))
	lon = math.atan2(y, x)
	h = (p / math.cos(lat)) - v
	return lat, lon, h
	
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
	def __init__(self, config, configType, fName, zone, hemisphere):
		self.config = config
		self.antennas = {}
		nhead = 2
		csvreader = csv.reader(open(fName, "rU"), dialect="excel")
		for row in csvreader:
			if nhead > 0:
				nhead -= 1
				continue
			if len(row) < 19 or len(row[0]) == 0:
				continue
			pad = int(row[0])
			name = row[2]
			if configType == "nominal":
				self.antennas[pad] = Antenna(pad, name, float(row[4]), float(row[5]), -1.0e9, zone, hemisphere)
			elif configType == "proposed":
				self.antennas[pad] = Antenna(pad, name, float(row[10]), float(row[11]), -1.0e9, zone, hemisphere)
			else:
				if len(row[17]) > 0:
					self.antennas[pad] = Antenna(pad, name, float(row[17]), float(row[18]), float(row[19]), zone, hemisphere)
				else:
					self.antennas[pad] = Antenna(pad, name, float(row[10]), float(row[11]), -1.0e9, zone, hemisphere)

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

class AntennaConfig:
	def __init__(self, name, antennas):
		self.name = name
		self.antennas = antennas

def usage(config):
	k = config.keys()
	k.sort()
	print "Usage:\n%s array type [format]\n" %(sys.argv[0])
	print "  array = %s" %(" | ".join(k))
	print "         This specifies the array or subset of the array to use.\n"
	print "  type = nominal | proposed | installed"
	print "         This specifies the type of antenna position to use. Nominal positions"
	print "         are base on the original ASKAP configuration document. Proposed"
	print "         positions include changes based on other considerations (geographical"
	print "         and rare species protection). Installed positions include the proposed"
	print "         positions updated with installed antenna postions where available.\n"
	print "  format = parset | calc | csv"

config = {}
config["A27CR3P6B"] = AntennaConfig("A27CR3P6B", range(1, 37))
config["A27CR3"] = AntennaConfig("A27CR3", range(1, 31))
config["BETA"] = AntennaConfig("BETA", (1, 3, 6, 8, 9, 29))
config["BETA2"] = AntennaConfig("BETA2", (1, 3, 6, 8, 9, 2))
config["BETA4"] = AntennaConfig("BETA4", (1, 3, 6, 8, 9, 4))
config["BETA5"] = AntennaConfig("BETA5", (1, 3, 6, 8, 9, 5))
config["BETA7"] = AntennaConfig("BETA7", (1, 3, 6, 8, 9, 7))
config["BETA10"] = AntennaConfig("BETA10", (1, 3, 6, 8, 9, 10))
config["BETA11"] = AntennaConfig("BETA11", (1, 3, 6, 8, 9, 11))
config["BETA12"] = AntennaConfig("BETA12", (1, 3, 6, 8, 9, 12))
config["BETA13"] = AntennaConfig("BETA13", (1, 3, 6, 8, 9, 13))
config["BETA14"] = AntennaConfig("BETA14", (1, 3, 6, 8, 9, 14))
config["BETA15"] = AntennaConfig("BETA15", (1, 3, 6, 8, 9, 15))
config["BETA16"] = AntennaConfig("BETA16", (1, 3, 6, 8, 9, 16))
config["BETA17"] = AntennaConfig("BETA17", (1, 3, 6, 8, 9, 17))
config["BETA18"] = AntennaConfig("BETA18", (1, 3, 6, 8, 9, 18))
config["BETA19"] = AntennaConfig("BETA19", (1, 3, 6, 8, 9, 19))
config["BETA20"] = AntennaConfig("BETA20", (1, 3, 6, 8, 9, 20))
config["BETA21"] = AntennaConfig("BETA21", (1, 3, 6, 8, 9, 21))
config["BETA22"] = AntennaConfig("BETA22", (1, 3, 6, 8, 9, 22))
config["BETA23"] = AntennaConfig("BETA23", (1, 3, 6, 8, 9, 23))
config["BETA24"] = AntennaConfig("BETA24", (1, 3, 6, 8, 9, 24))
config["BETA25"] = AntennaConfig("BETA25", (1, 3, 6, 8, 9, 25))
config["BETA26"] = AntennaConfig("BETA26", (1, 3, 6, 8, 9, 26))
config["BETA27"] = AntennaConfig("BETA27", (1, 3, 6, 8, 9, 27))
config["BETA28"] = AntennaConfig("BETA28", (1, 3, 6, 8, 9, 28))
config["BETA29"] = AntennaConfig("BETA29", (1, 3, 6, 8, 9, 29))
config["BETA29A"] = AntennaConfig("BETA29A", (1, 3, 6, 8, 11, 29))
config["BETA29B"] = AntennaConfig("BETA29B", (1, 3, 8, 9, 11, 29))
config["BETA30"] = AntennaConfig("BETA30", (1, 3, 6, 8, 9, 30))
config["BETA31"] = AntennaConfig("BETA31", (1, 3, 6, 8, 9, 31))
config["BETA32"] = AntennaConfig("BETA32", (1, 3, 6, 8, 9, 32))
config["BETA33"] = AntennaConfig("BETA33", (1, 3, 6, 8, 9, 33))
config["BETA34"] = AntennaConfig("BETA34", (1, 3, 6, 8, 9, 34))
config["BETA35"] = AntennaConfig("BETA35", (1, 3, 6, 8, 9, 35))
config["BETA36"] = AntennaConfig("BETA36", (1, 3, 6, 8, 9, 36))

configType = ("nominal", "proposed", "installed")
outputType = ("parset", "calc", "csv")
if len(sys.argv) < 3:
	usage(config)
	sys.exit(1)

at = string.upper(sys.argv[1])
ct = string.lower(sys.argv[2])
if not at in config or not ct in configType:
	usage(config)
	sys.exit(1)

ot = "parset"
if len(sys.argv) == 4:
	ot = string.lower(sys.argv[3])
	if not ot in outputType:
		usage(config)
		sys.exit(1)

antennas = AntennaList(config[string.upper(sys.argv[1])], string.lower(sys.argv[2]), "ASKAP Antenna Locations Master File.csv", 50, "south")

offxyz = [-2556743.707 - -2556745.438, 5097440.315 - 5097448.114, -2847749.657 - -2847753.833]

if ot == "parset":
	antennas.dump()
if ot == "calc":
	antennas.dumpcalc(3, 0.0, offxyz)
if ot == "csv":
	antennas.dumplatlong()
