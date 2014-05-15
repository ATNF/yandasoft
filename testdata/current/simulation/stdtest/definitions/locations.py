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
