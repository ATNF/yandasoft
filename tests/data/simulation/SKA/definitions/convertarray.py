#csvfile='aa.csv'
#outfile='SKA1AA2km.in'
csvfile='dish.csv'
outfile='SKA1Dish.in'
dmax=100000.0
lonref=40.019
latref=-30.7131

import math
import csv

sm_a = 6378137.0
invf = 298.257223563
f = 1.0 / invf

# Convert WGS-84 to ITRF
# lat and lon are the latitude and longitude in radians, h is the height in metres.
def WGS84ToITRF (lat, lon, h):
    SINK = math.sin(lat)
    COSK = math.cos(lat)
    e2 = 2.0 * f - f * f
    v = sm_a / math.sqrt(1.0 - e2 * SINK * SINK)
    x = (v + h) * COSK * math.cos(lon)
    y = (v + h) * COSK * math.sin(lon)
    z = ((1 - e2) * v + h) * SINK
    return x, y, z

def distance(lat, lon, latref, longref):
    return math.sqrt((lat-latref)*(lat-latref)+(lon-lonref)*(lon-lonref))*17687.5

ants=[]
lons=[]
lats=[]

fileid=file(csvfile, 'U')
ant=0
while(True):
    line=fileid.readline()
    if(line==''):
        break;
    line=line.split(',')
#    ant=int(line[0])
    ant=ant+1
    lon=float(line[0])
    lat=float(line[1])
    if(distance(lat, lon, latref, lonref)<dmax):
        ants.append(ant)
        lons.append(lon)
        lats.append(lat)
fileid.close()

print ants

s='antennas.SKA1.names = ['
for i in range(len(ants)):
    if(i<len(ants)-1):
        s=s+'Dish%d,'%ants[i]
    else:
        s=s+'Dish%d'%ants[i]

s=s+']\n'

print "Read %d lines" % len(ants)

outfileid=file(outfile, 'w')
outfileid.write('antennas.telescope = SKA1\n')
outfileid.write('antennas.SKA1.coordinates = global\n')
outfileid.write(s)
outfileid.write('antennas.SKA1.diameter = 80m\n')
outfileid.write('antennas.SKA1.scale = 1.0\n')
outfileid.write('antennas.SKA1.mount = alt-az\n')


for ant in range(len(ants)):
    lon=lons[ant]
    lat=lats[ant]
    el=300.0
    outfileid.write('# lat: %s; long: %s; el: %s\n' % (lat, lon, el))
    (xx, yy, zz) = WGS84ToITRF(lat*math.pi/180.0, lon*math.pi/180.0, el)
    outfileid.write('antennas.SKA1.Dish%d=[%s, %s, %s]\n' % (ants[ant], xx, yy, zz))

fileid.close()
outfileid.close()
