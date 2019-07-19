import os, sys


# helper method to do one simulation
def getResidualW(dec,ha):
  """
     dec - declination (degrees)
     ha - hour angle (hours)

     Return: A tuple with residual w-term in wavelengths and start elevation. 
  """

  print "Simulating dec=%f deg ha=%f hours" % (dec,ha)
  os.system("cp -f template_parset.in tmpparset.in")
  f = open("tmpparset.in","a")
  try:
    f.write("Csimulator.sources.field1.direction = [12h00m00.000, %fdeg, J2000]\n" % (dec,))
    f.write("tangent = [12h00m00.000, %fdeg, J2000]\n" % (dec,))
    half_dur = 2.5/60.
    f.write("Csimulator.observe.scan0 =[field1, Wide0, %fh, %fh]\n" % (ha-half_dur,ha+half_dur))
  finally:
    f.close()
  os.system("/Users/vor010/ASKAP/ASKAPsoft/Code/Components/Synthesis/synthesis/current/apps/csimulator.sh -inputs tmpparset.in > csim.out")
  elevbuf = []
  f = open("csim.out")
  try:
    for line in f:
        if 'el   =' in line:
           parts= line.split()
           if len(parts)!= 11:
              raise RuntimeError, "Expect 11 elements, you have %s" % (parts,)
           if parts[0]!='INFO' or parts[10]!='deg' or parts[7]!='el' or parts[8]!='=':
              raise RuntimeError, "Parse error for %s" % (parts,)
           elevbuf.append(float(parts[9]))
  finally:
    f.close()
  if len(elevbuf)!=2:
     raise RuntimeError, "Expect start and stop elevations in csim.out, you have %s" % (elevbuf,)
  elev = (elevbuf[0]+elevbuf[1])/2.
  # now execute advise utility
  os.system("/Users/vor010/ASKAP/ASKAPsoft/Code/Components/Synthesis/synthesis/current/apps/cadvise.sh -inputs tmpparset.in > cadv.out")
  residW = None
  f = open("cadv.out")
  try:
    for line in f:
        if 'Largest residual W:' in line:
           parts= line.split()
           if len(parts)!= 11:
              raise RuntimeError, "Expect 11 elements, you have %s" % (parts,)
           if parts[2]!='INFO' or parts[10]!='wavelengths' or parts[6]!='Largest' or parts[7]!='residual' or parts[8]!='W:':
              raise RuntimeError, "Parse error for %s" % (parts,)
           if residW != None:
              raise RuntimeError, "More than one residual W reported in the output cadv.out"
           residW = float(parts[9])
           
  finally:
    f.close()

  return (residW,elev)

#############

start_ha = -12
stop_ha = 12
start_dec = -90.
stop_dec = 60.

npoints = 100

f=file("result.dat","w")
try:
   for y in range(npoints):
     dec = start_dec + float(stop_dec - start_dec)/npoints * y
     for x in range(npoints):
        ha = start_ha + float(stop_ha - start_ha)/npoints * x
        res = getResidualW(dec,ha)
        if len(res)!=2:
           raise RuntimeError, "Expect 2-element tuple, you have %s" % (res,)
        if res[1]>0:
           f.write("%f %f %f %f\n" % ((dec,ha) + res))
finally:
  f.close()        

