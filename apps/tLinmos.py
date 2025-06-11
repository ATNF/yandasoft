#!/usr/bin/env python
# coding: utf-8
# # tLinmos - linmos test script
import os
from subprocess import call
import argparse as ap


def argInit():
  parser = ap.ArgumentParser(prog='tLinmos.py',formatter_class=ap.ArgumentDefaultsHelpFormatter,
                                 description='Test linmos speed')
  parser.add_argument('-i','--imsize',type=int,default=2048,help="Image size for beam images")
  parser.add_argument('-n','--nchan',type=int,default=16,help="Number of channels")
  parser.add_argument('-p','--nchanpercore',type=int,default=1, help="Number of channels to process per core")
  parser.add_argument('-N','--nodes',type=int,default=1, help="Number of nodes to use")
  parser.add_argument('-f','--footprint',default=None,help="Footprint file (use builtin if None)")
  parser.add_argument('-b','--sbatch',action='store_true', help="Submit as sbatch job")
  parser.add_argument('-c','--clean',action='store_true', help="Clean up work directory afterwards")
  parser.add_argument('-P','--parallel',action='store_true', help="Write output in parallel (default=serial)")
  parser.add_argument('-C','--collective',action='store_true', help="Write output using collective IO (default=individual)")
  parser.add_argument('-t','--time',default='1:0:0', help="sbatch max run time")

  return parser


try:
  args, unknown = argInit().parse_known_args()
except SystemExit:
  exit(1)

imsize = args.imsize
nchan = args.nchan
footprintFile = args.footprint
sbatch = args.sbatch
time = args.time
nchanpercore = args.nchanpercore
nodes = args.nodes
clean = args.clean
parallel = args.parallel
collective = args.collective
if (collective):
    parallel = True
    print("Setting parallel true, since collective is selected")
print("tLinmos called with imsize=",imsize,", nchan=",nchan,", nchanpercore=",nchanpercore,", nodes=",nodes,
 ", clean=",clean,", parallel=",parallel,", collective=",collective)

ncores = 6 # Mac
if sbatch:
    ncores = 128
else:
    nodes = 1
ntasks = nchan // nchanpercore
ntaskspernode = ntasks // nodes
ompnumthreads = ncores*nodes//ntasks # needs to be <= 128*nodes/ntasks
if ompnumthreads < 1:
    ompnumthreads = 1
if nchan != ntasks * nchanpercore:
    print("nchanpercore should divide evenly into nchan")
if ntaskspernode*ompnumthreads > ncores:
    print("not enough cores for requested nchanpercore and nchan")
    exit(1)

modules = "module load singularity/3.11.4-askap askapsoft/1.16.1-lustrefix askappy/2.7"
name = "tLinmos"
runid = "imsize%s_nchan%s_n%s_N%s_c%s_%s_%s" % (imsize,nchan,ntasks,nodes,ompnumthreads,"SP"[parallel],"IC"[collective])
workDir = "./tLinmos-"+runid # directory to use
batch = "#!/bin/bash\n"
if (sbatch):
    batch += """#SBATCH --time=%s
#SBATCH --partition=askaprt
#SBATCH --threads-per-core=1
#SBATCH --nodes=%i
#SBATCH --mem=0
#SBATCH --job-name=%s
#SBATCH --export=NONE
""" % (time,nodes,name)
    batch += "#SBATCH --output=tLinmos.%j.out\n\n"
    batch += modules+"\n"
    myRunner="srun --export=ALL --ntasks=%i --nodes=%i --cpus-per-task=%i "%(ntasks,nodes,ompnumthreads)
else:
    batch += "export TMPDIR=/tmp\n"
    myRunner="mpirun -np %i" % ntasks
if ompnumthreads>0:
    batch += "export OMP_NUM_THREADS=%d\n" % ompnumthreads

workDir = os.path.expandvars(workDir)
os.makedirs(workDir, exist_ok=True)
os.chdir(workDir)



# run tImageWrite to create one 'beam' image
with open("tImageWrite.in","w") as f:
    f.write("""imagetype=fits
name=image.beam00
spectral_first=false
stokes=[V]
size=%s
nchan=%s
imageaccess=collective
imageaccess.axis=3
imageaccess.write=parallel""" % (imsize,nchan))

batch += "for i in {1..35}; do\n"
batch += "   %s tImageWrite -c tImageWrite.in > tImageWrite_$i.log\n" % myRunner
batch += "   mv image.beam00.fits \"image.beam$(printf \"%02d\" $i).fits\"\n"
batch += "done\n"
batch += "%s tImageWrite -c tImageWrite.in > tImageWrite_0.log\n" % myRunner
# batch += """
# for i in {1..35}; do
#   %s tImageWrite -c tImageWrite.in > tImageWrite_$i.log
#   mv image.beam00.fits "image.beam$(printf "\%02d" $i).fits"
# done
# %s tImageWrite -c tImageWrite.in > tImageWrite_0.log
# """ % (myRunner, myRunner)

with open("createImages.py","w") as f:
    f.write("""#!/usr/bin/env python
# coding: utf-8
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from subprocess import call

footprintFile="%s"
imsize=%s
imageaccess="%s"
writemode="%s"
"""  % (footprintFile, imsize, ["individual","collective"][collective], ["serial","parallel"][parallel]))

    f.write("""
#get beam coordinates from footprint file
def decodeFootprint(footprintFile):
    with open(footprintFile) as f:
        beams=[]
        for line in f:
            radec = line[21:-1].split(",")
            beams.append(radec)
    return beams

if footprintFile == "None":
    footprintFile = "footprint.txt"
    with open(footprintFile,"w") as f:
        f.write(\"\"\" 0  (-2.467 -1.959)  05:06:21.014,-71:55:04.22
 1  (-1.567 -1.955)  05:17:55.442,-72:00:37.76
 2  (-0.667 -1.951)  05:29:34.872,-72:03:37.19
 3  ( 0.233 -1.948)  05:41:16.296,-72:04:01.31
 4  ( 1.133 -1.944)  05:52:56.645,-72:01:50.05
 5  ( 2.033 -1.940)  06:04:32.882,-71:57:04.21
 6  (-2.020 -1.177)  05:13:07.714,-71:11:37.90
 7  (-1.120 -1.174)  05:24:17.285,-71:15:44.21
 8  (-0.220 -1.170)  05:35:29.938,-71:17:21.91
 9  ( 0.680 -1.166)  05:46:42.991,-71:16:30.47
10  ( 1.580 -1.163)  05:57:53.760,-71:13:10.20
11  ( 2.480 -1.159)  06:08:59.614,-71:07:22.22
12  (-2.473 -0.400)  05:08:42.818,-70:22:13.48
13  (-1.573 -0.396)  05:19:24.665,-70:27:21.35
14  (-0.673 -0.392)  05:30:10.387,-70:30:06.84
15  ( 0.227 -0.389)  05:40:57.650,-70:30:29.09
16  ( 1.127 -0.385)  05:51:44.078,-70:28:28.02
17  ( 2.027 -0.381)  06:02:27.319,-70:24:04.21
18  (-2.027  0.381)  05:14:54.187,-69:38:31.52
19  (-1.127  0.385)  05:25:14.731,-69:42:19.66
20  (-0.227  0.389)  05:35:37.685,-69:43:50.12
21  ( 0.673  0.392)  05:46:00.948,-69:43:02.50
22  ( 1.573  0.396)  05:56:22.423,-69:39:57.02
23  ( 2.473  0.400)  06:06:40.044,-69:34:34.57
24  (-2.480  1.159)  05:10:43.411,-68:49:17.04
25  (-1.580  1.163)  05:20:40.421,-68:54:03.31
26  (-0.680  1.166)  05:30:40.495,-68:56:37.10
27  ( 0.220  1.170)  05:40:41.779,-68:56:57.77
28  ( 1.120  1.174)  05:50:42.406,-68:55:05.27
29  ( 2.020  1.177)  06:00:40.510,-68:51:00.04
30  (-2.033  1.940)  05:16:25.224,-68:05:20.18
31  (-1.133  1.944)  05:26:03.785,-68:08:52.94
32  (-0.233  1.948)  05:35:44.261,-68:10:17.29
33  ( 0.667  1.951)  05:45:24.984,-68:09:32.90
34  ( 1.567  1.955)  05:55:04.282,-68:06:39.96
35  ( 2.467  1.959)  06:04:40.502,-68:01:39.14
\"\"\")

beams = decodeFootprint(footprintFile)


i=0
files=[]
for beam in beams:
    file = "image.beam%2.2d.fits"%i
    print(file," ",)
    files.append(file[:-5])
    c = SkyCoord(beam[0]+' '+beam[1], unit=(u.hourangle, u.deg))
    #print("%.16E" %c.ra.deg,c.dec.deg)
    with fits.open(file, 'update') as hdul:
        header = hdul[0].header
        #print(repr(header))
        header['CRVAL1']=c.ra.deg
        header['CRVAL2']=c.dec.deg
        header['CDELT1']=header['CDELT1']/imsize*2048*5/3 # 6" pixels for 2048 image
        header['CDELT2']=header['CDELT2']/imsize*2048*5/3
    i = i + 1

# linmos files together
with open("linmos.in","w") as f:
    f.write(\"\"\"linmos.imagetype=fits
linmos.names=[%s]
linmos.outname=image.allbeams.linmos
linmos.outweight=weights.allbeams.linmos
linmos.weighttype=FromPrimaryBeamModel
linmos.weightstate=Inherent
linmos.cutoff=0.1
linmos.primarybeam=GaussianPB
linmos.imageaccess=%s
linmos.imageaccess.axis=3
linmos.imageaccess.write=%s
\"\"\" % (", ".join(files),imageaccess,writemode))
""")

batch += """
python createImages.py

%s linmos-mpi -c linmos.in > ../tLinmos_%s.log
""" % (myRunner,runid)

if clean:
    batch += """#cleanup
cd ..
rm -rf tLinmos-%s
""" % (runid)

#write out batch script
batch_name = name + ".sbatch"
with open(batch_name,"w") as b:
    b.write(batch)

if sbatch:
    # need to use native python module on setonix for sbatch command to work
    call(["sbatch",batch_name])
else:
    call(["/bin/bash","tLinmos.sbatch"])
