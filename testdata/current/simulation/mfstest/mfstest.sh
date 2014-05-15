#!/bin/bash -l
#export AIPSPATH=${ASKAP_ROOT}/Code/Components/Synthesis/testdata/current

HOSTNAME=`hostname -s`

echo "This is the ASKAPsoft mfstest. It will run for about 3 minutes."  | tee mfstest.$HOSTNAME.out
echo "Started " `date` " on " $HOSTNAME  | tee -a mfstest.$HOSTNAME.out

echo "Extracting Urvashi's test measurement set" | tee -a mfstest.$HOSTNAME.out
tar zxvf ${ASKAP_ROOT}/Code/Components/Synthesis/testdata/current/simulation/mfstest/mfsdata.tgz | tee -a mfstest.$HOSTNAME.out

cat > mfstest.clean.in <<EOF
Cimager.dataset                                 = ptest.ms
Cimager.dataset0                                 = ptest.ms

Cimager.visweights		= MFS
#Cimager.visweights.MFS.reffreq	= 1.35e9
Cimager.visweights.MFS.reffreq	= 1.2e9

Cimager.Images.reuse			= false
Cimager.Images.writeAtMajorCycle	= false

#Cimager.Images.Names                            = [image.taylor.0, image.taylor.1, image.taylor.2]
Cimager.Images.Names                            = [image]
Cimager.Images.shape				= [1024,1024]
Cimager.Images.cellsize	        		= [8.0arcsec, 8.0arcsec]


Cimager.Images.image.frequency	= [1.40e9,1.40e9]
Cimager.Images.image.nchan	= 1
Cimager.Images.image.nterms     = 2
Cimager.Images.image.polarisation     = ["RRLL"]
Cimager.Images.image.direction	= [19h59m28.5000, +40.44.01.5000, J2000]   

#
Cimager.gridder                          	= SphFunc
#
# Use a multiscale Clean solver
#
Cimager.solver                           	= Clean
Cimager.solver.Clean.algorithm                 	= MSMFS

Cimager.solver.Clean.scales			= [0]
Cimager.solver.Clean.niter                     	= 100
Cimager.solver.Clean.gain                      	= 0.5
Cimager.solver.Clean.tolerance                  = 0.1
Cimager.solver.Clean.verbose                   	= True
Cimager.threshold.minorcycle                    = [0.05mJy, 10%]
Cimager.threshold.majorcycle                    = 0.045mJy
# 
Cimager.ncycles                                 = 2
#
# Restore the image at the end
#
Cimager.restore                          	= True
Cimager.restore.beam                     	= [28arcsec, 28arcsec, 0deg]
#
# Use preconditioning for deconvolution
#
Cimager.preconditioner.Names			= [None]

EOF
echo "Running cimager to form MSMFS images" | tee -a  mfstest.$HOSTNAME.out
${ASKAP_ROOT}/Code/Components/Synthesis/synthesis/current/apps/cimager.sh -c mfstest.clean.in | tee -a mfstest.$HOSTNAME.out
echo "Ended " `date` " on " $HOSTNAME  | tee -a mfstest.$HOSTNAME.out
