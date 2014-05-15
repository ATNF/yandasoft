#!/bin/bash -l
export AIPSPATH=${ASKAP_ROOT}/Code/Base/accessors/current

HOSTNAME=`hostname -s`

echo "This is the ASKAPsoft stdtest. It will run for about 40-50 minutes."  | tee stdtest.$HOSTNAME.out
echo "Started " `date` " on " $HOSTNAME  | tee -a stdtest.$HOSTNAME.out

echo "Extracting 10uJy model" | tee -a stdtest.$HOSTNAME.out
tar zxvf ${ASKAP_ROOT}/Code/Components/Synthesis/testdata/current/simulation/models/10uJy.model.small.tgz | tee -a stdtest.$HOSTNAME.out
ERR=${PIPESTATUS[0]}
if [ $ERR -ne 0 ]; then
    exit $ERR
fi

cat > stdtest.simulator.in <<EOF 
Csimulator.dataset                              =       10uJy_stdtest.ms
#
# The name of the model source is 10uJy. Specify direction and model file
#
Csimulator.sources.names                        =       [10uJy]
Csimulator.sources.10uJy.direction              =       [12h30m00.000, -45.00.00.000, J2000]
Csimulator.sources.10uJy.model                  =       10uJy.model.small
#
# Define the antenna locations, feed locations, and spectral window definitions
#
Csimulator.antennas.definition                  =       definitions/ASKAP45.in
Csimulator.feeds.definition                     =       definitions/ASKAP1feeds.in
Csimulator.spws.definition                      =       definitions/ASKAPspws.in
#
# Standard settings for the simulaton step
#
Csimulator.simulation.blockage                  =       0.01
Csimulator.simulation.elevationlimit            =       8deg
Csimulator.simulation.autocorrwt                =       0.0
Csimulator.simulation.usehourangles             =       True
Csimulator.simulation.referencetime             =       [2007Mar07, UTC]
#
# Undersample in time by ~ 10 to make the processing run quickly
#
Csimulator.simulation.integrationtime           =       1500s
#
# Observe source 10uJy for 12 hours with a single channel spectral window
#
Csimulator.observe.number                       =       1
Csimulator.observe.scan0                        =       [10uJy, Continuum0, -6h, 6h]
#
# Use a gridder to apply primary beam during the W projection step.
#
Csimulator.gridder                              = AWProject
Csimulator.gridder.AWProject.wmax            	= 15000
Csimulator.gridder.AWProject.nwplanes        	= 129
Csimulator.gridder.AWProject.oversample     	= 8
Csimulator.gridder.AWProject.diameter		    = 12m
Csimulator.gridder.AWProject.blockage		    = 2m
Csimulator.gridder.AWProject.maxfeeds		    = 1
Csimulator.gridder.AWProject.maxsupport		    = 2048
Csimulator.gridder.AWProject.frequencydependent = false
Csimulator.gridder.AWProject.variablesupport    = true
Csimulator.gridder.AWProject.offsetsupport      = true
Csimulator.gridder.AWProject.tablename		    = AWProject.tab
EOF
echo "Running csimulator to create MeasurementSet for a single pointing" | tee -a  stdtest.$HOSTNAME.out
${ASKAP_ROOT}/Code/Components/Synthesis/synthesis/current/apps/csimulator.sh -c stdtest.simulator.in | tee -a stdtest.$HOSTNAME.out
ERR=${PIPESTATUS[0]}
if [ $ERR -ne 0 ]; then
    exit $ERR
fi

cat > stdtest.dirty.in <<EOF
Cimager.dataset                                 = 10uJy_stdtest.ms

Cimager.Images.Names                            = [image.i.10uJy_dirty_stdtest]
Cimager.Images.shape				= [2048,2048]
Cimager.Images.cellsize	        		= [6.0arcsec, 6.0arcsec]
Cimager.Images.image.i.10uJy_dirty_stdtest.frequency	= [1.420e9,1.420e9]
Cimager.Images.image.i.10uJy_dirty_stdtest.nchan	= 1
Cimager.Images.image.i.10uJy_dirty_stdtest.direction	= [12h30m00.00, -45.00.00.00, J2000]
#
Cimager.gridder                          	= AProjectWStack
Cimager.gridder.AProjectWStack.wmax            	= 15000
Cimager.gridder.AProjectWStack.nwplanes        	= 33
Cimager.gridder.AProjectWStack.oversample     	= 4
Cimager.gridder.AProjectWStack.diameter		= 12m
Cimager.gridder.AProjectWStack.blockage		= 2m
Cimager.gridder.AProjectWStack.maxfeeds		= 1
Cimager.gridder.AProjectWStack.maxsupport       = 512
Cimager.gridder.AProjectWStack.frequencydependent = false
Cimager.gridder.AProjectWStack.variablesupport    = true
Cimager.gridder.AProjectWStack.offsetsupport      = true
Cimager.gridder.AProjectWStack.tablename = AProjectWStack.tab
#
Cimager.solver                           	= Dirty
Cimager.solver.Dirty.tolerance                  = 0.1
Cimager.solver.Dirty.verbose                   	= True
Cimager.ncycles                                 = 0

Cimager.preconditioner.Names			= None

Cimager.restore                          	= True
Cimager.restore.beam                     	= [30arcsec, 30arcsec, 0deg]

EOF
echo "Running cimager to form Dirty image of single pointing" | tee -a  stdtest.$HOSTNAME.out
${ASKAP_ROOT}/Code/Components/Synthesis/synthesis/current/apps/cimager.sh -c stdtest.dirty.in | tee -a stdtest.$HOSTNAME.out
ERR=${PIPESTATUS[0]}
if [ $ERR -ne 0 ]; then
    exit $ERR
fi

cat > stdtest.clean.in <<EOF
Cimager.dataset                                 = 10uJy_stdtest.ms

Cimager.Images.Names                            = [image.i.10uJy_clean_stdtest]
Cimager.Images.shape				= [2048,2048]
Cimager.Images.cellsize	        		= [6.0arcsec, 6.0arcsec]
Cimager.Images.image.i.10uJy_clean_stdtest.frequency	= [1.420e9,1.420e9]
Cimager.Images.image.i.10uJy_clean_stdtest.nchan	= 1
Cimager.Images.image.i.10uJy_clean_stdtest.direction	= [12h30m00.00, -45.00.00.00, J2000]
#
Cimager.gridder                          	= WProject
Cimager.gridder.WProject.wmax            	= 15000
Cimager.gridder.WProject.nwplanes        	= 129
Cimager.gridder.WProject.oversample     	= 8
Cimager.gridder.WProject.maxsupport         = 1024
Cimager.gridder.WProject.variablesupport    = true
Cimager.gridder.WProject.offsetsupport      = true
Cimager.gridder.WProject.tablename          = WProject.tab
#
# Use a multiscale Clean solver
#
Cimager.solver                           	    = Clean
Cimager.solver.Clean.algorithm                 	= MultiScale
Cimager.solver.Clean.scales			            = [0, 3, 10]
Cimager.solver.Clean.niter                     	= 10000
Cimager.solver.Clean.gain                      	= 0.1
Cimager.solver.Clean.tolerance                  = 0.1
Cimager.solver.Clean.verbose                   	= True
Cimager.threshold.minorcycle                    = [0.27mJy, 10%]
Cimager.threshold.majorcycle                    = 0.3mJy
# 
Cimager.ncycles                                 = 10
#
# Restore the image at the end
#
Cimager.restore                          	= True
Cimager.restore.beam                     	= [30arcsec, 30arcsec, 0deg]
#
# Use preconditioning for deconvolution
#
Cimager.preconditioner.Names			    = [Wiener, GaussianTaper]
Cimager.preconditioner.Wiener.noisepower	= 100.0
Cimager.preconditioner.GaussianTaper		= [20arcsec, 20arcsec, 0deg]

EOF
echo "Running cimager to form Clean image of single pointing" | tee -a  stdtest.$HOSTNAME.out
${ASKAP_ROOT}/Code/Components/Synthesis/synthesis/current/apps/cimager.sh -c stdtest.clean.in | tee -a stdtest.$HOSTNAME.out
ERR=${PIPESTATUS[0]}
echo "Ended " `date` " on " $HOSTNAME  | tee -a stdtest.$HOSTNAME.out
if [ $ERR -ne 0 ]; then
        exit $ERR
fi
