#!/bin/bash

OUTPUT=output.txt

export AIPSPATH=${ASKAP_ROOT}/Code/Base/accessors/current
export IMAGER=../../../yanda/bin/imager
if [ ! -x ${IMAGER} ]; then
    echo $IMAGER does not exit
fi

#
IMAGE=image.IQUV.vela.fits
RESTORED=image.IQUV.vela.restored.fits
PSF=psf.IQUV.vela.fits
RESIDUAL=residual.IQUV.vela.fits
WEIGHTS=weights.IQUV.vela.fits

echo -n "Removing images..."
rm -f *.fits beamlog*

echo Done

echo -n Extracting measurement set...
tar -xvf ../vela.tar
echo Done

mpirun -np 2 $IMAGER -c polim.in | tee $OUTPUT
if [ $? -ne 0 ]; then
    echo Error: mpirun returned an error
    exit 1
fi

echo -n Removing measurement set...
rm -rf vela-60.ms
echo Done

# Check for instances of "Askap error"
grep -c "Askap error" $OUTPUT > /dev/null
if [ $? -ne 1 ]; then
    echo "Askap error reported in output.txt"
    exit 1
fi

# Check for instances of "Unexpected exception"
grep -c "Unexpected exception" $OUTPUT > /dev/null
if [ $? -ne 1 ]; then
    echo "Exception reported in output.txt"
    exit 1
fi
grep -c "BAD TERMINATION" $OUTPUT > /dev/null
if [ $? -ne 1 ]; then
    echo "MPI spotted bad termination (SEGV?)"
    exit 1
fi
# Check for the existance of the various image cubes
if [ ! -f ${IMAGE} ]; then
    echo "Error ${IMAGE} not created"
    exit 1
fi

if [ ! -f ${PSF}${IDX} ]; then
    echo "Error ${PSF} not created"
    exit 1
fi

if [ ! -f ${WEIGHTS}${IDX} ]; then
    echo "Error ${WEIGHTS} not created"
    exit 1
fi

# Check Peak values in dirty images are ok within 0.01
# Note: this only check the absolute values - the log has no info on the sign
PeakIQUV=`grep 'Iteration 0' ${OUTPUT} | head -4 | awk '{ print $12 }'`
PeakI=`echo ${PeakIQUV} | awk '{ print $1 - 0 }'` # - 0 gets rid of comma
PeakQ=`echo ${PeakIQUV} | awk '{ print $2 - 0 }'`
PeakU=`echo ${PeakIQUV} | awk '{ print $3 - 0 }'`
PeakV=`echo ${PeakIQUV} | awk '{ print $4 - 0 }'`
BadI=`echo ${PeakI} | awk '{ print $1 - 2.65 }' | sed 's/-//' | awk '{ print ($1 > 0.01) }'`
BadQ=`echo ${PeakQ} | awk '{ print $1 - 0.28 }' | sed 's/-//' | awk '{ print ($1 > 0.01) }'`
BadU=`echo ${PeakU} | awk '{ print $1 - 1.89 }' | sed 's/-//' | awk '{ print ($1 > 0.01) }'`
BadV=`echo ${PeakV} | awk '{ print $1 - 0.11 }' | sed 's/-//' | awk '{ print ($1 > 0.01) }'`
if [ $BadI -eq 1 ]; then
    echo "Error Stokes I peak value is wrong, $PeakI"
    exit 1
fi
if [ $BadQ -eq 1 ]; then
    echo "Error Stokes Q peak value is wrong, $PeakQ"
    exit 1
fi
if [ $BadU -eq 1 ]; then
    echo "Error Stokes U peak value is wrong, $PeakU"
    exit 1
fi
if [ $BadV -eq 1 ]; then
    echo "Error Stokes V peak value is wrong, $PeakV"
    exit 1
fi
echo "Peak values in dirty images are correct"

echo Done
