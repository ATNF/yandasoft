#!/bin/bash

cd `dirname $0`

INITIALDIR=`pwd`
echo Running test cases...

FAIL=0

# Non-cycling (ncycles 0) Spectral Line Imager Test
cd spectralline
./run.sh
if [ $? -eq 0 ]; then
    R1="spectralline (ncycles 0)     PASS"
else
    R1="spectralline (ncycles 0)     FAIL"
    FAIL=1
fi
rm *.fits

cd $INITIALDIR
cd continuum
./run.sh
if [ $? -eq 0 ]; then
    R2="continuum (ncycles 0)        PASS"
else
    R2="continuum (ncycles 0)        FAIL"
    FAIL=1
fi
rm *.fits

cd $INITIALDIR
cd polarisation
./run.sh
if [ $? -eq 0 ]; then
    R3="polarisation                 PASS"
else
    R3="polarisation                 FAIL"
    FAIL=1
fi
rm *.fits

cd $INITIALDIR
cd polarisation
./run-new.sh
if [ $? -eq 0 ]; then
    R4="polarisation 2               PASS"
else
    R4="polarisation 2               FAIL"
    FAIL=1
fi
rm *.fits

cd $INITIALDIR

# Print Results
echo
echo Result Summary:
echo ============================
echo $R1
echo $R2
echo $R3
echo $R4

if [ $FAIL -eq 0 ]; then
    exit 0
else
    exit 1
fi
