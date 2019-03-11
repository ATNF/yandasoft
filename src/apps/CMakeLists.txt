set (_progs
cadvise
cbpcalibrator
ccalapply
ccalibrator
ccontsubtract
cdeconvolver
cimager
cmpgains
correctdelay
csimulator
cubemerge
cubeslice
cubetest
gridderExperiments
illumextractor
imager
imgNaN2zero
imgscale
imgstat
linmos-mpi-old
linmos-mpi
linmos
msconcat
mslist
msmerge
mstabtest
noiseadder
#opcal ONLY BECAUSE IT IS THE SAME NAME AS A DIRECTORY UPSTREAM
randomgains
smear
tCFGeneration
tConvolve
tConvolveBLAS
tConvolveCASA
tConvolveResid
tGridding
tParallelIterator
tPreconditioning
tSaxpy
tWienerdbg
tcmult
testlinmos
tmerge
tstlayout
tstuvw
tstuvw2
)

foreach(prog ${_progs})
	add_executable (${prog} ${prog}.cc)
	target_link_libraries (${prog} yanda_synthesis ${GSL_LIBRARIES})
	install (TARGETS ${prog} DESTINATION bin)
endforeach(prog ${_progs})
