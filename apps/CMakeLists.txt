set (_progs
delaysolver
extractdata
genoffsets
makescans
bfweights
cadvise
cbpcalibrator
ccalapply
ccalibrator
ccontsubtract
cddcalibrator
cdeconvolver
cdeconvolver-mpi
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
imcontsub
imgNaN2zero
imgscale
imgstat
linmos-mpi
linmos
maskBadChannels
msconcat
mslist
msmerge
mstabtest
noiseadder
opcal 
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
tMultiScaleBasisFunction
tcmult
testlinmos
tmerge
tstlayout
tstuvw
tstuvw2
)



foreach(prog ${_progs})
	add_executable (${prog} ${prog}.cc)
        if(MPI_COMPILE_FLAGS)
  		set_target_properties(${prog} PROPERTIES
    		COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
	endif()

	if(MPI_LINK_FLAGS)
  		set_target_properties( ${prog} PROPERTIES
    		LINK_FLAGS "${MPI_LINK_FLAGS}")
	endif()


	target_link_libraries (${prog}
		askap::yandasoft
	)

	install (TARGETS ${prog} DESTINATION bin)

endforeach(prog ${_progs})
