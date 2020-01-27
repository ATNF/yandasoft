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
imcontsub
imgNaN2zero
imgscale
imgstat
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
	add_executable (${prog} 
            ${prog}.cc
            $<TARGET_OBJECTS:deconvolution>
  	    $<TARGET_OBJECTS:distributedimager>
	    $<TARGET_OBJECTS:gridding>
	    $<TARGET_OBJECTS:measurementequation>
	    $<TARGET_OBJECTS:messages>
	    $<TARGET_OBJECTS:opcal>
	    $<TARGET_OBJECTS:parallel>
	    $<TARGET_OBJECTS:simulation>
	    $<TARGET_OBJECTS:utils>
        )
        if(MPI_COMPILE_FLAGS)
  		set_target_properties(${prog} PROPERTIES
    		COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
	endif()

	if(MPI_LINK_FLAGS)
  		set_target_properties( ${prog} PROPERTIES
    		LINK_FLAGS "${MPI_LINK_FLAGS}")
	endif()

        target_link_libraries (${prog} 
        ${ASKAP_LIBRARY}
	${log4cxx_LIBRARY}
	${CASACORE_LIBRARIES}
	${FFTW_LIBRARIES}
	${SCIMATH_LIBRARY}
	${IMAGEMATH_LIBRARY}
	${LofarBlob_LIBRARY}
	${LofarCommon_LIBRARY}
	${Boost_LIBRARIES}
	${ACCESSORS_LIBRARY}
        ${XercesC_LIBRARY}
	${COMPONENTS_LIBRARY}
	${ASKAPPARALLEL_LIBRARY}
        ${GSL_LIBRARIES}
        ${CPPUNIT_LIBRARY})
if (MPI_FOUND)
	target_link_libraries(${prog}
		${MPI_LIBRARIES}
	)
endif (MPI_FOUND)


	install (TARGETS ${prog} DESTINATION bin)

endforeach(prog ${_progs})
