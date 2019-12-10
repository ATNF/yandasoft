add_library( utils OBJECT 
	CommandLineParser.cc
	LinmosUtils.cc
)

target_link_libraries(utils
	lofar::Common
	askap::scimath
	askap::accessors
)

set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)
if (MPI_COMPILE_FLAGS)
	set_property(TARGET utils PROPERTY COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
endif()


install (FILES
	CommandLineParser.h
	LinmosUtils.h
DESTINATION include/askap/utils
)
