add_library( simulation OBJECT 
Simulator.cc
)

set_property(TARGET simulation PROPERTY POSITION_INDEPENDENT_CODE ON)
if (MPI_COMPILE_FLAGS)
	set_property(TARGET simulation PROPERTY COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
endif()


install (FILES
Simulator.h
DESTINATION include/askap/simulation
)
