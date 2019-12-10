add_library( opcal OBJECT 
BaselineSolver.cc
GenericCalInfo.cc
ObservationDescription.cc
OpCalImpl.cc
PointingSolver.cc
ScanStats.cc
)

target_link_libraries(opcal
	askap::askap
	askap::accessors
)

set_property(TARGET opcal PROPERTY POSITION_INDEPENDENT_CODE ON)
if (MPI_COMPILE_FLAGS)
	set_property(TARGET opcal PROPERTY COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
endif()


install (FILES
BaselineSolver.h
GenericCalInfo.h
IGenericCalSolver.h
ObservationDescription.h
OpCalImpl.h
PointingSolver.h
ScanStats.h
DESTINATION include/askap/opcal
)
