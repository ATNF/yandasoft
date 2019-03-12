add_library( opcal OBJECT 
BaselineSolver.cc
GenericCalInfo.cc
ObservationDescription.cc
OpCalImpl.cc
PointingSolver.cc
ScanStats.cc
)

set_property(TARGET opcal PROPERTY POSITION_INDEPENDENT_CODE ON)

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
