add_sources_to_yandasoft(
BaselineSolver.cc
GenericCalInfo.cc
ObservationDescription.cc
OpCalImpl.cc
PointingSolver.cc
ScanStats.cc
)

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
