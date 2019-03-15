add_library( distributedimager OBJECT 
AdviseDI.cc
CalcCore.cc
ContinuumImager.cc
ContinuumMaster.cc
ContinuumWorker.cc
CubeBuilder.cc
CubeComms.cc
MSGroupInfo.cc
MSSplitter.cc
)

set_property(TARGET distributedimager PROPERTY POSITION_INDEPENDENT_CODE ON)
if (MPI_COMPILE_FLAGS)
	set_property(TARGET distributedimager  PROPERTY COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
endif()


install (FILES
AdviseDI.h
CalcCore.h
ContinuumImager.h
ContinuumMaster.h
ContinuumWorker.h
CubeBuilder.h
CubeComms.h
CubeManager.h
MSGroupInfo.h
MSSplitter.h
DESTINATION include/askap/distributedimager
)
