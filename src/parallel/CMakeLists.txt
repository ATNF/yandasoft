add_library( parallel OBJECT 
AdviseParallel.cc
BPCalibratorParallel.cc
CalibratorParallel.cc
ContSubtractParallel.cc
GroupVisAggregator.cc
ImagerParallel.cc
MEParallel.cc
MEParallelApp.cc
ParallelAccessor.cc
ParallelIteratorStatus.cc
ParallelWriteIterator.cc
SimParallel.cc
SynParallel.cc
)

set_property(TARGET parallel PROPERTY POSITION_INDEPENDENT_CODE ON)
if (MPI_COMPILE_FLAGS)
	set_property(TARGET parallel  PROPERTY COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
endif()


install (FILES
AdviseParallel.h
BPCalibratorParallel.h
CalibratorParallel.h
ContSubtractParallel.h
GroupVisAggregator.h
ImagerParallel.h
MEParallel.h
MEParallelApp.h
ParallelAccessor.h
ParallelIteratorStatus.h
ParallelWriteIterator.h
SimParallel.h
SynParallel.h
DESTINATION include/askap/parallel
)
