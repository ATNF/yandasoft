add_sources_to_yandasoft(
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
