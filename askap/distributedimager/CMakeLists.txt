add_sources_to_yandasoft(
AdviseDI.cc
CalcCore.cc
ContinuumImager.cc
ContinuumMaster.cc
ContinuumWorker.cc
CubeComms.cc
WorkUnitContainer.cc
DataSourceManager.cc
)

install (FILES
AdviseDI.h
CalcCore.h
ContinuumImager.h
ContinuumMaster.h
ContinuumWorker.h
CubeBuilder.h
CubeBuilder.tcc
CubeComms.h
CubeManager.h
WorkUnitContainer.h
DataSourceManager.h
DESTINATION include/askap/distributedimager
)
