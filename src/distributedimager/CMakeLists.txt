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
