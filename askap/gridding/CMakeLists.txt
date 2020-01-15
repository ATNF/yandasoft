add_library( gridding OBJECT 
AProjectGridderBase.cc
AProjectWStackVisGridder.cc
ATCAIllumination.cc
AWProjectVisGridder.cc
AltWProjectVisGridder.cc
BasicCompositeIllumination.cc
BoxVisGridder.cc
DiskIllumination.cc
FrequencyMapper.cc
GaussianWSampling.cc
GridKernel.cc
IBasicIllumination.cc
IVisGridder.cc
IVisWeights.cc
PowerWSampling.cc
SKA_LOWIllumination.cc
SmearingGridderAdapter.cc
SnapShotImagingGridderAdapter.cc
SphFuncVisGridder.cc
SupportSearcher.cc
TableVisGridder.cc
TestCFGenPerformance.cc
VisGridderFactory.cc
VisWeightsMultiFrequency.cc
WDependentGridderBase.cc
WProjectVisGridder.cc
WStackVisGridder.cc
)

set_property(TARGET gridding PROPERTY POSITION_INDEPENDENT_CODE ON)
if (MPI_COMPILE_FLAGS)
	set_property(TARGET gridding PROPERTY COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
endif()


install (FILES
AProjectGridderBase.h
AProjectGridderBase.tcc
AProjectWStackVisGridder.h
ATCAIllumination.h
AWProjectVisGridder.h
AltWProjectVisGridder.h
BasicCompositeIllumination.h
BoxVisGridder.h
DiskIllumination.h
FrequencyMapper.h
GaussianWSampling.h
GridKernel.h
IBasicIllumination.h
IVisGridder.h
IVisWeights.h
IWSampling.h
PowerWSampling.h
SKA_LOWIllumination.h
SmearingGridderAdapter.h
SnapShotImagingGridderAdapter.h
SphFuncVisGridder.h
SphFuncVisGridder.tcc
SupportSearcher.h
SupportSearcher.tcc
TableVisGridder.h
TestCFGenPerformance.h
UVPattern.h
VisGridderFactory.h
VisGridderWithPadding.h
VisWeightsMultiFrequency.h
WDependentGridderBase.h
WProjectVisGridder.h
WStackVisGridder.h
DESTINATION include/askap/gridding
)

