add_library( deconvolution OBJECT 
DeconvolverFactory.cc
DeconvolverHelpers.cc
)

set_property(TARGET deconvolution PROPERTY POSITION_INDEPENDENT_CODE ON)

if (MPI_COMPILE_FLAGS)
	set_property(TARGET deconvolution PROPERTY COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
endif()

install (FILES
BasisFunction.h
BasisFunction.tcc
DeconvolverBase.h
DeconvolverBase.tcc
DeconvolverBasisFunction.h
DeconvolverBasisFunction.tcc
DeconvolverControl.h
DeconvolverControl.tcc
DeconvolverEntropy.h
DeconvolverEntropy.tcc
DeconvolverFactory.h
DeconvolverFista.h
DeconvolverFista.tcc
DeconvolverHelpers.h
DeconvolverHogbom.h
DeconvolverHogbom.tcc
DeconvolverMonitor.h
DeconvolverMonitor.tcc
DeconvolverMultiTermBasisFunction.h
DeconvolverMultiTermBasisFunction.tcc
DeconvolverState.h
DeconvolverState.tcc
Emptiness.h
Emptiness.tcc
EntropyBase.h
EntropyBase.tcc
EntropyI.h
EntropyI.tcc
MultiScaleBasisFunction.h
MultiScaleBasisFunction.tcc
PointBasisFunction.h
PointBasisFunction.tcc

DESTINATION include/askap/deconvolution
)
