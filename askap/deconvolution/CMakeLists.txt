add_sources_to_yandasoft(
DeconvolverFactory.cc
DeconvolverHelpers.cc
)

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
