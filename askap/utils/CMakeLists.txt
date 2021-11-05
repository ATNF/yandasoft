add_sources_to_yandasoft(
	LinmosUtils.cc
    EigenSolve.cc
	IlluminationUtils.cc
	ImplCalWeightSolver.cc
	SkyCatalogTabWriter.cc
    DelaySolverImpl.cc
)

install (FILES
	LinmosUtils.h
	EigenSolve.h
	IlluminationUtils.h
	ImplCalWeightSolver.h
	SkyCatalogTabWriter.h
    DelaySolverImpl.h
    DESTINATION include/askap/utils
)
