add_sources_to_yandasoft(
	CommandLineParser.cc
	LinmosUtils.cc
    EigenSolve.cc
	IlluminationUtils.cc
	ImplCalWeightSolver.cc
	SkyCatalogTabWriter.cc
    DelaySolverImpl.cc
)

install (FILES
	CommandLineParser.h
	LinmosUtils.h
	EigenSolve.h
	IlluminationUtils.h
	ImplCalWeightSolver.h
	SkyCatalogTabWriter.h
    DelaySolverImpl.h
    DESTINATION include/askap/utils
)
