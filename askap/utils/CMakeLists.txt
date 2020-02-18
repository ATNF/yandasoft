add_sources_to_yandasoft(
	CommandLineParser.cc
	LinmosUtils.cc
)

install (FILES
	CommandLineParser.h
	LinmosUtils.h
DESTINATION include/askap/utils
)
