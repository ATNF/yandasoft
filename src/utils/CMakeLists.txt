add_library( utils OBJECT 
	CommandLineParser.cc
	LinmosUtils.cc
)
install (FILES
	CommandLineParser.h
	LinmosUtils.h
DESTINATION include/askap/utils
)
