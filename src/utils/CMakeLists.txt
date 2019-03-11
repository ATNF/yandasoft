add_library( utils OBJECT 
	CommandLineParser.cc
	LinmosUtils.cc
)

set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)

install (FILES
	CommandLineParser.h
	LinmosUtils.h
DESTINATION include/askap/utils
)
