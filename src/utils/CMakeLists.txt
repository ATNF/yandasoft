add_library( utils OBJECT 
LinmosUtils.cc
)
install (FILES
LinmosUtils.h
DESTINATION include/askap/utils
)
