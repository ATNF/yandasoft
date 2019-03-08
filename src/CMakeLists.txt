#
# Manage the synthesis build
#
project(yanda_synthesis)

set ( _modules  

deconvolution
#distributedimager
#gridding
#measurementequation
#messages
#opcal
#parallel
#simulation
#utils

)
# This is not a recommended way of doing this as if I add a file to the list
# cmake will not automatically know to rebuild. But it does save adding all the subdirectories
# and the individual libraries

foreach (module ${_modules})
	FILE(GLOB my_list "${module}/*.cc")
	list(APPEND all_sources ${my_list})
endforeach (module)


message(STATUS ${all_sources})

add_library(yanda_synthesis ${all_sources})


install (
TARGETS yanda_synthesis
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

