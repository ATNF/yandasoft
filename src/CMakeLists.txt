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


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
find_package(Askap REQUIRED)
find_package(Scimath REQUIRED)

find_package(log4cxx REQUIRED)
find_package(Casacore REQUIRED COMPONENTS  ms images mirlib coordinates fits lattices measures scimath scimath_f tables casa)
find_package(FFTW)
find_package(LofarCommon)
find_package(Boost REQUIRED COMPONENTS system)
find_package(Accessors REQUIRED)
find_package(ImageMath REQUIRED)
find_package(Scimath REQUIRED)

include_directories(${FFTW_INCLUDE_DIRS})
include_directories(${ASKAP_INCLUDE_DIRS})
include_directories(${SCIMATH_INCLUDE_DIRS})
include_directories(${LofarCommon_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ACCESSORS_INCLUDE_DIRS})
include_directories(${IMAGEMATH_INCLUDE_DIRS})
include_directories(${SCIMATH_INCLUDE_DIRS})

if (log4cxx_FOUND)
	include_directories(${log4cxx_INCLUDE_DIRS})
	add_definitions(-DHAVE_LOG4CXX)
endif()
if (CASACORE_FOUND)
	include_directories(${CASACORE_INCLUDE_DIRS})
	add_definitions(-DHAVE_AIPSPP)
endif()
add_library(yanda_synthesis ${all_sources})

target_link_libraries(yanda_synthesis 
	${ASKAP_LIBRARY}
	${log4cxx_LIBRARY}
	${CASACORE_LIBRARIES}
	${FFTW_LIBRARIES}
	${SCIMATH_LIBRARY}
	${IMAGEMATH_LIBRARY}
	${LofarCommon_LIBRARY}
	${Boost_LIBRARIES}
	${ACCESSORS_LIBRARY}
	)

install (
TARGETS yanda_synthesis
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

