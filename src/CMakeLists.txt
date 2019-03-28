#
# Manage the synthesis build
#
project(yanda_synthesis)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
find_package(Askap REQUIRED)
find_package(Scimath REQUIRED)
find_package(AskapParallel REQUIRED)
find_package(log4cxx REQUIRED)
find_package(Casacore REQUIRED COMPONENTS  ms images mirlib coordinates fits lattices measures scimath scimath_f tables casa)
find_package(FFTW REQUIRED)
find_package(GSL REQUIRED)
find_package(LofarCommon REQUIRED)
find_package(LofarBlob REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
find_package(Accessors REQUIRED)
find_package(Imagemath REQUIRED)
find_package(Scimath REQUIRED)
find_package(Components REQUIRED)
find_package(MPI REQUIRED)

include_directories(${FFTW_INCLUDE_DIRS})
include_directories(${ASKAP_INCLUDE_DIRS})
include_directories(${SCIMATH_INCLUDE_DIRS})
include_directories(${LofarCommon_INCLUDE_DIRS})
include_directories(${LofarBlob_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ACCESSORS_INCLUDE_DIRS})
include_directories(${IMAGEMATH_INCLUDE_DIRS})
include_directories(${SCIMATH_INCLUDE_DIRS})
include_directories(${COMPONENTS_INCLUDE_DIRS})

include_directories(${MPI_INCLUDE_PATH})

if (MPI_FOUND)
	add_definitions(-DHAVE_MPI)
endif()

include_directories(${log4cxx_INCLUDE_DIRS})
add_definitions(-DHAVE_LOG4CXX)
include_directories(${CASACORE_INCLUDE_DIRS})
add_definitions(-DHAVE_AIPSPP)
add_definitions(-Dcasa=casacore)

if (NOT CASACORE_VERSION STREQUAL trunk AND
    CASACORE_VERSION VERSION_GREATER 2.4.0)
	add_definitions(-DHAVE_CASACORE3)
endif()

add_subdirectory(deconvolution)
add_subdirectory(distributedimager)
add_subdirectory(gridding)
add_subdirectory(measurementequation)
add_subdirectory(messages)
add_subdirectory(opcal)
add_subdirectory(parallel)
add_subdirectory(simulation)
add_subdirectory(utils)


add_library(yanda_synthesis SHARED
	$<TARGET_OBJECTS:deconvolution>
	$<TARGET_OBJECTS:distributedimager>
	$<TARGET_OBJECTS:gridding>
	$<TARGET_OBJECTS:measurementequation>
	$<TARGET_OBJECTS:messages>
	$<TARGET_OBJECTS:opcal>
	$<TARGET_OBJECTS:parallel>
	$<TARGET_OBJECTS:simulation>
	$<TARGET_OBJECTS:utils>
	)

target_link_libraries(yanda_synthesis 
	${ASKAP_LIBRARY}
	${log4cxx_LIBRARY}
	${CASACORE_LIBRARIES}
	${FFTW_LIBRARIES}
	${SCIMATH_LIBRARY}
	${IMAGEMATH_LIBRARY}
	${LofarBlob_LIBRARY}
	${LofarCommon_LIBRARY}
	${Boost_LIBRARIES}
	${ACCESSORS_LIBRARY}
	${COMPONENTS_LIBRARY}
	${ASKAPPARALLEL_LIBRARY}
	)
if (MPI_FOUND)
	target_link_libraries(yanda_synthesis
		${MPI_LIBRARIES}
	)
endif (MPI_FOUND)


install (
TARGETS yanda_synthesis
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

add_subdirectory(apps)
