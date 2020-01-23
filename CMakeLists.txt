
cmake_minimum_required (VERSION 3.1.0)
include(CheckCXXCompilerFlag)
project(yandasoft VERSION 0.0.1 LANGUAGES CXX C)

if (EXISTS ${CMAKE_SOURCE_DIR}/askap-cmake)
    execute_process(COMMAND git pull 
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/askap-cmake
    TIMEOUT 10)
else ()	
    execute_process(COMMAND git clone https://bitbucket.csiro.au/scm/askapsdp/askap-cmake.git
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif ()


set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/askap-cmake")

option (CXX11 "Compile as C++11 if possible" YES)
# Detect if the compiler supports C++11 if we want to use it.
if (CXX11)
    check_cxx_compiler_flag(-std=c++11 HAS_CXX11)
    if (HAS_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      set(CXX11 ON)
    else()
	    message(STATUS "CASA_BUILD requires a c++11 compatible compiler")
      option(CXX11 "" NO)
    endif()
endif()

option (ENABLE_OPENMP "Build with OPENMP Support" NO)
option (ENABLE_SHARED "Build shared libraries" YES)
option (ENABLE_RPATH "Include rpath in executables and shared libraries" YES)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{CPPFLAGS}")
message(STATUS "CXXFLAGS are ${CMAKE_CXX_FLAGS}")

# uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

if (ENABLE_SHARED)
option (BUILD_SHARED_LIBS "" YES)
    if (ENABLE_RPATH)
        # Set RPATH to use for installed targets; append linker search path
        set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
        set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif (ENABLE_RPATH)
endif(ENABLE_SHARED)

# Yanda depends

# Yanda Packages
# YandaSoft dependencies .. in order
find_package(lofar-common REQUIRED)
find_package(lofar-blob REQUIRED)
find_package(askap-askap REQUIRED)
find_package(askap-imagemath REQUIRED)
find_package(askap-scimath REQUIRED)
find_package(askap-parallel REQUIRED)
find_package(askap-accessors REQUIRED)
find_package(log4cxx REQUIRED)
find_package(Casacore REQUIRED COMPONENTS  ms images mirlib coordinates fits lattices measures scimath scimath_f tables casa)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
find_package(Components REQUIRED)
find_package(MPI)
find_package(CPPUnit)
if (ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
endif()


if (CASACORE3 OR CXX11)
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

add_library(yandasoft SHARED)

macro(add_sources_to_yandasoft)
	foreach(arg IN ITEMS ${ARGN})
		target_sources(yandasoft PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${arg})
	endforeach()
endmacro()

add_subdirectory(askap/deconvolution)
add_subdirectory(askap/distributedimager)
add_subdirectory(askap/gridding)
add_subdirectory(askap/measurementequation)
add_subdirectory(askap/messages)
add_subdirectory(askap/opcal)
add_subdirectory(askap/parallel)
add_subdirectory(askap/simulation)
add_subdirectory(askap/utils)

set_target_properties(yandasoft PROPERTIES
	OUTPUT_NAME yanda_synthesis
)

add_library(askap::yandasoft ALIAS yandasoft)

target_link_libraries(yandasoft PUBLIC
	askap::askap
	askap::scimath
	askap::imagemath
	lofar::Common
	lofar::Blob
	askap::accessors
	askap::parallel
	${log4cxx_LIBRARY}
	${CASACORE_LIBRARIES}
	${COMPONENTS_LIBRARY}
	${Boost_LIBRARIES}
)

target_compile_definitions(yandasoft PUBLIC
    casa=casacore
    HAVE_AIPSPP
    HAVE_LOG4CXX
)

target_include_directories(yandasoft PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${Boost_INCLUDE_DIRS}
  ${COMPONENTS_INCLUDE_DIRS}
  ${log4cxx_INCLUDE_DIRS}
  ${CASACORE_INCLUDE_DIRS}
)

if (MPI_FOUND)
	target_link_libraries(yandasoft PUBLIC
		${MPI_LIBRARIES}
	)
	target_include_directories(yandasoft PUBLIC
		${MPI_INCLUDE_PATH}
	)
	target_compile_definitions(yandasoft PUBLIC
		HAVE_MPI
	)
   if (MPI_COMPILE_FLAGS)
		target_compile_options(yandasoft PROPERTY COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
	endif()
endif (MPI_FOUND)

if (ENABLE_OPENMP)
    target_link_libraries(yandasoft PRIVATE OpenMP::OpenMP_CXX)
endif (ENABLE_OPENMP)

if (NOT GSL_VERSION VERSION_LESS 2.0)
	target_compile_definitions(yandasoft PUBLIC
		HAVE_GSL2
	)
endif()

if (CASACORE3)
	target_compile_definitions(yandasoft PUBLIC
		HAVE_CASACORE3
	)
endif()


install (
TARGETS yandasoft
EXPORT yandasoft-targets
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

add_subdirectory(apps)

include(CTest)
if (BUILD_TESTING)
	if (NOT CPPUNIT_FOUND)
		message(FATAL_ERROR "CPPUnit is needed for building and running unit tests")
	endif()
	enable_testing()
	add_subdirectory(tests/deconvolution)
	add_subdirectory(tests/gridding)
	add_subdirectory(tests/measurementequation)
	add_subdirectory(tests/opcal)
endif ()


include(yanda_export)
yanda_export(askap)
