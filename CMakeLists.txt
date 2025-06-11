cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(VERSION "1.17.1")
project(yandasoft VERSION ${VERSION} LANGUAGES C CXX)


# ASKAP Integration (askapsoft)
#
# Creates a full build of the ASKAP Yanda system including associated packages suitable for deployment
# to the HPC system as a module.
#
# A typical dependency project is defined as follows:
#
# FetchContent_Declare(<project>
#   GIT_REPOSITORY    ssh://git@bitbucket.csiro.au:7999/askapsdp/<project>.git
#   GIT_TAG           ${<PROJECT_VERSION_TAG_VAR_from_dependencies.cmake>}
#   SOURCE_DIR        <project>
# )
#
# Author: Eric Bastholm (bas091)

find_package(Git REQUIRED)

option(BUILD_DOCKER       "Build for Docker" OFF)

include(CheckCXXCompilerFlag)
include (GNUInstallDirs)
# Load the version IDs of all the dependency sub repos
include(dependencies.cmake)


message (STATUS "Forming build label")
execute_process(
   COMMAND git describe --tags --always --dirty
   OUTPUT_VARIABLE BUILD_STR
   OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Build label is ${BUILD_STR}")

# Note: This has to be done prior to the cmake project command in case the
# CMAKE_CXX_COMPILER and CMAKE_C_COMPILER are set in cmake config overrides file.
# Otherwise, these two variable are not used if they are set after the project command
# is called.
# Load any developer overrides for FetchContent behaviour like external sources, feature branch inclusion, non default version inclusion, etc.
# DEV_OVERRIDES should be a full path to a cmake script.
message (STATUS "DEV_OVERRIDES = $ENV{DEV_OVERRIDES}")
if(DEFINED ENV{DEV_OVERRIDES})
  message (STATUS "cmake config overrides file specified, $ENV{DEV_OVERRIDES}")
  include($ENV{DEV_OVERRIDES} OPTIONAL)
else()
  message (STATUS "No cmake config overrides file specified, using defaults")
endif()

# Get a remote git repo url, either already set from config, or from current git dir.
message (STATUS "Determining remote git url prefix ...")
if (NOT ASKAPSDP_GIT_URL)
  message (STATUS "Remote git url prefix not specified with ASKAPSDP_GIT_URL value, using current git dir ...")
  execute_process(
     COMMAND git config --get remote.origin.url
     COMMAND xargs dirname
     OUTPUT_VARIABLE ASKAPSDP_GIT_URL
     OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if (NOT ASKAPSDP_GIT_URL)
    message(FATAL_ERROR
      "Could not determine a remote git url to retrieve dependencies. This can be because we are not in a git dir or you need to set ASKAPSDP_GIT_URL var in your dev-overrides or from the cmake command line with '-DASKAPSDP_GIT_URL=<repo-prefix>'.")
  endif()
endif()
message(STATUS "Remote git url prefix is ${ASKAPSDP_GIT_URL}")

message(STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}")


# Declare all of sub repos and tools
include(FetchContent)

# Put our askap cmake modules in the cmake path
# We need this before any other Fetch calls so we have the find_package routines available
# -------- askap cmake modules --------
FetchContent_Declare(askap-cmake
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/askap-cmake.git
  GIT_TAG           ${ASKAP_CMAKE_TAG}
  SOURCE_DIR        askap-cmake
)

message (STATUS "Fetching askap-cmake files")
FetchContent_GetProperties(askap-cmake)
if(NOT askap-cmake_POPULATED)
  FetchContent_Populate(askap-cmake)

endif()
list(APPEND CMAKE_MODULE_PATH "${askap-cmake_SOURCE_DIR}")
#message(STATUS "askap-cmake path is ${askap-cmake_SOURCE_DIR}")
message(STATUS "askap-cmake path is ${CMAKE_MODULE_PATH}")

include(version_utils)

# -------- lofar common --------
FetchContent_Declare(lofar-common
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/lofar-common.git
  GIT_TAG           ${LOFAR_COMMON_TAG}
  SOURCE_DIR        lofar-common
)

# -------- lofar blob --------
FetchContent_Declare(lofar-blob
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/lofar-blob.git
  GIT_TAG           ${LOFAR_BLOB_TAG}
  SOURCE_DIR        lofar-blob
)

# -------- base-askap --------
FetchContent_Declare(base-askap
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/base-askap.git
  GIT_TAG           ${BASE_ASKAP_TAG}
  SOURCE_DIR        base-askap
)

# -------- base-logfilters --------
FetchContent_Declare(base-logfilters
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/base-logfilters.git
  GIT_TAG           ${BASE_LOGFILTERS_TAG}
  SOURCE_DIR        base-logfilters
)

# -------- base-imagemath --------
FetchContent_Declare(base-imagemath
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/base-imagemath.git
  GIT_TAG           ${BASE_IMAGEMATH_TAG}
  SOURCE_DIR        base-imagemath
)

# -------- base-askapparallel --------
FetchContent_Declare(base-askapparallel
  GIT_REPOSITORY    git@github.com:alxndrwllmsn/base-askapparallel.git
  GIT_TAG           ${BASE_ASKAPPARALLEL_TAG}
  SOURCE_DIR        base-askapparallel
)

# -------- base-scimath --------
FetchContent_Declare(base-scimath
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/base-scimath.git
  GIT_TAG           ${BASE_SCIMATH_TAG}
  SOURCE_DIR        base-scimath
)

# -------- base-accessors --------
FetchContent_Declare(base-accessors
  GIT_REPOSITORY    git@github.com:alxndrwllmsn/base-accessors.git
  GIT_TAG           ${BASE_ACCESSORS_TAG}
  SOURCE_DIR        base-accessors
)

message (STATUS "Fetching sub repos")
FetchContent_MakeAvailable(lofar-common lofar-blob
                           base-askap base-logfilters base-imagemath base-askapparallel base-scimath
                           base-accessors
)
message (STATUS "Done - Fetching sub repos")

configure_file(askap_synthesis.cc.in askap/askap_synthesis.cc)

# Required Version of ASKAP dependencies is the MAJOR and MINOR version
# of this package. This allows the PATCH versions to change on everything
# but the version number cannot be less or more than the current version
#
set(REQUIRED_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}")

option (CXX11 "Compile as C++11 if possible" NO)
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

check_cxx_compiler_flag(-fcx-limited-range HAS_FCX)
if (HAS_FCX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcx-limited-range")
endif()

option (ENABLE_OPENMP "Build with OPENMP Support" NO)
option (ENABLE_SHARED "Build shared libraries" YES)
option (ENABLE_RPATH "Include rpath in executables and shared libraries" YES)
option (ENABLE_GRID_WITH_POINTERS "Grid with pointers" YES)
option (ENABLE_FLOAT_IMAGE_PARAMS "Use floats for image parameters" YES)

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
        if ( APPLE )
            # on OSX, the app/executable may not run because it can not find
            # the shared libraries (e.g: ICE and/or log4cxx libraries etc).
            # At least this is what I saw on my MacOS laptop when I tried to
            # get this repo to build and run. I could be wrong but I assume
            # this is because of the way I installed log4cxx on my laptop.
            # I got the old version of log4cxx (0.11.0) while brew install
            # will install log4cxx (0.12.0).
            list(APPEND CMAKE_INSTALL_RPATH "${log4cxx_LIBRARY}" "${ICE_LIB}")
        endif()
    endif (ENABLE_RPATH)
endif(ENABLE_SHARED)

# Yanda depends

# YandaSoft external dependencies. Note, dependencies provided via MakeAvailable mechanism have been removed
find_package(log4cxx REQUIRED)
find_package(Casacore REQUIRED COMPONENTS  ms images mirlib coordinates fits lattices measures scimath scimath_f tables casa)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
find_package(Components REQUIRED)
find_package(MPI)
find_package(CPPUnit)

if (ENABLE_OPENMP)
	find_package(OpenMP)
    if (OPENMP_FOUND)
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        else( )
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif( )
    else (OPENMP_FOUND)
        message(WARNING "Cannot fullfill USE_OPENMP, compiler does not support it")
    endif (OPENMP_FOUND)
endif()

if (CXX11)
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

add_library(yandasoft SHARED
	${CMAKE_CURRENT_BINARY_DIR}/askap/askap_synthesis.cc
)

macro(add_sources_to_yandasoft)
	foreach(arg IN ITEMS ${ARGN})
		target_sources(yandasoft PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${arg})
	endforeach()
endmacro()
add_subdirectory(askap/deconvolution)
add_subdirectory(askap/distributedimager)
add_subdirectory(askap/flagging)
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
  lofar::Common
  lofar::Blob
  askap::askap
  askap::scimath
  askap::imagemath
  askap::accessors
  askap::parallel
  ${log4cxx_LIBRARY}
  ${CASACORE_LIBRARIES}
  ${COMPONENTS_LIBRARY}
  ${Boost_LIBRARIES}
  ${GSL_LIBRARIES}
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
  ${GSL_INCLUDE_DIRS}
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
		target_compile_options(yandasoft PRIVATE ${MPI_COMPILE_FLAGS})
	endif()
endif (MPI_FOUND)

if (OPENMP_FOUND)
    target_link_libraries(yandasoft PRIVATE OpenMP::OpenMP_CXX)
endif (OPENMP_FOUND)

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
if (ENABLE_GRID_WITH_POINTERS)
    target_compile_definitions(yandasoft PUBLIC
        ASKAP_GRID_WITH_POINTERS
    )
endif()
if (ENABLE_FLOAT_IMAGE_PARAMS)
    target_compile_definitions(yandasoft PUBLIC
	    ASKAP_FLOAT_IMAGE_PARAMS
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
    add_subdirectory(tests/flagging)
    add_subdirectory(tests/gridding)
	add_subdirectory(tests/measurementequation)
	add_subdirectory(tests/opcal)
	add_subdirectory(tests/distributedimager)
	add_subdirectory(tests/utils)
endif ()
