
cmake_minimum_required (VERSION 3.1.0)
include(CheckCXXCompilerFlag)
project(yandasoft VERSION 0.0.1 LANGUAGES CXX C)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/askap-cmake)
    execute_process(COMMAND git pull 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/askap-cmake)
else ()	
    execute_process(COMMAND git clone https://bitbucket.csiro.au/scm/askapsdp/askap-cmake.git
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif ()


set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/askap-cmake")

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

if (ENABLE_OPENMP)
find_package(OpenMP REQUIRED) # Find the package
endif ()

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Manage the synthesis build

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Yanda depends

# Yanda Packages
# YandaSoft dependencies .. in order
include(install_yanda_depends)
find_package(LofarCommon)

set (EXTERNAL_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})

if (NOT LofarCommon_FOUND)
	install_yanda_depends(lofar-common)
	find_package(LofarCommon REQUIRED)
endif ()

find_package(LofarBlob)
if (NOT LofarBlob_FOUND)
	install_yanda_depends(lofar-blob)
	find_package(LofarBlob REQUIRED)
endif ()

find_package(Askap)
if (NOT ASKAP_FOUND)
	install_yanda_depends(base-askap)
	find_package(Askap REQUIRED)
endif ()

find_package(Imagemath)
if (NOT IMAGEMATH_FOUND)
	install_yanda_depends(base-imagemath)
	find_package(Imagemath REQUIRED)
endif ()

find_package(Scimath)
if (NOT SCIMATH_FOUND)
	install_yanda_depends(base-scimath)
	find_package(Scimath REQUIRED)
endif ()

find_package(AskapParallel)
if (NOT ASKAPPARALLEL_FOUND)
	install_yanda_depends(base-askapparallel)
	find_package(AskapParallel REQUIRED)
endif ()

find_package(Accessors)
if (NOT ACCESSORS_FOUND)
	install_yanda_depends(base-accessors)
	find_package(Accessors REQUIRED)
endif ()



find_package(log4cxx REQUIRED)
find_package(Casacore REQUIRED COMPONENTS  ms images mirlib coordinates fits lattices measures scimath scimath_f tables casa)
find_package(FFTW REQUIRED)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
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

if (NOT GSL_VERSION VERSION_LESS 2.0)
	add_definitions(-DHAVE_GSL2)
endif()

if (CASACORE3)
	add_definitions(-DHAVE_CASACORE3)
endif()

if (CASACORE3 OR CXX11)
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

add_subdirectory(askap/deconvolution)
add_subdirectory(askap/distributedimager)
add_subdirectory(askap/gridding)
add_subdirectory(askap/measurementequation)
add_subdirectory(askap/messages)
add_subdirectory(askap/opcal)
add_subdirectory(askap/parallel)
add_subdirectory(askap/simulation)
add_subdirectory(askap/utils)


add_library(yanda_synthesis SHARED
	askap/dummy.cpp
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
	${MPI_LIBRARIES}
	)
if (MPI_FOUND)
	target_link_libraries(yanda_synthesis
		${MPI_LIBRARIES}
	)

endif (MPI_FOUND)

if (ENABLE_OPENMP)
    target_link_libraries(yanda_synthesis
        ${OpenMP_CXX_LIBRARIES})
endif ()

install (
TARGETS yanda_synthesis
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

add_subdirectory(apps)

find_package(CPPUnit)
if (CPPUNIT_FOUND)
    include_directories(${CPPUNIT_INCLUDE_DIR})
    target_link_libraries(yanda_synthesis
       ${CPPUNIT_LIBRARY}
    )
	include(CTest)
	enable_testing()
	add_subdirectory(tests/deconvolution)
	add_subdirectory(tests/gridding)
endif ()

