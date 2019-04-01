cmake_minimum_required (VERSION 3.2.0)

project(yandasoft VERSION 0.0.1 LANGUAGES CXX)

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


add_subdirectory(src)


