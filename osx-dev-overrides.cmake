set ( CASACORE_ROOT_DIR     $ENV{CASACORE_ROOT_DIR}                     CACHE  PATH     "Path to the CASACORE root"         FORCE )
set ( COMPONENTS_ROOT_DIR   $ENV{COMPONENTS_ROOT_DIR}                   CACHE  PATH     "Path to the CASA components root"  FORCE )
set ( log4cxx_ROOT_DIR      $ENV{log4cxx_ROOT_DIR}                      CACHE  PATH     "Path to the log4cxx root"          FORCE )
set ( CPPUNIT_ROOT_DIR      $ENV{CPPUNIT_ROOT_DIR}                      CACHE  PATH     "Path to the CPPUNIT root"          FORCE )
set ( XercesC_ROOT_DIR      $ENV{XercesC_ROOT_DIR}                      CACHE  PATH     "Path to the XercesC includes"      FORCE )
set ( XercesC_INCLUDE_DIR   ${XercesC_ROOT_DIR}/include                 CACHE  PATH     "Path to the XercesC includes"      FORCE )

set ( CMAKE_CXX_COMPILER    "/usr/local/opt/llvm/bin/clang++"             CACHE  STRING   "Specify the C++ compiler"            FORCE )
set ( CMAKE_C_COMPILER      "/usr/local/opt/llvm/bin/clang"               CACHE  STRING   "Specify the C compiler"              FORCE )
