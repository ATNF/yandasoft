# A list of specific version for the dependencies of this build.
#
# The contents here will aid in 'pinning' the contents of the build. Changing this file constitutes changing the build
# and this means that the integration version should be promoted. The version numbers (tags) specified here will be
# injected into the CMakeLists.txt file in the External_project_add definitions and will, therefore, control the
# versions of various sub libs built and integrated. They need not be official versions and can be arbitrary branches or
# tags such as dev feature branches etc.
#
# You can affect the build content by manipulating the "tag/branch" column. Should be something that git expects for that 
# repo.
#     Identifier                   tag/branch       cache  type      description                        force it
set ( ASKAP_CMAKE_TAG              tags/1.6.0      CACHE  STRING    "askap.cmake tools"                FORCE )
set ( ASKAP_DEV_TAG                tags/2.13.0     CACHE  STRING    "askap dev tools"                  FORCE )
set ( LOFAR_COMMON_TAG             tags/1.3.0      CACHE  STRING    "lofar-common version"             FORCE )
set ( LOFAR_BLOB_TAG               tags/1.3.1      CACHE  STRING    "lofar-blob version"               FORCE )
set ( BASE_ASKAP_TAG               tags/1.8.0      CACHE  STRING    "base-askap version"               FORCE )
set ( BASE_LOGFILTERS_TAG          tags/1.5.0      CACHE  STRING    "base-logfilters version"          FORCE )
set ( BASE_IMAGEMATH_TAG           tags/1.15.0     CACHE  STRING    "base-imagemath version"           FORCE )
set ( BASE_ASKAPPARALLEL_TAG       adios_integration     CACHE  STRING    "base-askapparallel version"       FORCE )
set ( BASE_SCIMATH_TAG             tags/1.15.0     CACHE  STRING    "base-scimath version"             FORCE )
set ( BASE_ACCESSORS_TAG           adios_integration     CACHE  STRING    "base-accessors version"           FORCE )
set ( BASE_COMPONENTS_TAG          tags/1.14.0     CACHE  STRING    "base-components version"          FORCE )
