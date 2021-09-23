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
set ( ASKAP_CMAKE_TAG              tags/1.2.0                               CACHE  STRING    "askap.cmake tools"                FORCE )
set ( ASKAP_DEV_TAG                tags/2.3.1                               CACHE  STRING    "askap dev tools"                  FORCE )
set ( LOFAR_COMMON_TAG             tags/1.1.0                               CACHE  STRING    "lofar-common version"             FORCE )
set ( LOFAR_BLOB_TAG               tags/1.1.0                               CACHE  STRING    "lofar-blob version"               FORCE )
set ( BASE_ASKAP_TAG               tags/1.2.0                               CACHE  STRING    "base-askap version"               FORCE )
#set ( BASE_ASKAP_TAG               feature/AXA-1320-convert-base-askap-to-new-build-recipe                               CACHE  STRING    "base-askap version"               FORCE )
set ( BASE_LOGFILTERS_TAG          tags/1.1.0                               CACHE  STRING    "base-logfilters version"          FORCE )
set ( BASE_IMAGEMATH_TAG           tags/1.3.0                               CACHE  STRING    "base-imagemath version"           FORCE )
set ( BASE_ASKAPPARRALLEL_TAG      tags/1.1.0                               CACHE  STRING    "base-askapparrallel version"      FORCE )
set ( BASE_SCIMATH_TAG             tags/1.3.0                               CACHE  STRING    "base-scimath version"             FORCE )
set ( BASE_ACCESSORS_TAG           tags/1.3.0                               CACHE  STRING    "base-accessors version"           FORCE )
set ( BASE_COMPONENTS_TAG          tags/1.1.0                               CACHE  STRING    "base-components version"          FORCE )
set ( ASKAP_ANALYSIS_TAG           tags/1.2.1                               CACHE  STRING    "askap-analysis version"           FORCE )
set ( ASKAP_YANDASOFT_TAG          tags/1.3.0                               CACHE  STRING    "yandasoft version"                FORCE )
set ( ASKAP_PIPELINETASKS_TAG      tags/1.3.0                               CACHE  STRING    "askap-pipelinetasks version"      FORCE )
set ( ASKAP_INTERFACES_TAG         tags/1.2.0                               CACHE  STRING    "askap-interfaces version"         FORCE )
set ( ASKAP_SMS_TAG                tags/1.2.0                               CACHE  STRING    "askap-sms version" FORCE )
 
# TOS related repos are not versioned yet, so pinned this build with commit hash values.
set ( PYTHON_ASKAP_TAG             3c4871a07d7cdf8a71871074cd90e8b3bf8d16de CACHE  STRING    "python-askap (tos) version"       FORCE )
set ( PYTHON_ICEUTILS_TAG          8e34898ef30d530648b02cf81ca8e6e3c9d0781e CACHE  STRING    "python-iceutils (tos) version"    FORCE )
set ( PYTHON_LOGHANDLERS_TAG       e3b682ef6a79ab14943789ad5ef6e5a46c4d7bb7 CACHE  STRING    "python-loghandlers (tos) version" FORCE )
set ( PYTHON_PARSET_TAG            e6bfeb00e5b3a9e5a8aae20db48a1e6470aa1c35 CACHE  STRING    "python-parset (tos) version"      FORCE )
