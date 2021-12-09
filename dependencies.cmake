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
set ( LOFAR_COMMON_TAG             tags/1.1.0                               CACHE  STRING    "lofar-common version"             FORCE )
set ( LOFAR_BLOB_TAG               tags/1.1.0                               CACHE  STRING    "lofar-blob version"               FORCE )
# need to release new base-askap, this tag is the state after pull request #37
set ( BASE_ASKAP_TAG               30088d444dcdacc4d4c896a1cf6ad3cd0b1d6a97   CACHE  STRING    "base-askap version"               FORCE )
set ( BASE_LOGFILTERS_TAG          tags/1.1.0                               CACHE  STRING    "base-logfilters version"          FORCE )
# need to release new imagemath, this tag is the state after pull request #38 
set ( BASE_IMAGEMATH_TAG           483d77c96067f726e7078e7ec6799a39830d269a   CACHE  STRING    "base-imagemath version"           FORCE )
# need to release new base-parallel, this tag is the state after pull request #14
set ( BASE_ASKAPPARRALLEL_TAG      afd594d5906c24043dd10d4e8472b282fd10e289   CACHE  STRING    "base-askapparrallel version"      FORCE )
# need to release new scimath, this tag is the state after pull request #49
set ( BASE_SCIMATH_TAG             96404d4d75d2bc7dcaf7fa8ddf462f64f91e7488   CACHE  STRING    "base-scimath version"             FORCE )
# need to release new accessors, this tag is the state after pull request #39
set ( BASE_ACCESSORS_TAG           083682416cf7d1604688fa33f5f19e5ad6e70e98   CACHE  STRING    "base-accessors version"           FORCE )
set ( BASE_COMPONENTS_TAG          tags/1.1.0                               CACHE  STRING    "base-components version"          FORCE )
set ( ASKAP_SMS_TAG                tags/1.2.0                               CACHE  STRING    "askap-sms version" FORCE )
 
