# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

# From eigen

if (UMFPACK_INCLUDE_DIR AND UMFPACK_LIBRARIES)
  set(UMFPACK_FIND_QUIETLY TRUE)
endif (UMFPACK_INCLUDE_DIR AND UMFPACK_LIBRARIES)

find_package(PkgConfig)
include(LibFindMacros)

# Use pkg-config to get hints about paths
pkg_check_modules(UMFPACK QUIET umfpack)

IF(UMFPACK_FOUND)
  SET(UMFPACK_INCLUDE_DIR ${UMFPACK_INCLUDE_DIRS})
  IF(UMFPACK_LIBRARY_DIRS)
    SET(UMFPACK_LIBDIR ${UMFPACK_LIBRARY_DIRS})
  ELSE(UMFPACK_LIBRARY_DIRS)
    SET(UMFPACK_LIBDIR " ")
  ENDIF(UMFPACK_LIBRARY_DIRS)
ENDIF(UMFPACK_FOUND)

find_path(UMFPACK_INCLUDE_DIR
  NAMES
  umfpack.h
  PATHS
  $ENV{UMFPACKDIR}
  PATH_SUFFIXES
  suitesparse
  ufsparse
)

find_library(UMFPACK_LIBRARIES umfpack PATHS $ENV{UMFPACKDIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK DEFAULT_MSG
                                  UMFPACK_INCLUDE_DIR UMFPACK_LIBRARIES)

mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARIES)
