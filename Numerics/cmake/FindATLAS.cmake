#
# The Atlas find mecanism goes through FindBLAS and FindLAPACK macros
#


#
# Try to find Atlas with LAPACK (and then BLAS)
#
IF(ATLAS_FIND_REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
ELSE(ATLAS_FIND_REQUIRED)
  FIND_PACKAGE(LAPACK)
ENDIF(ATLAS_FIND_REQUIRED)

#
# Current state information
#
IF(NOT ATLAS_FIND_QUIETLY)
  MESSAGE(STATUS "BLAS_LINKER_FLAGS = ${BLAS_LINKER_FLAGS}")
  MESSAGE(STATUS "BLAS_LIBRARIES = ${BLAS_LIBRARIES}")

  MESSAGE(STATUS "LAPACK_LINKER_FLAGS = ${LAPACK_LINKER_FLAGS}")
  MESSAGE(STATUS "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
  MESSAGE(STATUS "COMPLETE_LAPACK_LIBRARIES = ${COMPLETE_LAPACK_LIBRARIES}")
ENDIF(NOT ATLAS_FIND_QUIETLY)

#
# Minimal config.h generation (cf autoheaders)
#
IF(BLAS_FOUND)
  SET(HAVE_BLAS 1)
ENDIF(BLAS_FOUND)

IF(LAPACK_FOUND)
  SET(HAVE_LAPACK 1)
ENDIF(LAPACK_FOUND)

# Try to find atlas/cblas.h and atlas/clapack.h
# On some systems (Debian, Ubuntu) they are misconfigured
# This can be a problem with some boost headers
IF(ATLAS_FOUND)
  MESSAGE(STATUS "ATLAS library found")
  FIND_PATH(ATLAS_INCLUDE_PATH 
    NAMES atlas/cblas.h atlas/clapack.h
    PATHS ${CMAKE_SOURCE_DIR}/src/utils/AtlasLocal)
  IF(NOT ATLAS_INCLUDE_PATH)
    MESSAGE(STATUS "Cannot find atlas/cblas.h and atlas/clapack.h")
    MESSAGE(STATUS "Failing back to cblas.h and clapack.h")
    FIND_PATH(CBLAS_INCLUDE_PATH 
      NAMES cblas.h clapack.h)
    IF(NOT CBLAS_INCLUDE_PATH)
      MESSAGE(STATUS "WARNING: cannot find cblas.h and clapack.h")
      MESSAGE(STATUS "Failing back to the fortran interface to Blas and Lapack")
      MESSAGE(STATUS "The ATLAS library may not be used correctly.")
    ELSE(NOT CBLAS_INCLUDE_PATH)
#      CHECK_SYMBOL_EXISTS(cblas_xerbla ${CBLAS_INCLUDE_PATH}/cblas.h HAVE_XERBLA)
      SET(HAVE_XERBLA 1)
      SET(HAVE_ATLAS 1)
      SET(HAVE_CBLAS_H 1)
      SET(HAVE_CLAPACK_H 1)
    ENDIF(NOT CBLAS_INCLUDE_PATH)
  ELSE(NOT ATLAS_INCLUDE_PATH)
    SET(ATLAS_INCLUDE_PATH ${ATLAS_INCLUDE_PATH}/atlas)
    SET(HAVE_ATLAS 1)
    SET(HAVE_CBLAS_H 1)
    SET(HAVE_CLAPACK_H 1)
  ENDIF(NOT ATLAS_INCLUDE_PATH)
ELSE(ATLAS_FOUND)
  IF(ATLAS_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "ATLAS library not found")
  ENDIF(ATLAS_FIND_REQUIRED)
ENDIF(ATLAS_FOUND)
