#
# Siconos Numerics : general configuration
#
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.4)
PROJECT(Siconos_Numerics Fortran)

SET(MAJOR_VERSION 2)
SET(MINOR_VERSION 2)
SET(PATCH_VERSION 1)

SET(VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")

SET( ENV{LC_ALL} C )

#
# SVN revision number
#
IF(EXISTS .svn)
  FIND_PACKAGE(Subversion)
  IF(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
    IF(Project_WC_REVISION)
      SET(SVN_REVISION ${Project_WC_REVISION})
      MESSAGE("-- Current svn revision is ${SVN_REVISION}")
    ENDIF(Project_WC_REVISION)
  ENDIF(Subversion_FOUND)
ENDIF(EXISTS .svn)
#
# Standard modules
#
INCLUDE(CPack)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckFunctionExists)
#
# Local modules
#
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

INCLUDE(fortran)
INCLUDE(LibraryProjectSetup)

#
# Required packages 
#
FIND_PACKAGE(LAPACK REQUIRED)

MESSAGE(STATUS "LAPACK_LINKER_FLAGS = ${LAPACK_LINKER_FLAGS}")
MESSAGE(STATUS "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")

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
  MESSAGE("-- ATLAS library found")
  FIND_PATH(ATLAS_INCLUDE_PATH 
    NAMES atlas/cblas.h atlas/clapack.h)
  IF(NOT ATLAS_INCLUDE_PATH)
    MESSAGE("-- Cannot find atlas/cblas.h and atlas/clapack.h")
    MESSAGE("-- Failing back to cblas.h and clapack.h")
    FIND_PATH(CBLAS_INCLUDE_PATH 
      NAMES cblas.h clapack.h)
    IF(NOT CBLAS_INCLUDE_PATH)
      MESSAGE("-- WARNING: cannot find cblas.h and clapack.h")
      MESSAGE("-- Failing back to the fortran interface to Blas and Lapack")
      MESSAGE("-- The ATLAS library may not be used correctly.")
    ELSE(NOT CBLAS_INCLUDE_PATH)
#      CHECK_SYMBOL_EXISTS(cblas_xerbla ${CBLAS_INCLUDE_PATH}/cblas.h HAVE_XERBLA)
      SET(HAVE_XERBLA 1)
      SET(HAVE_ATLAS 1)
      SET(HAVE_CBLAS_H 1)
      SET(HAVE_CLAPACK_H 1)
    ENDIF(NOT CBLAS_INCLUDE_PATH)
  ELSE(NOT ATLAS_INCLUDE_PATH)
    SET(HAVE_ATLAS 1)
    SET(HAVE_CBLAS_H 1)
    SET(HAVE_CLAPACK_H 1)
  ENDIF(NOT ATLAS_INCLUDE_PATH)
ENDIF(ATLAS_FOUND)

# config.h and include
CONFIGURE_FILE(config.h.cmake config.h)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

#
# Top level install
#
INSTALL(FILES AUTHORS ChangeLog COPYING INSTALL NEWS README 
  DESTINATION share/doc/siconos-${VERSION}/)

#
# Library sources
#
SUBDIRS(src)

INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Siconos/Numerics Package is dedicated to low-level algorithms")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

