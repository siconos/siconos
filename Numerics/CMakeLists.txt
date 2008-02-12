#
# Siconos Numerics : general configuration
#
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.4)

# Package name
PROJECT(siconos-numerics C CXX Fortran) 

#
# An encourage to out of source builds 
#
OPTION(IN_SOURCE_BUILD "if you really want a build in the source directory" OFF)

IF(NOT IN_SOURCE_BUILD)
  IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    MESSAGE(FATAL_ERROR "CMake generation for Siconos Numerics is not allowed within the source directory! 
Remove the CMakeCache.txt file and try again from another folder, e.g.: 

   rm CMakeCache.txt 
   cd <somewhere (preferably a local place on your computer and not a network folder)>
   cmake <source directory of Numerics>

If you really need an in source build, then run : cmake -DIN_SOURCE_BUILD=ON
")
  ENDIF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
ENDIF(NOT IN_SOURCE_BUILD)

#
# Version
#
SET(MAJOR_VERSION 2)
SET(MINOR_VERSION 2)
SET(PATCH_VERSION 1)

SET(VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")

#
# SVN revision number
#
SET( ENV{LC_ALL} C )
IF(EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  FIND_PACKAGE(Subversion)
  IF(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
    IF(Project_WC_REVISION)
      SET(SVN_REVISION ${Project_WC_REVISION})
      MESSAGE(STATUS "This is a build from sources under svn")
      MESSAGE(STATUS "Current svn revision is ${SVN_REVISION}")
    ENDIF(Project_WC_REVISION)
  ENDIF(Subversion_FOUND)
ELSE(EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  MESSAGE(STATUS "This is a build out of svn")
ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/.svn)

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
INCLUDE(FortranLibraries)
INCLUDE(LibraryProjectSetup)


#
# Documentation with Doxygen 
#

# Style sheet
FILE(GLOB CSS_FILES config/css/*.* config/css/*)
FOREACH(FILE ${CSS_FILES})
  IF(NOT IS_DIRECTORY ${FILE})
    FILE(RELATIVE_PATH RFILE ${CMAKE_SOURCE_DIR} ${FILE})
    CONFIGURE_FILE(${FILE} ${CMAKE_BINARY_DIR}/${RFILE} COPYONLY)
  ENDIF(NOT IS_DIRECTORY ${FILE})
ENDFOREACH(FILE ${CSS_FILES})

# If we are under svn the output should go in GeneratedDocs 
IF(SVN_REVISION)
  SET(DOXYGEN_OUTPUT ${CMAKE_SOURCE_DIR}/../GeneratedDocs)
ELSE(SVN_REVISION)
  SET(DOXYGEN_OUTPUT html)
ENDIF(SVN_REVISION)

SET(PROJECT_SHORT_NAME Numerics)

INCLUDE(TargetDoc OPTIONAL)


#
# Required packages 
#
FIND_PACKAGE(LAPACK REQUIRED)
#
# Current state information
#
MESSAGE(STATUS "BLAS_LINKER_FLAGS = ${BLAS_LINKER_FLAGS}")
MESSAGE(STATUS "BLAS_LIBRARIES = ${BLAS_LIBRARIES}")

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
  MESSAGE(STATUS "ATLAS library found")
  FIND_PATH(ATLAS_INCLUDE_PATH 
    NAMES atlas/cblas.h atlas/clapack.h)
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
    INCLUDE_DIRECTORIES(${ATLAS_INCLUDE_PATH}/atlas)
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
  DESTINATION share/doc/siconos-${VERSION}/Numerics)

#
# Library sources
#
ENABLE_TESTING()
INCLUDE(DartConfig)
INCLUDE(Dart)

SUBDIRS(src)


#
# Packaging
#
INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Siconos/Numerics Package is dedicated to low-level algorithms")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
