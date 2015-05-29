# Find the Siconos IO includes and libraries.
# The following variables are set if IO is found.  If IO is not
# found, SiconosIO_FOUND is set to false.
#  SiconosIO_FOUND        - True when the SiconosIO include directory is found.
#  SiconosIO_INCLUDE_DIRS - the path to where the Siconos IO include files are.
#  SiconosIO_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosIO_LIBRARIES    - The libraries to link against Siconos IO

# One may want to use a specific IO Library by setting
# SiconosIO_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosIO)
INCLUDE(FindPackageHandleStandardArgs)

IF(SiconosIO_LIBRARY_DIRECTORY)
  MESSAGE(STATUS "Looking for IO library in ${SiconosIO_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(SiconosIO_LIBRARY SiconosIO PATHS "${SiconosIO_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(SiconosIO_LIBRARY)
    MESSAGE(STATUS "Found : ${SiconosIO_LIBRARY}")
  ENDIF(SiconosIO_LIBRARY)
ELSE(SiconosIO_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosIO_LIBRARY SiconosIO)
ENDIF(SiconosIO_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SiconosIO
  REQUIRED_VARS SiconosIO_LIBRARY)

SET(SiconosIO_FOUND ${SICONOSIO_FOUND})

IF(SiconosIO_LIBRARY)
  GET_FILENAME_COMPONENT(SiconosIO_LIBRARY_DIRS ${SiconosIO_LIBRARY} PATH)
  SET(SiconosIO_LIBRARIES ${SiconosIO_LIBRARY})
  GET_FILENAME_COMPONENT(SiconosIO_LIBRARY_DIRS_DIR ${SiconosIO_LIBRARY_DIRS} PATH)

   FIND_PATH(SiconosIO_INCLUDE_DIRS MechanicsIO.hpp
    HINTS ${SiconosIO_LIBRARY_DIRS_DIR} ${SiconosIO_LIBRARY_DIRS_DIR_DIR} 
    ENV PATH
    PATH_SUFFIXES include/Siconos/IO)

  IF(NOT SiconosIO_INCLUDE_DIRS)
    IF(SiconosIO_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
        "Required Siconos IO headers not found. Please specify headers location in CMAKE_INCLUDE_PATH")
    ENDIF(SiconosIO_FIND_REQUIRED)
  ENDIF(NOT SiconosIO_INCLUDE_DIRS)

ELSE(SiconosIO_LIBRARY)
  IF(SiconosIO_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos IO library not found. Please specify library location in SiconosIO_LIBRARY_DIRECTORY")
  ENDIF(SiconosIO_FIND_REQUIRED)
ENDIF(SiconosIO_LIBRARY)
