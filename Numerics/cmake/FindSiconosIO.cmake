# Find the Siconos IO includes and libraries.
# The following variables are set if IO is found.  If IO is not
# found, SiconosIO_FOUND is set to false.
#  SiconosIO_FOUND        - True when the SiconosIO include directory is found.
#  SiconosIO_INCLUDE_DIRS - the path to where the Siconos IO include files are.
#  SiconosIO_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosIO_LIBRARIES    - The libraries to link against Siconos IO

# One may want to use a specific IO Library by setting
# SiconosIO_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosIO)

IF(SiconosIO_LIBRARY_DIRECTORY)
  MESSAGE(STATUS "Looking for IO library in ${SiconosIO_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(SiconosIO_FOUND SiconosIO PATHS "${SiconosIO_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(SiconosIO_FOUND)
    MESSAGE(STATUS "Found : ${SiconosIO_FOUND}")
  ENDIF(SiconosIO_FOUND)
ELSE(SiconosIO_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosIO_FOUND SiconosIO)
ENDIF(SiconosIO_LIBRARY_DIRECTORY)

IF(SiconosIO_FOUND)
  GET_FILENAME_COMPONENT(SiconosIO_LIBRARY_DIRS ${SiconosIO_FOUND} PATH)
  SET(SiconosIO_LIBRARIES ${SiconosIO_FOUND})
  GET_FILENAME_COMPONENT(SiconosIO_LIBRARY_DIRS_DIR ${SiconosIO_LIBRARY_DIRS} PATH)
  SET(SiconosIO_INCLUDE_DIRS ${SiconosIO_LIBRARY_DIRS_DIR}/include/Siconos/IO)
ELSE(SiconosIO_FOUND)
  IF(SiconosIO_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos IO library not found. Please specify library location in SiconosIO_LIBRARY_DIRECTORY")
  ENDIF(SiconosIO_FIND_REQUIRED)
ENDIF(SiconosIO_FOUND)
