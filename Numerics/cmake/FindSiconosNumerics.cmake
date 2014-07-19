# Find the Siconos Numerics includes and libraries.
# The following variables are set if Numerics is found.  If Numerics is not
# found, SiconosNumerics_FOUND is set to false.
#  SiconosNumerics_FOUND        - True when the SiconosNumerics include directory is found.
#  SiconosNumerics_INCLUDE_DIRS - the path to where the Siconos Numerics include files are.
#  SiconosNumerics_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosNumerics_LIBRARIES    - The libraries to link against Siconos Numerics

# One may want to use a specific Numerics Library by setting
# SiconosNumerics_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosNumerics)
INCLUDE(FindPackageHandleStandardArgs)

IF(SiconosNumerics_LIBRARY_DIRECTORY)
  MESSAGE(STATUS "Looking for Numerics library in ${SiconosNumerics_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(SiconosNumerics_LIBRARY SiconosNumerics PATHS "${SiconosNumerics_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(SiconosNumerics_LIBRARY)
    MESSAGE(STATUS "Found : ${SiconosNumerics_LIBRARY}")
  ENDIF(SiconosNumerics_LIBRARY)
ELSE(SiconosNumerics_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosNumerics_LIBRARY SiconosNumerics PATHS ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH)
ENDIF(SiconosNumerics_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SiconosNumerics 
  REQUIRED_VARS SiconosNumerics_LIBRARY)

SET(SiconosNumerics_FOUND ${SICONOSNUMERICS_FOUND})

IF(SiconosNumerics_LIBRARY)
  SET(SiconosNumerics_LIBRARIES ${SiconosNumerics_LIBRARY})
  GET_FILENAME_COMPONENT(SiconosNumerics_LIBRARY_NAME ${SiconosNumerics_LIBRARY} NAME)
  GET_FILENAME_COMPONENT(SiconosNumerics_LIBRARY_DIRS ${SiconosNumerics_LIBRARY} PATH)
  GET_FILENAME_COMPONENT(SiconosNumerics_LIBRARY_DIRS_DIR ${SiconosNumerics_LIBRARY_DIRS} PATH)
  GET_FILENAME_COMPONENT(SiconosNumerics_LIBRARY_DIRS_DIR_DIR ${SiconosNumerics_LIBRARY_DIRS_DIR} PATH)

  FIND_PATH(SiconosNumerics_INCLUDE_DIRS SiconosNumerics.h 
    HINTS ${SiconosNumerics_LIBRARY_DIRS_DIR} ${SiconosNumerics_LIBRARY_DIRS_DIR_DIR} 
    ENV PATH
    PATH_SUFFIXES include/Siconos/Numerics)
  
  IF(NOT SiconosNumerics_INCLUDE_DIRS)
    IF(SiconosNumerics_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
        "Required Siconos Numerics headers not found. Please specify headers location in CMAKE_INCLUDE_PATH")
    ENDIF(SiconosNumerics_FIND_REQUIRED)
  ENDIF(NOT SiconosNumerics_INCLUDE_DIRS)

ELSE(SiconosNumerics_LIBRARY)
  IF(SiconosNumerics_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos Numerics library not found. Please specify library location in SiconosNumerics_LIBRARY_DIRECTORY")
  ENDIF(SiconosNumerics_FIND_REQUIRED)
ENDIF(SiconosNumerics_LIBRARY)
