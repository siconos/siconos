# Find the Siconos Numerics includes and libraries.
# The following variables are set if Numerics is found.  If Numerics is not
# found, SiconosNumerics_FOUND is set to false.
#  SiconosNumerics_FOUND        - True when the SiconosNumerics include directory is found.
#  SiconosNumerics_INCLUDE_DIRS - the path to where the Siconos Numerics include files are.
#  SiconosNumerics_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosNumerics_LIBRARIES    - The libraries to link against Siconos Numerics

# One may want to use a specific Numerics Library by setting
# SiconosNumerics_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosNumerics)

IF(SiconosNumerics_LIBRARY_DIRECTORY)
  MESSAGE(STATUS "Looking for Numerics library in ${SiconosNumerics_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(SiconosNumerics_FOUND SiconosNumerics PATHS "${SiconosNumerics_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(SiconosNumerics_FOUND)
    MESSAGE(STATUS "Found : ${SiconosNumerics_FOUND}")
  ENDIF(SiconosNumerics_FOUND)
ELSE(SiconosNumerics_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosNumerics_FOUND SiconosNumerics ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH)
ENDIF(SiconosNumerics_LIBRARY_DIRECTORY)

IF(SiconosNumerics_FOUND)
  GET_FILENAME_COMPONENT(SiconosNumerics_LIBRARY_DIRS ${SiconosNumerics_FOUND} PATH)
  SET(SiconosNumerics_LIBRARIES ${SiconosNumerics_FOUND})
  GET_FILENAME_COMPONENT(SiconosNumerics_LIBRARY_DIRS_DIR ${SiconosNumerics_LIBRARY_DIRS} PATH)
  SET(SiconosNumerics_INCLUDE_DIRS ${SiconosNumerics_LIBRARY_DIRS_DIR}/include/Siconos/Numerics)
ELSE(SiconosNumerics_FOUND)
  IF(SiconosNumerics_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos Numerics library not found. Please specify library location in SiconosNumerics_LIBRARY_DIRECTORY")
  ENDIF(SiconosNumerics_FIND_REQUIRED)
ENDIF(SiconosNumerics_FOUND)
