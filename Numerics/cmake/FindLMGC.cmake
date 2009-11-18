# Find the LMGC includes and libraries.
# The following variables are set if LMGC is found.  If LMGC is not
# found, LMGC_FOUND is set to false.
#  LMGC_FOUND        - True when the LMGC include directory is found.
#  LMGC_INCLUDE_DIRS - the path to where the Siconos LMGC include files are.
#  LMGC_LIBRARY_DIRS - The path to where the Siconos library files are.
#  LMGC_LIBRARIES    - The libraries to link against Siconos LMGC

# One may want to use a specific LMGC Library by setting
# LMGC_LIBRARY_DIRECTORY before FIND_PACKAGE(LMGC)

IF(LMGC_LIBRARY_DIRECTORY)
  FIND_LIBRARY(LMGC_FOUND lmgc PATHS "${LMGC_LIBRARY_DIRECTORY}")
ELSE(LMGC_LIBRARY_DIRECTORY)
  FIND_LIBRARY(LMGC_FOUND lmgc)
ENDIF(LMGC_LIBRARY_DIRECTORY)

IF(LMGC_FOUND)
  GET_FILENAME_COMPONENT(LMGC_LIBRARY_DIRS ${LMGC_FOUND} PATH)
  SET(LMGC_LIBRARIES ${LMGC_FOUND})
  GET_FILENAME_COMPONENT(LMGC_LIBRARY_DIRS_DIR ${LMGC_LIBRARY_DIRS} PATH)
  SET(LMGC_INCLUDE_DIRS ${LMGC_LIBRARY_DIRS_DIR}/include/LMGC)
ELSE(LMGC_FOUND)
  IF(LMGC_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required LMGC library not found. Please specify library location in LMGC_LIBRARY_DIRECTORY")
  ENDIF(LMGC_FIND_REQUIRED)
ENDIF(LMGC_FOUND)
