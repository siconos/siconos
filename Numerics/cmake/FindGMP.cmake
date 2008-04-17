# Find the Siconos Numerics includes and libraries.
# The following variables are set if gmp library is found.  If gmp is not
# found, GMP_FOUND is set to false.
#  GMP_FOUND        - True when the GMP include directory is found.
#  GMP_INCLUDE_DIRS - the path to where the gmp include files are.
#  GMP_LIBRARY_DIRS - The path to where the gmp library files are.
#  GMP_LIBRARIES    - The libraries to link against gmp

# One may want to use a specific Numerics Library by setting
# GMP_LIBRARY_DIRECTORY before FIND_PACKAGE(GMP)

IF(GMP_LIBRARY_DIRECTORY)
  FIND_LIBRARY(GMP_FOUND gmp PATHS "${GMP_LIBRARY_DIRECTORY}")
ELSE(GMP_LIBRARY_DIRECTORY)
  FIND_LIBRARY(GMP_FOUND gmp)
ENDIF(GMP_LIBRARY_DIRECTORY)

IF(GMP_FOUND)
  GET_FILENAME_COMPONENT(GMP_LIBRARY_DIRS ${GMP_FOUND} PATH)
  SET(GMP_LIBRARIES ${GMP_FOUND})
  FIND_PATH(GMP_INCLUDE_DIRS NAMES gmp.h PATHS ${GMP_INCLUDE_DIR})
ELSE(GMP_FOUND)
  IF(GMP_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required gmp library not found. Please specify library location in GMP_LIBRARY_DIRECTORY")
  ENDIF(GMP_FIND_REQUIRED)
ENDIF(GMP_FOUND)
