# Find the external_mlcp_simplex includes and libraries.
# The following variables are set if libMlcpSimplex is found.  If libMlcpSimplex is not
# found, MlcpSimplex_FOUND is set to false.
#  MlcpSimplex_FOUND        - True when the MlcpSimplex include directory is found.
#  MlcpSimplex_INCLUDE_DIRS - the path to where the Path Ferris include files are.
#  MlcpSimplex_LIBRARY_DIRS - The path to where the Path library files are.
#  MlcpSimplex_LIBRARIES    - The libraries to link against Path Ferris

# One may want to use a specific Ferris Library by setting
# MlcpSimplex_LIBRARY_DIRECTORY before FIND_PACKAGE(MlcpSimplex)
INCLUDE(FindPackageHandleStandardArgs)

IF(MlcpSimplex_LIBRARY_DIRECTORY)
  FIND_LIBRARY(MlcpSimplex_LIBRARY MlcpSimplex PATHS "${MlcpSimplex_LIBRARY_DIRECTORY}")
ELSE(MlcpSimplex_LIBRARY_DIRECTORY)
  FIND_LIBRARY(MlcpSimplex_LIBRARY MlcpSimplex)
ENDIF(MlcpSimplex_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MlcpSimplex
  REQUIRED_VARS MlcpSimplex_LIBRARY)

SET(MlcpSimplex_FOUND ${MLCPSIMPLEX_FOUND})

IF(MlcpSimplex_LIBRARY)
  GET_FILENAME_COMPONENT(MlcpSimplex_LIBRARY_DIRS ${MlcpSimplex_LIBRARY} PATH)
  GET_FILENAME_COMPONENT(MlcpSimplex_LIBRARIES ${MlcpSimplex_LIBRARY} NAME)
  GET_FILENAME_COMPONENT(MlcpSimplex_LIBRARY_DIRS_DIR ${MlcpSimplex_LIBRARY_DIRS} PATH)
  SET(MlcpSimplex_INCLUDE_DIRS ${MlcpSimplex_LIBRARY_DIRS_DIR}/include)
  SET(HAVE_MLCPSIMPLEX 1)
ELSE(MlcpSimplex_LIBRARY)
  IF(MlcpSimplex_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required MlcpSimplex library not found. Please specify library location in MlcpSimplex_LIBRARY_DIRECTORY")
  ENDIF(MlcpSimplex_FIND_REQUIRED)
  
ENDIF(MlcpSimplex_LIBRARY)
