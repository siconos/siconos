# Find external cplex includes and libraries.
# The following variables are set if libcplex is found.  If libcplex is not
# found, Cplex_FOUND is set to false.
#  Cplex_FOUND        - True when the Cplex include directory is found.
#  Cplex_INCLUDE_DIRS - the path to where the Cplex include files are.
#  Cplex_LIBRARY_DIRS - The path to where the Cplex library files are.
#  Cplex_LIBRARIES    - The libraries to link against Cplex

# One may want to use a specific Ferris Library by setting
# Cplex_LIBRARY_DIRECTORY before FIND_PACKAGE(Cplex)

IF(Cplex_LIBRARY_DIRECTORY)
  FIND_LIBRARY(Cplex_FOUND cplex PATHS "${Cplex_LIBRARY_DIRECTORY}")
ELSE(Cplex_LIBRARY_DIRECTORY)
  FIND_LIBRARY(Cplex_FOUND cplex)
ENDIF(Cplex_LIBRARY_DIRECTORY)

IF(Cplex_FOUND)
  GET_FILENAME_COMPONENT(Cplex_LIBRARY_DIRS ${Cplex_FOUND} PATH)
  GET_FILENAME_COMPONENT(Cplex_LIBRARIES ${Cplex_FOUND} NAME)
  GET_FILENAME_COMPONENT(Cplex_LIBRARY_DIRS_DIR ${Cplex_LIBRARY_DIRS} PATH)
  SET(Cplex_INCLUDE_DIRS ${Cplex_LIBRARY_DIRS_DIR}/include)
ELSE(Cplex_FOUND)
  IF(Cplex_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Cplex library not found. Please specify library location in Cplex_LIBRARY_DIRECTORY")
  ENDIF(Cplex_FIND_REQUIRED)

ENDIF(Cplex_FOUND)
