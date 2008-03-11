# Find the Path Ferris includes and libraries.
# The following variables are set if Ferris is found.  If Ferris is not
# found, PathFerris_FOUND is set to false.
#  PathFerris_FOUND        - True when the PathFerris include directory is found.
#  PathFerris_INCLUDE_DIRS - the path to where the Path Ferris include files are.
#  PathFerris_LIBRARY_DIRS - The path to where the Path library files are.
#  PathFerris_LIBRARIES    - The libraries to link against Path Ferris

# One may want to use a specific Ferris Library by setting
# PathFerris_LIBRARY_DIRECTORY before FIND_PACKAGE(PathFerris)

IF(PathFerris_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PathFerris_FOUND PathFerris PATHS "${PathFerris_LIBRARY_DIRECTORY}")
ELSE(PathFerris_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PathFerris_FOUND PathFerris)
ENDIF(PathFerris_LIBRARY_DIRECTORY)

IF(PathFerris_FOUND)
  GET_FILENAME_COMPONENT(PathFerris_LIBRARY_DIRS ${PathFerris_FOUND} PATH)
  GET_FILENAME_COMPONENT(PathFerris_LIBRARIES ${PathFerris_FOUND} NAME)
  GET_FILENAME_COMPONENT(PathFerris_LIBRARY_DIRS_DIR ${PathFerris_LIBRARY_DIRS} PATH)
  SET(PathFerris_INCLUDE_DIRS ${PathFerris_LIBRARY_DIRS_DIR}/include)
  SET(HAVE_PATHFERRIS 1)
ELSE(PathFerris_FOUND)
  IF(PathFerris_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos Numerics library not found. Please specify library location in PathFerris_LIBRARY_DIRECTORY")
  ENDIF(PathFerris_REQUIRED)
  

ENDIF(PathFerris_FOUND)
