# Find the Path Ferris includes and libraries.
# The following variables are set if Ferris is found.  If Ferris is not
# found, PathFerris_FOUND is set to false.
#  PathFerris_FOUND        - True when the PathFerris include directory is found.
#  PathFerris_INCLUDE_DIRS - the path to where the Path Ferris include files are.
#  PathFerris_LIBRARY_DIRS - The path to where the Path library files are.
#  PathFerris_LIBRARIES    - The libraries to link against Path Ferris

# One may want to use a specific Ferris Library by setting
# PathFerris_LIBRARY_DIRECTORY before FIND_PACKAGE(PathFerris)
INCLUDE(FindPackageHandleStandardArgs)

IF(NOT PathFerris_VERSION)
  SET(PathFerris_VERSION "47")
ENDIF(NOT PathFerris_VERSION)

IF(NOT PathFerris_LIBRARY_DIRECTORY AND GAMS_FOUND)
  SET(PathFerris_LIBRARY_DIRECTORY ${GAMS_DIR})
ENDIF(NOT PathFerris_LIBRARY_DIRECTORY AND GAMS_FOUND)

IF(PathFerris_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PathFerris_LIBRARY path${PathFerris_VERSION} PATHS "${PathFerris_LIBRARY_DIRECTORY}")
ELSE(PathFerris_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PathFerris_LIBRARY path${PathFerris_VERSION})
ENDIF(PathFerris_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(PathFerris
  REQUIRED_VARS PathFerris_LIBRARY)

SET(PathFerris_FOUND ${PATHFERRIS_FOUND})

IF(PathFerris_LIBRARY)
  GET_FILENAME_COMPONENT(PathFerris_LIBRARY_DIRS ${PathFerris_LIBRARY} PATH)
  GET_FILENAME_COMPONENT(PathFerris_LIBRARIES ${PathFerris_LIBRARY} NAME)
  GET_FILENAME_COMPONENT(PathFerris_LIBRARY_DIRS_DIR ${PathFerris_LIBRARY_DIRS} PATH)
  SET(PathFerris_INCLUDE_DIRS ${PathFerris_LIBRARY_DIRS_DIR}/include)
  SET(HAVE_PATHFERRIS 1)
ELSE(PathFerris_LIBRARY)
  IF(PathFerris_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required path${PathFerris_VERSION} library not found. Please specify library location in PathFerris_LIBRARY_DIRECTORY")
  ENDIF(PathFerris_FIND_REQUIRED)

ENDIF(PathFerris_LIBRARY)
